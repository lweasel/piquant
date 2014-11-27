import piquant.parameters as parameters
import pytest
import schema
import tempfile


def _get_test_parameter(
        name="name", title="The Name",
        option_name="--option-name", option_validator=int,
        is_numeric=False, value_namer=None, file_namer=None):

    param = parameters._Parameter(
        name, title, option_name, option_validator, is_numeric,
        value_namer, file_namer)

    parameters._PARAMETERS.remove(param)
    parameters._RUN_PARAMETERS.remove(param)

    return param


def _get_ignore_params():
    return [
        parameters.READ_DEPTH,
        parameters.PAIRED_END,
        parameters.ERRORS,
        parameters.BIAS,
        parameters.STRANDED,
        parameters.NOISE_DEPTH_PERCENT,
        parameters.TRANSCRIPT_GTF,
        parameters.NOISE_TRANSCRIPT_GTF,
        parameters.GENOME_FASTA_DIR,
        parameters.NUM_MOLECULES,
        parameters.NUM_NOISE_MOLECULES,
        parameters.NUM_THREADS
    ]


def test_get_run_parameters_returns_parameters_instances():
    params = parameters.get_run_parameters()
    assert all([isinstance(p, parameters._Parameter) for p in params])


def test_parameter_name_is_correct():
    name = "parameter name"
    p = _get_test_parameter(name=name)
    assert p.name == name


def test_parameter_title_is_correct():
    title = "A Title"
    p = _get_test_parameter(title=title)
    assert p.title == title


def test_parameter_option_name_is_correct():
    option_name = "--polyA"
    p = _get_test_parameter(option_name=option_name)
    assert p.option_name == option_name


def test_parameter_option_validator_is_correct():
    option_validator = lambda x: x
    p = _get_test_parameter(option_validator=option_validator)
    assert p.option_validator == option_validator


def test_parameter_is_numeric_is_correct():
    is_numeric = True
    p = _get_test_parameter(is_numeric=is_numeric)
    assert p.is_numeric == is_numeric


def test_get_value_name_returns_correct_value_when_no_value_namer_supplied():
    value = 30
    p = _get_test_parameter()
    assert p.get_value_name(value) == value


def test_get_value_name_returns_correct_value_when_value_namer_supplied():
    value = 30
    p = _get_test_parameter(value_namer=lambda x: "VAL{v}".format(v=x))
    assert p.get_value_name(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_no_value_or_file_namer_supplied():
    value = 30
    p = _get_test_parameter()
    assert p.get_file_name_part(value) == value


def test_get_file_name_part_returns_correct_value_when_no_file_namer_supplied():
    value = 30
    p = _get_test_parameter(value_namer=lambda x: "VAL{v}".format(v=x))
    assert p.get_file_name_part(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_file_namer_supplied():
    value = 30
    p = _get_test_parameter(file_namer=lambda x: "{v}bp".format(v=x))
    assert p.get_file_name_part(value) == str(value) + "bp"


def test_validate_command_line_parameter_sets_returns_correct_number_of_param_sets():
    options = {
        "--quant-method": "Cufflinks",
        "--read-length": "10,20",
    }

    ignore_params = _get_ignore_params()

    param_vals = parameters.validate_command_line_parameter_sets(
        None, options, ignore_params)
    assert len(param_vals) == 2


def test_validate_command_line_parameter_sets_returns_correct_number_of_transformed_values():
    num_values = 5
    options = {
        "--read-length": ",".join([str(i) for i in range(0, num_values)])
    }

    ignore_params = _get_ignore_params()
    ignore_params.append(parameters.QUANT_METHOD)

    param_vals = parameters.validate_command_line_parameter_sets(
        None, options, ignore_params)
    assert len(param_vals[parameters.READ_LENGTH.name]) == num_values


def test_validate_command_line_parameter_sets_returns_correct_transformed_values():
    len1 = 10
    len2 = 20
    options = {
        "--read-length": str(len1) + "," + str(len2)
    }

    ignore_params = _get_ignore_params()
    ignore_params.append(parameters.QUANT_METHOD)

    param_vals = parameters.validate_command_line_parameter_sets(
        None, options, ignore_params)
    assert len1 in param_vals[parameters.READ_LENGTH.name]
    assert len2 in param_vals[parameters.READ_LENGTH.name]


def test_validate_command_line_parameter_sets_raises_exception_for_invalid_param_values():
    options = {
        "--read-length": "abc"
    }
    with pytest.raises(schema.SchemaError):
        parameters.validate_command_line_parameter_sets(None, options)


def test_validate_command_line_parameter_sets_raises_exception_if_param_values_not_supplied():
    options = {
        "--read-length": "10,20",
    }

    ignore_params = _get_ignore_params()

    with pytest.raises(schema.SchemaError):
        parameters.validate_command_line_parameter_sets(
            None, options, ignore_params)


def test_validate_command_line_parameter_sets_reads_parameters_from_file():
    len1 = 10
    len2 = 20
    with tempfile.NamedTemporaryFile() as f:
        f.write("--quant-method Cufflinks\n")
        f.write("--read-length " + str(len1) + "," + str(len2) + "\n")
        f.flush()

        ignore_params = _get_ignore_params()

        param_vals = parameters.validate_command_line_parameter_sets(
            f.name, {}, ignore_params)

        assert len1 in param_vals[parameters.READ_LENGTH.name]
        assert len2 in param_vals[parameters.READ_LENGTH.name]


def test_validate_command_line_parameter_sets_overrides_file_parameters_with_cl_options():
    len1 = 10
    len2 = 20
    len3 = 30
    with tempfile.NamedTemporaryFile() as f:
        f.write("--read-length " + str(len1) + "," + str(len2))
        f.flush()

        options = {
            "--read-length": str(len3)
        }

        ignore_params = _get_ignore_params()
        ignore_params.append(parameters.QUANT_METHOD)

        param_vals = parameters.validate_command_line_parameter_sets(
            f.name, options, ignore_params)

        assert len1 not in param_vals[parameters.READ_LENGTH.name]
        assert len2 not in param_vals[parameters.READ_LENGTH.name]
        assert len3 in param_vals[parameters.READ_LENGTH.name]


def test_get_file_name_returns_correct_name():
    assert parameters.get_file_name(
        read_depth=30, read_length=50,
        paired_end=True, bias=False) == "30x_50b_pe_no_bias"


def test_execute_for_param_sets_executes_for_correct_number_of_parameter_sets():
    len1 = 3
    len2 = 5
    len3 = 7
    params_values = {
        "read_length": [1]*len1,
        "read_depth": [2]*len2,
        "errors": [True]*len3
    }

    execute_counter = []

    def count_incrementer(logger, options, **params):
        execute_counter.append(1)

    parameters.execute_for_param_sets(
        [count_incrementer], None, None, **params_values)

    assert len(execute_counter) == len1 * len2 * len3


def test_execute_for_param_sets_executes_all_callables():
    execute_record = []

    def callable1(logger, options, **params):
        execute_record.append(1)

    def callable2(logger, options, **params):
        execute_record.append(2)

    parameters.execute_for_param_sets(
        [callable1, callable2], None, None, param=["a"])

    assert 1 in execute_record
    assert 2 in execute_record
    assert len(execute_record) == 2


def test_execute_for_param_sets_executes_for_correct_sets_of_parameters():
    params1 = [75, 100]
    params2 = [30, 50]
    params_values = {
        "read_length": params1,
        "read_depth": params2
    }

    execute_record = []

    def callable1(logger, options, **params):
        execute_record.append([v for v in params.values()])

    parameters.execute_for_param_sets([callable1], None, None, **params_values)

    execute_record = [set(params) for params in execute_record]
    assert set([params1[0], params2[0]]) in execute_record
    assert set([params1[0], params2[1]]) in execute_record
    assert set([params1[1], params2[0]]) in execute_record
    assert set([params1[1], params2[1]]) in execute_record
