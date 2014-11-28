import piquant.piquant_options as po
import pytest
import schema
import tempfile


def _get_test_multi_quant_run_option(
        name="name", title="The Name",
        description="description for the option",
        option_validator=int,
        is_numeric=False, value_namer=None, file_namer=None):

    option = po._PiquantOption(
        name, description, title=title,
        option_validator=option_validator, is_numeric=is_numeric,
        value_namer=value_namer, file_namer=file_namer,
        option_type=po._PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

    po._PiquantOption._QUANT_RUN_OPTIONS.remove(option)
    po._PiquantOption._MULTI_QUANT_RUN_OPTIONS.remove(option)

    return option


def _get_ignore_options():
    return [
        po.READ_DEPTH,
        po.PAIRED_END,
        po.ERRORS,
        po.BIAS,
        po.STRANDED,
        po.NOISE_DEPTH_PERCENT,
        po.TRANSCRIPT_GTF,
        po.NOISE_TRANSCRIPT_GTF,
        po.GENOME_FASTA_DIR,
        po.NUM_MOLECULES,
        po.NUM_NOISE_MOLECULES,
        po.NUM_THREADS
    ]


def test_get_multi_quant_run_options_returns_piquant_option_instances():
    opts = po.get_multiple_quant_run_options()
    assert all([isinstance(o, po._PiquantOption) for o in opts])


def test_mqr_option_name_is_correct():
    name = "mqr_option_name"
    o = _get_test_multi_quant_run_option(name=name)
    assert o.name == name


def test_mqr_option_title_is_correct():
    title = "A Title"
    o = _get_test_multi_quant_run_option(title=title)
    assert o.title == title


def test_mqr_option_option_name_is_correct():
    name = "mqr_option_name"
    o = _get_test_multi_quant_run_option(name=name)
    assert o.get_option_name() == "--" + name.replace("_", "-")


def test_mqr_option_option_validator_is_correct():
    option_validator = lambda x: x
    o = _get_test_multi_quant_run_option(option_validator=option_validator)
    assert o.option_validator == option_validator


def test_mqr_option_is_numeric_is_correct():
    is_numeric = True
    o = _get_test_multi_quant_run_option(is_numeric=is_numeric)
    assert o.is_numeric == is_numeric


def test_get_value_name_returns_correct_value_when_no_value_namer_supplied():
    value = 30
    o = _get_test_multi_quant_run_option()
    assert o.get_value_name(value) == value


def test_get_value_name_returns_correct_value_when_value_namer_supplied():
    value = 30
    o = _get_test_multi_quant_run_option(value_namer=lambda x: "VAL{v}".format(v=x))
    assert o.get_value_name(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_no_value_or_file_namer_supplied():
    value = 30
    o = _get_test_multi_quant_run_option()
    assert o.get_file_name_part(value) == value


def test_get_file_name_part_returns_correct_value_when_no_file_namer_supplied():
    value = 30
    o = _get_test_multi_quant_run_option(value_namer=lambda x: "VAL{v}".format(v=x))
    assert o.get_file_name_part(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_file_namer_supplied():
    value = 30
    o = _get_test_multi_quant_run_option(file_namer=lambda x: "{v}bp".format(v=x))
    assert o.get_file_name_part(value) == str(value) + "bp"


#def test_validate_command_line_parameter_sets_returns_correct_number_of_param_sets():
    #options = {
        #"--quant-method": "Cufflinks",
        #"--read-length": "10,20",
    #}
#
    #ignore_options = _get_ignore_options()
#
    #param_vals = parameters.validate_command_line_parameter_sets(
        #None, options, ignore_options)
    #assert len(param_vals) == 2


#def test_validate_command_line_parameter_sets_returns_correct_number_of_transformed_values():
    #num_values = 5
    #options = {
        #"--read-length": ",".join([str(i) for i in range(0, num_values)])
    #}

    #ignore_options = _get_ignore_options()
    #ignore_options.append(parameters.QUANT_METHOD)
#
    #param_vals = parameters.validate_command_line_parameter_sets(
        #None, options, ignore_options)
    #assert len(param_vals[parameters.READ_LENGTH.name]) == num_values


#def test_validate_command_line_parameter_sets_returns_correct_transformed_values():
    #len1 = 10
    #len2 = 20
    #options = {
        #"--read-length": str(len1) + "," + str(len2)
    #}
#
    #ignore_options = _get_ignore_options()
    #ignore_options.append(parameters.QUANT_METHOD)
#
    #param_vals = parameters.validate_command_line_parameter_sets(
        #None, options, ignore_options)
    #assert len1 in param_vals[parameters.READ_LENGTH.name]
    #assert len2 in param_vals[parameters.READ_LENGTH.name]


#def test_validate_command_line_parameter_sets_raises_exception_for_invalid_param_values():
    #options = {
        #"--read-length": "abc"
    #}
    #with pytest.raises(schema.SchemaError):
        #parameters.validate_command_line_parameter_sets(None, options)

#def test_validate_command_line_parameter_sets_raises_exception_if_param_values_not_supplied():
    #options = {
        #"--read-length": "10,20",
    #}
#
    #ignore_options = _get_ignore_options()
#
    #with pytest.raises(schema.SchemaError):
        #parameters.validate_command_line_parameter_sets(
            #None, options, ignore_options)


#def test_validate_command_line_parameter_sets_reads_parameters_from_file():
    #len1 = 10
    #len2 = 20
    #with tempfile.NamedTemporaryFile() as f:
        #f.write("--quant-method Cufflinks\n")
        #f.write("--read-length " + str(len1) + "," + str(len2) + "\n")
        #f.flush()
#
        #ignore_options = _get_ignore_options()
#
        #param_vals = parameters.validate_command_line_parameter_sets(
            #f.name, {}, ignore_options)
#
        #assert len1 in param_vals[parameters.READ_LENGTH.name]
        #assert len2 in param_vals[parameters.READ_LENGTH.name]


#def test_validate_command_line_parameter_sets_overrides_file_parameters_with_cl_options():
    #len1 = 10
    #len2 = 20
    #len3 = 30
    #with tempfile.NamedTemporaryFile() as f:
        #f.write("--read-length " + str(len1) + "," + str(len2))
        #f.flush()
#
        #options = {
            #"--read-length": str(len3)
        #}
#
        #ignore_options = _get_ignore_options()
        #ignore_options.append(parameters.QUANT_METHOD)
#
        #param_vals = parameters.validate_command_line_parameter_sets(
            #f.name, options, ignore_options)
#
        #assert len1 not in param_vals[parameters.READ_LENGTH.name]
        #assert len2 not in param_vals[parameters.READ_LENGTH.name]
        #assert len3 in param_vals[parameters.READ_LENGTH.name]


def test_get_file_name_returns_correct_name():
    assert po.get_file_name(
        read_depth=30, read_length=50,
        paired_end=True, bias=False) == "30x_50b_pe_no_bias"


def test_execute_for_mqr_option_sets_executes_for_correct_number_of_option_sets():
    len1 = 3
    len2 = 5
    len3 = 7
    mqr_options_values = {
        "read_length": [1] * len1,
        "read_depth": [2] * len2,
        "errors": [True] * len3
    }

    execute_counter = []

    def count_incrementer(logger, options, **mqr_options):
        execute_counter.append(1)

    po.execute_for_mqr_option_sets(
        [count_incrementer], None, None, **mqr_options_values)

    assert len(execute_counter) == len1 * len2 * len3


def test_execute_for_mqr_option_sets_executes_all_callables():
    execute_record = []

    def callable1(logger, options, **mqr_options):
        execute_record.append(1)

    def callable2(logger, options, **mqr_options):
        execute_record.append(2)

    po.execute_for_mqr_option_sets(
        [callable1, callable2], None, None, mqr_option=["a"])

    assert 1 in execute_record
    assert 2 in execute_record
    assert len(execute_record) == 2


def test_execute_for_mqr_option_sets_executes_for_correct_sets_of_mqr_options():
    mqr_options1 = [75, 100]
    mqr_options2 = [30, 50]
    mqr_options_values = {
        "read_length": mqr_options1,
        "read_depth": mqr_options2
    }

    execute_record = []

    def callable1(logger, options, **mqr_options):
        execute_record.append([v for v in mqr_options.values()])

    po.execute_for_mqr_option_sets([callable1], None, None, **mqr_options_values)

    execute_record = [set(mqr_options) for mqr_options in execute_record]
    assert set([mqr_options1[0], mqr_options2[0]]) in execute_record
    assert set([mqr_options1[0], mqr_options2[1]]) in execute_record
    assert set([mqr_options1[1], mqr_options2[0]]) in execute_record
    assert set([mqr_options1[1], mqr_options2[1]]) in execute_record
