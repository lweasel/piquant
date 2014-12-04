import piquant.log as log
import piquant.piquant_options as po
import piquant.piquant_commands as pc
#import pytest
#import schema
import tempfile
import sys


def _get_logger():
    return log.get_logger(sys.stderr, "critical")


def _get_test_option_value(
        default_value=None, option_validator=int,
        value_namer=lambda x: x, file_namer=None):

    return po._OptionValue(default_value, option_validator,
                           value_namer, file_namer)


def _get_test_option(
        name="name", title="The Name",
        description="description for the option",
        is_numeric=False, option_value=_get_test_option_value()):

    option = po._PiquantOption(
        name, description, title=title, option_value=option_value,
        is_numeric=is_numeric)

    po._PiquantOption._OPTIONS.remove(option)
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


def test_piquant_option_name_is_correct():
    name = "piquant_option_name"
    opt = _get_test_option(name=name)
    assert opt.name == name


def test_piquant_option_description_is_correct():
    description = "This is an option"
    opt = _get_test_option(description=description)
    assert opt.description == description


def test_piquant_option_default_value_is_correct():
    default_value = "DEF"
    opt = _get_test_option(
        option_value=_get_test_option_value(default_value=default_value))
    assert opt.default_value() == default_value


def test_piquant_option_option_validator_is_correct():
    option_validator = lambda x: x
    opt = _get_test_option(
        option_value=_get_test_option_value(option_validator=option_validator))
    assert opt.option_value.validator == option_validator


def test_piquant_option_has_value_is_correct():
    opt = _get_test_option(option_value=None)
    assert not opt.has_value()
    opt = _get_test_option(option_value=_get_test_option_value())
    assert opt.has_value()


def test_piquant_option_title_is_correct():
    title = "A Title"
    opt = _get_test_option(title=title)
    assert opt.title == title


def test_piquant_option_is_numeric_is_correct():
    is_numeric = True
    opt = _get_test_option(is_numeric=is_numeric)
    assert opt.is_numeric == is_numeric


def test_piquant_option_value_namer_is_correct():
    value_namer = lambda x: x + "bp"
    opt = _get_test_option(
        option_value=_get_test_option_value(value_namer=value_namer))
    assert opt.option_value.value_namer == value_namer


def test_piquant_option_file_namer_is_correct():
    file_namer = lambda x: x + "bp"
    opt = _get_test_option(
        option_value=_get_test_option_value(file_namer=file_namer))
    assert opt.option_value.file_namer == file_namer


def test_get_usage_string_returns_correct_string_when_option_has_value():
    opt = _get_test_option(name="the_opt")
    assert opt.get_usage_string() == "--the-opt=<the-opt>"


def test_get_usage_string_returns_correct_string_when_option_does_not_have_value():
    opt = _get_test_option(name="the_opt", option_value=None)
    assert opt.get_usage_string() == "--the-opt"


def test_get_option_name_returns_correct_value_when_dashes_included():
    opt = _get_test_option(name="the_opt")
    assert opt.get_option_name() == "--the-opt"


def test_get_option_name_returns_correct_value_when_dashes_not_included():
    opt = _get_test_option(name="the_opt")
    assert opt.get_option_name(include_dashes=False) == "the-opt"


def test_get_option_description_includes_usage_string():
    opt = _get_test_option()
    assert opt.get_option_description().find(opt.get_usage_string()) > -1


def test_get_option_description_includes_description():
    opt = _get_test_option()
    assert opt.get_option_description().find(opt.description) > -1


def test_get_option_description_includes_default_value_specification_when_default_specified():
    default_value = "DEF"
    opt = _get_test_option(
        option_value=_get_test_option_value(default_value=default_value))
    assert opt.get_option_description().find("default:") > -1
    assert opt.get_option_description().find(default_value) > -1


def test_get_option_description_does_not_include_default_value_when_default_value_not_specified():
    opt = _get_test_option()
    assert opt.get_option_description().find("default") == -1


def test_get_value_name_returns_correct_value_when_no_value_namer_supplied():
    value = 30
    opt = _get_test_option()
    assert opt.get_value_name(value) == value


def test_get_value_name_returns_correct_value_when_value_namer_supplied():
    value = 30
    opt = _get_test_option(
        option_value=_get_test_option_value(
            value_namer=lambda x: "VAL{v}".format(v=x)))
    assert opt.get_value_name(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_no_value_or_file_namer_supplied():
    value = 30
    opt = _get_test_option()
    assert opt.get_file_name_part(value) == value


def test_get_file_name_part_returns_correct_value_when_no_file_namer_supplied():
    value = 30
    opt = _get_test_option(
        option_value=_get_test_option_value(
            value_namer=lambda x: "VAL{v}".format(v=x)))
    assert opt.get_file_name_part(value) == "VAL" + str(value)


def test_get_file_name_part_returns_correct_value_when_file_namer_supplied():
    value = 30
    opt = _get_test_option(
        option_value=_get_test_option_value(
            file_namer=lambda x: "{v}bp".format(v=x)))
    assert opt.get_file_name_part(value) == str(value) + "bp"


def test_get_multi_quant_run_options_returns_piquant_option_instances():
    opts = po.get_multiple_quant_run_options()
    assert all([isinstance(o, po._PiquantOption) for o in opts])


def test_validate_option_values_returns_correct_number_of_options():
    options = {
        po.PLOT_FORMAT: "svg",
        po.QUANT_METHOD: "Cufflinks",
        po.READ_LENGTH: "10,20",
    }

    opt_vals, qr_opt_vals = po._validate_option_values(
        _get_logger(), {o.get_option_name(): v for o, v in options.items()},
        {}, options.keys())
    assert len(opt_vals) == 1
    assert len(qr_opt_vals) == 2


def test_validate_option_values_returns_correct_number_of_transformed_values():
    num_values = 5
    options = {
        po.READ_LENGTH: ",".join([str(i) for i in range(1, num_values + 1)])
    }

    opt_vals, qr_opt_vals = po._validate_option_values(
        _get_logger(), {o.get_option_name(): v for o, v in options.items()},
        {}, options.keys())
    assert len(qr_opt_vals[po.READ_LENGTH.name]) == num_values


def test_validate_option_values_returns_correct_transformed_values():
    len1 = 10
    len2 = 20
    options = {
        po.READ_LENGTH: str(len1) + "," + str(len2)
    }

    opt_vals, qr_opt_vals = po._validate_option_values(
        _get_logger(), {o.get_option_name(): v for o, v in options.items()},
        {}, options.keys())

    for val in [len1, len2]:
        assert val in qr_opt_vals[po.READ_LENGTH.name]


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


def test_read_file_options_reads_options_from_file():
    quant_method = "Cufflinks"
    read_length = "10,20"
    with tempfile.NamedTemporaryFile() as f:
        f.write(po.QUANT_METHOD.get_option_name() + " " + quant_method + "\n")
        f.write(po.READ_LENGTH.get_option_name() + " " + read_length + "\n")
        f.flush()

        options = po._read_file_options(f.name)

        assert options[po.QUANT_METHOD.get_option_name()] == quant_method
        assert options[po.READ_LENGTH.get_option_name()] == read_length


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
    piquant_options_values = {
        "read_length": [1] * len1,
        "read_depth": [2] * len2,
        "errors": [True] * len3
    }

    execute_counter = []

    def count_incrementer(logger, options, **piquant_options):
        execute_counter.append(1)

    command = pc._PiquantCommand("dummy", [])
    command.executables = [count_incrementer]
    po.execute_for_mqr_option_sets(
        command, None, None, **piquant_options_values)

    assert len(execute_counter) == len1 * len2 * len3


def test_execute_for_mqr_option_sets_executes_all_callables():
    execute_record = []

    def callable1(logger, options, **piquant_options):
        execute_record.append(1)

    def callable2(logger, options, **piquant_options):
        execute_record.append(2)

    command = pc._PiquantCommand("dummy", [])
    command.executables = [callable1, callable2]
    po.execute_for_mqr_option_sets(
        command, None, None, piquant_option=["a"])

    assert 1 in execute_record
    assert 2 in execute_record
    assert len(execute_record) == 2


def test_execute_for_mqr_option_sets_executes_for_correct_sets_of_piquant_options():
    piquant_options1 = [75, 100]
    piquant_options2 = [30, 50]
    piquant_options_values = {
        "read_length": piquant_options1,
        "read_depth": piquant_options2
    }

    execute_record = []

    def callable1(logger, options, **piquant_options):
        execute_record.append([v for v in piquant_options.values()])

    command = pc._PiquantCommand("dummy", [])
    command.executables = [callable1]
    po.execute_for_mqr_option_sets(
        command, None, None, **piquant_options_values)

    execute_record = \
        [set(piquant_options) for piquant_options in execute_record]
    assert set([piquant_options1[0], piquant_options2[0]]) in execute_record
    assert set([piquant_options1[0], piquant_options2[1]]) in execute_record
    assert set([piquant_options1[1], piquant_options2[0]]) in execute_record
    assert set([piquant_options1[1], piquant_options2[1]]) in execute_record
