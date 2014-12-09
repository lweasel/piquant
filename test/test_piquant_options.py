import piquant.log as log
import piquant.piquant_options as po
import piquant.piquant_commands as pc
import pytest
import schema
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

    po._PiquantOption.OPTIONS.remove(option)
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
        po.READ_LENGTH.get_option_name():
        ",".join([str(i) for i in range(1, num_values + 1)])
    }

    opt_vals, qr_opt_vals = po._validate_option_values(
        _get_logger(), options, {}, [po.READ_LENGTH])

    assert len(qr_opt_vals[po.READ_LENGTH.name]) == num_values


def test_validate_option_values_returns_correct_transformed_values():
    len1 = 10
    len2 = 20
    options = {
        po.READ_LENGTH.get_option_name(): str(len1) + "," + str(len2)
    }

    opt_vals, qr_opt_vals = po._validate_option_values(
        _get_logger(), options, {}, [po.READ_LENGTH])

    for val in [len1, len2]:
        assert val in qr_opt_vals[po.READ_LENGTH.name]


def test_validate_option_values_raises_exception_for_invalid_option_values():
    options = {
        po.READ_LENGTH.get_option_name(): "abc"
    }

    with pytest.raises(schema.SchemaError):
        po._validate_option_values(
            _get_logger(), options, {}, [po.READ_LENGTH])


def test_validate_option_values_raises_exception_if_quant_run_option_values_not_supplied():
    options = {
        po.READ_LENGTH.get_option_name(): "10,20",
    }

    with pytest.raises(schema.SchemaError):
        po._validate_option_values(
            _get_logger(), options, {}, [po.READ_LENGTH, po.READ_DEPTH])


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


def test_validate_option_values_overrides_file_parameters_with_cl_options():
    len1 = 10
    len2 = 20
    len3 = 30

    file_options = {
        po.READ_LENGTH.get_option_name(): str(len1) + "," + str(len2)
    }

    cl_options = {
        po.READ_LENGTH.get_option_name(): str(len3)
    }

    opts, qr_opts = po._validate_option_values(
        _get_logger(), cl_options, file_options, [po.READ_LENGTH])

    assert len1 not in qr_opts[po.READ_LENGTH.name]
    assert len2 not in qr_opts[po.READ_LENGTH.name]
    assert len3 in qr_opts[po.READ_LENGTH.name]


def test_get_file_name_returns_correct_name():
    assert po.get_file_name(
        read_depth=30, read_length=50,
        paired_end=True, bias=False) == "30x_50b_pe_no_bias"


def test_get_value_names_returns_correct_translated_values():
    mqr_options = {
        po.READ_LENGTH: 50,
        po.READ_DEPTH: 10,
        po.STRANDED: False,
        po.PAIRED_END: True,
        po.ERRORS: False,
        po.BIAS: True,
        po.NOISE_DEPTH_PERCENT: 5,
        po.QUANT_METHOD: "Cufflinks",
    }

    assert po.get_value_names(mqr_options) == \
        ["Cufflinks", "10x", "50b", "paired-end", "no errors",
         "unstranded", "with bias", "noise-5x"]


def test_execute_for_mqr_option_sets_executes_for_correct_number_of_option_sets():
    len1 = 3
    len2 = 5
    len3 = 7
    qr_options = {
        po.READ_LENGTH.name: [1] * len1,
        po.READ_DEPTH.name: [2] * len2,
        po.ERRORS.name: [True] * len3
    }

    execute_counter = []

    def count_incrementer(logger, script_dir, options, **qr_options):
        del logger
        del script_dir
        del options
        del qr_options
        execute_counter.append(1)

    command = pc._PiquantCommand("dummy", [])
    command.executables = [count_incrementer]
    po.execute_for_mqr_option_sets(
        command, None, None, None, qr_options)

    assert len(execute_counter) == len1 * len2 * len3


def test_execute_for_mqr_option_sets_executes_all_callables():
    qr_options = {
        po.READ_LENGTH.name: [50]
    }

    execute_record = []

    def callable1(logger, script_dir, options, **qr_options):
        del logger
        del script_dir
        del options
        del qr_options
        execute_record.append(1)

    def callable2(logger, script_dir, options, **qr_options):
        del logger
        del script_dir
        del options
        del qr_options
        execute_record.append(2)

    command = pc._PiquantCommand("dummy", [])
    command.executables = [callable1, callable2]

    po.execute_for_mqr_option_sets(
        command, None, None, None, qr_options)

    assert 1 in execute_record
    assert 2 in execute_record
    assert len(execute_record) == 2


def test_execute_for_mqr_option_sets_executes_for_correct_sets_of_piquant_options():
    piquant_options1 = [75, 100]
    piquant_options2 = [30, 50]
    qr_options = {
        po.READ_LENGTH.name: piquant_options1,
        po.READ_DEPTH.name: piquant_options2
    }

    execute_record = []

    def callable1(logger, script_dir, options, **qr_option):
        del logger
        del script_dir
        del options
        execute_record.append([v for v in qr_option.values()])

    command = pc._PiquantCommand("dummy", [])
    command.executables = [callable1]
    po.execute_for_mqr_option_sets(
        command, None, None, None, qr_options)

    execute_record = \
        [set(piquant_options) for piquant_options in execute_record]
    assert set([piquant_options1[0], piquant_options2[0]]) in execute_record
    assert set([piquant_options1[0], piquant_options2[1]]) in execute_record
    assert set([piquant_options1[1], piquant_options2[0]]) in execute_record
    assert set([piquant_options1[1], piquant_options2[1]]) in execute_record
