import os
import os.path
import piquant.log as log
import piquant.piquant as piq
import piquant.piquant_options as po
import piquant.quantifiers as quant
import pytest
import sys
import time
import utils


def _get_test_options(output_dir):
    return {
        po.READS_OUTPUT_DIR.name: output_dir,
        po.QUANT_OUTPUT_DIR.name: output_dir,
        po.NO_CLEANUP.name: True,
        po.PLOT_FORMAT.name: "pdf",
        po.GROUPED_THRESHOLD.name: 3000,
        po.ERROR_FRACTION_THRESHOLD.name: 10,
        po.NOT_PRESENT_CUTOFF.name: 0.1
    }


def get_test_qr_options(quant_method=None):
    qr_options = {
        po.READ_DEPTH.name: 30,
        po.READ_LENGTH.name: 50,
        po.PAIRED_END.name: True,
        po.BIAS.name: False
    }

    if quant_method:
        qr_options[po.QUANT_METHOD.name] = quant_method

    return qr_options


def _get_logger():
    return log.get_logger(sys.stderr, "critical")


def _check_file_exists(parent_dir, file_name):
    assert os.path.exists(parent_dir + os.path.sep + file_name)


def test_get_options_dir_returns_correct_path_for_reads_directory():
    output_dir = "dummy"
    assert piq._get_options_dir(
        False, _get_test_options(output_dir),
        **get_test_qr_options(quant_method="quant")) == \
        output_dir + os.path.sep + "30x_50b_pe_no_bias"


def test_get_options_dir_returns_correct_path_for_run_directory():
    output_dir = "dummy"
    quant_method = "quant"
    assert piq._get_options_dir(
        True, _get_test_options(output_dir),
        **get_test_qr_options(quant_method=quant_method)) == \
        output_dir + os.path.sep + quant_method + "_30x_50b_pe_no_bias"


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        options_dir = piq._get_options_dir(
            False, test_options, **get_test_qr_options())
        os.mkdir(options_dir)

        directory_checker = piq._reads_directory_checker(True)
        directory_checker(_get_logger(), test_options, **get_test_qr_options())


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(True)
        with pytest.raises(SystemExit):
            directory_checker(_get_logger(), test_options,
                              **get_test_qr_options())


def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(False)
        directory_checker(_get_logger(), test_options,
                          **get_test_qr_options())


def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        options_dir = piq._get_options_dir(
            False, test_options, **get_test_qr_options())
        os.mkdir(options_dir)

        directory_checker = piq._reads_directory_checker(False)
        with pytest.raises(SystemExit):
            directory_checker(_get_logger(), test_options,
                              **get_test_qr_options())


def test_prepare_read_simulation_creates_correct_files():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        qr_options = get_test_qr_options()
        piq._prepare_read_simulation(_get_logger(), options, **qr_options)

        reads_dir = piq._get_options_dir(False, options, **qr_options)
        _check_file_exists(reads_dir, "run_simulation.sh")
        _check_file_exists(reads_dir, "flux_simulator_main_expression.par")
        _check_file_exists(reads_dir, "flux_simulator_main_simulation.par")


def test_create_reads_executes_run_simulation_script():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        qr_options = get_test_qr_options()

        reads_dir = piq._get_options_dir(False, options, **qr_options)
        os.mkdir(reads_dir)

        test_filename = "test"
        utils.write_executable_script(
            reads_dir, "run_simulation.sh", "touch " + test_filename)

        piq._create_reads(_get_logger(), options, **qr_options)
        time.sleep(0.1)

        assert os.path.exists(reads_dir + os.path.sep + test_filename)


def test_prepare_quantification_creates_correct_file():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        qr_options = get_test_qr_options(quant_method=quant._Cufflinks())
        piq._prepare_quantification(
            _get_logger(), options, **qr_options)

        quant_dir = piq._get_options_dir(True, options, **qr_options)
        _check_file_exists(quant_dir, "run_quantification.sh")
