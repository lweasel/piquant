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
        po.OUTPUT_DIRECTORY: output_dir,
        po.NO_CLEANUP: True,
        po.PLOT_FORMAT: "pdf",
        po.GROUPED_THRESHOLD: 3000,
        po.ERROR_FRACTION_THRESHOLD: 10
    }


def _get_test_params(quant_method=None):
    params = {
        "read_depth": 30,
        "read_length": 50,
        "paired_end": True,
        "bias": False
    }

    if quant_method:
        params["quant_method"] = quant_method

    return params


def _get_logger():
    return log.get_logger(sys.stderr, "critical")


def _check_file_exists(parent_dir, file_name):
    assert os.path.exists(parent_dir + os.path.sep + file_name)


def test_get_parameters_dir_returns_correct_path():
    output_dir = "dummy"
    assert piq._get_parameters_dir(
        _get_test_options(output_dir), **_get_test_params()) == \
        output_dir + os.path.sep + "30x_50b_pe_no_bias"


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        params_dir = piq._get_parameters_dir(
            test_options, **_get_test_params())
        os.mkdir(params_dir)

        directory_checker = piq._reads_directory_checker(True)
        directory_checker(_get_logger(), test_options, **_get_test_params())


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(True)
        with pytest.raises(SystemExit):
            directory_checker(_get_logger(), test_options, **_get_test_params())


def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(False)
        directory_checker(_get_logger(), test_options, **_get_test_params())


def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        params_dir = piq._get_parameters_dir(
            test_options, **_get_test_params())
        os.mkdir(params_dir)

        directory_checker = piq._reads_directory_checker(False)
        with pytest.raises(SystemExit):
            directory_checker(_get_logger(), test_options, **_get_test_params())


def test_prepare_read_simulation_creates_correct_files():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        params = _get_test_params()
        piq._prepare_read_simulation(_get_logger(), options, **params)

        reads_dir = piq._get_parameters_dir(options, **params)
        _check_file_exists(reads_dir, "run_simulation.sh")
        _check_file_exists(reads_dir, "flux_simulator_expression.par")
        _check_file_exists(reads_dir, "flux_simulator_simulation.par")


def test_create_reads_executes_run_simulation_script():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        params = _get_test_params()

        reads_dir = piq._get_parameters_dir(options, **params)
        os.mkdir(reads_dir)

        test_filename = "test"
        utils.write_executable_script(
            reads_dir, "run_simulation.sh", "touch " + test_filename)

        piq._create_reads(_get_logger(), options, **params)
        time.sleep(0.1)

        assert os.path.exists(reads_dir + os.path.sep + test_filename)


def test_prepare_quantification_creates_correct_file():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        params = _get_test_params(quant_method=quant._Cufflinks())
        piq._prepare_quantification(_get_logger(), options, **params)

        quant_dir = piq._get_parameters_dir(options, **params)
        _check_file_exists(quant_dir, "run_quantification.sh")
