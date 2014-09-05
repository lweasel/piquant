import os
import os.path
import piquant.piquant as piq
import piquant.piquant_options as po
import pytest
import utils

TEST_OUTPUT_DIR = "dummy"

TEST_PARAMS = {
    "read_depth": 30,
    "read_length": 50,
    "paired_end": True,
    "bias": False
}


def _get_test_options(test_output_dir=TEST_OUTPUT_DIR):
    return {
        po.OUTPUT_DIRECTORY: test_output_dir,
        po.NO_CLEANUP: True
    }


def _check_file_exists(parent_dir, file_name):
    assert os.path.exists(parent_dir + os.path.sep + file_name)


def test_get_parameters_dir_returns_correct_path():
    assert piq._get_parameters_dir(_get_test_options(), **TEST_PARAMS) == \
        TEST_OUTPUT_DIR + os.path.sep + "30x_50b_pe_no_bias"


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        params_dir = piq._get_parameters_dir(test_options, **TEST_PARAMS)
        os.mkdir(params_dir)

        directory_checker = piq._reads_directory_checker(True)
        directory_checker(None, test_options, **TEST_PARAMS)


def test_read_directory_checker_returns_correct_checker_if_directory_should_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(True)
        with pytest.raises(SystemExit):
            directory_checker(None, test_options, **TEST_PARAMS)

def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_doesnt_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)

        directory_checker = piq._reads_directory_checker(False)
        directory_checker(None, test_options, **TEST_PARAMS)


def test_read_directory_checker_returns_correct_checker_if_directory_shouldnt_exist_and_does_exist():
    with utils.temp_dir_created() as temp_dir:
        test_options = _get_test_options(temp_dir)
        params_dir = piq._get_parameters_dir(test_options, **TEST_PARAMS)
        os.mkdir(params_dir)

        directory_checker = piq._reads_directory_checker(False)
        with pytest.raises(SystemExit):
            directory_checker(None, test_options, **TEST_PARAMS)


def test_prepare_read_simulation_creates_correct_files():
    with utils.temp_dir_created() as dir_path:
        options = _get_test_options(dir_path)
        piq._prepare_read_simulation(None, options, **TEST_PARAMS)
        reads_dir = dir_path + os.path.sep + "30x_50b_pe_no_bias"
        _check_file_exists(reads_dir, "run_simulation.sh")
        _check_file_exists(reads_dir, "flux_simulator_expression.par")
        _check_file_exists(reads_dir, "flux_simulator_simulation.par")
