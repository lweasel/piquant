import piquant.process as ps
import os.path
import time
import utils

SCRIPT_NAME = "./script.sh"


def _write_and_run_script(dirname, command):
    utils.write_executable_script(dirname, SCRIPT_NAME, command)
    ps.run_in_directory(dirname, SCRIPT_NAME)
    time.sleep(0.1)


def test_run_in_directory_executes_command_in_directory():
    with utils.temp_dir_created() as dirname:
        _write_and_run_script(dirname, "pwd > out.txt")

        with open(dirname + os.path.sep + 'out.txt') as f:
            path = f.readlines()[0].strip()
            assert path == os.path.realpath(dirname)


def test_run_in_directory_include_command_line_args():
    with utils.temp_dir_created() as dirname:
        ps.run_in_directory(dirname, "touch", [SCRIPT_NAME])
        time.sleep(0.1)
        assert os.path.exists(dirname + os.path.sep + SCRIPT_NAME)
