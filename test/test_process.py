import piquant.process as ps
import os.path
import stat
import time

from utils import temp_dir_created

SCRIPT_NAME = "./script.sh"


def _write_and_run_script(dirname, command):
        filename = dirname + os.path.sep + SCRIPT_NAME
        with open(filename, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(command + "\n")

        os.chmod(filename,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                 stat.S_IRGRP | stat.S_IROTH)

       #writer = fw.BashScriptWriter()
       #writer.add_line(command)
       #writer.write_to_file(dirname, SCRIPT_NAME)
        ps.run_in_directory(dirname, SCRIPT_NAME)
        time.sleep(0.1)


def test_run_in_directory_executes_command_in_directory():
    with temp_dir_created() as dirname:
        _write_and_run_script(dirname, "pwd > out.txt")

        with open(dirname + os.path.sep + 'out.txt') as f:
            path = f.readlines()[0].strip()
            assert path == dirname


def test_run_in_directory_include_command_line_args():
    with temp_dir_created() as dirname:
        ps.run_in_directory(dirname, "touch", [SCRIPT_NAME])
        time.sleep(0.1)
        assert os.path.exists(dirname + os.path.sep + SCRIPT_NAME)
