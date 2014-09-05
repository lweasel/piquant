import contextlib
import shutil
import os.path
import stat
import tempfile


@contextlib.contextmanager
def temp_dir_created():
    dirname = tempfile.mkdtemp()
    try:
        yield dirname
    finally:
        shutil.rmtree(dirname)


def write_executable_script(dir_name, script_name, command):
    file_name = dir_name + os.path.sep + script_name
    with open(file_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(command + "\n")

    os.chmod(file_name,
             stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
             stat.S_IRGRP | stat.S_IROTH)
