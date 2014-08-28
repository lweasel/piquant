import contextlib
import shutil
import tempfile


@contextlib.contextmanager
def temp_dir_created():
    dirname = tempfile.mkdtemp()
    try:
        yield dirname
    finally:
        shutil.rmtree(dirname)
