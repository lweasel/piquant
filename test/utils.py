import shutil
import tempfile


class TempDir:
    def __enter__(self):
        self.dirname = tempfile.mkdtemp()
        return self.dirname

    def __exit__(self, type, value, traceback):
        shutil.rmtree(self.dirname)
