# Liberally adapted from:
# "http://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/"

import codecs
import os
import re
import sys

from distutils.core import setup
from setuptools.command.test import test as TestCommand

here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

long_description = read('README.rst')


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='piquant-pipeline',
    version=find_version('piquant', '__init__.py'),
    packages=['piquant'],
    license='MIT License',
    description='Pipeline for investigating the quantification of transcripts',
    long_description=long_description,
    author="Owen Dando",
    author_email='owen.dando@ed.ac.uk',
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    extras_require={
        'testing': ['pytest']
    }
)
