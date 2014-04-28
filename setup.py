# Liberally adapted from:
# "http://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/"

import assess_isoform_quantification
import sys

from distutils.core import setup
from setuptools.command.test import test as TestCommand


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
    name='piquant',
    version=assess_isoform_quantification.__version__,
    description='Pipeline for investigating the quantification of transcripts',
    author="Owen Dando",
    author_email='owen.dando@ed.ac.uk',
    packages=['piquant'],
    install_requires=['docopt>=0.6.1',
                      'ez-setup>=0.9',
                      'gtf-to-genes>=1.09',
                      'py>=1.4.20',
                      'pytest>=2.5.2',
                      'schema>=0.2.0'],
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    extras_require={
        'testing': ['pytest']
    }
)
