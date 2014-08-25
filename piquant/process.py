"""
Utility functions for running scripts. Exports:

run_in_directory: Run a command in a directory.
"""

import os
import subprocess


def run_in_directory(run_dir, command, cl_args=None, nohup=True):
    """
    Run a command in the specified directory.

    Run a command in the specified directory. Unlike subprocess.Popen(), the
    command's path can be specified relative to the run directory.
    run_dir: the directory in which to run the command.
    command: the command or script to run.
    cl_args: a list of command line cl_args for the command.
    nohup: If true, invoke the command immune to hangups.
    """
    cwd = os.getcwd()
    os.chdir(run_dir)
    args = [command]
    if cl_args is not None:
        args = args + cl_args
    if nohup:
        args = ['nohup'] + args
    subprocess.Popen(args)
    os.chdir(cwd)
