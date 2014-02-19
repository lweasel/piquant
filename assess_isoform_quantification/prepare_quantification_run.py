#!/usr/bin/python

"""Usage:
    prepare_quantification_run [--log-level=<log-level>] [--run-directory <run-directory>]

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
-d --run-directory=<run-directory>  Directory to create to which run files will be written [default: out].
"""

from docopt import docopt
from schema import SchemaError

import log
import options as opt
import os
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
RUN_DIRECTORY = "--run-directory"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="prepare_quantification_run v0.1")

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    opt.validate_dir_option(
        options[RUN_DIRECTORY],
        "Run directory should not already exist",
        should_exist=False)
except SchemaError as exc:
    exit("Exiting. " + exc.code)

logger = log.getLogger(sys.stderr, options[LOG_LEVEL])

logger.info("Creating run directory {dir}...".
            format(dir=options[RUN_DIRECTORY]))

os.mkdir(options[RUN_DIRECTORY])
