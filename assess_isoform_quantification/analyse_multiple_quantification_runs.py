#!/usr/bin/python

"""Usage:
    analyse_multiple_quantification_runs [--log-level=<log-level>] <stats-prefix>

-h --help                 Show this message.
-v --version              Show version.
--log-level=<log-level>   Set logging level (one of {log_level_vals}) [default: info].
<stats-prefix>            Path prefix for files containing assembled overall statistics.
"""

import classifiers
import docopt
import ordutils.log as log
import ordutils.options as opt
import schema
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
STATS_PREFIX = "<stats-prefix>"

QUANTIFICATION_METHOD_PARAM = "quant-method"
READ_LENGTH_PARAM = "read-length"
READ_DEPTH_PARAM = "read-depth"
PAIRED_END_PARAM = "paired-end"
ERRORS_PARAM = "errors"

NUMERICAL_PARAMS = [
    READ_LENGTH_PARAM,
    READ_DEPTH_PARAM
]

GROUPING_PARAMS = [
    QUANTIFICATION_METHOD_PARAM,
    PAIRED_END_PARAM,
    ERRORS_PARAM
]

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="analyse_multiple_quantification_runs v0.1")

clsfrs = classifiers.get_classifiers()

# TODO: duplicate from analyse_run_quantification.py
space_to_underscore = lambda x: x.replace(' ', '_')

STATS_FILES = [""] + ["_by_" + space_to_underscore(c.get_column_name())
                      for c in clsfrs if c.produces_grouped_stats()]
STATS_FILES = [options[STATS_PREFIX] + t + ".csv" for t in STATS_FILES]

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    for stats_file in STATS_FILES:
        opt.validate_file_option(
            stats_file, "Statistics file should exist")
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger

logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Read in overall statistics
