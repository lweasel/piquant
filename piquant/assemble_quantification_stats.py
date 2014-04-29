#!/usr/bin/python

# TODO: tidy
# TODO: add logging
# TODO: get STATS_TYPES from statistics module.

"""Usage:
    assemble_quantification_stats [--log-level=<log-level>] [--out-dir=<out-dir>] [--run-dir=<run-dir>] --quant-methods=<quant-methods> --read-lengths=<read-lengths> --read-depths=<read-depths> --paired-ends=<paired-ends> --errors=<errors> --biases=<biases>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
--run-dir=<out-dir>                 Parent directory for quantification run directories [default: output].
--out-dir=<out-dir>                 Directory to output assembled stats files to [default: output/overall_stats].
-q --quant-methods=<quant-methods>  Comma-separated list of quantification methods to assemble stats for.
-l --read-lengths=<read-lengths>    Comma-separated list of read-lengths to assemble stats for.
-d --read-depths=<read-depths>      Comma-separated list of read-depths to assemble stats for.
-p --paired-ends=<paired-ends>      Comma-separated list of True/False strings indicating whether stats should be assembled for single or paired-end reads.
-e --errors=<errors>                Comma-separated list of True/False strings indicating whether stats should be assembled with or without read errors.
-b --biases=<biases>                Comma-separated list of True/False strings indicating whether stats should be assembled with or without read sequence bias.
"""

import docopt
import itertools
import ordutils.log as log
import ordutils.options as opt
import os.path
import pandas as pd
import piquant_options as popt
import quantification_run as qr
import schema
import statistics as stats
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
OUTPUT_DIRECTORY = "--out-dir"
RUN_DIRECTORY = "--run-dir"
QUANT_METHODS = "--quant-methods"
READ_LENGTHS = "--read-lengths"
READ_DEPTHS = "--read-depths"
PAIRED_ENDS = "--paired-ends"
ERRORS = "--errors"
BIASES = "--biases"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="assemble_quantification_stats v0.1")

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    options[QUANT_METHODS] = opt.validate_list_option(
        options[QUANT_METHODS], popt.check_quantification_method)
    options[READ_LENGTHS] = set(opt.validate_list_option(
        options[READ_LENGTHS], int))
    options[READ_DEPTHS] = set(opt.validate_list_option(
        options[READ_DEPTHS], int))
    options[PAIRED_ENDS] = set(opt.validate_list_option(
        options[PAIRED_ENDS], opt.check_boolean_value))
    options[ERRORS] = set(opt.validate_list_option(
        options[ERRORS], opt.check_boolean_value))
    options[BIASES] = set(opt.validate_list_option(
        options[BIASES], opt.check_boolean_value))
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

if not os.path.exists(options[OUTPUT_DIRECTORY]):
    os.mkdir(options[OUTPUT_DIRECTORY])

for quant_method, length, depth, paired_end, error, bias in \
    itertools.product(
        options[QUANT_METHODS], options[READ_LENGTHS],
        options[READ_DEPTHS], options[PAIRED_ENDS],
        options[ERRORS], options[BIASES]):

    # Check all run directories exist
    run_dir = qr.get_run_dir(
        options[RUN_DIRECTORY], quant_method,
        length, depth, paired_end, error, bias)
    if not os.path.exists(run_dir):
        sys.exit("Run directory '{d}' should exist.".format(d=run_dir))

for pset in stats.get_stats_param_sets():
    overall_stats_df = pd.DataFrame()

    for quant_method, length, depth, paired_end, error, bias in \
        itertools.product(
            options[QUANT_METHODS], options[READ_LENGTHS],
            options[READ_DEPTHS], options[PAIRED_ENDS],
            options[ERRORS], options[BIASES]):

        run_name = qr.get_run_name(
            quant_method, length, depth, paired_end, error, bias)
        run_dir = qr.get_run_dir(
            options[RUN_DIRECTORY], quant_method,
            length, depth, paired_end, error, bias)

        stats_file = stats.get_stats_file(run_dir, run_name, **pset)
        stats_df = pd.read_csv(stats_file)
        overall_stats_df = overall_stats_df.append(stats_df)

    overall_stats_file = stats.get_stats_file(
        options[OUTPUT_DIRECTORY], stats.OVERALL_STATS_PREFIX, **pset)
    stats.write_stats_data(overall_stats_file, overall_stats_df, index=False)
