#!/usr/bin/python

# TODO: tidy
# TODO: add logging
# TODO: get STATS_TYPES from statistics module.

"""Usage:
    assemble_quantification_stats [--log-level=<log-level>] [--out-dir=<out-dir>] [--run-dir=<run-dir>] --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --polya=<polya>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
--run-dir=<out-dir>                 Parent directory for quantification run directories [default: output].
--out-dir=<out-dir>                 Directory to output assembled stats files to [default: output/overall_stats].
-q --quant-method=<quant-methods>  Comma-separated list of quantification methods to assemble stats for.
-l --read-length=<read-lengths>    Comma-separated list of read-lengths to assemble stats for.
-d --read-depth=<read-depths>      Comma-separated list of read-depths to assemble stats for.
-p --paired-end=<paired-ends>      Comma-separated list of True/False strings indicating whether stats should be assembled for single or paired-end reads.
-e --error=<errors>                Comma-separated list of True/False strings indicating whether stats should be assembled with or without read errors.
-b --bias=<biases>                Comma-separated list of True/False strings indicating whether stats should be assembled with or without read sequence bias.
-a --polya=<polya>                  Comma-separated list of True/False strings indicating whether stats should be assembled assuming transcripts do or do not have polyA tails.
"""

import docopt
import ordutils.log as log
import ordutils.options as opt
import os.path
import pandas as pd
import parameters
import schema
import statistics as stats
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
OUTPUT_DIRECTORY = "--out-dir"
RUN_DIRECTORY = "--run-dir"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="assemble_quantification_stats v0.1")

# Validate command-line options
param_values = None

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    param_values = parameters.validate_command_line_parameter_sets(options)
except schema.SchemaError as exc:
    exit(exc.code)


def check_run_directory(**params):
    run_dir = options[RUN_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)
    if not os.path.exists(run_dir):
        sys.exit("Run directory '{d}' should exist.".format(d=run_dir))


class StatsAccumulator:
    def __init__(self):
        self.overall_stats_df = pd.DataFrame()

    def __call__(self, **params):
        run_name = parameters.get_file_name(**params)
        run_dir = options[RUN_DIRECTORY] + os.path.sep + run_name

        stats_file = stats.get_stats_file(run_dir, run_name, **pset)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

if not os.path.exists(options[OUTPUT_DIRECTORY]):
    os.mkdir(options[OUTPUT_DIRECTORY])

parameters.execute_for_param_sets(
    [check_run_directory],
    **param_values)

for pset in stats.get_stats_param_sets():
    stats_acc = StatsAccumulator()
    parameters.execute_for_param_sets([stats_acc], **param_values)

    overall_stats_file = stats.get_stats_file(
        options[OUTPUT_DIRECTORY], stats.OVERALL_STATS_PREFIX, **pset)
    stats.write_stats_data(
        overall_stats_file, stats_acc.overall_stats_df, index=False)
