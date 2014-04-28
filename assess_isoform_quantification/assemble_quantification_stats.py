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
import quantifiers as qs
import schema
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

STATS_TYPES = ["_stats", "_stats_by_gene_transcript_number", "_stats_by_log10_real_FPKM", "_stats_by_transcript_length", "_stats_by_unique_sequence_percentage", "_distribution_stats_asc_by_absolute_percent_error", "_distribution_stats_desc_by_absolute_percent_error"]

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="assemble_quantification_stats v0.1")

# Validate command-line options


def validate_list_option(option, item_validator, separator=','):
    items = option.split(separator)
    return [schema.Schema(
        schema.Use(item_validator)).validate(i) for i in items]


def check_boolean_value(data):
    data = data.lower()
    if data in ["true", "t", "yes", "y"]:
        return True
    elif data in ["false", "f", "no", "n"]:
        return False
    else:
        raise Exception("Can't convert '{d}' to bool.".format(d=data))


def check_quantification_method(data):
    available_methods = qs.get_quantification_methods()
    return available_methods[data]
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    options[QUANT_METHODS] = validate_list_option(
        options[QUANT_METHODS], check_quantification_method)
    options[READ_LENGTHS] = set(validate_list_option(
        options[READ_LENGTHS], int))
    options[READ_DEPTHS] = set(validate_list_option(
        options[READ_DEPTHS], int))
    options[PAIRED_ENDS] = set(validate_list_option(
        options[PAIRED_ENDS], check_boolean_value))
    options[ERRORS] = set(validate_list_option(
        options[ERRORS], check_boolean_value))
    options[BIASES] = set(validate_list_option(
        options[BIASES], check_boolean_value))
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger

logger = log.get_logger(sys.stderr, options[LOG_LEVEL])


def _get_reads_name(length, depth, paired_end, errors, bias):
    return "{d}x_{l}b_{p}_{e}_{b}".format(
        d=depth, l=length,
        p="pe" if paired_end else "se",
        e="errors" if errors else "no_errors",
        b="bias" if bias else "no_bias")


def _get_run_name(quant_method, length, depth, paired_end, errors, bias):
    return "{qm}_{r}".format(
        qm=quant_method.get_name(),
        r=_get_reads_name(length, depth, paired_end, errors, bias))


def _get_run_dir(output_dir, quant_method, length, depth,
                 paired_end, error, bias):
    return output_dir + os.path.sep + \
        _get_run_name(quant_method, length, depth, paired_end, error, bias)

if not os.path.exists(options[OUTPUT_DIRECTORY]):
    os.mkdir(options[OUTPUT_DIRECTORY])

for quant_method, length, depth, paired_end, error, bias in \
    itertools.product(
        options[QUANT_METHODS], options[READ_LENGTHS],
        options[READ_DEPTHS], options[PAIRED_ENDS],
        options[ERRORS], options[BIASES]):

    # Check all run directories exist
    run_dir = _get_run_dir(
        options[RUN_DIRECTORY], quant_method,
        length, depth, paired_end, error, bias)
    if not os.path.exists(run_dir):
        sys.exit("Run directory '{d}' should exist.".format(d=run_dir))

for stats_type in STATS_TYPES:
    overall_stats_df = pd.DataFrame()

    for quant_method, length, depth, paired_end, error, bias in \
        itertools.product(
            options[QUANT_METHODS], options[READ_LENGTHS],
            options[READ_DEPTHS], options[PAIRED_ENDS],
            options[ERRORS], options[BIASES]):

        run_name = _get_run_name(
            quant_method, length, depth, paired_end, error, bias)
        run_dir = _get_run_dir(
            options[RUN_DIRECTORY], quant_method,
            length, depth, paired_end, error, bias)

        stats_file = run_dir + os.path.sep + run_name + stats_type + ".csv"
        stats_df = pd.read_csv(stats_file)
        overall_stats_df = overall_stats_df.append(stats_df)

    overall_stats_file = options[OUTPUT_DIRECTORY] + os.path.sep + \
        "overall" + stats_type + ".csv"
    overall_stats_df.to_csv(
        overall_stats_file, float_format="%.5f", index=False)
