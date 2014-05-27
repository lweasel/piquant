#!/usr/bin/python

"""Usage:
    analyse_multiple_quantification_runs [--log-level=<log-level>] <overall-stats-dir>

-h --help                 Show this message.
-v --version              Show version.
--log-level=<log-level>   Set logging level (one of {log_level_vals}) [default: info].
<overall_stats_dir>       Parent directory for files containing assembled overall statistics.
"""

import classifiers
import docopt
import tpms_plotting as plot
import itertools
import ordutils.log as log
import ordutils.options as opt
import os.path
import pandas as pd
import parameters
import schema
import statistics
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
STATS_DIR = "<overall-stats-dir>"

ORDER_VALUES = [True, False]

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(
    __doc__, version="analyse_multiple_quantification_runs v0.1")

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    stats_files = [statistics.get_stats_file(
        options[STATS_DIR], statistics.OVERALL_STATS_PREFIX, **pset)
        for pset in statistics.get_stats_param_sets()]
    for stats_file in stats_files:
        opt.validate_file_option(
            stats_file, "Statistics file should exist")
except schema.SchemaError as exc:
    exit(exc.code)

# Utility functions for manipulating sets of parameters


def degenerate_param(param, param_values):
    return len(param_values[param]) <= 1


def get_non_degenerate_params(params, param_values):
    return [p for p in params if not degenerate_param(p, param_values)]


def remove_from(params, to_remove):
    get_pset = lambda x: x if isinstance(x, set) \
        else (set(x) if isinstance(x, list) else set([x]))
    return get_pset(params) - get_pset(to_remove)


def get_fixed_params(all_params, non_fixed):
    fixed_params = [p for p in remove_from(all_params, non_fixed)
                    if not degenerate_param(p, param_values)]
    value_sets = [v for v in itertools.product(
                  *[param_values[p] for p in fixed_params])]
    return fixed_params, value_sets


def get_stats_for_fixed_params(stats_df, fixed_params, fp_values_set):
    fixed_param_values = {}
    for i, fp in enumerate(fixed_params):
        fp_value = fp_values_set[i]
        stats_df = stats_df[stats_df[fp.name] == fp_value]
        fixed_param_values[fp] = fp_value
    return stats_df, fixed_param_values

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Read in overall statistics
logger.info("Reading in statistics calculated for the whole set of TPMs...")

overall_stats_file = statistics.get_stats_file(
    options[STATS_DIR], statistics.OVERALL_STATS_PREFIX)
overall_stats = pd.read_csv(overall_stats_file)

param_values = {}
for param in parameters.get_parameters():
    param_series = overall_stats[param.name]
    param_values[param] = param_series.value_counts().index.tolist()

# Draw graphs derived from statistics calculated for the whole set of TPMs.
# e.g. the Spearman correlation of calculated and real TPMs graphed as
# read-depth varies, for each quantification method, in the case of paired-end
# reads with errors and bias.
logger.info("Drawing graphs derived from statistics calculated for the " +
            "whole set of TPMs...")

graph_file_basename = \
    options[STATS_DIR] + os.path.sep + statistics.OVERALL_STATS_PREFIX

numerical_params = [p for p in parameters.get_parameters() if p.is_numeric]

for param in get_non_degenerate_params(
        parameters.get_parameters(), param_values):
    for num_p in get_non_degenerate_params(
            remove_from(numerical_params, param), param_values):
        fixed_params, fp_values_sets = \
            get_fixed_params(parameters.get_parameters(), [param, num_p])

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                overall_stats, fixed_params, fp_values_set)

            for stat in statistics.get_graphable_statistics():
                plot.plot_statistic_by_parameter_values(
                    stats_df, graph_file_basename,
                    stat, param, num_p, fixed_param_values)

# Draw graphs derived from statistics calculated on groups of TPMs that have
# been stratified into sets based on some classifier of transcripts. e.g. the
# median percentage error of calculated vs real TPMs graphed as the percentage
# of unique sequence per-transcript varies, for single and paired-end reads, in
# the case of reads with errors and bias, and a particular quantification
# method.
logger.info("Drawing graphs derived from statistics calculated on " +
            "subsets of TPMs...")

more_than_100_filter = lambda x: x[statistics.NUM_TPMS] > 100

clsfrs = classifiers.get_classifiers()
grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]
dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]

for clsfr in grp_clsfrs:
    stats_file = statistics.get_stats_file(
        options[STATS_DIR], statistics.OVERALL_STATS_PREFIX, clsfr)
    clsfr_stats = pd.read_csv(stats_file)

    for param in get_non_degenerate_params(
            parameters.get_parameters(), param_values):
        fixed_params, fp_values_sets = \
            get_fixed_params(parameters.get_parameters(), param)

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                clsfr_stats, fixed_params, fp_values_set)

            for stat in statistics.get_graphable_by_classifier_statistics():
                filtered_stats_df = stats_df[more_than_100_filter(stats_df)]
                plot.plot_statistic_by_transcript_classifier_values(
                    filtered_stats_df, graph_file_basename, stat, param,
                    clsfr, fixed_param_values)

# Draw distributions illustrating the percentage of TPMs above or below some
# threshold as that threshold changes. e.g. the percentage of TPMs whose
# absolute percentage error in calculated TPM, as compared to real TPM, is
# below a particular threshold.
logger.info("Drawing distribution plots...")

for clsfr, asc in itertools.product(dist_clsfrs, ORDER_VALUES):
    stats_file = statistics.get_stats_file(
        options[STATS_DIR], statistics.OVERALL_STATS_PREFIX, clsfr, asc)
    clsfr_stats = pd.read_csv(stats_file)

    for param in get_non_degenerate_params(
            parameters.get_parameters(), param_values):
        fixed_params, fp_values_sets = \
            get_fixed_params(parameters.get_parameters(), param)

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                clsfr_stats, fixed_params, fp_values_set)

            plot.plot_cumulative_transcript_distribution_grouped(
                stats_df, graph_file_basename, param,
                clsfr, asc, fixed_param_values)
