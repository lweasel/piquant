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
import fpkms_plotting as plot
import itertools
import ordutils.log as log
import ordutils.options as opt
import os.path
import pandas as pd
import parameters as params
import schema
import statistics as stats
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

    stats_files = [stats.get_stats_file(
        options[STATS_DIR], stats.OVERALL_STATS_PREFIX, **pset)
        for pset in stats.get_stats_param_sets()]
    for stats_file in stats_files:
        opt.validate_file_option(
            stats_file, "Statistics file should exist")
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Read in overall statistics
overall_stats_file = stats.get_stats_file(
    options[STATS_DIR], stats.OVERALL_STATS_PREFIX)
overall_stats = pd.read_csv(overall_stats_file)

param_values = {}
for param in params.PARAMETERS:
    param_series = overall_stats[param.name]
    param_values[param] = param_series.value_counts().index.tolist()

# Utility functions for manipulating sets of parameters
degenerate_param = lambda x: len(param_values[x]) <= 1
non_degenerate_params = lambda x: [p for p in x if not degenerate_param(p)]


def remove_from(parameters, to_remove):
    get_pset = lambda x: x if isinstance(x, set) \
        else (set(x) if isinstance(x, list) else set([x]))
    return get_pset(parameters) - get_pset(to_remove)


def get_fixed_params(all_params, non_fixed):
    fixed_params = [p for p in remove_from(all_params, non_fixed)
                    if not degenerate_param(p)]
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

# Create graphs based on overall statistics
# TODO: better description!
graph_file_basename = \
    options[STATS_DIR] + os.path.sep + stats.OVERALL_STATS_PREFIX

numerical_params = [p for p in params.PARAMETERS if p.is_numeric]

for param in non_degenerate_params(params.PARAMETERS):
    for num_p in non_degenerate_params(remove_from(numerical_params, param)):
        fixed_params, fp_values_sets = \
            get_fixed_params(params.PARAMETERS, [param, num_p])

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                overall_stats, fixed_params, fp_values_set)

            for stat in stats.get_graphable_statistics():
                plot.plot_statistic(
                    stats_df, graph_file_basename,
                    stat, param, num_p, fixed_param_values)

# Create graphs based on statistics stratified by classifier
# TODO: better description!
more_than_100_filter = lambda x: x[stats.NUM_FPKMS] > 100

clsfrs = classifiers.get_classifiers()
grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]
dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]

for clsfr in grp_clsfrs:
    stats_file = stats.get_stats_file(
        options[STATS_DIR], stats.OVERALL_STATS_PREFIX, clsfr)
    clsfr_stats = pd.read_csv(stats_file)

    for param in non_degenerate_params(params.PARAMETERS):
        fixed_params, fp_values_sets = \
            get_fixed_params(params.PARAMETERS, param)

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                clsfr_stats, fixed_params, fp_values_set)

            for stat in stats.get_graphable_by_classifier_statistics():
                plot.plot_statistic_by_classifier(
                    stats_df, graph_file_basename, stat, param,
                    clsfr, more_than_100_filter, fixed_param_values)

# Create distribution plots
# TODO: better description!
for clsfr, asc in itertools.product(dist_clsfrs, ORDER_VALUES):
    stats_file = stats.get_stats_file(
        options[STATS_DIR], stats.OVERALL_STATS_PREFIX, clsfr, asc)
    clsfr_stats = pd.read_csv(stats_file)

    for param in non_degenerate_params(params.PARAMETERS):
        fixed_params, fp_values_sets = \
            get_fixed_params(params.PARAMETERS, param)

        for fp_values_set in fp_values_sets:
            stats_df, fixed_param_values = get_stats_for_fixed_params(
                clsfr_stats, fixed_params, fp_values_set)

            plot.plot_cumulative_transcript_distribution_grouped(
                stats_df, graph_file_basename, param,
                clsfr, asc, fixed_param_values)
