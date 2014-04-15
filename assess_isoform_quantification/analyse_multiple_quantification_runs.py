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
import fpkms_plotting as plot
import itertools
import ordutils.log as log
import ordutils.options as opt
import pandas as pd
import parameters as params
import schema
import statistics as stats
import sys
import utils

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
STATS_PREFIX = "<stats-prefix>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="analyse_multiple_quantification_runs v0.1")

# Validate command-line options
clsfrs = classifiers.get_grouped_stats_classifiers()
dist_clsfrs = classifiers.get_distribution_plot_classifiers()


def get_stats_file(stats_type):
    ret = options[STATS_PREFIX] + "_stats"
    if stats_type is not None:
        ret += "_by_" + utils.spaces_to_underscore(stats_type)
    return ret + ".csv"


def get_distribution_stats_file(stats_type, ascending):
    ret = options[STATS_PREFIX] + "_distribution_stats_"
    ret += utils.get_order_string(ascending) + "_by_"
    ret += utils.spaces_to_underscore(stats_type)
    return ret + ".csv"

OVERALL_STATS_FILE = get_stats_file(None)

try:
    opt.validate_dict_option(options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    stats_files = [OVERALL_STATS_FILE] + \
        [get_stats_file(c.get_column_name()) for c in clsfrs] + \
        [get_distribution_stats_file(c.get_column_name(), asc)
            for c, asc in itertools.product(dist_clsfrs, [True, False])]
    for stats_file in stats_files:
        opt.validate_file_option(
            stats_file, "Statistics file should exist")
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger

logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Read in overall statistics

overall_stats = pd.read_csv(OVERALL_STATS_FILE)

# Create graphs based on overall statistics
# TODO: better description!

param_values = {}

NUMERICAL_PARAMS = set([p for p in params.PARAMETERS if p.is_numeric])

for param in params.PARAMETERS:
    param_values[param] = overall_stats[param.name].value_counts().index.tolist()

for param in []:
#for param in params.PARAMETERS:
    if len(param_values[param]) <= 1:
        continue

    for numerical_param in NUMERICAL_PARAMS - set([param]):
        if len(param_values[numerical_param]) <= 1:
            continue

        fixed_params = [p for p in (set(params.PARAMETERS) - set([param, numerical_param])) if len(param_values[p]) > 1]
        fixed_param_values_sets = [v for v in itertools.product(*[param_values[p] for p in fixed_params])]

        for fixed_param_values_set in fixed_param_values_sets:
            stats_df = overall_stats
            fixed_param_values = {}
            for i, fp in enumerate(fixed_params):
                fp_value = fixed_param_values_set[i]
                stats_df = stats_df[stats_df[fp.name] == fp_value]
                fixed_param_values[fp] = fp_value

            opts = plot.PlotOptions("dummy", "this is the label", options[STATS_PREFIX])
            for stat in [s for s in stats.get_statistics() if s.graphable]:
                plot.plot_statistic(stats_df, opts, stat, param, numerical_param, fixed_param_values)

# Create graphs based on statistics stratified by classifier

more_than_100_filter = lambda x: x["count"] > 100

for clsfr in []:
#for clsfr in clsfrs:
    stats_file = get_stats_file(clsfr.get_column_name())
    clsfr_stats = pd.read_csv(stats_file)

    for param in params.PARAMETERS:
        if len(param_values[param]) <= 1:
            continue

        fixed_params = [p for p in (set(params.PARAMETERS) - set([param])) if len(param_values[p]) > 1]
        fixed_param_values_sets = [v for v in itertools.product(*[param_values[p] for p in fixed_params])]

        for fixed_param_values_set in fixed_param_values_sets:
            stats_df = clsfr_stats
            fixed_param_values = {}
            for i, fp in enumerate(fixed_params):
                fp_value = fixed_param_values_set[i]
                stats_df = stats_df[stats_df[fp.name] == fp_value]
                fixed_param_values[fp] = fp_value

            opts = plot.PlotOptions("dummy", "this is the label", options[STATS_PREFIX])

            for stat in [s for s in stats.get_statistics() if s.graphable_by_classifier]:
                plot.plot_statistic_by_classifier(
                    stats_df, opts, stat, param,
                    clsfr, more_than_100_filter, fixed_param_values)

# Create distribution plots
for clsfr, asc in itertools.product(dist_clsfrs, [True, False]):
    stats_file = get_distribution_stats_file(clsfr.get_column_name(), asc)
    clsfr_stats = pd.read_csv(stats_file)

    for param in params.PARAMETERS:
        if len(param_values[param]) <= 1:
            continue

        fixed_params = [p for p in (set(params.PARAMETERS) - set([param])) if len(param_values[p]) > 1]
        fixed_param_values_sets = [v for v in itertools.product(*[param_values[p] for p in fixed_params])]

        for fixed_param_values_set in fixed_param_values_sets:
            stats_df = clsfr_stats
            fixed_param_values = {}
            for i, fp in enumerate(fixed_params):
                fp_value = fixed_param_values_set[i]
                stats_df = stats_df[stats_df[fp.name] == fp_value]
                fixed_param_values[fp] = fp_value

            opts = plot.PlotOptions("dummy", "this is the label", options[STATS_PREFIX])

            plot.plot_cumulative_transcript_distribution_grouped(
                stats_df, opts, param, clsfr, asc, fixed_param_values)
