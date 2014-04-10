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
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
STATS_PREFIX = "<stats-prefix>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="analyse_multiple_quantification_runs v0.1")

clsfrs = classifiers.get_classifiers()

# TODO: duplicate from analyse_run_quantification.py
space_to_underscore = lambda x: x.replace(' ', '_')

STATS_FILES = [""] + ["_by_" + space_to_underscore(c.get_column_name())
                      for c in clsfrs if c.produces_grouped_stats()]
STATS_FILES = [options[STATS_PREFIX] + t + ".csv" for t in STATS_FILES]

OVERALL_STATS_FILE = STATS_FILES[0]

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

overall_stats = pd.read_csv(OVERALL_STATS_FILE)

# Create graphs based on overall statistics
# TODO: better description!

param_values = {}

NUMERICAL_PARAMS = set([p for p in params.PARAMETERS if p.is_numeric])

for param in params.PARAMETERS:
    param_values[param] = overall_stats[param.name].value_counts().index.tolist()

for param in params.PARAMETERS:
    if len(param_values[param]) <= 1:
        continue

    for numerical_param in NUMERICAL_PARAMS - set([param]):
        if len(param_values[numerical_param]) <= 1:
            continue

        fixed_params = [p for p in (set(params.PARAMETERS) - set([param, numerical_param])) if len(param_values[p]) > 1]
        fixed_param_values_sets = [v for v in itertools.product(*[param_values[p] for p in fixed_params])]

        for fixed_param_values_set in fixed_param_values_sets:
            stats = overall_stats
            fixed_param_values = {}
            for i, fp in enumerate(fixed_params):
                fp_value = fixed_param_values_set[i]
                stats = stats[stats[fp.name] == fp_value]
                fixed_param_values[fp] = fp_value

            opts = plot.PlotOptions("dummy", "this is the label", options[STATS_PREFIX])
            plot.plot_statistic(stats, opts, "sensitivity", param, numerical_param, fixed_param_values)
