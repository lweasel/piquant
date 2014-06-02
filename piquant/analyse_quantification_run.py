#!/usr/bin/python

"""Usage:
    analyse_quantification_run [--log-level=<log-level> --scatter-max=<scatter-max-val> --log10-scatter-min=<log10-scatter-min-val> --log10-scatter-max=<log10-scatter-max-val>] --quant-method=<quant-method> --read-length=<read-length> --read-depth=<read-depth> --paired-end=<paired-end> --error=<errors> --bias=<bias> <tpm-file> [<out-file>]

-h --help                                    Show this message.
-v --version                                 Show version.
--log-level=<log-level>                      Set logging level (one of {log_level_vals}) [default: info].
--scatter-max=<scatter-max-val>              Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log10-scatter-min=<log10-scatter-min-val>  Minimum x and y values for log10 scatter plot; a value of 0 means do not impose a minimum [default: 0].
--log10-scatter-max=<log10-scatter-max-val>  Maximum x and y values for log10 scatter plot; a value of 0 means do not impose a maximum [default: 0].
--quant-method=<quant-method>                Method used to quantify transcript abundances.
--read-length=<read-length>                  The length of sequence reads.
--read-depth=<read-depth>                    The depth of reads sequenced across the transcriptome.
--paired-end=<paired-end>                    Whether paired-end sequence reads were used.
--error=<errors>                             Whether the reads contain sequencing errors.
--bias=<bias>                                Whether the reads contain sequence bias.
<tpm-file>                                   File containing real and calculated TPMs.
<out-file>                                   Basename for output graph and data files.
"""

import classifiers
import docopt
import itertools
import ordutils.log as log
import ordutils.options as opt
import pandas as pd
import parameters
import statistics
import tpms as t
import plot
import schema
import sys

TRANSCRIPT_COUNT_LABEL = "No. transcripts per gene"
TRUE_POSITIVES_LABEL = "true positives"
NON_ZERO_LABEL = "non-zero"
ALL_LABEL = "all"

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
TPM_FILE = "<tpm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG10_SCATTER_MIN = "--log10-scatter-min"
LOG10_SCATTER_MAX = "--log10-scatter-max"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    opt.validate_file_option(options[TPM_FILE], "Could not open TPM file")
    options[SCATTER_MAX] = opt.validate_int_option(
        options[SCATTER_MAX], "Invalid maximum value for scatter plot axes")
    options[LOG10_SCATTER_MIN] = opt.validate_int_option(
        options[LOG10_SCATTER_MIN],
        "Invalid minimum value for log10 scatter plot axes")
    options[LOG10_SCATTER_MAX] = opt.validate_int_option(
        options[LOG10_SCATTER_MAX],
        "Invalid maximum value for log10 scatter plot axes")
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Read TPMs into a data frame
logger.info("Reading TPMs...")
tpms = pd.read_csv(options[TPM_FILE])

# Determine whether each TPM measurement is a true/false positive/negative.
# For our purposes, marking an TPM as positive or negative is determined
# by whether it is greater or less than an "isoform not present" cutoff value.
logger.info("Marking positives and negatives...")
t.mark_positives_and_negatives(tpms)

# Calculate the percent error (positive or negative) of the calculated TPM
# values from the real values
logger.info("Calculating TPM percentage errors...")
t.calculate_percent_error(tpms)

# Calculate log ratio of calculated and real TPMs
logger.info("Calculating TPM log ratios...")
t.calculate_log_ratios(tpms)

# Apply various classification measures to the TPM data
logger.info("Applying classifiers...")
clsfrs = classifiers.get_classifiers()
t.apply_classifiers(tpms, clsfrs)

# Get a data frame containing only true positive TPMs
tp_tpms = t.get_true_positives(tpms)

# Write statistics pertaining to the set of quantified transcripts as a whole.
logger.info("Writing overall statistics...")


def add_overall_stats(stats):
    for param in parameters.get_parameters():
        stats[param.name] = options[param.option_name]

stats_types = statistics.get_statistics()

if options[OUT_FILE_BASENAME]:
    stats = t.get_stats(tpms, tp_tpms, stats_types)
    add_overall_stats(stats)

    stats_file_name = \
        statistics.get_stats_file(".", options[OUT_FILE_BASENAME])
    statistics.write_stats_data(stats_file_name, stats, index=False)

# Write statistics for TPMS stratified by various classification measures
logger.info("Writing statistics for stratified TPMs")

non_zero = tpms[(tpms[t.REAL_TPM] > 0) & (tpms[t.CALCULATED_TPM] > 0)]

if options[OUT_FILE_BASENAME]:
    for classifier in clsfrs:
        if classifier.produces_grouped_stats():
            column_name = classifier.get_column_name()
            stats = t.get_grouped_stats(tpms, tp_tpms,
                                        column_name, stats_types)
            add_overall_stats(stats)

            stats_file_name = statistics.get_stats_file(
                ".", options[OUT_FILE_BASENAME], classifier)
            statistics.write_stats_data(stats_file_name, stats)

        elif classifier.produces_distribution_plots():
            for ascending in [True, False]:
                stats = t.get_distribution_stats(
                    non_zero, tp_tpms, classifier, ascending)
                add_overall_stats(stats)

                stats_file_name = statistics.get_stats_file(
                    ".", options[OUT_FILE_BASENAME], classifier, ascending)
                statistics.write_stats_data(stats_file_name, stats)

# Make a scatter plot of log transformed calculated vs real TPMs
logger.info("Plotting graphs...")

if options[OUT_FILE_BASENAME]:
    plot.log_tpm_scatter_plot(
        tp_tpms, options[OUT_FILE_BASENAME],
        options[parameters.QUANT_METHOD.option_name], TRUE_POSITIVES_LABEL)

# Make boxplots of log ratios stratified by various classification measures
# (e.g. the number of transcripts per-originating gene of each transcript)
# TODO: need something in graph titles to indicate filter being applied
more_than_100_filter = lambda x: len(x[t.REAL_TPM]) > 100
tpm_infos = [(non_zero, NON_ZERO_LABEL), (tp_tpms, TRUE_POSITIVES_LABEL)]

if options[OUT_FILE_BASENAME]:
    classifiers = [c for c in clsfrs if c.produces_grouped_stats()]
    filters = [None, more_than_100_filter]
    for c, f, ti in itertools.product(classifiers, filters, tpm_infos):
        plot.log_ratio_boxplot(
            ti[0], options[OUT_FILE_BASENAME],
            parameters.QUANT_METHOD.option_name, ti[1], c, f)

# Make plots showing the percentage of isoforms above or below threshold values
# according to various classification measures
if options[OUT_FILE_BASENAME]:
    classifiers = [c for c in clsfrs if c.produces_distribution_plots()]
    ascending = [True, False]
    for c, asc, ti in itertools.product(classifiers, ascending, tpm_infos):
        plot.plot_cumulative_transcript_distribution(
            ti[0], options[OUT_FILE_BASENAME],
            options[parameters.QUANT_METHOD.option_name], ti[1], c, asc)
