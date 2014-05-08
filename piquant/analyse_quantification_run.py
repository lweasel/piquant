#!/usr/bin/python

# TODO: add logging
# TODO: method of specifying options here (e.g. <bias>) is horrible

"""Usage:
    analyse_quantification_run [--scatter-max=<scatter-max-val>] [--log10-scatter-min=<log10-scatter-min-val>] [--log10-scatter-max=<log10-scatter-max-val>] <quant-method> <read-length> <read-depth> <paired-end> <errors> <bias> <tpm-file> [<out-file>]

-h --help                                  Show this message.
-v --version                               Show version.
--scatter-max=<scatter-max-val>            Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log10-scatter-min=<log10-scatter-min-val>  Minimum x and y values for log10 scatter plot; a value of 0 means do not impose a minimum [default: 0].
--log10-scatter-max=<log10-scatter-max-val>  Maximum x and y values for log10 scatter plot; a value of 0 means do not impose a maximum [default: 0].
<quant-method>                             Method used to quantify transcript abundances.
<read-length>                              The length of sequence reads.
<read-depth>                               The depth of reads sequenced across the transcriptome.
<paired-end>                               Whether paired-end sequence reads were used.
<errors>                                   Whether the reads contain sequencing errors.
<bias>                                     Whether the reads contain sequence bias.
<tpm-file>                                File containing real and calculated TPMs.
<out-file>                                 Basename for output graph and data files.
"""

from docopt import docopt
from schema import SchemaError

import classifiers
import ordutils.options as opt
import pandas as pd
import statistics
import tpms as t
import tpms_plotting as plot

QUANT_METHOD = "quant-method"
READ_DEPTH = "read-depth"
READ_LENGTH = "read-length"
PAIRED_END = "paired-end"
ERRORS = "errors"
BIAS = "bias"

TRANSCRIPT_COUNT_LABEL = "No. transcripts per gene"
TRUE_POSITIVES_LABEL = "true positives"
NON_ZERO_LABEL = "non-zero"
ALL_LABEL = "all"

opt_string = lambda x: "<{x}>".format(x=x)

TPM_FILE = "<tpm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG10_SCATTER_MIN = "--log10-scatter-min"
LOG10_SCATTER_MAX = "--log10-scatter-max"
QUANT_METHOD_OPT = opt_string(QUANT_METHOD)
READ_DEPTH_OPT = opt_string(READ_DEPTH)
READ_LENGTH_OPT = opt_string(READ_LENGTH)
PAIRED_END_OPT = opt_string(PAIRED_END)
ERRORS_OPT = opt_string(ERRORS)
BIAS_OPT = opt_string(BIAS)

# Read in command-line options
options = docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    opt.validate_file_option(options[TPM_FILE], "Could not open TPM file")
    options[SCATTER_MAX] = opt.validate_int_option(
        options[SCATTER_MAX], "Invalid maximum value for scatter plot axes")
    options[LOG10_SCATTER_MIN] = opt.validate_int_option(
        options[LOG10_SCATTER_MIN],
        "Invalid minimum value for log10 scatter plot axes")
    options[LOG10_SCATTER_MAX] = opt.validate_int_option(
        options[LOG10_SCATTER_MAX],
        "Invalid maximum value for log10 scatter plot axes")
except SchemaError as exc:
    exit(exc.code)

# Read TPMs into a data frame
tpms = pd.read_csv(options[TPM_FILE])

# Determine whether each TPM measurement is a true/false positive/negative.
# For our purposes, marking an TPM as positive or negative is determined
# by whether it is greater or less than an "isoform not present" cutoff value.
t.mark_positives_and_negatives(tpms)

# Calculate the percent error (positive or negative) of the calculated TPM
# values from the real values
t.calculate_percent_error(tpms)

# Calculate log ratio of calculated and real TPMs
t.calculate_log_ratios(tpms)

# Apply various classification measures to the TPM data
clsfrs = classifiers.get_classifiers()
t.apply_classifiers(tpms, clsfrs)

# Get a data frame containing only true positive TPMs
tp_tpms = t.get_true_positives(tpms)

# Write statistics pertaining to the set of quantified transcripts as a whole.


def add_overall_stats(stats):
    stats[QUANT_METHOD] = options[QUANT_METHOD_OPT]
    stats[READ_LENGTH] = options[READ_LENGTH_OPT]
    stats[READ_DEPTH] = options[READ_DEPTH_OPT]
    stats[PAIRED_END] = options[PAIRED_END_OPT]
    stats[ERRORS] = options[ERRORS_OPT]
    stats[BIAS] = options[BIAS_OPT]

if options[OUT_FILE_BASENAME]:
    stats = t.get_stats(tpms, tp_tpms)
    add_overall_stats(stats)

    stats_file_name = \
        statistics.get_stats_file(".", options[OUT_FILE_BASENAME])
    statistics.write_stats_data(stats_file_name, stats, index=False)

# Write statistics for TPMS stratified by various classification measures
non_zero = tpms[(tpms[t.REAL_TPM] > 0) & (tpms[t.CALCULATED_TPM] > 0)]

if options[OUT_FILE_BASENAME]:
    for classifier in clsfrs:
        if classifier.produces_grouped_stats():
            column_name = classifier.get_column_name()
            stats = t.get_grouped_stats(tpms, tp_tpms, column_name)
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
TP_PLOT_OPTIONS = plot.PlotOptions(
    options[QUANT_METHOD_OPT], TRUE_POSITIVES_LABEL,
    options[OUT_FILE_BASENAME])

if options[OUT_FILE_BASENAME]:
    plot.log_tpm_scatter_plot(tp_tpms, TP_PLOT_OPTIONS)

# Make boxplots of log ratios stratified by various classification measures
# (e.g. the number of transcripts per-originating gene of each transcript)
# TODO: need something in graph titles to indicate filter being applied
more_than_100_filter = lambda x: len(x[t.REAL_TPM]) > 100

NON_ZERO_PLOT_OPTIONS = plot.PlotOptions(
    options[QUANT_METHOD_OPT], NON_ZERO_LABEL,
    options[OUT_FILE_BASENAME])

if options[OUT_FILE_BASENAME]:
    for classifier in [c for c in clsfrs if c.produces_grouped_stats()]:
        plot.log_ratio_boxplot(
            non_zero, NON_ZERO_PLOT_OPTIONS, classifier)
        plot.log_ratio_boxplot(
            non_zero, NON_ZERO_PLOT_OPTIONS, classifier,
            filter=more_than_100_filter)

        plot.log_ratio_boxplot(
            tp_tpms, TP_PLOT_OPTIONS, classifier)
        plot.log_ratio_boxplot(
            tp_tpms, TP_PLOT_OPTIONS, classifier,
            filter=more_than_100_filter)

# Make plots showing the percentage of isoforms above or below threshold values
# according to various classification measures
if options[OUT_FILE_BASENAME]:
    for classifier in [c for c in clsfrs if c.produces_distribution_plots()]:
        for ascending in [True, False]:
            plot.plot_cumulative_transcript_distribution(
                non_zero, NON_ZERO_PLOT_OPTIONS, classifier, ascending)
            plot.plot_cumulative_transcript_distribution(
                tp_tpms, TP_PLOT_OPTIONS, classifier, ascending)
