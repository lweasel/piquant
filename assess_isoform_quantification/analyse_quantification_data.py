#!/usr/bin/python

"""Usage:
    analyse_quantification_data [--scatter-max=<scatter-max-val>] [--log10-scatter-min=<log10-scatter-min-val>] [--log10-scatter-max=<log10-scatter-max-val>] <quant-method> <read-length> <read-depth> <paired-end> <errors> <fpkm-file> [<out-file>]

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
<fpkm-file>                                File containing real and calculated FPKMs.
<out-file>                                 Basename for output graph and data files.
"""

from docopt import docopt
from schema import SchemaError

import fpkms as f
import matplotlib.pyplot as plt
import numpy as np
import options as opt
import pandas as pd
import seaborn as sb
import classifiers

QUANT_METHOD = "quant-method"
READ_DEPTH = "read-depth"
READ_LENGTH = "read-length"
PAIRED_END = "paired-end"
ERRORS = "errors"
NUM_FPKMS = "num-fpkms"
TP_NUM_FPKMS = "tp-num-fpkms"
TP_COUNT = "tp-count"
TP_LOG10_RATIO_MEAN = "tp-log-ratio-mean"
TP_LOG10_RATIO_STD = "tp-log-ratio-std"
TP_LOG10_RATIO_MEDIAN = "tp-log-ratio-med"
TP_LOG10_FPKM_RHO = "tp-log-fpkm-rho"
TP_MEDIAN_PERCENT_ERROR = "tp-median-percent-error"
FALSE_POSITIVE = "false-pos"
FALSE_NEGATIVE = "false-neg"
TRUE_POSITIVE = "true-pos"
TRUE_NEGATIVE = "true-neg"
SENSITIVITY = "sensitivity"
SPECIFICITY = "specificity"
TP_ERROR_FRACTION = "tp-error-frac"

TRANSCRIPT_COUNT_LABEL = "No. transcripts per gene"
TRUE_POSITIVES_LABEL = "true positives"
NON_ZERO_LABEL = "non-zero"
ALL_LABEL = "all"
NO_FILTER_LABEL = "no filter"

NOT_PRESENT_CUTOFF = 0.1
CUMULATIVE_DISTRIBUTION_POINTS = 20

opt_string = lambda x: "<{x}>".format(x=x)

FPKM_FILE = "<fpkm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG10_SCATTER_MIN = "--log10-scatter-min"
LOG10_SCATTER_MAX = "--log10-scatter-max"
QUANT_METHOD_OPT = opt_string(QUANT_METHOD)
READ_DEPTH_OPT = opt_string(READ_DEPTH)
READ_LENGTH_OPT = opt_string(READ_LENGTH)
PAIRED_END_OPT = opt_string(PAIRED_END)
ERRORS_OPT = opt_string(ERRORS)

# Read in command-line options
options = docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    opt.validate_file_option(options[FPKM_FILE], "Could not open FPKM file")
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

fpkms = pd.read_csv(options[FPKM_FILE])

# Determine whether each FPKM measurement is a true/false positive/negative.
# For our purposes, marking an FPKM as positive or negative is determined
# by whether it is greater or less than an "isoform not present" cutoff value.

fpkms[FALSE_NEGATIVE] = (fpkms[f.REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
                        (fpkms[f.CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[FALSE_POSITIVE] = (fpkms[f.CALCULATED_FPKM] > NOT_PRESENT_CUTOFF) & \
                        (fpkms[f.REAL_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[TRUE_NEGATIVE] = (fpkms[f.REAL_FPKM] <= NOT_PRESENT_CUTOFF) & \
                       (fpkms[f.CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[TRUE_POSITIVE] = (fpkms[f.REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
                       (fpkms[f.CALCULATED_FPKM] > NOT_PRESENT_CUTOFF)

# Calculate the percent error (positive or negative) of the calculated FPKM
# values from the real values
fpkms[f.PERCENT_ERROR] = \
    100 * (fpkms[f.CALCULATED_FPKM] - fpkms[f.REAL_FPKM]) / fpkms[f.REAL_FPKM]

# Calculate log ratio of calculated and real FPKMs
fpkms[f.LOG10_REAL_FPKM] = np.log10(fpkms[f.REAL_FPKM])
fpkms[f.LOG10_CALCULATED_FPKM] = np.log10(fpkms[f.CALCULATED_FPKM])
fpkms[f.LOG10_RATIO] = \
    fpkms[f.LOG10_CALCULATED_FPKM] - fpkms[f.LOG10_REAL_FPKM]

cfrs = [s() for s in classifiers.CLASSIFIERS]

for classifier in cfrs:
    column_name = classifier.get_column_name()
    fpkms[column_name] = fpkms.apply(
        classifier.get_classification_value, axis=1)

tp_fpkms = fpkms[fpkms[TRUE_POSITIVE]]


def add_overall_stats(stats, fpkms, tp_fpkms):
    stats[QUANT_METHOD] = options[QUANT_METHOD_OPT]
    stats[READ_LENGTH] = options[READ_LENGTH_OPT]
    stats[READ_DEPTH] = options[READ_DEPTH_OPT]
    stats[PAIRED_END] = options[PAIRED_END_OPT]
    stats[ERRORS] = options[ERRORS_OPT]
    stats[NUM_FPKMS] = len(fpkms)
    stats[TP_NUM_FPKMS] = len(tp_fpkms)


# Write statistics pertaining to the set of quantified transcripts as a whole.
if options[OUT_FILE_BASENAME]:
    stats = f.get_stats(fpkms, tp_fpkms)
    add_overall_stats(stats, fpkms, tp_fpkms)

    with open(options[OUT_FILE_BASENAME] + "_stats.csv", "w") as out_file:
        stats.to_csv(out_file, float_format="%.5f", index=False)

space_to_underscore = lambda x: x.replace(' ', '_')

# Write statistics for FPKMS stratified by various classification measures

if options[OUT_FILE_BASENAME]:
    for classifier in [c for c in cfrs if c.produces_grouped_stats()]:
        column_name = classifier.get_column_name()
        stats = f.get_grouped_stats(fpkms, tp_fpkms, column_name)
        add_overall_stats(stats, fpkms, tp_fpkms)

        stats_file_name = options[OUT_FILE_BASENAME] + "_stats_by_" + \
            space_to_underscore(column_name) + ".csv"

        with open(stats_file_name, "w") as out_file:
            stats.to_csv(out_file, float_format="%.5f")

# Make a scatter plot of log transformed calculated vs real FPKMs


def log_fpkm_scatter_plot(name, fpkms, min_val, max_val):
    plt.figure()
    plt.scatter(fpkms[f.LOG10_REAL_FPKM].values,
                fpkms[f.LOG10_CALCULATED_FPKM].values,
                c="lightblue", alpha=0.4)

    plt.suptitle("Scatter plot of log calculated vs real FPKMs: " +
                 options[QUANT_METHOD_OPT] + " " + name)
    plt.xlabel("Log10 real FPKM")
    plt.ylabel("Log10 calculated FPKM")

    if min_val == 0:
        min_val = np.log10(NOT_PRESENT_CUTOFF) - 0.2
    plt.xlim(xmin=min_val)
    plt.ylim(ymin=min_val)
    if max_val != 0:
        plt.xlim(xmax=max_val)
        plt.ylim(ymax=max_val)

    plt.savefig(options[OUT_FILE_BASENAME] + "_" + space_to_underscore(name) +
                "_log10_scatter.pdf", format="pdf")
    plt.close()

if options[OUT_FILE_BASENAME]:
    log_fpkm_scatter_plot(TRUE_POSITIVES_LABEL, tp_fpkms,
                          options[LOG10_SCATTER_MIN],
                          options[LOG10_SCATTER_MAX])

# Make boxplots of log ratios stratified by various classification measures
# (e.g. the number of transcripts per-originating gene of each transcript)


def log_ratio_boxplot(name_elements, classifier, fpkms, filter=None):
    grouping_column = classifier.get_column_name()
    if filter:
        grouped_fpkms = fpkms.groupby(grouping_column)
        fpkms = grouped_fpkms.filter(filter)

    plt.figure()
    sb.boxplot(fpkms[f.LOG10_RATIO], groupby=fpkms[grouping_column],
               sym='', color='lightblue')

    plt.suptitle("Log ratios of calculated to real FPKMs: " +
                 ", ".join([options[QUANT_METHOD_OPT]] + name_elements))

    plt.xlabel(grouping_column[:1].upper() + grouping_column[1:])
    plt.ylabel("Log ratio (calculated/real FPKM)")

    locs, labels = plt.xticks()
    plt.xticks(locs, classifier.get_value_labels(len(labels)))

    filename = options[OUT_FILE_BASENAME] + "_" \
        + space_to_underscore(" ".join([grouping_column] + name_elements)) \
        + "_boxplot.pdf"
    plt.savefig(filename, format="pdf")
    plt.close()

more_than_100_filter = lambda x: len(x[f.REAL_FPKM]) > 100

non_zero = fpkms[(fpkms[f.REAL_FPKM] > 0) & (fpkms[f.CALCULATED_FPKM] > 0)]

if options[OUT_FILE_BASENAME]:
    for classifier in [c for c in cfrs if c.produces_grouped_stats()]:
        log_ratio_boxplot([NON_ZERO_LABEL, NO_FILTER_LABEL],
                          classifier, non_zero)
        log_ratio_boxplot([NON_ZERO_LABEL], classifier, non_zero,
                          filter=more_than_100_filter)

        log_ratio_boxplot([TRUE_POSITIVES_LABEL, NO_FILTER_LABEL],
                          classifier, tp_fpkms)
        log_ratio_boxplot([TRUE_POSITIVES_LABEL],
                          classifier, tp_fpkms,
                          filter=more_than_100_filter)

# Make plots showing the percentage of isoforms above or below threshold values
# according to various classification measures


def plot_cumulative_transcript_distribution(
        name_elements, fpkms, classifier, ascending):

    values = fpkms.apply(classifier.get_value, axis=1)
    values.sort(ascending=ascending)

    xbounds = classifier.get_distribution_plot_range()
    if xbounds is None:
        xbounds = (values.min(), values.max())

    xvals = np.linspace(xbounds[0], xbounds[1], CUMULATIVE_DISTRIBUTION_POINTS)

    size = float(len(values))
    yvals = [100 * len(values[values < x if ascending else values > x]) / size
             for x in xvals]

    plt.figure()
    plt.plot(xvals, yvals, '-o')

    plt.ylim(ymin=-2.5, ymax=102.5)

    xmargin = (xbounds[1] - xbounds[0]) / 40.0
    plt.xlim(xmin=xbounds[0]-xmargin, xmax=xbounds[1]+xmargin)

    grouping_column = classifier.get_column_name()
    capitalized = grouping_column[:1].upper() + grouping_column[1:]
    plt.xlabel(capitalized)
    plt.ylabel("Percentage of isoforms " +
               ("less" if ascending else "greater") + " than threshold")

    plt.suptitle(capitalized + " threshold: " +
                 ", ".join([options[QUANT_METHOD_OPT]] + name_elements))

    filename = options[OUT_FILE_BASENAME] + "_" \
        + space_to_underscore(" ".join([grouping_column] + name_elements)) \
        + "_" + ("asc" if ascending else "desc") + "_distribution.pdf"
    plt.savefig(filename, format="pdf")
    plt.close()

if options[OUT_FILE_BASENAME]:
    for classifier in [s for s in cfrs if s.produces_distribution_plots()]:
        for ascending in [True, False]:
            plot_cumulative_transcript_distribution(
                [NON_ZERO_LABEL], non_zero, classifier, ascending)
            plot_cumulative_transcript_distribution(
                [TRUE_POSITIVES_LABEL], tp_fpkms, classifier, ascending)
