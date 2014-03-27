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
import stratifiers

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


def get_sensitivity(fpkms):
    num_tp = len(fpkms[fpkms[TRUE_POSITIVE]])
    num_fn = len(fpkms[fpkms[FALSE_NEGATIVE]])
    if num_tp + num_fn == 0:
        return 1
    return float(num_tp) / (num_tp + num_fn)


def get_specificity(fpkms):
    num_fp = len(fpkms[fpkms[FALSE_POSITIVE]])
    num_tn = len(fpkms[fpkms[TRUE_NEGATIVE]])
    if num_fp + num_tn == 0:
        return 1
    return float(num_tn) / (num_tn + num_fp)

# Calculate the percent error (positive or negative) of the calculated FPKM
# values from the real values
fpkms[f.PERCENT_ERROR] = \
    100 * (fpkms[f.CALCULATED_FPKM] - fpkms[f.REAL_FPKM]) / fpkms[f.REAL_FPKM]


def get_error_fraction(fpkms, error_percent):
    num_errors = len(fpkms[abs(fpkms[f.PERCENT_ERROR]) > error_percent])
    return float(num_errors) / len(fpkms)

# Calculate log ratio of calculated and real FPKMs
fpkms[f.LOG10_REAL_FPKM] = np.log10(fpkms[f.REAL_FPKM])
fpkms[f.LOG10_CALCULATED_FPKM] = np.log10(fpkms[f.CALCULATED_FPKM])
fpkms[f.LOG10_RATIO] = \
    fpkms[f.LOG10_CALCULATED_FPKM] - fpkms[f.LOG10_REAL_FPKM]

strats = [s() for s in stratifiers.STRATIFIERS]

for stratifier in strats:
    column_name = stratifier.get_column_name()
    fpkms[column_name] = fpkms.apply(
        stratifier.get_stratification_value, axis=1)

tp_fpkms = fpkms[fpkms[TRUE_POSITIVE]]

# Write statistics pertaining to the set of quantified transcripts as a whole.
if options[OUT_FILE_BASENAME]:
    with open(options[OUT_FILE_BASENAME] + "_stats.csv", "w") as out_file:
        # Spearman correlation coefficient between real and calculated FPKMs.
        rho = tp_fpkms[f.LOG10_CALCULATED_FPKM].corr(
            tp_fpkms[f.LOG10_REAL_FPKM], method='spearman')

        # The median percent error - i.e. the median of the percent errors of
        # the calculated values from the real ones
        tp_mpe = tp_fpkms[f.PERCENT_ERROR].median()

        stats_dict = {
            QUANT_METHOD: options[QUANT_METHOD_OPT],
            READ_LENGTH: options[READ_LENGTH_OPT],
            READ_DEPTH: options[READ_DEPTH_OPT],
            PAIRED_END: options[PAIRED_END_OPT],
            ERRORS: options[ERRORS_OPT],
            NUM_FPKMS: len(fpkms),
            TP_NUM_FPKMS: len(tp_fpkms),
            TP_LOG10_FPKM_RHO: rho,
            TP_MEDIAN_PERCENT_ERROR: tp_mpe,
            SENSITIVITY: get_sensitivity(fpkms),
            SPECIFICITY: get_specificity(fpkms),
            TP_ERROR_FRACTION: get_error_fraction(tp_fpkms, 10)
        }

        stats = pd.DataFrame([stats_dict])
        stats.to_csv(out_file, float_format="%.5f", index=False)

space_to_underscore = lambda x: x.replace(' ', '_')

# Write statistics for FPKMS stratified by various stratification measures

if options[OUT_FILE_BASENAME]:
    for stratifier in [s for s in strats if s.produces_grouped_stats()]:
        column_name = stratifier.get_column_name()
        grouped = fpkms.groupby(column_name)
        tp_grouped = tp_fpkms.groupby(column_name)

        stats_file_name = options[OUT_FILE_BASENAME] + "_stats_by_" + \
            space_to_underscore(column_name) + ".csv"

        with open(stats_file_name, "w") as out_file:
            summary = grouped.describe()
            main_stats = summary[f.REAL_FPKM].unstack()
            main_stats = main_stats.drop(
                ["mean", "std", "min", "25%", "50%", "75%", "max"], 1)

            main_stats[SENSITIVITY] = grouped.apply(get_sensitivity)
            main_stats[SPECIFICITY] = grouped.apply(get_specificity)

            tp_summary = tp_grouped.describe()
            tp_stats = tp_summary[f.LOG10_RATIO].unstack()
            tp_stats = tp_stats.drop(["min", "25%", "75%", "max"], 1)
            tp_stats = tp_stats.rename(columns={
                "count": TP_COUNT,
                "mean": TP_LOG10_RATIO_MEAN,
                "std": TP_LOG10_RATIO_STD,
                "50%": TP_LOG10_RATIO_MEDIAN})

            tp_stats[TP_LOG10_FPKM_RHO] = tp_grouped.apply(
                lambda x: x[f.LOG10_CALCULATED_FPKM].corr(
                    x[f.LOG10_REAL_FPKM], method="spearman"))

            pe_stats = tp_summary[f.PERCENT_ERROR].unstack()
            tp_stats[TP_MEDIAN_PERCENT_ERROR] = pe_stats["50%"]

            tp_stats[TP_ERROR_FRACTION] = tp_grouped.apply(
                get_error_fraction, 10)

            output_stats = pd.concat([main_stats, tp_stats], axis=1)
            output_stats[QUANT_METHOD] = options[QUANT_METHOD_OPT]
            output_stats[READ_LENGTH] = options[READ_LENGTH_OPT]
            output_stats[READ_DEPTH] = options[READ_DEPTH_OPT]
            output_stats[PAIRED_END] = options[PAIRED_END_OPT]
            output_stats[ERRORS] = options[ERRORS_OPT]
            output_stats[TP_NUM_FPKMS] = len(fpkms)
            output_stats.to_csv(out_file, float_format="%.5f")

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

# Make boxplots of log ratios stratified by various stratification measures
# (e.g. the number of transcripts per-originating gene of each transcript)


def log_ratio_boxplot(name_elements, stratifier, fpkms, filter=None):
    grouping_column = stratifier.get_column_name()
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
    plt.xticks(locs, stratifier.get_value_labels(len(labels)))

    filename = options[OUT_FILE_BASENAME] + "_" \
        + space_to_underscore(" ".join([grouping_column] + name_elements)) \
        + "_boxplot.pdf"
    plt.savefig(filename, format="pdf")
    plt.close()

more_than_100_filter = lambda x: len(x[f.REAL_FPKM]) > 100

non_zero = fpkms[(fpkms[f.REAL_FPKM] > 0) & (fpkms[f.CALCULATED_FPKM] > 0)]

if options[OUT_FILE_BASENAME]:
    for stratifier in [s for s in strats if s.produces_grouped_stats()]:
        log_ratio_boxplot([NON_ZERO_LABEL, NO_FILTER_LABEL],
                          stratifier, non_zero)
        log_ratio_boxplot([NON_ZERO_LABEL], stratifier, non_zero,
                          filter=more_than_100_filter)

        log_ratio_boxplot([TRUE_POSITIVES_LABEL, NO_FILTER_LABEL],
                          stratifier, tp_fpkms)
        log_ratio_boxplot([TRUE_POSITIVES_LABEL],
                          stratifier, tp_fpkms,
                          filter=more_than_100_filter)

# Make plots showing the percentage of isoforms above or below threshold values
# according to various stratification measures


def plot_cumulative_transcript_distribution(
        name_elements, fpkms, stratifier, ascending):

    values = fpkms.apply(stratifier.get_value, axis=1)
    values.sort(ascending=ascending)

    xbounds = stratifier.get_distribution_plot_range()
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

    grouping_column = stratifier.get_column_name()
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
    for stratifier in [s for s in strats if s.produces_distribution_plots()]:
        for ascending in [True, False]:
            plot_cumulative_transcript_distribution(
                [NON_ZERO_LABEL], non_zero, stratifier, ascending)
            plot_cumulative_transcript_distribution(
                [TRUE_POSITIVES_LABEL], tp_fpkms, stratifier, ascending)
