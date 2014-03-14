#!/usr/bin/python

"""Usage:
    analyse_quantification_data [--scatter-max=<scatter-max-val>] [--log2-scatter-min=<log2-scatter-min-val>] [--log2-scatter-max=<log2-scatter-max-val>] <quant-method> <read-length> <read-depth> <paired-end> <errors> <fpkm-file> <out-file>

-h --help                                  Show this message.
-v --version                               Show version.
--scatter-max=<scatter-max-val>            Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log2-scatter-min=<log2-scatter-min-val>  Minimum x and y values for log2 scatter plot; a value of 0 means do not impose a minimum [default: 0].
--log2-scatter-max=<log2-scatter-max-val>  Maximum x and y values for log2 scatter plot; a value of 0 means do not impose a maximum [default: 0].
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

import matplotlib.pyplot as plt
import numpy as np
import options as opt
import pandas as pd

QUANT_METHOD = "quant-method"
READ_DEPTH = "read-depth"
READ_LENGTH = "read-length"
PAIRED_END = "paired-end"
ERRORS = "errors"
NUM_FPKMS = "num-fpkms"
TP_NUM_FPKMS = "tp-num-fpkms"
TP_COUNT = "tp-count"
REAL_FPKM = "real-fpkm"
LOG2_REAL_FPKM = "log2-real-fpkm"
CALCULATED_FPKM = "calc-fpkm"
LOG2_CALCULATED_FPKM = "log2-calc-fpkm"
LOG2_RATIO = "log-ratio"
TP_LOG2_RATIO_MEAN = "tp-log-ratio-mean"
TP_LOG2_RATIO_STD = "tp-log-ratio-std"
TP_LOG2_RATIO_MEDIAN = "tp-log-ratio-med"
TP_LOG2_FPKM_RHO = "tp-log-fpkm-rho"
TRANSCRIPT_COUNT = "num-transcripts"
PERCENT_ERROR = "percent-error"
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
NON_ZERO_LABEL = "non zero"
ALL_LABEL = "all"
NO_FILTER_LABEL = "no filter"

NOT_PRESENT_CUTOFF = 0.1

opt_string = lambda x: "<{x}>".format(x=x)

FPKM_FILE = "<fpkm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG2_SCATTER_MIN = "--log2-scatter-min"
LOG2_SCATTER_MAX = "--log2-scatter-max"
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
    options[LOG2_SCATTER_MIN] = opt.validate_int_option(
        options[LOG2_SCATTER_MIN],
        "Invalid minimum value for log2 scatter plot axes")
    options[LOG2_SCATTER_MAX] = opt.validate_int_option(
        options[LOG2_SCATTER_MAX],
        "Invalid maximum value for log2 scatter plot axes")
except SchemaError as exc:
    exit(exc.code)

fpkms = pd.read_csv(options[FPKM_FILE])

# Determine whether each FPKM measurement is a true/false positive/negative.
# For our purposes, marking an FPKM as positive or negative is determined
# by whether it is greater or less than an "isoform not present" cutoff value.

fpkms[FALSE_NEGATIVE] = (fpkms[REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
                        (fpkms[CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[FALSE_POSITIVE] = (fpkms[CALCULATED_FPKM] > NOT_PRESENT_CUTOFF) & \
                        (fpkms[REAL_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[TRUE_NEGATIVE] = (fpkms[REAL_FPKM] <= NOT_PRESENT_CUTOFF) & \
                       (fpkms[CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
fpkms[TRUE_POSITIVE] = (fpkms[REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
                       (fpkms[CALCULATED_FPKM] > NOT_PRESENT_CUTOFF)


def get_sensitivity(fpkms):
    num_tp = len(fpkms[fpkms[TRUE_POSITIVE]])
    num_fn = len(fpkms[fpkms[FALSE_NEGATIVE]])
    return float(num_tp) / (num_tp + num_fn)


def get_specificity(fpkms):
    num_fp = len(fpkms[fpkms[FALSE_POSITIVE]])
    num_tn = len(fpkms[fpkms[TRUE_NEGATIVE]])
    return float(num_tn) / (num_tn + num_fp)

# Calculate the percent error (positive or negative) of the calculated FPKM
# values from the real values
fpkms[PERCENT_ERROR] = \
    100 * (fpkms[CALCULATED_FPKM] - fpkms[REAL_FPKM]) / fpkms[REAL_FPKM]


def get_error_fraction(fpkms, error_percent):
    num_errors = len(fpkms[abs(fpkms[PERCENT_ERROR]) > error_percent])
    return float(num_errors) / len(fpkms)

# Calculate log ratio of calculated and real FPKMs
fpkms[LOG2_REAL_FPKM] = np.log2(fpkms[REAL_FPKM])
fpkms[LOG2_CALCULATED_FPKM] = np.log2(fpkms[CALCULATED_FPKM])
fpkms[LOG2_RATIO] = \
    fpkms[LOG2_CALCULATED_FPKM] - fpkms[LOG2_REAL_FPKM]

tp_fpkms = fpkms[fpkms[TRUE_POSITIVE]]

# Write log ratio summary statistics stratified by the number of transcripts
# per originating gene
with open(options[OUT_FILE_BASENAME] + "_stats.csv", "w") as out_file:
    # Spearman correlation coefficient between real and calculated FPKMs.
    rho = tp_fpkms[LOG2_CALCULATED_FPKM].corr(
        tp_fpkms[LOG2_REAL_FPKM], method='spearman')

    # The median percent error - i.e. the median of the percent errors of
    # the calculated values from the real ones
    tp_mpe = tp_fpkms[PERCENT_ERROR].median()

    stats_dict = {
        QUANT_METHOD: options[QUANT_METHOD_OPT],
        READ_LENGTH: options[READ_LENGTH_OPT],
        READ_DEPTH: options[READ_DEPTH_OPT],
        PAIRED_END: options[PAIRED_END_OPT],
        ERRORS: options[ERRORS_OPT],
        NUM_FPKMS: len(fpkms),
        TP_NUM_FPKMS: len(tp_fpkms),
        TP_LOG2_FPKM_RHO: rho,
        TP_MEDIAN_PERCENT_ERROR: tp_mpe,
        SENSITIVITY: get_sensitivity(fpkms),
        SPECIFICITY: get_specificity(fpkms),
        TP_ERROR_FRACTION: get_error_fraction(tp_fpkms, 10)
    }

    stats = pd.DataFrame([stats_dict])
    stats.to_csv(out_file, float_format="%.5f", index=False)

# Group transcripts by the number of transcripts for their originating gene,
# and discard those groups with fewer than 100 members
grouped = fpkms.groupby(TRANSCRIPT_COUNT)
tp_grouped = tp_fpkms.groupby(TRANSCRIPT_COUNT)

with open(options[OUT_FILE_BASENAME] +
          "_stats_by_num_transcripts_per_gene.csv", "w") as out_file:

    summary = grouped.describe()
    main_stats = summary[REAL_FPKM].unstack()
    main_stats = main_stats.drop(
        ["mean", "std", "min", "25%", "50%", "75%", "max"], 1)

    main_stats[SENSITIVITY] = grouped.apply(get_sensitivity)
    main_stats[SPECIFICITY] = grouped.apply(get_specificity)

    tp_summary = tp_grouped.describe()
    tp_stats = tp_summary[LOG2_RATIO].unstack()
    tp_stats = tp_stats.drop(["min", "25%", "75%", "max"], 1)
    tp_stats = tp_stats.rename(columns={
        "count": TP_COUNT,
        "mean": TP_LOG2_RATIO_MEAN,
        "std": TP_LOG2_RATIO_STD,
        "50%": TP_LOG2_RATIO_MEDIAN})

    tp_stats[TP_LOG2_FPKM_RHO] = tp_grouped.apply(
        lambda x: x[LOG2_CALCULATED_FPKM].corr(x[LOG2_REAL_FPKM],
                                               method="spearman"))

    pe_stats = tp_summary[PERCENT_ERROR].unstack()
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

# Make scatter plots of calculated vs real FPKMs

space_to_underscore = lambda x: x.replace(' ', '_')


def fpkm_scatter_plot(name, fpkms, max_val):
    plt.figure()
    plt.scatter(fpkms[REAL_FPKM].values, fpkms[CALCULATED_FPKM].values)
    name = options[QUANT_METHOD_OPT] + " " + name
    plt.suptitle("Scatter plot of calculated vs real FPKMs: " + name)
    plt.xlabel("Real FPKM")
    plt.ylabel("Calculated FPKM")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    if max_val != 0:
        plt.xlim(xmax=max_val)
        plt.ylim(ymax=max_val)
    plt.savefig(options[OUT_FILE_BASENAME] + "_" + space_to_underscore(name) +
                "_scatter.pdf", format="pdf")


fpkm_scatter_plot(ALL_LABEL, fpkms, options[SCATTER_MAX])
fpkm_scatter_plot(TRUE_POSITIVES_LABEL, tp_fpkms, options[SCATTER_MAX])

# Make a scatter plot of log transformed calculated vs real FPKMs


def log_fpkm_scatter_plot(name, fpkms, min_val, max_val):
    plt.figure()
    plt.scatter(fpkms[LOG2_REAL_FPKM].values,
                fpkms[LOG2_CALCULATED_FPKM].values)
    name = options[QUANT_METHOD_OPT] + " " + name
    plt.suptitle("Scatter plot of log calculated vs real FPKMs: " + name)
    plt.xlabel("Log2 real FPKM")
    plt.ylabel("Log2 calculated FPKM")
    if min_val == 0:
        min_val = np.log2(NOT_PRESENT_CUTOFF)
    plt.xlim(xmin=min_val)
    plt.ylim(ymin=min_val)
    if max_val != 0:
        plt.xlim(xmax=max_val)
        plt.ylim(ymax=max_val)
    plt.savefig(options[OUT_FILE_BASENAME] + "_" + space_to_underscore(name) +
                "_log2_scatter.pdf", format="pdf")

log_fpkm_scatter_plot(TRUE_POSITIVES_LABEL, tp_fpkms,
                      options[LOG2_SCATTER_MIN], options[LOG2_SCATTER_MAX])

# Make boxplots of log ratios stratified by the number of transcripts per
# originating gene.


def log_ratio_boxplot(name_elements, grouping_label, fpkms, filter=None):
    if filter:
        grouped_fpkms = fpkms.groupby(filter[0])
        fpkms = grouped_fpkms.filter(filter[1])
    bp = fpkms.boxplot(column=LOG2_RATIO, by=TRANSCRIPT_COUNT, sym="")
    bp.set_title("")

    name = " ".join([options[QUANT_METHOD_OPT]] + name_elements)
    plt.suptitle("Log ratios of calculated to real FPKMs:" + name)

    plt.xlabel(grouping_label)
    plt.ylabel("Log ratio (calculated/real FPKM)")
    plt.savefig(options[OUT_FILE_BASENAME] + "_" + space_to_underscore(name) +
                "_boxplot.pdf", format="pdf")

more_than_100_filter = lambda x: len(x[REAL_FPKM]) > 100

non_zero = fpkms[(fpkms[REAL_FPKM] > 0) & (fpkms[CALCULATED_FPKM] > 0)]
log_ratio_boxplot([NON_ZERO_LABEL, NO_FILTER_LABEL],
                  TRANSCRIPT_COUNT_LABEL, non_zero)
log_ratio_boxplot([NON_ZERO_LABEL], TRANSCRIPT_COUNT_LABEL, non_zero,
                  filter=(TRANSCRIPT_COUNT, more_than_100_filter))

log_ratio_boxplot([TRUE_POSITIVES_LABEL, NO_FILTER_LABEL],
                  TRANSCRIPT_COUNT_LABEL, tp_fpkms)
log_ratio_boxplot([TRUE_POSITIVES_LABEL],
                  TRANSCRIPT_COUNT_LABEL, tp_fpkms,
                  filter=(TRANSCRIPT_COUNT, more_than_100_filter))
