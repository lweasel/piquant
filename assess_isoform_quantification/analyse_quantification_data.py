#!/usr/bin/python

"""Usage:
    analyse_quantification_data [--scatter-max=<scatter-max-val>] [--log2-scatter-min=<log2-scatter-min-val>] [--log2-scatter-max=<log2-scatter-max-val>] <quant-method> <read-length> <read-depth> <paired-end> <errors> <fpkm-file> <out-file>

-h --help                                  Show this message.
-v --version                               Show version.
--scatter-max=<scatter-max-val>            Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log2-scatter-min=<log2-scatter-min-val>  Minimum x and y values for log2 scatter plot; a value of 0 means do not impose a minimum [default: -10].
--log2-scatter-max=<log2-scatter-max-val>  Maximum x and y values for log2 scatter plot; a value of 0 means do not impose a maximum [default: 15].
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
SPEARMAN_RHO = "spearman-rho"
REAL_FPKM = "real-fpkm"
LOG2_REAL_FPKM = "log2-real-fpkm"
CALCULATED_FPKM = "calc-fpkm"
LOG2_CALCULATED_FPKM = "log2-calc-fpkm"
LOG2_RATIO = "log-ratio"
TRANSCRIPT_COUNT = "num-transcripts"

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

# Keep only transcripts for which both real and calculated FPKM are non-zero
fpkms = fpkms[fpkms[REAL_FPKM] > 0]
fpkms = fpkms[fpkms[CALCULATED_FPKM] > 0]

# Calculate log ratio of calculated and real FPKMs
fpkms[LOG2_REAL_FPKM] = np.log2(fpkms[REAL_FPKM])
fpkms[LOG2_CALCULATED_FPKM] = np.log2(fpkms[CALCULATED_FPKM])
fpkms[LOG2_RATIO] = \
    fpkms[LOG2_CALCULATED_FPKM] - fpkms[LOG2_REAL_FPKM]

# Group transcripts by the number of transcripts for their originating gene,
# and discard those groups with fewer than 100 members
grouped = fpkms.groupby(TRANSCRIPT_COUNT)
filtered = grouped.filter(lambda x: len(x[REAL_FPKM]) > 100)

# Write log ratio summary statistics stratified by the number of transcripts
# per originating gene
with open(options[OUT_FILE_BASENAME] + "_stats.csv", "w") as out_file:
    out_file.write(",".join([
        QUANT_METHOD, READ_LENGTH, READ_DEPTH, PAIRED_END,
        ERRORS, NUM_FPKMS, SPEARMAN_RHO
    ]) + "\n")

    # Spearman correlation coefficient between real and calculated FPKMs.
    rho = fpkms[LOG2_CALCULATED_FPKM].corr(
        fpkms[LOG2_REAL_FPKM], method='spearman')

    out_file.write(','.join([
        options[QUANT_METHOD_OPT], options[READ_LENGTH_OPT],
        options[READ_DEPTH_OPT], options[PAIRED_END_OPT],
        options[ERRORS_OPT], len(fpkms), rho]) + "\n")

with open(options[OUT_FILE_BASENAME] +
          "_stats_by_num_transcripts_per_gene.csv", "w") as out_file:

    grouped = filtered.groupby(TRANSCRIPT_COUNT)
    summary = grouped.describe()
    log_ratio_stats = summary[LOG2_RATIO].unstack()
    log_ratio_stats = log_ratio_stats.drop(["min", "25%", "75%", "max"], 1)
    log_ratio_stats = log_ratio_stats.rename(columns={"50%": "median"})

    log_ratio_stats[SPEARMAN_RHO] = grouped.apply(
        lambda x: x[LOG2_CALCULATED_FPKM].corr(x[LOG2_REAL_FPKM],
                                                   method="spearman"))

    log_ratio_stats[QUANT_METHOD] = options[QUANT_METHOD_OPT]
    log_ratio_stats[READ_LENGTH] = options[READ_LENGTH_OPT]
    log_ratio_stats[READ_DEPTH] = options[READ_DEPTH_OPT]
    log_ratio_stats[PAIRED_END] = options[PAIRED_END_OPT]
    log_ratio_stats[ERRORS] = options[ERRORS_OPT]
    log_ratio_stats[NUM_FPKMS] = len(fpkms)

    log_ratio_stats.to_csv(out_file)

# Make a scatter plot of calculated vs real FPKMs
plt.figure()
scatter = plt.scatter(
    fpkms[REAL_FPKM].values,
    fpkms[CALCULATED_FPKM].values)
plt.suptitle("Scatter plot of calculated vs real FPKMs")
plt.xlabel("Real FPKM")
plt.ylabel("Calculated FPKM")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
if options[SCATTER_MAX] != 0:
    plt.xlim(xmax=options[SCATTER_MAX])
    plt.ylim(ymax=options[SCATTER_MAX])

plt.savefig(options[OUT_FILE_BASENAME] + "_scatter.pdf", format="pdf")

# Make a scatter plot of log transformed calculated vs real FPKMs
plt.figure()
scatter = plt.scatter(
    fpkms[LOG2_REAL_FPKM].values,
    fpkms[LOG2_CALCULATED_FPKM].values)
plt.suptitle("Scatter plot of log transformed calculated vs real FPKMs")
plt.xlabel("Log2 real FPKM")
plt.ylabel("Log2 calculated FPKM")
if options[LOG2_SCATTER_MIN] != 0:
    plt.xlim(xmin=options[LOG2_SCATTER_MIN])
    plt.ylim(ymin=options[LOG2_SCATTER_MIN])
if options[LOG2_SCATTER_MAX] != 0:
    plt.xlim(xmax=options[LOG2_SCATTER_MAX])
    plt.ylim(ymax=options[LOG2_SCATTER_MAX])

plt.savefig(options[OUT_FILE_BASENAME] + "_log2_scatter.pdf", format="pdf")

# Make a boxplot of log ratios stratified by the number of transcripts per
# originating gene.
bp = filtered.boxplot(column=LOG2_RATIO, by=TRANSCRIPT_COUNT, sym="")
bp.set_title("")
plt.suptitle("Distribution of log ratios of calculated to real FPKMs")
plt.xlabel("No. transcripts per gene")
plt.ylabel("Log ratio (calculated/real FPKM)")

plt.savefig(options[OUT_FILE_BASENAME] + "_boxplot.pdf", format="pdf")
