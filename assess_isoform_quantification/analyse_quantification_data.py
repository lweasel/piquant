#!/usr/bin/python

# TODO: validate_file_option does not return file

"""Usage:
    analyse_quantification_data [--scatter-max=<scatter-max-val>] [--log2-scatter-min=<log2-scatter-min-val>] [--log2-scatter-max=<log2-scatter-max-val>] <fpkm-file> <out-file>

-h --help                                  Show this message.
-v --version                               Show version.
--scatter-max=<scatter-max-val>            Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log2-scatter-min=<log2-scatter-min-val>  Minimum x and y values for log2 scatter plot; a value of 0 means do not impose a minimum [default: -10].
--log2-scatter-max=<log2-scatter-max-val>  Maximum x and y values for log2 scatter plot; a value of 0 means do not impose a maximum [default: 15].
<fpkm-file>                                File containing real and calculated FPKMs.
<out-file>                                 Basename for output graph and data files.
"""

from docopt import docopt
from schema import SchemaError

import matplotlib.pyplot as plt
import numpy as np
import options as opt
import pandas as pd

import sys

FPKM_FILE = "<fpkm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG2_SCATTER_MIN = "--log2-scatter-min"
LOG2_SCATTER_MAX = "--log2-scatter-max"

REAL_FPKM_COL = "Real FPKM"
LOG2_REAL_FPKM_COL = "Log2 real FPKM"
CALCULATED_FPKM_COL = "Calculated FPKM"
LOG2_CALCULATED_FPKM_COL = "Log2 calculated FPKM"
LOG2_RATIO_COL = "Log ratio"
TRANSCRIPT_COUNT_COL = "Num transcripts for gene"

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
fpkms = fpkms[fpkms[REAL_FPKM_COL] > 0]
fpkms = fpkms[fpkms[CALCULATED_FPKM_COL] > 0]

# Calculate log ratio of calculated and real FPKMs
fpkms[LOG2_REAL_FPKM_COL] = np.log2(fpkms[REAL_FPKM_COL])
fpkms[LOG2_CALCULATED_FPKM_COL] = np.log2(fpkms[CALCULATED_FPKM_COL])
fpkms[LOG2_RATIO_COL] = \
    fpkms[LOG2_CALCULATED_FPKM_COL] - fpkms[LOG2_REAL_FPKM_COL]

# Group transcripts by the number of transcripts for their originating gene,
# and discard those groups with fewer than 100 members
grouped = fpkms.groupby(TRANSCRIPT_COUNT_COL)
filtered = grouped.filter(lambda x: len(x[REAL_FPKM_COL]) > 100)

# Write log ratio summary statistics stratified by the number of transcripts
# per originating gene
with open(options[OUT_FILE_BASENAME] + "_stats.txt", "w") as out_file:
    num_fpkms = len(fpkms)
    out_file.write("{n} transcripts with real and calculated FPKMS > 0.\n".
                   format(n=num_fpkms))

    rho = fpkms[LOG2_CALCULATED_FPKM_COL].corr(
        fpkms[LOG2_REAL_FPKM_COL], method='spearman')
    out_file.write("Rho for log real and calculated FPKMs: {rho}\n".
                   format(rho=rho))

    out_file.write("\n" + ("=" * 30) + "\n\n")

    grouped = filtered.groupby(TRANSCRIPT_COUNT_COL)
    summary = grouped.describe()
    log_ratio_stats = summary["Log ratio"].unstack()
    log_ratio_stats = log_ratio_stats.drop(["min", "25%", "75%", "max"], 1)
    log_ratio_stats = log_ratio_stats.rename(columns={"50%": "median"})

    log_ratio_stats["rho"] = grouped.apply(
        lambda x: x[LOG2_CALCULATED_FPKM_COL].corr(x[LOG2_REAL_FPKM_COL],
                                                   method="spearman"))

    out_file.write("Stats stratified by number of transcripts " +
                   "per originating gene.\n")
    out_file.write(log_ratio_stats.to_string())

# Make a scatter plot of calculated vs real FPKMs
plt.figure()
scatter = plt.scatter(
    fpkms[REAL_FPKM_COL].values,
    fpkms[CALCULATED_FPKM_COL].values)
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
    fpkms[LOG2_REAL_FPKM_COL].values,
    fpkms[LOG2_CALCULATED_FPKM_COL].values)
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
bp = filtered.boxplot(column=LOG2_RATIO_COL, by=TRANSCRIPT_COUNT_COL, sym="")
bp.set_title("")
plt.suptitle("Distribution of log ratios of calculated to real FPKMs")
plt.xlabel("No. transcripts per gene")
plt.ylabel("Log ratio (calculated/real FPKM)")

plt.savefig(options[OUT_FILE_BASENAME] + "_boxplot.pdf", format="pdf")
