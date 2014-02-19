#!/usr/bin/python

"""Usage:
    analyse_quantification_data [--scatter-max=<scatter-max-val>] <fpkm-file> <out-file>

-h --help                      Show this message.
-v --version                   Show version.
--scatter-max=scatter-max-val  Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
<fpkm-file>                    File containing real and calculated FPKMs.
<out-file>                     Basename for output graph and data files.
"""

from docopt import docopt
from schema import SchemaError

import matplotlib.pyplot as plt
import numpy as np
import options as opt
import pandas as pd

FPKM_FILE = "<fpkm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"

REAL_FPKM_COL = "Real FPKM"
CALCULATED_FPKM_COl = "Calculated FPKM"
TRANSCRIPT_COUNT_COL = "Num transcripts for gene"
LOG_RATIO_COL = "Log ratio"

# Read in command-line options
options = docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    options[FPKM_FILE] = opt.validate_file_option(
        options[FPKM_FILE], "Could not open FPKM file")
    options[SCATTER_MAX] = opt.validate_int_option(
        options[SCATTER_MAX],
        "Invalid maximum value for scatter plot axes",
        nullable=True)
except SchemaError as exc:
    exit(exc.code)

fpkms = pd.read_csv(options[FPKM_FILE])

# Keep only transcripts for which both real and calculated FPKM are non-zero
fpkms = fpkms[fpkms[REAL_FPKM_COL] > 0]
fpkms = fpkms[fpkms[CALCULATED_FPKM_COl] > 0]

# Calculate log ratio of calculated and real FPKMs
fpkms[LOG_RATIO_COL] = np.log2(fpkms[CALCULATED_FPKM_COl]/fpkms[REAL_FPKM_COL])

# Group transcripts by the number of transcripts for their originating gene,
# and discard those groups with fewer than 100 members
grouped = fpkms.groupby(TRANSCRIPT_COUNT_COL)
filtered = grouped.filter(lambda x: len(x[REAL_FPKM_COL]) > 100)

# Write log ratio summary statistics stratified by the number of transcripts
# per originating gene
with open(options[OUT_FILE_BASENAME] + "_stats.txt", "w") as out_file:
    summary = filtered.groupby(TRANSCRIPT_COUNT_COL).describe()
    summary = summary.drop(TRANSCRIPT_COUNT_COL, 1)
    out_file.write(str(summary))

# Make a scatter plot of calculated vs real FPKMs
scatter = plt.scatter(
    filtered[CALCULATED_FPKM_COl].values,
    filtered[REAL_FPKM_COL].values)
plt.suptitle("Scatter plot of log transformed calculated vs real FPKMs")
plt.xlabel("Real FPKM")
plt.ylabel("Calculated FPKM")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
if options[SCATTER_MAX] > 0:
    plt.xlim(xmax=options[SCATTER_MAX])
    plt.ylim(ymax=options[SCATTER_MAX])

plt.savefig(options[OUT_FILE_BASENAME] + "_scatter.pdf", format="pdf")

# Make a boxplot of log ratios stratified by the number of transcripts per
# originating gene.
bp = filtered.boxplot(column=LOG_RATIO_COL, by=TRANSCRIPT_COUNT_COL, sym="")
bp.set_title("")
plt.suptitle("Distribution of log ratios of calculated to real FPKMs")
plt.xlabel("No. transcripts per gene")
plt.ylabel("Log ratio (calculated/real FPKM)")

plt.savefig(options[OUT_FILE_BASENAME] + "_boxplot.pdf", format="pdf")
