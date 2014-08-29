#!/usr/bin/env python

"""Usage:
    assemble_quantification_data [{log_option_spec} --method=<quantification-method> --out=<output-file> <pro-file> <quantification-file> <transcript-count-file> <unique-sequence-file>

-h --help                    Show this message.
-v --version                 Show version.
{log_option_spec}      {log_option_description}
-m --method=<quant-method>   Method used to quantify transcript abundances.
-o <output-file> --out=<output-file>
                             Output file for real and calculated TPMs.
<pro-file>                   Flux Simulator gene expression profile file.
<quantification-file>        Output file from transcript quantifier from which abundances will be read.
<transcript-count-file>      File containing per-gene transcript counts.
<unique-sequence-file>       File containing unique sequence lengths per-transcript.
"""

from docopt import docopt
from schema import SchemaError

import flux_simulator as fs
import options as opt
import pandas
import quantifiers as qs
import tpms

QUANT_METHOD = "--method"
OUT_FILE = "--out"
PRO_FILE = "<pro-file>"
QUANT_FILE = "<quantification-file>"
COUNT_FILE = "<transcript-count-file>"
UNIQUE_SEQ_FILE = "<unique-sequence-file>"

TRANSCRIPT_COL = "transcript"
COUNT_COL = "transcript_count"
UNIQUE_SEQ_LENGTH_COL = "unique-length"

SORTED_PREFIX = "sorted"
SORTED_BAM_FILE = SORTED_PREFIX + ".bam"

# Read in command-line options
__doc__ = opt.substitute_into_usage(__doc__)
options = docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    opt.validate_log_level(options)

    opt.validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    opt.validate_file_option(
        options[QUANT_FILE], "Could not open transcript abundance file")
    opt.validate_file_option(
        options[COUNT_FILE], "Could not open transcript count file")
    opt.validate_file_option(
        options[UNIQUE_SEQ_FILE],
        "Could not open unique sequence lengths file")
    options[QUANT_METHOD] = opt.validate_dict_option(
        options[QUANT_METHOD], qs.get_quantification_methods(),
        "Unknown quantification method")
except SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = opt.get_logger_for_options(options)

# Read in the expression profile file, and calculate the true TPM
# for each transcript
logger.info("Reading expression profiles...")

profiles = fs.read_expression_profiles(options[PRO_FILE])
profiles[tpms.REAL_TPM] = \
    profiles[fs.PRO_FILE_FRAC_COL].map(lambda x: 1000000 * x)

# Read calculated TPM values for each transcript produced by a particular
# quanfication method
logger.info("Reading calculated TPMs...")

options[QUANT_METHOD].calculate_transcript_abundances(options[QUANT_FILE])

profiles[tpms.CALCULATED_TPM] = profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].\
    map(options[QUANT_METHOD].get_transcript_abundance)

# Read per-gene transcript counts
logger.info("Reading per-gene transcript counts...")

transcript_counts = pandas.read_csv(
    options[COUNT_FILE], index_col=TRANSCRIPT_COL)


set_transcript_count = lambda t_id: \
    transcript_counts.ix[t_id][COUNT_COL] \
    if t_id in transcript_counts.index else 0

profiles[tpms.TRANSCRIPT_COUNT] = \
    profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].map(set_transcript_count)

# Read unique sequence lengths per-transcript
logger.info("Reading unique sequence lengths per-transcript")

unique_seqs = pandas.read_csv(
    options[UNIQUE_SEQ_FILE], index_col=TRANSCRIPT_COL)

set_unique_length = lambda t_id: \
    unique_seqs.ix[t_id][UNIQUE_SEQ_LENGTH_COL] \
    if t_id in unique_seqs.index else 0

profiles[tpms.UNIQUE_SEQ_LENGTH] = \
    profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].map(set_unique_length)

# Write TPMs and other relevant data to output file
logger.info("Writing TPMs to file {out}".format(out=options[OUT_FILE]))

profiles.rename(
    columns={
        fs.PRO_FILE_TRANSCRIPT_ID_COL: TRANSCRIPT_COL,
        fs.PRO_FILE_LENGTH_COL: tpms.LENGTH
    },
    inplace=True)

profiles.to_csv(options[OUT_FILE], index=False,
                cols=[TRANSCRIPT_COL, tpms.LENGTH, tpms.UNIQUE_SEQ_LENGTH,
                      tpms.TRANSCRIPT_COUNT, tpms.REAL_TPM,
                      tpms.CALCULATED_TPM])
