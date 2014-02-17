#!/usr/bin/python

# TODO: Think I should be using DataFrame.map, not DataFrame.apply.

"""Usage:
    assess_isoform_quantification [--log-level=<log-level>] --method=<quantification-method> --out=<output-file> <pro-file> <read-file> <quantification-file> <transcript-count-file>

-h --help                 Show this message.
-v --version              Show version.
--log-level=<log-level>   Set logging level (one of {log_level_vals}) [default: info].
-m <quantification-method> --method=<quantification-method>
                          Method used to quantify transcript abundances.
-o <output-file> --out=<output-file>
                          Output file for real and calculated FPKMs.
<pro-file>                Flux Simulator gene expression profile file.
<read-file>               BAM file containing mapped reads used as input to transcript quantifier
<quantification-file>     Output file from transcript quantifier from which abundances will be read.
<transcript-count-file>   File containing per-gene transcript counts.
"""

from docopt import docopt
from schema import SchemaError

import flux_simulator as fs
import log
import options as opt
import pandas
import pysam
import quantifiers as qs
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
QUANT_METHOD = "--method"
OUT_FILE = "--out"
PRO_FILE = "<pro-file>"
READ_FILE = "<read-file>"
QUANT_FILE = "<quantification-file>"
COUNT_FILE = "<transcript-count-file>"

COUNTS_TRANSCRIPT_COL = "transcript"
COUNTS_COUNT_COL = "transcript_count"

TRANSCRIPT_COL = 'Transcript'
REAL_FPKM_COL = 'Real FPKM'
CALCULATED_FPKM_COL = 'Calculated FPKM'
TRANSCRIPT_COUNT_COL = 'Num transcripts for gene'

CUFFLINKS_METHOD = "cufflinks"

QUANT_METHODS = {
    CUFFLINKS_METHOD: qs.Cufflinks
}

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="assemble_quantification_data v0.1")

# Validate command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    options[PRO_FILE] = opt.validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    options[QUANT_FILE] = opt.validate_file_option(
        options[QUANT_FILE], "Could not open transcript abundance file")
    opt.validate_file_option(
        options[READ_FILE], "Could not open BAM file containing reads")
    opt.validate_file_option(
        options[COUNT_FILE], "Could not open transcript count file")
    options[QUANT_METHOD] = opt.validate_dict_option(
        options[QUANT_METHOD], QUANT_METHODS, "Unknown quantification method")
except SchemaError as exc:
    exit(exc.code)

logger = log.getLogger(sys.stderr, options[LOG_LEVEL])

# Read in mapped reads and record transcript length and number of
# mapped reads per transcript, and also the total number of millions
# of mapped reads

logger.info("Reading mapped reads...")

transcript_info = {}

pysam.index(options[READ_FILE])
read_bam = pysam.Samfile(options[READ_FILE], "rb")
count = 0

for read in read_bam.fetch():
    transcript_id = fs.get_transcript_id(read.qname)
    if transcript_id in transcript_info:
        info = transcript_info[transcript_id]
        info[1] += 1
    else:
        info = [fs.get_transcript_length(read.qname), 1]
        transcript_info[transcript_id] = info
    count += 1
    if count % 100000 == 0:
        logger.debug("...read {num}".format(num=count))

mapped_reads = sum([info[1] for info in transcript_info.values()])
m_r_millons = float(mapped_reads)/1000000

logger.info("{num} mapped reads read.".format(num=mapped_reads))

# Read in the expression profile file, and calculate the true FPKM
# for each transcript

logger.info("Reading expression profiles...")


def fpkm_calc(profile_row):
    t_id = profile_row[fs.PRO_FILE_TRANSCRIPT_ID_COL]
    if t_id not in transcript_info:
        return 0
    t_info = transcript_info[t_id]
    t_kb = float(t_info[0])/1000
    return t_info[1] / (t_kb * m_r_millons)

profiles = fs.read_expression_profiles(options[PRO_FILE])
profiles[REAL_FPKM_COL] = profiles.apply(fpkm_calc, axis=1)

# Read calculated FPKM values for each transcript produced by a particular
# quanfication method

logger.info("Reading calculated FPKMs...")

quant_method = options[QUANT_METHOD]()
quant_method.calculate_transcript_abundances(options[QUANT_FILE])

set_calculated_fpkm = lambda row: \
    quant_method.get_transcript_abundance(row[fs.PRO_FILE_TRANSCRIPT_ID_COL])
profiles[CALCULATED_FPKM_COL] = profiles.apply(set_calculated_fpkm, axis=1)

# Read per-gene transcript counts

logger.info("Reading per-gene transcript counts...")

transcript_counts = pandas.read_csv(
    options[COUNT_FILE], index_col=COUNTS_TRANSCRIPT_COL)


def set_transcript_count(profile_row):
    t_id = profile_row[fs.PRO_FILE_TRANSCRIPT_ID_COL]
    if t_id not in transcript_counts.index:
        return 0
    return transcript_counts.ix[t_id][COUNTS_COUNT_COL]

profiles[TRANSCRIPT_COUNT_COL] = profiles.apply(set_transcript_count, axis=1)

# Write FPKMs and other relevant data to output file

logger.info("Writing FPKMs to file {out}".format(out=options[OUT_FILE]))

profiles.rename(columns={1: TRANSCRIPT_COL}, inplace=True)

profiles.to_csv(options[OUT_FILE], index=False,
                cols=[TRANSCRIPT_COL, TRANSCRIPT_COUNT_COL,
                      REAL_FPKM_COL, CALCULATED_FPKM_COL])
