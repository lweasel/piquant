#!/usr/bin/python

"""Usage:
    assess_isoform_quantification {quant_method}={quant_method_val} {pro_file} {read_file} {quant_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{quant_method_short} {quant_method_val} {quant_method}={quant_method_val}       Method used to quantify transcript abundances.
{pro_file}          Flux Simulator gene expression profile file.
{read_file}         BAM file containing mapped reads used as input to transcript quantifier
{quant_file}        Output file from transcript quantifier from which abundances will be read.
"""

from docopt import docopt
from schema import SchemaError

import flux_simulator as fs
import options as opt
import pysam
import quantifiers as qs

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"
READ_FILE = "<read-file>"
QUANT_METHOD_SHORT = "-m"
QUANT_METHOD = "--method"
QUANT_METHOD_VAL = "<quantification-method>"
QUANT_FILE = "<quantification-file>"

REAL_FPKM_COL = 'Real FPKM'
CALCULATED_FPKM_COL = 'Calculated FPKM'

CUFFLINKS_METHOD = "cufflinks"

QUANT_METHODS = {
    CUFFLINKS_METHOD: qs.Cufflinks
}

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE,
    read_file=READ_FILE,
    quant_method_short=QUANT_METHOD_SHORT,
    quant_method=QUANT_METHOD,
    quant_method_val=QUANT_METHOD_VAL,
    quant_file=QUANT_FILE)

# Read in command-line options
options = docopt(__doc__, version="assess_isoform_quantification v0.1")

# Validate command-line options
try:
    options[PRO_FILE] = opt.validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    options[QUANT_FILE] = opt.validate_file_option(
        options[QUANT_FILE], "Could not open transcript abundance file")
    opt.validate_file_option(
        options[READ_FILE], "Could not open BAM file containing reads")
    options[QUANT_METHOD] = opt.validate_dict_option(
        options[QUANT_METHOD], QUANT_METHODS, "Unknown quantification method")
except SchemaError as exc:
    exit(exc.code)

# Read in mapped reads and record transcript length and number of
# mapped reads per transcript, and also the total number of millions
# of mapped reads

transcript_info = {}

pysam.index(options[READ_FILE])
read_bam = pysam.Samfile(options[READ_FILE], "rb")
for read in read_bam.fetch():
    transcript_id = fs.get_transcript_id(read.qname)
    if transcript_id in transcript_info:
        info = transcript_info[transcript_id]
        info[1] += 1
    else:
        info = [fs.get_transcript_length(read.qname), 1]
        transcript_info[transcript_id] = info

mapped_reads = sum([info[1] for info in transcript_info.values()])
m_r_millons = float(mapped_reads)/1000000

# Read in the expression profile file, and calculate the true FPKM
# for each transcript


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
quant_method = options[QUANT_METHOD]()
quant_method.calculate_transcript_abundances(options[QUANT_FILE])

set_calculated_fpkm = lambda row: \
    quant_method.get_transcript_abundance(row[fs.PRO_FILE_TRANSCRIPT_ID_COL])
profiles[CALCULATED_FPKM_COL] = profiles.apply(set_calculated_fpkm, axis=1)

# Temporarily just print info for those transcripts with non-zero real FPKM
non_zero = profiles[profiles[REAL_FPKM_COL] != 0]
print(non_zero)
