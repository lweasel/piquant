#!/usr/bin/python

"""Usage:
    assess_isoform_quantification {pro_file} {read_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{pro_file}          Flux Simulator gene expression profile file.
{read_file}         BAM file containing mapped reads used as input to
                    transcript quantifier
"""

from docopt import docopt
from options import validate_file_option
from schema import SchemaError

import flux_simulator as fs
import pysam

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"
READ_FILE = "<read-file>"

FPKM_COL = 'FPKM'

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE,
    read_file=READ_FILE)

# Read in command-line options
options = docopt(__doc__, version="assess_isoform_quantification v0.1")

# Validate command-line options
try:
    options[PRO_FILE] = validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    validate_file_option(
        options[READ_FILE], "Could not open BAM file containing reads")
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
profiles[FPKM_COL] = profiles.apply(fpkm_calc, axis=1)

# Temporarily just print info for those transcripts with non-zero FPKM
print(profiles[profiles[FPKM_COL] != 0])
