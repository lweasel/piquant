#!/usr/bin/python

"""Usage:
    clean_mapped_reads {in_bam_file} {out_bam_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{in_bam_file}             Input BAM file containing mapped reads to be cleaned.
{out_bam_file}            Output BAM file to write cleaned reads to.
"""

from docopt import docopt
from flux_simulator import get_fragment_bounds
from options import validate_file_option
from pysam import Samfile
from schema import SchemaError

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
IN_BAM_FILE = "<in-bam-file>"
OUT_BAM_FILE = "<out-bam-file>"

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    in_bam_file=IN_BAM_FILE,
    out_bam_file=OUT_BAM_FILE)

# Read in command-line options
options = docopt(__doc__, version="clean_mapped_reads v0.1")

# Validate command-line options
try:
    validate_file_option(options[IN_BAM_FILE], "Could not open input BAM file")
except SchemaError as exc:
    exit(exc.code)

input_bam = Samfile(options[IN_BAM_FILE], "rb")
output_bam = Samfile(options[OUT_BAM_FILE], "wb", template=input_bam)

for read in input_bam.fetch():
    # Find region and start, end positions within which the originating
    # fragment should lie
    f_region, f_start, f_end = get_fragment_bounds(read.qname)

    # Find region and start, end positions of mapped read
    r_region = input_bam.getrname(read.tid)
    r_start = read.pos
    r_end = read.aend

    # Check read region and position are consistent with originating
    # transcript - otherwise the incorrectly mapped read is suppressed.
    # NB. this does not eliminate reads incorrectly mapped within the
    # bounds of the originating transcript.
    if f_region == r_region and r_start >= f_start and r_end <= f_end:
        output_bam.write(read)
