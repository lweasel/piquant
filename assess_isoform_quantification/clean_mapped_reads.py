#!/usr/bin/python

# TODO: validate_file_option does not return file

"""Usage:
    clean_mapped_reads [{help}] [{version}] [{write_rejected}] {in_bam_file}
                       {out_bam_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{write_rejected}            If specfied, write rejected reads to the file
                            "rejected.bam".
{in_bam_file}             Input BAM file containing mapped reads to be cleaned.
{out_bam_file}            Output BAM file to write cleaned reads to.
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
WRITE_REJECTED = "--write-rejected"
IN_BAM_FILE = "<in-bam-file>"
OUT_BAM_FILE = "<out-bam-file>"

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    write_rejected=WRITE_REJECTED,
    in_bam_file=IN_BAM_FILE,
    out_bam_file=OUT_BAM_FILE)

# Read in command-line options
options = docopt(__doc__, version="clean_mapped_reads v0.1")

# Validate command-line options
try:
    validate_file_option(options[IN_BAM_FILE], "Could not open input BAM file")
except SchemaError as exc:
    exit(exc.code)

pysam.index(options[IN_BAM_FILE])

input_bam = pysam.Samfile(options[IN_BAM_FILE], "rb")
output_bam = pysam.Samfile(options[OUT_BAM_FILE], "wb", template=input_bam)
rejected_bam = pysam.Samfile("rejected.bam", "wb", template=input_bam) \
    if options[WRITE_REJECTED] else None

read_ids = {}

for read in input_bam.fetch():
    bam = rejected_bam

    if read.qname not in read_ids:
        # Now check read region and position are consistent with originating
        # transcript - otherwise the incorrectly mapped read is suppressed.
        # We also suppress correctly mapped reads which lie wholly outside
        # the bounds of the original transcript (e.g. ones lying entirely
        # within the polyA tail.
        # NB. this procedure does not eliminate reads incorrectly mapped
        # within the bounds of the originating transcript.

        # Find region and start, end positions of mapped read
        r_region = input_bam.getrname(read.tid)
        r_start = read.pos
        r_end = read.aend

        # Find region and start, end positions within which the originating
        # fragment should lie
        f_region, f_start, f_end = fs.get_fragment_bounds(read.qname)

        if f_region == r_region and r_start >= f_start and r_end <= f_end:
            t_region, t_start, t_end = fs.get_transcript_bounds(read.qname)
            if r_end >= t_start and r_start <= t_end:
                # We retain only one alignment per read, and set all retained
                # alignments to be primary
                read_ids[read.qname] = True
                read.is_secondary = False
                bam = output_bam

    if bam is not None:
        bam.write(read)
