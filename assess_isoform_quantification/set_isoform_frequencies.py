#!/usr/bin/python

"""Usage:
    set_isoform_frequencies [{help}] [{version}] {pro_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{pro_file}  Flux simulator gene expression profile file.
"""

from docopt import docopt
from options import validate_file_option
from pandas import read_csv
from schema import SchemaError

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"

LOCUS_COL = 'loc'
TRANSCRIPT_ID_COL = 't_id'
CODING_COL = 'c'
LENGTH_COL = 'len'
FRACTION_COL = 'f'
NUM_TRANSCRIPTS_COL = 'n'
UNKNOWN_COL_1 = 'u1'
UNKNOWN_COL_2 = 'u2'

PRO_FILE_COLUMNS = [
    LOCUS_COL,
    TRANSCRIPT_ID_COL,
    CODING_COL,
    LENGTH_COL,
    FRACTION_COL,
    NUM_TRANSCRIPTS_COL,
    UNKNOWN_COL_1,
    UNKNOWN_COL_2]

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE)

# Read in command-line options
options = docopt(__doc__, version="set_isoform_frequencies v0.1")

# Validate command-line options

try:
    options[PRO_FILE] = validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
except SchemaError as exc:
    exit(exc.code)

df = read_csv(options[PRO_FILE], sep='\s*',
              names=PRO_FILE_COLUMNS, index_col=TRANSCRIPT_ID_COL)
