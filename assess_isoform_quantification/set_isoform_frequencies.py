#!/usr/bin/python

"""Usage:
    set_isoform_frequencies [{help}] [{version}] [{log_level}={log_level_val}] {pro_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{log_level}={log_level_val}   Set logging level (one of {log_level_vals}) [default: 'info'].
{pro_file}                Flux simulator gene expression profile file.
"""

from docopt import docopt
from log import getLogger, LEVELS
from options import validate_file_option, validate_dict_option
from pandas import read_csv
from schema import SchemaError
from sys import stderr

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"
LOG_LEVEL = "--log-level"
LOG_LEVEL_VAL = "<log-level>"
LOG_LEVEL_VALS = str(LEVELS.keys())

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE,
    log_level=LOG_LEVEL,
    log_level_val=LOG_LEVEL_VAL,
    log_level_vals=LOG_LEVEL_VALS)

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

# Read in command-line options
options = docopt(__doc__, version="set_isoform_frequencies v0.1")

# Validate command-line options
try:
    options[PRO_FILE] = validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    validate_dict_option(
        options[LOG_LEVEL], LEVELS, "Invalid log level")
except SchemaError as exc:
    exit(exc.code)

logger = getLogger(stderr, options[LOG_LEVEL])
logger.info("Reading transcript abundances from expression profile file.")

transcript_abundances = read_csv(options[PRO_FILE], sep='\s*',
                                 names=PRO_FILE_COLUMNS,
                                 index_col=TRANSCRIPT_ID_COL)

num_molecules = transcript_abundances[NUM_TRANSCRIPTS_COL].sum()
