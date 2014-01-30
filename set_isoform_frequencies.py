#!/usr/bin/python

"""Usage:
    set_isoform_frequencies [{help}] [{version}] {pro_file}

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
{pro_file}  Flux simulator gene expression profile file.
"""

from docopt import docopt
from schema import Schema, SchemaError

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE)

# Read in command-line options
options = docopt(__doc__, version="set_isoform_frequencies v0.1")


# Validate command-line options

def validate_file_option(file_option, msg):
    msg = "{msg} '{file}'.".format(msg=msg, file=file_option)
    return Schema(open, error=msg).validate(file_option)

try:
    options[PRO_FILE] = validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
except SchemaError as exc:
    exit(exc.code)

print("Now do all the things.")
