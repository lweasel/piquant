#!/usr/bin/python

"""Usage:
    set_isoform_frequencies [{help}] [{version}]

{help_short} {help}                 Show this message.
{version_short} {version}              Show version.
"""

from docopt import docopt
from schema import Schema, SchemaError

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION)

# Read in command-line options
options = docopt(__doc__, version="set_isoform_frequencies v0.1")

# Validate command-line options
try:
    pass
except SchemaError as exc:
    exit(exc.code)

print("Now do all the things.")
