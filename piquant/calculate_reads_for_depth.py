# pylint: disable=E1103

"""
Usage:
    calculate_reads_for_depth [{log_option_spec}]
        <pro-file> <read-length> <read-depth>

Options:
{help_option_spec}
   {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<pro-file>
    Flux Simulator gene expression profile file.
<read-length>
    The length of simulated reads.
<read-depth>
    The (approximate) read depth that is required given the specified read
    length.

Calculate the approximate number of reads required to be simulated for a set of
transcripts in order to provide the specified sequencing depth, given a particular
length of read.
"""

from __future__ import print_function

import docopt
import schema

from . import flux_simulator as fs
from . import options as opt
from .__init__ import __version__

PRO_FILE = "<pro-file>"
READ_LENGTH = "<read-length>"
READ_DEPTH = "<read-depth>"


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_file_option(
            options[PRO_FILE], "Could not open expression profile file")
        options[READ_LENGTH] = opt.validate_int_option(
            options[READ_LENGTH], "Read length must be a positive integer",
            min_val=1)
        options[READ_DEPTH] = opt.validate_float_option(
            options[READ_DEPTH], "Read depth must be a positive float",
            min_val=0)
    except schema.SchemaError as exc:
        exit("Exiting. " + exc.code)


def _read_expression_profiles(logger, pro_file):
    logger.info("Reading expression profiles from '{f}'...".format(f=pro_file))
    profiles = fs.read_expression_profiles(pro_file)
    profiles.set_index(fs.PRO_FILE_TRANSCRIPT_ID_COL, inplace=True)
    logger.info("...read profiles for {n} transcripts.".
                format(n=len(profiles)))

    # Filter those transcripts with zero expression in the profile file
    profiles = profiles[profiles[fs.PRO_FILE_NUM_COL] > 0]
    logger.info("Retained {n} transcripts with non-zero expression.".
                format(n=len(profiles)))

    return profiles


def _calculate_reads_for_depth(profiles, read_length, required_depth):
    total_transcript_length = profiles[fs.PRO_FILE_LENGTH_COL].sum()
    bases_to_sequence = total_transcript_length * required_depth
    return int(bases_to_sequence // read_length)


def calculate_reads_for_depth(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        docstring, argv=args,
        version="calculate_reads_for_depth v" + __version__)

    # Validate and process command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Read in Flux Simulator expression profiles
    profiles = _read_expression_profiles(logger, options[PRO_FILE])

    # Calculate the number of reads required to approximately give the
    # specified overall average depth of coverage
    print(_calculate_reads_for_depth(
        profiles, options[READ_LENGTH], options[READ_DEPTH]))
