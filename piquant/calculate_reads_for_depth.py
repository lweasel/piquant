#!/usr/bin/env python

"""Usage:
    calculate_reads_for_depth [--log-level=<log-level>] <pro-file> <read-length> <read-depth>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
<pro-file>                          Flux Simulator gene expression profile file.
<read-length>                       The length of simulated reads.
<read-depth>                        The (approximate) read depth that is required given the specified read length.
"""

import docopt
import flux_simulator as fs
import log
import options as opt
import schema
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())

PRO_FILE = "<pro-file>"
READ_LENGTH = "<read-length>"
READ_DEPTH = "<read-depth>"


def _validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[PRO_FILE], "Could not open expression profile file")
        options[READ_LENGTH] = opt.validate_int_option(
            options[READ_LENGTH], "Read length must be a positive integer",
            nonneg=True)
        options[READ_DEPTH] = opt.validate_int_option(
            options[READ_DEPTH], "Read depth must be a positive integer",
            nonneg=True)
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
    return bases_to_sequence // read_length


if __name__ == "__main__":
    # Read in command-line options
    __doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(__doc__, version="calculate_reads_for_depth v0.1")

    # Validate and process command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    # Read in Flux Simulator expression profiles
    profiles = _read_expression_profiles(logger, options[PRO_FILE])

    # Calculate the number of reads required to approximately give the
    # specified overall average depth of coverage
    print(_calculate_reads_for_depth(
        profiles, options[READ_LENGTH], options[READ_DEPTH]))
