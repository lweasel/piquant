#!/usr/bin/python

"""Usage:
    prepare_quantification_run [--log-level=<log-level>] [--run-directory <run-directory>] <transcript-gtf-file> <genome-dir>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
-d --run-directory=<run-directory>  Directory to create to which run files will be written [default: out].
<transcript-gtf-file>               GTF formatted file describing transcripts to be simulated.
<genome-dir>                        Directory containing per-chromosome sequences as FASTA files.
"""

from docopt import docopt
from schema import SchemaError

import log
import options as opt
import os
import os.path
import stat
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
RUN_DIRECTORY = "--run-directory"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_DIRECTORY = "<genome-dir>"

FLUX_SIMULATOR_PARAMS_FILE = "flux_simulator.par"
FLUX_SIMULATOR_PRO_FILE = "flux_simulator.pro"
RUN_SCRIPT = "run_quantification.sh"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="prepare_quantification_run v0.1")

# Validate and process command-line options
try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    opt.validate_dir_option(
        options[RUN_DIRECTORY],
        "Run directory should not already exist",
        should_exist=False)
    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcripts GTF file does not exist")
    opt.validate_dir_option(
        options[GENOME_DIRECTORY], "Genome sequence directory does not exist")
except SchemaError as exc:
    exit("Exiting. " + exc.code)

options[TRANSCRIPT_GTF_FILE] = os.path.abspath(options[TRANSCRIPT_GTF_FILE])
options[GENOME_DIRECTORY] = os.path.abspath(options[GENOME_DIRECTORY])

# Create directory for run files

logger = log.getLogger(sys.stderr, options[LOG_LEVEL])

logger.info("Creating run directory {dir}.".
            format(dir=options[RUN_DIRECTORY]))

os.mkdir(options[RUN_DIRECTORY])

# Write Flux Simulator parameters file

logger.info("Creating Flux Simulator parameters file.")


def get_output_file(filename):
    return open(options[RUN_DIRECTORY] + os.path.sep + filename, "w")


def write_lines(f, lines):
    f.write("\n".join(lines))
    f.write("\n")

with get_output_file(FLUX_SIMULATOR_PARAMS_FILE) as params:
    lines = [
        "REF_FILE_NAME {f}".format(f=options[TRANSCRIPT_GTF_FILE]),
        "GEN_DIR {d}".format(d=options[GENOME_DIRECTORY]),
        "PCR_DISTRIBUTION none",
        "LIB_FILE_NAME flux_simulator.lib",
        "SEQ_FILE_NAME reads.bed",
        "FASTA YES"
    ]
    write_lines(params, lines)

# Write shell script to run quantification
script_path = None

with get_output_file(RUN_SCRIPT) as script:
    lines = [
        "#!/bin/bash",
        "",
        "# Run Flux Simulator to create expression profiles then simulate reads",
        "flux-simulator -t simulator -x -p {f}".format(f=FLUX_SIMULATOR_PARAMS_FILE),
        "",
        "# (this is a hack - Flux Simulator seems to sometimes incorrectly",
        "# output transcripts with zero length)",
        "ZERO_LENGTH_COUNT=$(awk 'BEGIN {i=0} $4 == 0 {i++;} END{print i}'" +
        " {f})".format(f=FLUX_SIMULATOR_PRO_FILE),
        "echo",
        "echo Removing $ZERO_LENGTH_COUNT transcripts with zero length...",
        "echo",
        "awk '$4 > 0' {f} > tmp; mv tmp {f}".format(f=FLUX_SIMULATOR_PRO_FILE),
        "",
        "flux-simulator -t simulator -l -s -p {f}".format(f=FLUX_SIMULATOR_PARAMS_FILE)
    ]
    write_lines(script, lines)

    script_path = os.path.abspath(script.name)

os.chmod(script_path,
         stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
         stat.S_IRGRP | stat.S_IROTH)
