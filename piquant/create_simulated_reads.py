#!/usr/bin/python

# TODO: add logging

"""Usage:
    create_simulated_reads [--log-level=<log-level>] [--out-dir=<out_dir>] [--num-fragments=<num-fragments>] [--prepare-only|--run-only] --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --polya=<polya> <transcript-gtf-file> <genome-fasta-dir>

-h --help                        Show this message.
-v --version                     Show version.
--log-level=<log-level>          Set logging level (one of {log_level_vals}) [default: info].
--out-dir=<out-dir>              Parent output directory to which read directories will be written [default: output].
--num-fragments=<num-fragments>  Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--prepare-only                   Read simulation scripts will be created, but not run (they must not already exist).
--run-only                       Read simulation scripts will be run, but not created (they must already exist).
-l --read-length=<read-lengths>  Comma-separated list of read-lengths to create reads for.
-d --read-depth=<read-depths>    Comma-separated list of read-depths to create reads for
-p --paired-end=<paired-ends>    Comma-separated list of True/False strings indicating whether paired-end reads should be created.
-e --error=<errors>              Comma-separated list of True/False strings indicating whether reads should be created with errors.
-b --bias=<biases>               Comma-separated list of True/False strings indicating whether reads should be created with sequence bias.
-a --polya=<polya>               Comma-separated list of True/False strings indicating whether reads should be created with polyA tails.
<transcript-gtf-file>            GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>               Directory containing per-chromosome sequences as FASTA files.
"""

import docopt
import ordutils.log as log
import ordutils.options as opt
import os.path
import parameters
import prepare_read_simulation as prs
import schema
import subprocess
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
OUTPUT_DIRECTORY = "--out-dir"
NUM_FRAGMENTS = "--num-fragments"
PREPARE_ONLY = "--prepare-only"
RUN_ONLY = "--run-only"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="create_simulated_reads v0.1")

# Validate command-line options
param_values = None

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    opt.validate_dir_option(
        options[OUTPUT_DIRECTORY], "Output parent directory does not exist.")
    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")
    opt.validate_dir_option(
        options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")
    options[NUM_FRAGMENTS] = opt.validate_int_option(
        options[NUM_FRAGMENTS],
        "Number of fragments must be a positive integer.",
        nonneg=True)

    param_values = parameters.validate_command_line_parameter_sets(options)
except schema.SchemaError as exc:
    exit(exc.code)


def check_reads_directory(**params):
    reads_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)
    if options[RUN_ONLY] != os.path.exists(reads_dir):
        sys.exit("Reads directory '{d}' should ".format(d=reads_dir) +
                 ("" if options[RUN_ONLY] else "not ") + "already exist.")


def create_and_run_reads_simulation(**params):
    reads_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)

    if not options[RUN_ONLY]:
        prs.create_simulation_files(
            reads_dir, options[TRANSCRIPT_GTF_FILE], options[GENOME_FASTA_DIR],
            options[NUM_FRAGMENTS], **params)

    if not options[PREPARE_ONLY]:
        cwd = os.getcwd()
        os.chdir(reads_dir)
        args = ['nohup', './run_simulation.sh']
        subprocess.Popen(args)
        os.chdir(cwd)


# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

parameters.execute_for_param_sets(
    [check_reads_directory, create_and_run_reads_simulation],
    **param_values)
