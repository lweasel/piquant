#!/usr/bin/python

# TODO: add logging

"""Usage:
    create_simulated_reads [--log-level=<log-level>] [--out-dir=<out_dir>] [--num-fragments=<num-fragments>] [--prepare-only|--run-only] --read-lengths=<read-lengths> --read-depths=<read-depths> --paired-ends=<paired-ends> --errors=<errors> --biases=<biases> <transcript-gtf-file> <genome-fasta-dir>

-h --help                Show this message.
-v --version             Show version.
--log-level=<log-level>  Set logging level (one of {log_level_vals}) [default: info].
--out-dir=<out-dir>      Parent output directory to which read directories will be written [default: output].
--num-fragments=<num-fragments>  Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--prepare-only           Read simulation scripts will be created, but not run (they must not already exist).
--run-only               Read simulation scripts will be run, but not created (they must already exist).
-l --read-lengths=<read-lengths>  Comma-separated list of read-lengths to create reads for.
-d --read-depths=<read-depths>    Comma-separated list of read-depths to create reads for
-p --paired-ends=<paired-ends>    Comma-separated list of True/False strings indicating whether paired-end reads should be created.
-e --errors=<errors>              Comma-separated list of True/False strings indicating whether reads should be created with errors.
-b --biases=<biases>              Comma-separated list of True/False strings indicating whether reads should be created with sequence bias.
<transcript-gtf-file>    GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>       Directory containing per-chromosome sequences as FASTA files.
"""

import docopt
import itertools
import ordutils.log as log
import ordutils.options as opt
import os.path
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
READ_LENGTHS = "--read-lengths"
READ_DEPTHS = "--read-depths"
PAIRED_ENDS = "--paired-ends"
ERRORS = "--errors"
BIASES = "--biases"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="create_simulated_reads v0.1")

# Validate command-line options


def validate_list_option(option, item_validator, separator=','):
    items = option.split(separator)
    return [schema.Schema(
        schema.Use(item_validator)).validate(i) for i in items]


def check_boolean_value(data):
    data = data.lower()
    if data in ["true", "t", "yes", "y"]:
        return True
    elif data in ["false", "f", "no", "n"]:
        return False
    else:
        raise Exception("Can't convert '{d}' to bool.".format(d=data))


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
    options[READ_LENGTHS] = set(validate_list_option(
        options[READ_LENGTHS], int))
    options[READ_DEPTHS] = set(validate_list_option(
        options[READ_DEPTHS], int))
    options[PAIRED_ENDS] = set(validate_list_option(
        options[PAIRED_ENDS], check_boolean_value))
    options[ERRORS] = set(validate_list_option(
        options[ERRORS], check_boolean_value))
    options[BIASES] = set(validate_list_option(
        options[BIASES], check_boolean_value))
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])


def _get_reads_name(length, depth, paired_end, errors, bias):
    return "{d}x_{l}b_{p}_{e}_{b}".format(
        d=depth, l=length,
        p="pe" if paired_end else "se",
        e="errors" if errors else "no_errors",
        b="bias" if bias else "no_bias")


def _get_reads_dir(output_dir, length, depth, paired_end, error, bias):
    return output_dir + os.path.sep + \
        _get_reads_name(length, depth, paired_end, error, bias)


for length, depth, paired_end, error, bias in \
        itertools.product(
            options[READ_LENGTHS], options[READ_DEPTHS],
            options[PAIRED_ENDS], options[ERRORS], options[BIASES]):

    reads_dir = _get_reads_dir(
        options[OUTPUT_DIRECTORY], length, depth, paired_end, error, bias)

    if options[RUN_ONLY] != os.path.exists(reads_dir):
        sys.exit("Reads directory '{d}' should ".format(d=reads_dir) +
                 ("" if options[RUN_ONLY] else "not ") + "already exist.")

for length, depth, paired_end, error, bias in \
        itertools.product(
            options[READ_LENGTHS], options[READ_DEPTHS],
            options[PAIRED_ENDS], options[ERRORS], options[BIASES]):

    reads_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        _get_reads_name(length, depth, paired_end, error, bias)

    if not options[RUN_ONLY]:
        prs.create_simulation_files(
            reads_dir,
            options[TRANSCRIPT_GTF_FILE],
            options[GENOME_FASTA_DIR], options[NUM_FRAGMENTS],
            length, depth, paired_end, error, bias)

    if not options[PREPARE_ONLY]:
        os.chdir(reads_dir)
        args = ['nohup', './run_simulation.sh']
        subprocess.Popen(args)
