#!/usr/bin/python

# TODO: should be able to separately specify parent directory for reads directories.

"""Usage:
    run_quantifiers prepare [--log-level=<log-level>] [--out-dir=<out-dir>] [--num-fragments=<num-fragments>] --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --polya=<polya> <transcript-gtf-file> <genome-fasta-dir>
    run_quantifiers prequantify [--log-level=<log-level>] [--out-dir=<out-dir>] [--num-fragments=<num-fragments>] --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --polya=<polya> <transcript-gtf-file> <genome-fasta-dir>
    run_quantifiers quantify [--log-level=<log-level>] [--out-dir=<out-dir>] [--num-fragments=<num-fragments>] --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --polya=<polya> <transcript-gtf-file> <genome-fasta-dir>

-h --help                          Show this message.
-v --version                       Show version.
--log-level=<log-level>            Set logging level (one of {log_level_vals}) [default: info].
--out-dir=<out-dir>                Parent output directory to which quantification run directories will be written [default: output].
--num-fragments=<num-fragments>    Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--prepare-only                     Quantification run scripts will be created, but not run.
--run-only                         Quantification run scripts will be run, but not created (they must already exist).
-q --quant-method=<quant-methods>  Comma-separated list of quantification methods to run.
-l --read-length=<read-lengths>    Comma-separated list of read-lengths to perform quantification for.
-d --read-depth=<read-depths>      Comma-separated list of read-depths to perform quantification for.
-p --paired-end=<paired-ends>      Comma-separated list of True/False strings indicating whether quantification should be performed for single or paired-end reads.
-e --error=<errors>                Comma-separated list of True/False strings indicating whether quantification should be performed with or without read errors.
-b --bias=<biases>                 Comma-separated list of True/False strings indicating whether quantification should be performed with or without read sequence bias.
-a --polya=<polya>                 Comma-separated list of True/False strings indicating whether quantification should be performed assuming transcripts do or do not have polyA tails.
<transcript-gtf-file>              GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>                 Directory containing per-chromosome sequences as FASTA files.
"""

from docopt import docopt
from schema import SchemaError

import ordutils.log as log
import ordutils.options as opt
import os
import os.path
import parameters
import prepare_quantification_run as prq
import quantifiers as qs
import subprocess
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
OUTPUT_DIRECTORY = "--out-dir"
NUM_FRAGMENTS = "--num-fragments"
PREPARE = "prepare"
PREQUANTIFY = "prequantify"
QUANTIFY = "quantify"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="prepare_quantification_run v0.1")

# Validate and process command-line options
param_values = None

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    param_values = parameters.validate_command_line_parameter_sets(options)

    options[NUM_FRAGMENTS] = opt.validate_int_option(
        options[NUM_FRAGMENTS],
        "Number of fragments must be a positive integer.",
        nonneg=True)

    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")

    opt.validate_dir_option(
        options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")
except SchemaError as exc:
    exit("Exiting. " + exc.code)

quantifiers_used = []

if False in param_values[parameters.PAIRED_END]:
    for qm in param_values[parameters.QUANT_METHOD]:
        if qm.requires_paired_end_reads():
            exit("Exiting. Quantification method {m} ".format(m=qm.get_name())
                 + "does not support single end reads.")


def check_reads_directory(**params):
    params = dict(params)
    del params[parameters.QUANT_METHOD]
    reads_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)
    if not os.path.exists(reads_dir):
        sys.exit("Reads directory '{d}' should exist.".format(d=reads_dir))


def check_run_directory(**params):
    run_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)
    should_exist = options[PREQUANTIFY] or options[QUANTIFY]
    if should_exist != os.path.exists(run_dir):
        sys.exit("Run directory '{d}' should ".format(d=run_dir) +
                 ("" if should_exist else "not ") + "already exist.")


def execute_quantification_script(run_dir, cl_opts):
    cwd = os.getcwd()
    os.chdir(run_dir)
    args = ['nohup', './run_quantification.sh', cl_opts]
    subprocess.Popen(args)
    os.chdir(cwd)


def create_and_run_quantification(**params):
    run_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)

    if options[PREPARE]:
        # Create the run quantification script
        reads_params = dict(params)
        del reads_params[parameters.QUANT_METHOD]
        reads_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
            parameters.get_file_name(**reads_params)
        quantifier_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
            "quantifier_scratch"

        params_spec = {
            qs.TRANSCRIPT_GTF_FILE: options[TRANSCRIPT_GTF_FILE],
            qs.GENOME_FASTA_DIR: options[GENOME_FASTA_DIR]
        }

        prq.write_run_quantification_script(
            reads_dir, run_dir, quantifier_dir, options[TRANSCRIPT_GTF_FILE],
            params_spec, **params)
    elif options[PREQUANTIFY]:
        # Execute prequantification commands, but only once per quantifier
        quant_method = param_values[parameters.QUANT_METHOD]
        if quant_method not in quantifiers_used:
            quantifiers_used.append(quant_method)
            logger.info("Executing prequantification for " +
                        quant_method.get_name())
            execute_quantification_script(run_dir, "-p")
    elif options[QUANTIFY]:
        # Execute the run quantification script
        logger.info("Executing shell script to run quantification analysis.")
        execute_quantification_script(run_dir, "-qa")


# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

options[OUTPUT_DIRECTORY] = os.path.abspath(options[OUTPUT_DIRECTORY])

parameters.execute_for_param_sets(
    [check_reads_directory, check_run_directory,
     create_and_run_quantification],
    **param_values)
