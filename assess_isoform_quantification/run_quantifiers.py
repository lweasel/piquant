#!/usr/bin/python

"""Usage:
    run_quantifiers [--log-level=<log-level>] --method=<quant-method> --params=<param-values> [--run-dir=<run-dir] [--input-dir=<input-dir] [--num-fragments=<num-fragments>] [--read-depth=<read-depth>] [--read-length=<read-length>] [--errors] [--paired-end] [--bias] [--prepare-only|--run-only] <transcript-gtf-file> <genome-fasta-dir>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
-d --run-dir=<run-dir>              Directory to create to which run files will be written [default: out].
-m --method=<quant-method>          Method used to quantify transcript abundances.
-p --params=<param-values>          Comma-separated list of key=value parameters required by the specified quantification method.
--input-dir=<input-dir>             Directory containing pre-created Flux Simulator parameters files and simulated reads.
--num-fragments=<num-fragments>     Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--read-depth=<read-depth>           The approximate depth of reads required across the expressed transcriptome [default: 30].
--read-length=<read-length>         The length of sequence reads [default: 50].
--paired-end                        Create and use paired-end sequence reads (note: option must still be specified even if pre-created reads are being used).
--errors                            Flux Simulator will use a position-dependent error model to simulate sequencing errors.
--bias                              Flux Simulator will introduce sequence bias into the RNA fragmentation process.
--prepare-only                      Quantification run scripts will be created, but not run.
--run-only                          Quantification run scripts will be run, but not created (they must already exist).
<transcript-gtf-file>               GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>                  Directory containing per-chromosome sequences as FASTA files.
"""

from docopt import docopt
from schema import Schema, SchemaError

import flux_simulator as fs
import ordutils.log as log
import ordutils.options as opt
import os
import os.path
import prepare_quantification_run as prq
import quantifiers as qs
import subprocess
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
RUN_DIRECTORY = "--run-dir"
INPUT_DIRECTORY = "--input-dir"
QUANT_METHOD = "--method"
PARAMS_SPEC = "--params"
NUM_FRAGMENTS = "--num-fragments"
READ_DEPTH = "--read-depth"
READ_LENGTH = "--read-length"
PAIRED_END = "--paired-end"
ERRORS = "--errors"
BIAS = "--bias"
PREPARE_ONLY = "--prepare-only"
RUN_ONLY = "--run-only"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

TOPHAT_OUTPUT_DIR = "tho"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="prepare_quantification_run v0.1")

# Validate and process command-line options

quant_method_name = options[QUANT_METHOD]

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    opt.validate_dir_option(
        options[RUN_DIRECTORY],
        "Run directory should " + ("" if options[RUN_ONLY] else "not ") +
        "already exist",
        should_exist=options[RUN_ONLY])

    opt.validate_dir_option(
        options[INPUT_DIRECTORY],
        "Input directory should exist if specified",
        nullable=True)

    options[QUANT_METHOD] = opt.validate_dict_option(
        options[QUANT_METHOD], qs.get_quantification_methods(),
        "Unknown quantification method")

    options[PARAMS_SPEC] = Schema(
        options[QUANT_METHOD].get_params_validator()).\
        validate(options[PARAMS_SPEC])
    options[PARAMS_SPEC][qs.TRANSCRIPT_GTF_FILE] = options[TRANSCRIPT_GTF_FILE]
    options[PARAMS_SPEC][qs.GENOME_FASTA_DIR] = options[GENOME_FASTA_DIR]

    options[NUM_FRAGMENTS] = opt.validate_int_option(
        options[NUM_FRAGMENTS],
        "Number of fragments must be a positive integer.",
        nonneg=True)

    options[READ_DEPTH] = opt.validate_int_option(
        options[READ_DEPTH], "Read depth must be a positive integer",
        nonneg=True)

    options[READ_LENGTH] = opt.validate_int_option(
        options[READ_LENGTH], "Read length must be a positive integer",
        nonneg=True)

    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")

    opt.validate_dir_option(
        options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")
except SchemaError as exc:
    exit("Exiting. " + exc.code)

if options[QUANT_METHOD].requires_paired_end_reads() \
        and not options[PAIRED_END]:
    exit("Exiting. Quantification method {m} ".format(m=quant_method_name)
         + "does not support single end reads.")

# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

# Create directory for run files
if not options[RUN_ONLY]:
    logger.info("Creating run directory '{dir}'.".
                format(dir=options[RUN_DIRECTORY]))
    os.mkdir(options[RUN_DIRECTORY])

# Write Flux Simulator parameters files
if not options[RUN_ONLY]:
    logger.info("Creating Flux Simulator parameters files.")

    fs_pro_file = fs.EXPRESSION_PARAMS_FILE.replace("par", "pro")

    if options[INPUT_DIRECTORY]:
        options[INPUT_DIRECTORY] = os.path.abspath(options[INPUT_DIRECTORY])
        fs_pro_file = options[INPUT_DIRECTORY] + os.path.sep + fs_pro_file
    else:
        fs.write_flux_simulator_params_files(
            options[TRANSCRIPT_GTF_FILE],
            options[PARAMS_SPEC][qs.GENOME_FASTA_DIR],
            options[NUM_FRAGMENTS], options[READ_LENGTH], options[PAIRED_END],
            options[ERRORS], options[BIAS], fs_pro_file,
            options[RUN_DIRECTORY])

# Write shell script to run quantification analysis
if not options[RUN_ONLY]:
    logger.info("Creating shell script to run quantification analysis.")

    prq.write_run_quantification_script(
        options[RUN_DIRECTORY], options[INPUT_DIRECTORY],
        options[TRANSCRIPT_GTF_FILE], fs_pro_file,
        options[QUANT_METHOD], options[READ_LENGTH], options[READ_DEPTH],
        options[PAIRED_END], options[ERRORS], options[BIAS],
        dict(options[PARAMS_SPEC]))

# Execute the run quantification script
if not options[PREPARE_ONLY]:
    logger.info("Executing shell script to run quantification analysis.")
    os.chdir(options[RUN_DIRECTORY])

    run_params = "-" + ("" if options[INPUT_DIRECTORY] else "r") + "qa"
    args = ['nohup', './run_quantification.sh', run_params]
    subprocess.Popen(['nohup', './run_quantification.sh', run_params])
