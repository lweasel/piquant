#!/usr/bin/python

# TODO: should be able to separately specify parent directory for reads directories.
# TODO: rename out-dir

"""Usage:
    piquant prepare_read_dirs [--log-level=<log-level> --out-dir=<out_dir> --num-fragments=<num-fragments> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> <transcript-gtf-file> <genome-fasta-dir>]
    piquant create_reads [--log-level=<log-level> --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_reads [--log-level=<log-level> --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant prepare_quant_dirs [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>] <transcript-gtf-file> <genome-fasta-dir>
    piquant prequantify [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant quantify [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_quant [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant analyse_runs [--log-level=<log-level> --out-dir=<out-dir> --stats-dir=<stats-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]

Options:
-h --help                          Show this message.
-v --version                       Show version.
--log-level=<log-level>            Set logging level (one of {log_level_vals}) [default: info].
--out-dir=<out-dir>                Parent output directory to which quantification run directories will be written [default: output].
--stats-dir=<stats-dir>            Directory to output assembled stats and graphs to [default: output/overall_stats].
--num-fragments=<num-fragments>    Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
-f --params-file=<params-file>     File containing specification of quantification methods, read-lengths, read-depths and end, error and bias parameter values to create reads for.
-q --quant-method=<quant-methods>  Comma-separated list of quantification methods to run.
-l --read-length=<read-lengths>    Comma-separated list of read-lengths to perform quantification for.
-d --read-depth=<read-depths>      Comma-separated list of read-depths to perform quantification for.
-p --paired-end=<paired-ends>      Comma-separated list of True/False strings indicating whether quantification should be performed for single or paired-end reads.
-e --error=<errors>                Comma-separated list of True/False strings indicating whether quantification should be performed with or without read errors.
-b --bias=<biases>                 Comma-separated list of True/False strings indicating whether quantification should be performed with or without read sequence bias.
<transcript-gtf-file>              GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>                 Directory containing per-chromosome sequences as FASTA files.
"""

import docopt
import flux_simulator as fs
import ordutils.log as log
import ordutils.options as opt
import os
import os.path
import pandas as pd
import parameters
import prepare_quantification_run as prq
import prepare_read_simulation as prs
import process
import quantifiers as qs
import schema
import statistics
import sys
import time
import plot

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
OUTPUT_DIRECTORY = "--out-dir"
STATS_DIRECTORY = "--stats-dir"
NUM_FRAGMENTS = "--num-fragments"
PARAMS_FILE = "--params-file"
PREPARE_READ_DIRS = "prepare_read_dirs"
CREATE_READS = "create_reads"
CHECK_READS = "check_reads"
PREPARE_QUANT_DIRS = "prepare_quant_dirs"
PREQUANTIFY = "prequantify"
QUANTIFY = "quantify"
CHECK_QUANTIFICATION = "check_quant"
ANALYSE_RUNS = "analyse_runs"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt.docopt(__doc__, version="piquant v0.1")

# Validate and process command-line options
param_values = None

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    opt.validate_dir_option(
        options[OUTPUT_DIRECTORY], "Output parent directory does not exist")
    options[NUM_FRAGMENTS] = opt.validate_int_option(
        options[NUM_FRAGMENTS],
        "Number of fragments must be a positive integer",
        nonneg=True, nullable=True)

    if options[PREPARE_READ_DIRS] or options[PREPARE_QUANT_DIRS]:
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")
        opt.validate_dir_option(
            options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")

    opt.validate_file_option(
        options[PARAMS_FILE],
        "Parameter specification file should exist",
        nullable=True)

    processing_reads = options[PREPARE_READ_DIRS] or \
        options[CREATE_READS] or options[CHECK_READS]

    ignore_params = [parameters.QUANT_METHOD] if processing_reads else []
    param_values = parameters.validate_command_line_parameter_sets(
        options[PARAMS_FILE], options, ignore_params=ignore_params)

    if not processing_reads and False in param_values[parameters.PAIRED_END]:
        for qm in param_values[parameters.QUANT_METHOD]:
            if qm.requires_paired_end_reads():
                raise schema.SchemaError(
                    None, "Quantification method " + qm.get_name() +
                    " does not support single-end reads.")
except schema.SchemaError as exc:
    exit("Exiting. " + exc.code)


def get_parameters_dir(**params):
    return options[OUTPUT_DIRECTORY] + os.path.sep + \
        parameters.get_file_name(**params)


def reads_directory_checker(should_exist):
    def check_reads_directory(**params):
        params = dict(params)
        if parameters.QUANT_METHOD in params:
            del params[parameters.QUANT_METHOD]

        reads_dir = get_parameters_dir(**params)
        if should_exist != os.path.exists(reads_dir):
            sys.exit("Reads directory '{d}' should {n}already exist.".
                     format(d=reads_dir, n=("" if should_exist else "not ")))

    return check_reads_directory


def prepare_read_simulation(**params):
    reads_dir = get_parameters_dir(**params)
    prs.create_simulation_files(
        reads_dir, options[TRANSCRIPT_GTF_FILE], options[GENOME_FASTA_DIR],
        options[NUM_FRAGMENTS], **params)


def create_reads(**params):
    run_dir = get_parameters_dir(**params)
    process.run_in_directory(run_dir, './run_simulation.sh')


def check_reads_created(**params):
    reads_dir = get_parameters_dir(**params)
    reads_file = fs.get_reads_file(
        params[parameters.ERRORS],
        'l' if params[parameters.PAIRED_END] else None)
    reads_file_path = reads_dir + os.path.sep + reads_file
    if not os.path.exists(reads_file_path):
        run_name = os.path.basename(reads_dir)
        logger.error("Run " + run_name + " did not complete.")


def run_directory_checker(should_exist):
    def check_run_directory(**params):
        run_dir = get_parameters_dir(**params)
        if should_exist != os.path.exists(run_dir):
            sys.exit("Run directory '{d}' should {n}already exist.".
                     format(d=run_dir, n=("" if should_exist else "not ")))

    return check_run_directory


def execute_quantification_script(run_dir, cl_opts):
    process.run_in_directory(run_dir, './run_quantification.sh', cl_opts)


def prepare_quantification(**params):
    run_dir = get_parameters_dir(**params)

    reads_params = dict(params)
    del reads_params[parameters.QUANT_METHOD]
    reads_dir = get_parameters_dir(**params)

    quantifier_dir = options[OUTPUT_DIRECTORY] + os.path.sep + \
        "quantifier_scratch"

    params_spec = {
        qs.TRANSCRIPT_GTF_FILE: options[TRANSCRIPT_GTF_FILE],
        qs.GENOME_FASTA_DIR: options[GENOME_FASTA_DIR]
    }

    prq.write_run_quantification_script(
        reads_dir, run_dir, quantifier_dir, options[TRANSCRIPT_GTF_FILE],
        params_spec, **params)

quantifiers_used = []


def prequantify(**params):
    run_dir = get_parameters_dir(**params)

    quant_method = params[parameters.QUANT_METHOD]
    if quant_method not in quantifiers_used:
        quantifiers_used.append(quant_method)
        logger.info("Executing prequantification for " +
                    quant_method.get_name())
        execute_quantification_script(run_dir, "-p")
        time.sleep(1)


def quantify(**params):
    run_dir = get_parameters_dir(**params)

    logger.info("Executing shell script to run quantification analysis.")
    execute_quantification_script(run_dir, "-qa")


def check_quantification_completed(**params):
    run_dir = get_parameters_dir(**params)

    main_stats_file = statistics.get_stats_file(
        run_dir, os.path.basename(run_dir))
    if not os.path.exists(main_stats_file):
        run_name = parameters.get_file_name(**params)
        logger.error("Run " + run_name + " did not complete")


class StatsAccumulator:
    def __init__(self, pset):
        self.overall_stats_df = pd.DataFrame()
        self.pset = pset

    def __call__(self, **params):
        # TODO: get rid of need to pass run_name into get_stats_file()
        run_name = parameters.get_file_name(**params)
        run_dir = get_parameters_dir(**params)

        stats_file = statistics.get_stats_file(run_dir, run_name, **self.pset)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

    def write_stats(self):
        overall_stats_file = statistics.get_stats_file(
            options[STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX,
            **self.pset)
        statistics.write_stats_data(
            overall_stats_file, self.overall_stats_df, index=False)


# Set up logger
logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

options[OUTPUT_DIRECTORY] = os.path.abspath(options[OUTPUT_DIRECTORY])

stats_accumulators = []

if options[PREPARE_READ_DIRS]:
    to_execute = [reads_directory_checker(False),
                  prepare_read_simulation]
elif options[CREATE_READS]:
    to_execute = [reads_directory_checker(True),
                  create_reads]
elif options[CHECK_READS]:
    to_execute = [reads_directory_checker(True),
                  check_reads_created]
elif options[PREPARE_QUANT_DIRS]:
    to_execute = [run_directory_checker(False),
                  prepare_quantification]
elif options[PREQUANTIFY]:
    to_execute = [prequantify]
elif options[QUANTIFY]:
    to_execute = [reads_directory_checker(True),
                  run_directory_checker(True),
                  quantify]
elif options[CHECK_QUANTIFICATION]:
    to_execute = [run_directory_checker(True),
                  check_quantification_completed]
elif options[ANALYSE_RUNS]:
    to_execute = [run_directory_checker(True)]
    for pset in statistics.get_stats_param_sets():
        stats_acc = StatsAccumulator(pset)
        to_execute.append(stats_acc)
        stats_accumulators.append(stats_acc)

parameters.execute_for_param_sets(to_execute, **param_values)

if options[ANALYSE_RUNS]:
    if not os.path.exists(options[STATS_DIRECTORY]):
        os.mkdir(options[STATS_DIRECTORY])
    for stats_acc in stats_accumulators:
        stats_acc.write_stats()

    overall_stats_file = statistics.get_stats_file(
        options[STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX)
    overall_stats = pd.read_csv(overall_stats_file)

    param_values = {}
    for param in parameters.get_parameters():
        param_series = overall_stats[param.name]
        param_values[param] = param_series.value_counts().index.tolist()

    logger.info("Drawing graphs derived from statistics calculated for the " +
                "whole set of TPMs...")
    plot.draw_overall_stats_graphs(
        options[STATS_DIRECTORY], overall_stats, param_values)

    logger.info("Drawing graphs derived from statistics calculated on " +
                "subsets of TPMs...")
    plot.draw_grouped_stats_graphs(options[STATS_DIRECTORY], param_values)

    logger.info("Drawing distribution plots...")
    plot.draw_distribution_graphs(options[STATS_DIRECTORY], param_values)
