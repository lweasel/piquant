#!/usr/bin/env python

# TODO: should be able to separately specify parent directory for reads directories.
# TODO: rename out-dir

"""Usage:
    piquant prepare_read_dirs [--log-level=<log-level> --out-dir=<out_dir> --num-fragments=<num-fragments> --nocleanup --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --transcript-gtf=<transcript-gtf-file> --genome-fasta=<genome-fasta-dir>]
    piquant create_reads [--log-level=<log-level> --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_reads [--log-level=<log-level> --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant prepare_quant_dirs [--log-level=<log-level> --out-dir=<out-dir> --nocleanup --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --transcript-gtf=<transcript-gtf-file> --genome-fasta=<genome-fasta-dir> --plot-format=<plot-format> --boxplot-threshold=<boxplot-threshold>]
    piquant prequantify [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant quantify [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_quant [--log-level=<log-level> --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant analyse_runs [--log-level=<log-level> --out-dir=<out-dir> --stats-dir=<stats-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --plot-format=<plot-format>]

Options:
-h --help                                Show this message.
-v --version                             Show version.
--log-level=<log-level>                  Set logging level (one of {log_level_vals}) [default: info].
--out-dir=<out-dir>                      Parent output directory to which quantification run directories will be written [default: output].
--stats-dir=<stats-dir>                  Directory to output assembled stats and graphs to [default: output/analysis].
--num-fragments=<num-fragments>          Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--nocleanup                              If not specified, files non-essential for subsequent quantification (when creating reads) and assessing quantification accuracy (when quantifying) will be deleted.
-f --params-file=<params-file>           File containing specification of quantification methods, read-lengths, read-depths and end, error and bias parameter values to create reads for.
-q --quant-method=<quant-methods>        Comma-separated list of quantification methods to run.
-l --read-length=<read-lengths>          Comma-separated list of read-lengths to perform quantification for.
-d --read-depth=<read-depths>            Comma-separated list of read-depths to perform quantification for.
-p --paired-end=<paired-ends>            Comma-separated list of True/False strings indicating whether quantification should be performed for single or paired-end reads.
-e --error=<errors>                      Comma-separated list of True/False strings indicating whether quantification should be performed with or without read errors.
-b --bias=<biases>                       Comma-separated list of True/False strings indicating whether quantification should be performed with or without read sequence bias.
--transcript-gtf=<transcript-gtf-file>   GTF formatted file describing the transcripts to be simulated.
--genome-fasta=<genome-fasta-dir>        Directory containing per-chromosome sequences as FASTA files.
--plot-format=<plot-format>              Output format for graphs (one of {plot_formats}) [default: pdf].
--boxplot-threshold=<boxplot-threshold>  Minimum number of data points required for a group of transcripts to be shown on a boxplot [default: 300].
"""

import docopt
import flux_simulator as fs
import log
import os
import os.path
import pandas as pd
import parameters
import piquant_options as po
import plot
import prepare_quantification_run as prq
import prepare_read_simulation as prs
import process
import schema
import statistics
import sys
import time


def get_parameters_dir(**params):
    return os.path.join(options[po.OUTPUT_DIRECTORY],
                        parameters.get_file_name(**params))


def reads_directory_checker(should_exist):
    def check_reads_directory(**params):
        params = dict(params)
        if parameters.QUANT_METHOD.name in params:
            del params[parameters.QUANT_METHOD.name]

        reads_dir = get_parameters_dir(**params)
        if should_exist != os.path.exists(reads_dir):
            sys.exit("Reads directory '{d}' should {n}already exist.".
                     format(d=reads_dir, n=("" if should_exist else "not ")))

    return check_reads_directory


def prepare_read_simulation(**params):
    reads_dir = get_parameters_dir(**params)
    cleanup = not options[po.NO_CLEANUP]
    prs.create_simulation_files(reads_dir, cleanup, **params)


def create_reads(**params):
    run_dir = get_parameters_dir(**params)
    process.run_in_directory(run_dir, './run_simulation.sh')


def check_reads_created(**params):
    reads_dir = get_parameters_dir(**params)
    reads_file = fs.get_reads_file(
        params[parameters.ERRORS.name],
        fs.LEFT_READS if params[parameters.PAIRED_END.name] else None)

    if not os.path.exists(os.path.join(reads_dir, reads_file)):
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
    del reads_params[parameters.QUANT_METHOD.name]
    reads_dir = get_parameters_dir(**reads_params)

    prq.write_run_quantification_script(reads_dir, run_dir, options, **params)

quantifiers_used = []


def prequantify(**params):
    run_dir = get_parameters_dir(**params)

    quant_method = params[parameters.QUANT_METHOD.name]
    if quant_method not in quantifiers_used:
        quantifiers_used.append(quant_method)
        logger.info("Executing prequantification for " +
                    quant_method.get_name())
        execute_quantification_script(run_dir, ["-p"])
        time.sleep(1)


def quantify(**params):
    run_dir = get_parameters_dir(**params)

    logger.info("Executing shell script to run quantification analysis.")
    execute_quantification_script(run_dir, ["-qa"])


def check_quantification_completed(**params):
    run_dir = get_parameters_dir(**params)

    main_stats_file = statistics.get_stats_file(
        run_dir, os.path.basename(run_dir))
    if not os.path.exists(main_stats_file):
        run_name = parameters.get_file_name(**params)
        logger.error("Run " + run_name + " did not complete")


class StatsAccumulator:
    ACCUMULATORS = []

    def __init__(self, stratified_stats_type):
        self.overall_stats_df = pd.DataFrame()
        self.stratified_stats_type = stratified_stats_type
        StatsAccumulator.ACCUMULATORS.append(self)

    def __call__(self, **params):
        # TODO: get rid of need to pass run_name into get_stats_file()
        run_name = parameters.get_file_name(**params)
        run_dir = get_parameters_dir(**params)

        stats_file = statistics.get_stats_file(
            run_dir, run_name, **self.stratified_stats_type)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

    def write_stats(self):
        overall_stats_file = statistics.get_stats_file(
            options[po.STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX,
            **self.stratified_stats_type)
        statistics.write_stats_data(
            overall_stats_file, self.overall_stats_df, index=False)


def get_executables_for_commands():
    execs = {}
    execs[po.PREPARE_READ_DIRS] = \
        [reads_directory_checker(False), prepare_read_simulation]
    execs[po.CREATE_READS] = \
        [reads_directory_checker(True), create_reads]
    execs[po.CHECK_READS] = \
        [reads_directory_checker(True), check_reads_created]
    execs[po.PREPARE_QUANT_DIRS] = \
        [run_directory_checker(False), prepare_quantification]
    execs[po.PREQUANTIFY] = \
        [prequantify]
    execs[po.QUANTIFY] = \
        [reads_directory_checker(True), run_directory_checker(True), quantify]
    execs[po.CHECK_QUANTIFICATION] = \
        [run_directory_checker(True), check_quantification_completed]
    execs[po.ANALYSE_RUNS] = \
        [run_directory_checker(True)] + \
        [StatsAccumulator(t) for t in statistics.get_stratified_stats_types()]
    return execs


def get_piquant_command(options):
    return [opt for opt, val in options.items()
            if (val and opt in get_executables_for_commands())][0]


# Read in command-line options
__doc__ = __doc__.format(
    log_level_vals=str(log.LEVELS.keys()),
    plot_formats=plot.PLOT_FORMATS)
options = docopt.docopt(__doc__, version="piquant v0.1")

# Validate and process command-line options
param_values = None
try:
    options, param_values = po.validate_command_line_options(options)
except schema.SchemaError as exc:
    exit("Exiting. " + exc.code)

# Set up logger
logger = log.get_logger(sys.stderr, options[po.LOG_LEVEL])

piquant_command = get_piquant_command(options)
parameters.execute_for_param_sets(
    get_executables_for_commands()[piquant_command], **param_values)

if piquant_command == po.ANALYSE_RUNS:
    if not os.path.exists(options[po.STATS_DIRECTORY]):
        os.mkdir(options[po.STATS_DIRECTORY])
    for stats_acc in StatsAccumulator.ACCUMULATORS:
        stats_acc.write_stats()

    overall_stats_file = statistics.get_stats_file(
        options[po.STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX)
    overall_stats = pd.read_csv(overall_stats_file)

    stats_param_values = \
        {p: overall_stats[p.name].value_counts().index.tolist()
            for p in parameters.get_run_parameters()}

    logger.info("Drawing graphs derived from statistics calculated for the " +
                "whole set of TPMs...")
    plot.draw_overall_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        overall_stats, stats_param_values)

    logger.info("Drawing graphs derived from statistics calculated on " +
                "subsets of TPMs...")
    plot.draw_grouped_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values)

    logger.info("Drawing distribution plots...")
    plot.draw_distribution_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values)
