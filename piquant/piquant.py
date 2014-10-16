#!/usr/bin/env python

"""Usage:
    piquant prepare_read_dirs [{log_option_spec} --out-dir=<out_dir> --num-molecules=<num-molecules> --nocleanup --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --transcript-gtf=<transcript-gtf-file> --genome-fasta=<genome-fasta-dir>]
    piquant create_reads [{log_option_spec} --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_reads [{log_option_spec} --out-dir=<out_dir> --params-file=<params-file> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant prepare_quant_dirs [{log_option_spec} --out-dir=<out-dir> --nocleanup --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --transcript-gtf=<transcript-gtf-file> --genome-fasta=<genome-fasta-dir> --plot-format=<plot-format> --grouped-threshold=<threshold>]
    piquant prequantify [{log_option_spec} --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant quantify [{log_option_spec} --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant check_quant [{log_option_spec} --out-dir=<out-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases>]
    piquant analyse_runs [{log_option_spec} --out-dir=<out-dir> --stats-dir=<stats-dir> --params-file=<params-file> --quant-method=<quant-methods> --read-length=<read-lengths> --read-depth=<read-depths> --paired-end=<paired-ends> --error=<errors> --bias=<biases> --plot-format=<plot-format>]

Options:
{help_option_spec}                                {help_option_description}
{ver_option_spec}                             {ver_option_description}
{log_option_spec}                  {log_option_description}
--out-dir=<out-dir>                      Parent output directory to which quantification run directories will be written [default: output].
--stats-dir=<stats-dir>                  Directory to output assembled stats and graphs to [default: output/analysis].
--num-molecules=<num-molecules>          Flux Simulator parameters will be set for simulation to start with this number of transcript molecules in the initial population [default: 30000000].
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
--grouped-threshold=<threshold>          Minimum number of data points required for a group of transcripts to be shown on a plot [default: 300].
"""

import docopt
import flux_simulator as fs
import options as opt
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


def _get_parameters_dir(options, **params):
    """
    Return the path of a reads or quantification directory.

    Return the path of a reads or quantification directory given a dictionary
    of run parameters (e.g. quantification method, read depth etc.).

    options: A dictionary mapping from piquant command line option names to
    option values.
    params: A dictionary mapping from parameters._Parameter instances to
    parameter values.
    """
    return os.path.join(options[po.OUTPUT_DIRECTORY],
                        parameters.get_file_name(**params))


def _reads_directory_checker(should_exist):
    """
    Return a function checking the existence of a reads directory.

    Return a function which, when called, will exit the Python interpreter if
    the specified reads directory does or doesn't exists.

    should_exist: If True, the returned function will exit if the specified
    reads or quantification directory does not exist. If False, the function
    will exit if the reads or quantification directory does already exist.
    """
    def check_reads_directory(logger, options, **params):
        params = dict(params)
        if parameters.QUANT_METHOD.name in params:
            del params[parameters.QUANT_METHOD.name]

        reads_dir = _get_parameters_dir(options, **params)
        if should_exist != os.path.exists(reads_dir):
            sys.exit("Reads directory '{d}' should {n}already exist.".
                     format(d=reads_dir, n=("" if should_exist else "not ")))

    return check_reads_directory


def _prepare_read_simulation(logger, options, **params):
    """
    Write bash script and support files to perform RNA-seq read simulation.

    Write a bash script and support files such that when the bash script is
    executed, Flux Simulator will be used to simulate RNA-seq reads for the
    simulation parameters encapsulated by 'params'.

    logger: Logs messages to standard error.
    options: A dictionary mapping from piquant command line option names to
    option values.
    params: A dictionary mapping from parameters._Parameter instances to
    parameter values, describing properties of the read simulation to be
    performed.
    """
    reads_dir = _get_parameters_dir(options, **params)
    cleanup = not options[po.NO_CLEANUP]
    prs.create_simulation_files(reads_dir, cleanup, **params)


def _create_reads(logger, options, **params):
    run_dir = _get_parameters_dir(options, **params)
    process.run_in_directory(run_dir, './run_simulation.sh')


def _check_reads_created(logger, options, **params):
    reads_dir = _get_parameters_dir(options, **params)
    reads_file = fs.get_reads_file(
        params[parameters.ERRORS.name],
        fs.LEFT_READS if params[parameters.PAIRED_END.name] else None)

    if not os.path.exists(os.path.join(reads_dir, reads_file)):
        run_name = os.path.basename(reads_dir)
        logger.error("Run " + run_name + " did not complete.")


def _run_directory_checker(should_exist):
    def check_run_directory(logger, options, **params):
        run_dir = _get_parameters_dir(options, **params)
        if should_exist != os.path.exists(run_dir):
            sys.exit("Run directory '{d}' should {n}already exist.".
                     format(d=run_dir, n=("" if should_exist else "not ")))

    return check_run_directory


def _execute_quantification_script(run_dir, cl_opts):
    process.run_in_directory(run_dir, './run_quantification.sh', cl_opts)


def _prepare_quantification(logger, options, **params):
    """
    Write bash script to perform transcriptome quantification.

    Write a bash script which, when executed, with use a specified
    transcriptome quantification tool to estimate transcript abundances from a
    set of simulated RNA-seq reads.

    logger: Logs messages to standard error.
    options: A dictionary mapping from piquant command line option names to
    option values.
    params: A dictionary mapping from parameters._Parameter instances to
    parameter values, describing properties of the quantification run to be
    performed.
    """
    run_dir = _get_parameters_dir(options, **params)

    reads_params = dict(params)
    del reads_params[parameters.QUANT_METHOD.name]
    reads_dir = _get_parameters_dir(options, **reads_params)

    prq.write_run_quantification_script(reads_dir, run_dir, options, **params)

quantifiers_used = []


def _prequantify(logger, options, **params):
    run_dir = _get_parameters_dir(options, **params)

    quant_method = params[parameters.QUANT_METHOD.name]
    if quant_method not in quantifiers_used:
        quantifiers_used.append(quant_method)
        logger.info("Executing prequantification for " + str(quant_method))
        _execute_quantification_script(run_dir, ["-p"])
        time.sleep(1)


def _quantify(logger, options, **params):
    run_dir = _get_parameters_dir(options, **params)

    logger.info("Executing shell script to run quantification analysis.")
    _execute_quantification_script(run_dir, ["-qa"])


def _check_quantification_completed(logger, options, **params):
    run_dir = _get_parameters_dir(options, **params)

    main_stats_file = statistics.get_stats_file(
        run_dir, os.path.basename(run_dir))
    if not os.path.exists(main_stats_file):
        run_name = parameters.get_file_name(**params)
        logger.error("Run " + run_name + " did not complete")


class _StatsAccumulator:
    ACCUMULATORS = []

    def __init__(self, stratified_stats_type):
        self.overall_stats_df = pd.DataFrame()
        self.stratified_stats_type = stratified_stats_type
        _StatsAccumulator.ACCUMULATORS.append(self)

    def __call__(self, logger, options, **params):
        run_name = parameters.get_file_name(**params)
        run_dir = _get_parameters_dir(options, **params)

        stats_file = statistics.get_stats_file(
            run_dir, run_name, **self.stratified_stats_type)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

    def _write_stats(self):
        overall_stats_file = statistics.get_stats_file(
            options[po.STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX,
            **self.stratified_stats_type)
        statistics.write_stats_data(
            overall_stats_file, self.overall_stats_df, index=False)


def _get_executables_for_commands():
    execs = {}
    execs[po.PREPARE_READ_DIRS] = \
        [_reads_directory_checker(False), _prepare_read_simulation]
    execs[po.CREATE_READS] = \
        [_reads_directory_checker(True), _create_reads]
    execs[po.CHECK_READS] = \
        [_reads_directory_checker(True), _check_reads_created]
    execs[po.PREPARE_QUANT_DIRS] = \
        [_run_directory_checker(False), _prepare_quantification]
    execs[po.PREQUANTIFY] = \
        [_prequantify]
    execs[po.QUANTIFY] = \
        [_reads_directory_checker(True),
         _run_directory_checker(True), _quantify]
    execs[po.CHECK_QUANTIFICATION] = \
        [_run_directory_checker(True), _check_quantification_completed]
    execs[po.ANALYSE_RUNS] = \
        [_run_directory_checker(True)] + \
        [_StatsAccumulator(t) for t in statistics.get_stratified_stats_types()]
    return execs


def _get_piquant_command(options):
    return [opt for opt, val in options.items()
            if (val and opt in _get_executables_for_commands())][0]


def _write_accumulated_stats(options):
    if not os.path.exists(options[po.STATS_DIRECTORY]):
        os.mkdir(options[po.STATS_DIRECTORY])
    for stats_acc in _StatsAccumulator.ACCUMULATORS:
        stats_acc._write_stats()


def _get_overall_stats(options):
    overall_stats_file = statistics.get_stats_file(
        options[po.STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX)
    return pd.read_csv(overall_stats_file)


def _get_stats_param_values(overall_stats):
    return {p: overall_stats[p.name].value_counts().index.tolist()
            for p in parameters.get_run_parameters()}


def _draw_overall_stats_graphs(options, overall_stats, stats_param_values):
    logger.info("Drawing graphs derived from statistics calculated for the " +
                "whole set of TPMs...")
    plot.draw_overall_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        overall_stats, stats_param_values)


def _draw_grouped_stats_graphs(options, stats_param_values):
    logger.info("Drawing graphs derived from statistics calculated on " +
                "subsets of TPMs...")
    plot.draw_grouped_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values, options[po.GROUPED_THRESHOLD])


def _draw_distribution_graphs(options, stats_param_values):
    logger.info("Drawing distribution plots...")
    plot.draw_distribution_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values)


def _analyse_runs():
    _write_accumulated_stats(options)

    overall_stats = _get_overall_stats(options)
    stats_param_values = _get_stats_param_values(overall_stats)

    _draw_overall_stats_graphs(options, overall_stats, stats_param_values)
    _draw_grouped_stats_graphs(options, stats_param_values)
    _draw_distribution_graphs(options, stats_param_values)


def _run_piquant_command(logger, options):
    piquant_command = _get_piquant_command(options)

    parameters.execute_for_param_sets(
        _get_executables_for_commands()[piquant_command],
        logger, options, **param_values)

    if piquant_command == po.ANALYSE_RUNS:
        _analyse_runs()


if __name__ == "__main__":
    # Read in command-line options
    __doc__ = opt.substitute_common_options_into_usage(
        __doc__, plot_formats=plot.PLOT_FORMATS)
    options = docopt.docopt(__doc__, version="piquant v0.1")

    # Validate and process command-line options
    param_values = None
    try:
        options, param_values = po.validate_command_line_options(options)
    except schema.SchemaError as exc:
        exit("Exiting. " + exc.code)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Run the specified piquant command
    _run_piquant_command(logger, options)
