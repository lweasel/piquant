#!/usr/bin/env python

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
import tpms

from __init__ import __version__


def _get_parameters_dir(run_dir, options, **params):
    """
    Return the path of a reads or quantification directory.

    Return the path of a reads or quantification directory given a dictionary
    of run parameters (e.g. quantification method, read depth etc.).

    run_dir: True if a quantification directory path is to be returned, False
    if the associated reads directory is to be returned.
    options: A dictionary mapping from piquant command line option names to
    option values.
    params: A dictionary mapping from parameters._Parameter instances to
    parameter values.
    """
    params = dict(params)
    if not run_dir and parameters.QUANT_METHOD.name in params:
        del params[parameters.QUANT_METHOD.name]

    dir_option = po.QUANT_OUTPUT_DIR if run_dir else po.READS_OUTPUT_DIR

    return os.path.join(
        options[dir_option], parameters.get_file_name(**params))


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
        reads_dir = _get_parameters_dir(False, options, **params)
        logger.debug("Checking reads directory {d}, should_exist={s}".
                     format(d=reads_dir, s=should_exist))

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
    if not os.path.exists(options[po.READS_OUTPUT_DIR]):
        os.mkdir(options[po.READS_OUTPUT_DIR])

    reads_dir = _get_parameters_dir(False, options, **params)
    cleanup = not options[po.NO_CLEANUP]
    logger.debug("Creating simulation files in " + reads_dir)

    prs.create_simulation_files(reads_dir, cleanup, **params)


def _create_reads(logger, options, **params):
    run_dir = _get_parameters_dir(False, options, **params)
    logger.debug("Creating reads in " + run_dir)

    process.run_in_directory(run_dir, './run_simulation.sh')


def _check_reads_created(logger, options, **params):
    reads_dir = _get_parameters_dir(False, options, **params)
    reads_file = fs.get_reads_file(
        params[parameters.ERRORS.name],
        paired_end=(fs.LEFT_READS if params[parameters.PAIRED_END.name]
                    else None))

    if not os.path.exists(os.path.join(reads_dir, reads_file)):
        run_name = os.path.basename(reads_dir)
        logger.error("Run " + run_name + " did not complete.")


def _run_directory_checker(should_exist):
    def check_run_directory(logger, options, **params):
        run_dir = _get_parameters_dir(True, options, **params)
        logger.debug("Checking run directory {d}, should_exist={s}".
                     format(d=run_dir, s=should_exist))

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
    if not os.path.exists(options[po.QUANT_OUTPUT_DIR]):
        os.mkdir(options[po.QUANT_OUTPUT_DIR])

    reads_dir = _get_parameters_dir(False, options, **params)
    run_dir = _get_parameters_dir(True, options, **params)
    logger.debug("Creating quantification files in " + run_dir)

    prq.write_run_quantification_script(reads_dir, run_dir, options, **params)


def _prequantifier():
    quantifiers_used = []

    def prequantify(logger, options, **params):
        run_dir = _get_parameters_dir(True, options, **params)

        quant_method = params[parameters.QUANT_METHOD.name]
        if quant_method not in quantifiers_used:
            quantifiers_used.append(quant_method)
            logger.info("Executing prequantification for " + str(quant_method))
            _execute_quantification_script(run_dir, ["-p"])
            time.sleep(1)

    return prequantify


def _quantify(logger, options, **params):
    run_dir = _get_parameters_dir(True, options, **params)

    logger.info("Executing shell script to run quantification analysis.")
    _execute_quantification_script(run_dir, ["-qa"])


def _check_quantification_completed(logger, options, **params):
    run_dir = _get_parameters_dir(True, options, **params)

    main_stats_file = statistics.get_stats_file(
        run_dir, os.path.basename(run_dir), tpms.TRANSCRIPT)
    if not os.path.exists(main_stats_file):
        run_name = parameters.get_file_name(**params)
        logger.error("Run " + run_name + " did not complete")


class _StatsAccumulator:
    ACCUMULATORS = []

    def __init__(self, tpm_level, stratified_stats_type):
        self.overall_stats_df = pd.DataFrame()
        self.tpm_level = tpm_level
        self.stratified_stats_type = stratified_stats_type
        _StatsAccumulator.ACCUMULATORS.append(self)

    def __call__(self, logger, options, **params):
        run_name = parameters.get_file_name(**params)
        run_dir = _get_parameters_dir(True, options, **params)

        stats_file = statistics.get_stats_file(
            run_dir, run_name, self.tpm_level, **self.stratified_stats_type)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

    def _write_stats(self, stats_dir):
        overall_stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX,
            self.tpm_level, **self.stratified_stats_type)
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
    execs[po.PREQUANTIFY] = [_prequantifier()]
    execs[po.QUANTIFY] = \
        [_reads_directory_checker(True),
         _run_directory_checker(True), _quantify]
    execs[po.CHECK_QUANTIFICATION] = \
        [_run_directory_checker(True), _check_quantification_completed]
    execs[po.ANALYSE_RUNS] = \
        [_run_directory_checker(True)] + \
        [_StatsAccumulator(tpms.TRANSCRIPT, {})] + \
        [_StatsAccumulator(tpms.GENE, {})] + \
        [_StatsAccumulator(tpms.TRANSCRIPT, t) for t
            in statistics.get_stratified_stats_types()]
    return execs


def _write_accumulated_stats(options):
    if not os.path.exists(options[po.STATS_DIRECTORY]):
        os.mkdir(options[po.STATS_DIRECTORY])
    for stats_acc in _StatsAccumulator.ACCUMULATORS:
        stats_acc._write_stats(options[po.STATS_DIRECTORY])


def _get_overall_stats(options, tpm_level):
    overall_stats_file = statistics.get_stats_file(
        options[po.STATS_DIRECTORY], statistics.OVERALL_STATS_PREFIX, tpm_level)
    return pd.read_csv(overall_stats_file)


def _get_stats_param_values(overall_stats):
    return {p: overall_stats[p.name].value_counts().index.tolist()
            for p in parameters.get_run_parameters()}


def _draw_overall_stats_graphs(
        logger, options, overall_stats, stats_param_values, tpm_level):
    logger.info("Drawing graphs derived from statistics calculated for the " +
                "whole set of TPMs...")
    plot.draw_overall_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        overall_stats, stats_param_values, tpm_level)


def _draw_grouped_stats_graphs(logger, options, stats_param_values):
    logger.info("Drawing graphs derived from statistics calculated on " +
                "subsets of TPMs...")
    plot.draw_grouped_stats_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values, options[po.GROUPED_THRESHOLD])


def _draw_distribution_graphs(logger, options, stats_param_values):
    logger.info("Drawing distribution plots...")
    plot.draw_distribution_graphs(
        options[po.PLOT_FORMAT], options[po.STATS_DIRECTORY],
        stats_param_values)


def _analyse_runs(logger, options):
    _write_accumulated_stats(options)

    overall_transcript_stats = _get_overall_stats(options, tpms.TRANSCRIPT)
    transcript_stats_param_values = \
        _get_stats_param_values(overall_transcript_stats)

    overall_gene_stats = _get_overall_stats(options, tpms.GENE)
    gene_stats_param_values = _get_stats_param_values(overall_gene_stats)

    _draw_overall_stats_graphs(
        logger, options, overall_transcript_stats,
        transcript_stats_param_values, tpms.TRANSCRIPT)
    _draw_overall_stats_graphs(
        logger, options, overall_gene_stats, gene_stats_param_values, tpms.GENE)
    _draw_grouped_stats_graphs(
        logger, options, transcript_stats_param_values)
    _draw_distribution_graphs(
        logger, options, transcript_stats_param_values)


def _run_piquant_command(logger, piquant_command, options, param_values):
    parameters.execute_for_param_sets(
        _get_executables_for_commands()[piquant_command],
        logger, options, **param_values)

    if piquant_command == po.ANALYSE_RUNS:
        _analyse_runs(logger, options)


def piquant(args):
    # Check that a piquant command has been specified
    if len(args) < 2:
        exit(("Exiting - one of the following piquant commands must be " +
              "specified:\n    {cs}\nTo get help on a particular command, " +
              "run that command with the --help option.\n\nFor further " +
              "documentation on piquant, please see " +
              "http://piquant.readthedocs.org/en/latest/.").format(
             cs="\n    ".join(po.get_command_names())))

    command_name, args = sys.argv[1], sys.argv[2:]

    # Check the specified command is valid
    command = po.get_command(command_name)

    # Read command-line options
    usage = po.get_usage_message(command)
    options = docopt.docopt(usage, version="piquant v" + __version__)

    # Validate and process command-line options
    param_values = None
    try:
        options, param_values = \
            po.validate_command_line_options(command, options)
    except schema.SchemaError as exc:
        exit("Exiting. " + exc.code)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Run the specified piquant command
    _run_piquant_command(logger, command, options, param_values)

if __name__ == "__main__":
    piquant(sys.argv)
