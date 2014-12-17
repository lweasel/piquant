import docopt
import os
import os.path
import pandas as pd
import schema
import sys
import time

from . import flux_simulator as fs
from . import options as opt
from . import piquant_commands as pc
from . import piquant_options as po
from . import plot
from . import prepare_quantification_run as prq
from . import prepare_read_simulation as prs
from . import process
from . import resource_usage as ru
from . import statistics
from . import tpms
from .__init__ import __version__


def _get_options_dir(run_dir, options, **qr_options):
    """
    Return the path of a reads or quantification directory.

    Return the path of a reads or quantification directory given a dictionary
    of quantification run options (e.g. quantification method, read depth
    etc.).

    run_dir: True if a quantification directory path is to be returned, False
    if the associated reads directory is to be returned.
    options: A dictionary mapping from piquant command line option names to
    option values.
    qr_options: A dictionary mapping from piquant_options._QuantRunOption
    instances to option values.
    """
    qr_options = dict(qr_options)
    if not run_dir and po.QUANT_METHOD.name in qr_options:
        del qr_options[po.QUANT_METHOD.name]

    dir_option = po.QUANT_OUTPUT_DIR if run_dir else po.READS_OUTPUT_DIR

    return os.path.join(
        options[dir_option.name], po.get_run_name(qr_options))


def _reads_directory_checker(should_exist):
    """
    Return a function checking the existence of a reads directory.

    Return a function which, when called, will exit the Python interpreter if
    the specified reads directory does or doesn't exists.

    should_exist: If True, the returned function will exit if the specified
    reads or quantification directory does not exist. If False, the function
    will exit if the reads or quantification directory does already exist.
    """
    def check_reads_directory(logger, options, **qr_options):
        reads_dir = _get_options_dir(False, options, **qr_options)
        logger.debug("Checking reads directory {d}, should_exist={s}".
                     format(d=reads_dir, s=should_exist))

        if should_exist != os.path.exists(reads_dir):
            sys.exit("Reads directory '{d}' should {n}already exist.".
                     format(d=reads_dir, n=("" if should_exist else "not ")))

    return check_reads_directory


def _prepare_read_simulation(logger, options, **qr_options):
    """
    Write bash script and support files to perform RNA-seq read simulation.

    Write a bash script and support files such that when the bash script is
    executed, Flux Simulator will be used to simulate RNA-seq reads for the
    quantification run options encapsulated by 'qr_options'.

    logger: Logs messages to standard error.
    options: A dictionary mapping from piquant command line option names to
    option values.
    qr_options: A dictionary mapping from piquant_options._QuantRunOption
    instances to option values, describing properties of the read simulation
    to be performed.
    """
    if not os.path.exists(options[po.READS_OUTPUT_DIR.name]):
        os.mkdir(options[po.READS_OUTPUT_DIR.name])

    reads_dir = _get_options_dir(False, options, **qr_options)
    cleanup = not options[po.NO_CLEANUP.name]
    logger.debug("Creating simulation files in " + reads_dir)

    prs.create_simulation_files(reads_dir, cleanup, **qr_options)


def _create_reads(logger, options, **qr_options):
    run_dir = _get_options_dir(False, options, **qr_options)
    logger.debug("Creating reads in " + run_dir)

    process.run_in_directory(run_dir, './run_simulation.sh')


def _check_reads_created(logger, options, **qr_options):
    reads_dir = _get_options_dir(False, options, **qr_options)
    reads_file = fs.get_reads_file(
        qr_options[po.ERRORS.name],
        paired_end=(fs.LEFT_READS if qr_options[po.PAIRED_END.name]
                    else None))

    if not os.path.exists(os.path.join(reads_dir, reads_file)):
        run_name = os.path.basename(reads_dir)
        logger.error("Run " + run_name + " did not complete.")


def _run_directory_checker(should_exist):
    def check_run_directory(logger, options, **qr_options):
        run_dir = _get_options_dir(True, options, **qr_options)
        logger.debug("Checking run directory {d}, should_exist={s}".
                     format(d=run_dir, s=should_exist))

        if should_exist != os.path.exists(run_dir):
            sys.exit("Run directory '{d}' should {n}already exist.".
                     format(d=run_dir, n=("" if should_exist else "not ")))

    return check_run_directory


def _execute_quantification_script(run_dir, cl_opts):
    process.run_in_directory(run_dir, './run_quantification.sh', cl_opts)


def _prepare_quantification(logger, options, **qr_options):
    """
    Write bash script to perform transcriptome quantification.

    Write a bash script which, when executed, with use a specified
    transcriptome quantification tool to estimate transcript abundances from a
    set of simulated RNA-seq reads.

    logger: Logs messages to standard error.
    options: A dictionary mapping from piquant command line option names to
    option values.
    qr_options: A dictionary mapping from piquant_options._QuantRunOption
    instances to option values, describing properties of the quantification run
    to be performed.
    """
    if not os.path.exists(options[po.QUANT_OUTPUT_DIR.name]):
        os.mkdir(options[po.QUANT_OUTPUT_DIR.name])

    reads_dir = _get_options_dir(False, options, **qr_options)
    run_dir = _get_options_dir(True, options, **qr_options)
    logger.debug("Creating quantification files in " + run_dir)

    prq.write_script(reads_dir, run_dir, options, **qr_options)


def _prequantifier():
    quantifiers_used = []

    def prequantify(logger, options, **qr_options):
        run_dir = _get_options_dir(True, options, **qr_options)

        quant_method = qr_options[po.QUANT_METHOD.name]
        if quant_method not in quantifiers_used:
            quantifiers_used.append(quant_method)
            logger.info("Executing prequantification for " + str(quant_method))
            _execute_quantification_script(run_dir, ["-p"])
            time.sleep(1)

    return prequantify


def _quantify(logger, options, **qr_options):
    run_dir = _get_options_dir(True, options, **qr_options)

    logger.info("Executing shell script to run quantification analysis.")
    _execute_quantification_script(run_dir, ["-qa"])


def _check_quantification_completed(logger, options, **qr_options):
    run_dir = _get_options_dir(True, options, **qr_options)

    main_stats_file = statistics.get_stats_file(
        run_dir, os.path.basename(run_dir), tpms.TRANSCRIPT)
    if not os.path.exists(main_stats_file):
        run_name = po.get_run_name(qr_options)
        logger.error("Run " + run_name + " did not complete")


class _StatsAccumulator(object):
    ACCUMULATORS = []

    def __init__(self, tpm_level, classifier=None, ascending=False):
        self.overall_stats_df = pd.DataFrame()
        self.tpm_level = tpm_level
        self.classifier = classifier
        self.ascending = ascending
        _StatsAccumulator.ACCUMULATORS.append(self)

    def __call__(self, logger, options, **qr_options):
        run_name = po.get_run_name(qr_options)
        run_dir = _get_options_dir(True, options, **qr_options)

        stats_file = statistics.get_stats_file(
            run_dir, run_name, self.tpm_level, self.classifier, self.ascending)
        stats_df = pd.read_csv(stats_file)
        self.overall_stats_df = self.overall_stats_df.append(stats_df)

    def write_accumulated_data(self, stats_dir):
        overall_stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX,
            self.tpm_level, self.classifier, self.ascending)
        statistics.write_stats_data(
            overall_stats_file, self.overall_stats_df, index=False)


class _ResourceUsageAccumulator(object):
    ACCUMULATORS = []

    def __init__(self, resource_type):
        self.resource_usage_df = pd.DataFrame()
        self.resource_type = resource_type

    def __call__(self, logger, options, **qr_options):
        run_name = po.get_run_name(qr_options)
        run_dir = _get_options_dir(True, options, **qr_options)

        usage_file = ru.get_resource_usage_file(
            self.resource_type, prefix=run_name, directory=run_dir)

        if os.path.exists(usage_file):
            usage_df = pd.read_csv(usage_file)
            self.resource_usage_df = self.resource_usage_df.append(usage_df)

    def write_accumulated_data(self, stats_dir):
        usage_file_name = ru.get_resource_usage_file(
            self.resource_type, prefix=ru.OVERALL_USAGE_PREFIX,
            directory=stats_dir)
        ru.write_usage_summary(usage_file_name, self.resource_usage_df)


def _set_executables_for_commands(record_usage):
    pc.PREPARE_READ_DIRS.executables = [
        _reads_directory_checker(False),
        _prepare_read_simulation]
    pc.CREATE_READS.executables = [
        _reads_directory_checker(True),
        _create_reads]
    pc.CHECK_READS.executables = [
        _reads_directory_checker(True),
        _check_reads_created]
    pc.PREPARE_QUANT_DIRS.executables = [
        _run_directory_checker(False),
        _prepare_quantification]
    pc.PREQUANTIFY.executables = [
        _run_directory_checker(True),
        _prequantifier()]
    pc.QUANTIFY.executables = [
        _reads_directory_checker(True),
        _run_directory_checker(True),
        _quantify]
    pc.CHECK_QUANTIFICATION.executables = [
        _run_directory_checker(True),
        _check_quantification_completed]

    pc.ANALYSE_RUNS.executables = [
        _run_directory_checker(True),
        _StatsAccumulator(tpms.TRANSCRIPT),
        _StatsAccumulator(tpms.GENE)] + \
        [_StatsAccumulator(tpms.TRANSCRIPT, classifier=clsfr, ascending=asc)
            for clsfr, asc in statistics.get_stratified_stats_types()]

    if record_usage:
        pc.ANALYSE_RUNS.executables += [
            _ResourceUsageAccumulator(ru.PREQUANT_RESOURCE_TYPE),
            _ResourceUsageAccumulator(ru.QUANT_RESOURCE_TYPE),
        ]


def _write_accumulated_stats_and_usage(options):
    stats_dir = options[po.STATS_DIRECTORY.name]
    if not os.path.exists(stats_dir):
        os.mkdir(stats_dir)
    for acc in _StatsAccumulator.ACCUMULATORS + \
            _ResourceUsageAccumulator.ACCUMULATORS:
        acc.write_accumulated_data(stats_dir)


def _get_overall_stats(options, tpm_level):
    overall_stats_file = statistics.get_stats_file(
        options[po.STATS_DIRECTORY.name],
        statistics.OVERALL_STATS_PREFIX, tpm_level)
    return pd.read_csv(overall_stats_file)


def _get_overall_usage(options, resource_type):
    overall_usage_file = ru.get_resource_usage_file(
        resource_type, prefix=ru.OVERALL_USAGE_PREFIX,
        directory=options[po.STATS_DIRECTORY.name])
    return pd.read_csv(overall_usage_file)


def _get_stats_option_values(overall_stats):
    return {p: overall_stats[p.name].value_counts().index.tolist()
            for p in po.get_multiple_quant_run_options()}


def _draw_overall_stats_graphs(
        logger, plot_format, stats_dir, overall_stats,
        stats_option_values, tpm_level):

    logger.info("Drawing graphs derived from statistics calculated for the " +
                "whole set of {tpm_level} TPMs...".format(tpm_level=tpm_level))
    plot.draw_overall_stats_graphs(
        plot_format, stats_dir, overall_stats, stats_option_values, tpm_level)


def _draw_usage_graphs(
        logger, plot_format, stats_dir, usage_quant, usage_prequant,
        stats_option_values):

    logger.info("Draw graphs of time and memory resource usage...")
    plot.draw_quant_res_usage_graphs(
        plot_format, stats_dir, usage_quant, stats_option_values)
    plot.draw_prequant_usage_barplot(
        plot_format, stats_dir, usage_prequant)


def _draw_grouped_stats_graphs(
        logger, plot_format, stats_dir, grouped_threshold, stats_option_values):

    logger.info("Drawing graphs derived from statistics calculated on " +
                "subsets of transcript TPMs...")
    plot.draw_grouped_stats_graphs(
        plot_format, stats_dir, stats_option_values, grouped_threshold)


def _draw_distribution_graphs(
        logger, plot_format, stats_dir, stats_option_values):

    logger.info("Drawing distribution plots...")
    plot.draw_distribution_graphs(
        plot_format, stats_dir, stats_option_values)


def _analyse_runs(logger, record_usage, options):
    _write_accumulated_stats_and_usage(options)

    overall_transcript_stats = _get_overall_stats(options, tpms.TRANSCRIPT)
    overall_gene_stats = _get_overall_stats(options, tpms.GENE)

    stats_option_values = _get_stats_option_values(overall_transcript_stats)

    plot_format = options[po.PLOT_FORMAT.name]
    stats_dir = options[po.STATS_DIRECTORY.name]

    _draw_overall_stats_graphs(
        logger, plot_format, stats_dir, overall_transcript_stats,
        stats_option_values, tpms.TRANSCRIPT)
    _draw_overall_stats_graphs(
        logger, plot_format, stats_dir, overall_gene_stats,
        stats_option_values, tpms.GENE)
    _draw_grouped_stats_graphs(
        logger, plot_format, stats_dir, options[po.GROUPED_THRESHOLD.name],
        stats_option_values)
    _draw_distribution_graphs(
        logger, plot_format, stats_dir, stats_option_values)

    if record_usage:
        usage_prequant = _get_overall_usage(options, ru.PREQUANT_RESOURCE_TYPE)
        usage_quant = _get_overall_usage(options, ru.QUANT_RESOURCE_TYPE)
        _draw_usage_graphs(
            logger, plot_format, stats_dir, usage_quant, usage_prequant,
            stats_option_values)


def _run_piquant_command(logger, piquant_command, options, qr_options):
    record_usage = (po.NO_USAGE.name not in options) or \
        (not options[po.NO_USAGE.name])
    _set_executables_for_commands(record_usage)

    po.execute_for_mqr_option_sets(piquant_command, logger, options, qr_options)

    if piquant_command == pc.ANALYSE_RUNS:
        _analyse_runs(logger, record_usage, options)


def piquant(args):
    # Check that a piquant command has been specified
    if len(args) < 2:
        exit(("Exiting - one of the following piquant commands must be " +
              "specified:\n    {cs}\nTo get help on a particular command, " +
              "run that command with the --help option.\n\nFor further " +
              "documentation on piquant, please see " +
              "http://piquant.readthedocs.org/en/latest/.").format(
             cs="\n    ".join(pc.get_command_names())))

    # Check the specified command is valid
    command_name = args[0]
    command = pc.get_command(command_name)

    # Read command-line options
    usage = command.get_usage_message()
    options = docopt.docopt(
        usage, argv=args, version="piquant v" + __version__)

    # Set up logger
    opt.validate_log_level(options)
    logger = opt.get_logger_for_options(options)

    # Validate and process command-line options
    qr_options = None
    try:
        options, qr_options = \
            po.validate_options(logger, command, options)
    except schema.SchemaError as exc:
        exit("Exiting. " + exc.code)

    # Run the specified piquant command
    _run_piquant_command(logger, command, options, qr_options)
