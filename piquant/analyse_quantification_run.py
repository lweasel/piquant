"""
Usage:
    analyse_quantification_run [{log_option_spec} --plot-format=<plot-format> --grouped-threshold=<grouped-threshold> --error-fraction-threshold=<ef-threshold> --not-present-cutoff=<cutoff> --prequant-usage-file=<prequant-usage-file> --quant-usage-file=<quant-usage-file>] --quant-method=<quant-method> --read-length=<read-length> --read-depth=<read-depth> --paired-end=<paired-end> --errors=<errors> --bias=<bias> --stranded=<stranded> --noise-perc=<noise-depth-percentage> <tpm-file> <out-file>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
--plot-format=<plot-format>
    Output format for graphs (one of {plot_formats}) [default: pdf].
--grouped-threshold=<grouped-threshold>
    Minimum number of data points required for a group of transcripts to be
    shown on a plot [default: 300].
--error-fraction-threshold=<ef-threshold>
    Transcripts whose estimated TPM is greater than this percentage higher or
    lower than their real TPM are considered above threshold for the "error
    fraction" statistic [default: 10].
--not-present-cutoff=<cutoff>
    Cut-off value for the number of transcripts per-million below which a
    transcript is considered to be "not present" [default: 0.1].
--prequant-usage-file=<prequant-usage-file>
    CSV file recording time and memory usage of prequantification steps.
--quant-usage-file=<quant-usage-file>
    CSV file recording time and memory usage of quantification steps.
--quant-method=<quant-method>
    Method used to quantify transcript abundances.
--read-length=<read-length>
    The length of sequence reads.
--read-depth=<read-depth>
    The depth of reads sequenced across the transcriptome.
--paired-end=<paired-end>
    Whether paired-end sequence reads were used.
--errors=<errors>
    Whether the reads contain sequencing errors.
--bias=<bias>
    Whether the reads contain sequence bias.
--stranded=<stranded>
    Whether the reads are strand-specific or unstranded.
--noise-perc=<noise-depth-percentage>
    The depth of reads sequenced across a set of transcripts defined as "noise"
    (as a percentage of the depth of reads for the set of transcripts to be
    quantified).
<tpm-file>
    CSV file containing real and calculated TPMs.
<out-file>
    Basename for output graph and data files.

analyse_quantification_run reads a tpms.csv file produced by the
assemble_quantification_data script, then calculates statistics and plots
graphs to assess the accuracy of transcript abundance estimates produced in a
single quantification run.
"""

import collections
import docopt
import itertools
import numpy as np
import os.path
import pandas as pd
import schema

from . import classifiers
from . import options as opt
from . import piquant_options as po
from . import resource_usage as ru
from . import statistics
from . import tpms as t
from . import plot
from .__init__ import __version__

TRANSCRIPT_COUNT_LABEL = "No. transcripts per gene"
TRUE_POSITIVES_LABEL = "true positive TPMs"
TPM_FILE = "<tpm-file>"
PREQUANT_USAGE_FILE = "--prequant-usage-file"
QUANT_USAGE_FILE = "--quant-usage-file"
OUT_FILE_BASENAME = "<out-file>"

PIQUANT_OPTIONS = [
    po.PLOT_FORMAT,
    po.GROUPED_THRESHOLD,
    po.ERROR_FRACTION_THRESHOLD,
    po.NOT_PRESENT_CUTOFF
]

TpmInfo = collections.namedtuple("TpmInfo", ["tpms", "label"])


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)
        opt.validate_file_option(
            options[TPM_FILE], "Could not open TPM file")

        for option in PIQUANT_OPTIONS:
            opt_name = option.get_option_name()
            options[option.name] = option.validator()(options[opt_name])
            del options[opt_name]
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_tpm_infos(non_zero, tp_tpms):
    return [TpmInfo(non_zero, "non-zero real TPMs"),
            TpmInfo(tp_tpms, TRUE_POSITIVES_LABEL)]


def _add_mqr_option_values(stats, options):
    for option in po.get_multiple_quant_run_options():
        stats[option.name] = options[option.get_option_name()]


def _write_overall_stats(transcript_tpms, tp_transcript_tpms,
                         gene_tpms, tp_gene_tpms, options):

    for tpms, tp_tpms, tpm_level in \
            [(transcript_tpms, tp_transcript_tpms, t.TRANSCRIPT),
             (gene_tpms, tp_gene_tpms, t.GENE)]:
        stats = t.get_stats(tpms, tp_tpms, statistics.get_statistics())
        _add_mqr_option_values(stats, options)

        stats_file_name = statistics.get_stats_file(
            ".", options[OUT_FILE_BASENAME], tpm_level)
        statistics.write_stats_data(stats_file_name, stats, index=False)


def _write_stratified_stats(tpms, tp_tpms, non_zero, options):
    clsfr_stats = {}

    for classifier in classifiers.get_classifiers():
        if classifier.produces_grouped_stats():
            column_name = classifier.get_column_name()
            stats = t.get_grouped_stats(
                tpms, tp_tpms, column_name, statistics.get_statistics())
            _add_mqr_option_values(stats, options)
            clsfr_stats[classifier] = stats

            stats_file_name = statistics.get_stats_file(
                ".", options[OUT_FILE_BASENAME], t.TRANSCRIPT, classifier)
            statistics.write_stats_data(stats_file_name, stats)

        elif classifier.produces_distribution_plots():
            for ascending in [True, False]:
                stats = t.get_distribution_stats(
                    non_zero, tp_tpms, classifier, ascending)
                _add_mqr_option_values(stats, options)

                stats_file_name = statistics.get_stats_file(
                    ".", options[OUT_FILE_BASENAME], t.TRANSCRIPT,
                    classifier, ascending)
                statistics.write_stats_data(
                    stats_file_name, stats, index=False)

    return clsfr_stats


def _draw_tpm_scatter_plots(tp_transcript_tpms, tp_gene_tpms, plot_format,
                            basename, not_present_cutoff):

    plot.log_tpm_scatter_plot(
        plot_format, tp_transcript_tpms, basename + "_transcript",
        TRUE_POSITIVES_LABEL, not_present_cutoff)
    plot.log_tpm_scatter_plot(
        plot_format, tp_gene_tpms, basename + "_gene",
        TRUE_POSITIVES_LABEL, not_present_cutoff)


def _draw_log_ratio_boxplots(non_zero, tp_tpms, options):
    tpm_infos = _get_tpm_infos(non_zero, tp_tpms)
    clsfrs = [c for c in classifiers.get_classifiers()
              if c.produces_grouped_stats()]

    for clsfr, tpm_info in itertools.product(clsfrs, tpm_infos):
        plot.log_ratio_boxplot(
            options[po.PLOT_FORMAT.name], tpm_info.tpms,
            options[OUT_FILE_BASENAME], tpm_info.label, clsfr,
            options[po.GROUPED_THRESHOLD.name])


def _draw_stats_vs_clsfr_plots(clsfr_stats, options):
    for classifier, stats in clsfr_stats.items():
        stats.reset_index(level=0, inplace=True)
        for statistic in statistics.get_graphable_statistics():
            plot.plot_statistic_vs_classifier(
                options[po.PLOT_FORMAT.name], stats,
                options[OUT_FILE_BASENAME], statistic, classifier,
                options[po.GROUPED_THRESHOLD.name])


def _draw_cumulative_dist_plots(
        tp_tpms, non_zero, options):

    tpm_infos = _get_tpm_infos(non_zero, tp_tpms)
    clsfrs = [c for c in classifiers.get_classifiers()
              if c.produces_distribution_plots()]
    ascending = [True, False]

    for clsfr, asc, tpm_info in itertools.product(clsfrs, ascending, tpm_infos):
        plot.plot_transcript_cumul_dist(
            options[po.PLOT_FORMAT.name], tpm_info.tpms,
            options[OUT_FILE_BASENAME], tpm_info.label, clsfr, asc)


def _prepare_data(transcript_tpms, gene_tpms, not_present_cutoff, logger):
    # Determine whether each TPM measurement is a true/false positive/negative.
    # For our purposes, marking an TPM as positive or negative is determined by
    # whether it is greater or less than an "isoform not present" cutoff value.
    logger.info("Marking positives and negatives...")
    t.mark_positives_and_negatives(
        not_present_cutoff, transcript_tpms, gene_tpms)

    # Calculate the percent error (positive or negative) of the calculated TPM
    # values from the real values
    logger.info("Calculating TPM percentage errors...")
    t.calculate_percent_error(transcript_tpms, gene_tpms)

    # Calculate log ratio of calculated and real TPMs
    logger.info("Calculating TPM log ratios...")
    t.calculate_log_ratios(transcript_tpms, gene_tpms)

    # Apply various classification measures to the transcript TPM data
    logger.info("Applying classifiers...")
    t.apply_classifiers(transcript_tpms, classifiers.get_classifiers())


def _write_statistics(
        options, logger, transcript_tpms, tp_transcript_tpms,
        non_zero_transcript_tpms, gene_tpms, tp_gene_tpms):

    for stat in statistics.get_statistics():
        stat.set_options(options)

    # Write statistics pertaining to the set of quantified transcripts and
    # genes as a whole.
    logger.info("Writing overall statistics...")
    _write_overall_stats(transcript_tpms, tp_transcript_tpms,
                         gene_tpms, tp_gene_tpms, options)

    # Write statistics for TPMS stratified by various classification measures
    logger.info("Writing statistics for stratified TPMs")
    return _write_stratified_stats(
        transcript_tpms, tp_transcript_tpms, non_zero_transcript_tpms, options)


def _draw_graphs(options, tp_transcript_tpms, non_zero_transcript_tpms,
                 tp_gene_tpms, clsfr_stats):

    # Make a scatter plot of log transformed calculated vs real TPMs
    _draw_tpm_scatter_plots(
        tp_transcript_tpms, tp_gene_tpms,
        options[po.PLOT_FORMAT.name],
        options[OUT_FILE_BASENAME],
        options[po.NOT_PRESENT_CUTOFF.name])

    # Make boxplots of log ratios stratified by various classification measures
    # (e.g. the number of transcripts per-originating gene of each transcript)
    _draw_log_ratio_boxplots(
        non_zero_transcript_tpms, tp_transcript_tpms, options)

    # Make plots of statistics calculated on groups of transcripts stratified
    # by classification measures
    _draw_stats_vs_clsfr_plots(clsfr_stats, options)

    # Make plots showing the percentage of isoforms above or below threshold
    # values according to various classification measures
    _draw_cumulative_dist_plots(
        tp_transcript_tpms, non_zero_transcript_tpms, options)


def _analyse_run(logger, options):
    # Read TPMs into a data frame
    logger.info("Reading TPMs from " + options[TPM_FILE])

    transcript_tpms = pd.read_csv(options[TPM_FILE])
    gene_tpms = transcript_tpms[[t.GENE, t.REAL_TPM, t.CALCULATED_TPM]].\
        groupby(t.GENE).aggregate(np.sum)

    _prepare_data(transcript_tpms, gene_tpms,
                  options[po.NOT_PRESENT_CUTOFF.name], logger)

    # Get data frames containing only true positive TPMs and TPMs with non-zero
    # real and estimated abundances
    logger.info("Getting filtered TPMs...")
    tp_transcript_tpms = t.get_true_positives(transcript_tpms)
    non_zero_transcript_tpms = t.get_non_zero_tpms(transcript_tpms)
    tp_gene_tpms = t.get_true_positives(gene_tpms)

    clsfr_stats = _write_statistics(
        options, logger, transcript_tpms, tp_transcript_tpms,
        non_zero_transcript_tpms, gene_tpms, tp_gene_tpms)

    # Draw graphs
    logger.info("Plotting graphs...")
    _draw_graphs(options, tp_transcript_tpms, non_zero_transcript_tpms,
                 tp_gene_tpms, clsfr_stats)


def _summarise_resource_usage(logger, usage_file, resource_type, options):
    logger.info("Reading timing and memory usage from " + usage_file)

    # The resource usage file may not have been specified (e.g. if the user is
    # not interested in resource usage statistics), or it may not exist even if
    # specified (i.e. for prequantification, which is executed in only one
    # quantification directory per quantifier) - in either case, no further
    # action needs to be taken.
    if not (usage_file and os.path.exists(usage_file)):
        return

    # Read timing and memory usage info, then calculate sums of times over
    # multiple command invocations, and maximum memory usage by any command
    usage_summary = ru.get_usage_summary(usage_file)

    # Add run identification data and write resource usage summary to file
    _add_mqr_option_values(usage_summary, options)

    usage_file_name = ru.get_resource_usage_file(
        resource_type, prefix=options[OUT_FILE_BASENAME], directory=".")
    ru.write_usage_summary(usage_file_name, usage_summary)


def _analyse_resource_usage(logger, options):
    _summarise_resource_usage(
        logger, options[PREQUANT_USAGE_FILE],
        ru.PREQUANT_RESOURCE_TYPE, options)
    _summarise_resource_usage(
        logger, options[QUANT_USAGE_FILE],
        ru.QUANT_RESOURCE_TYPE, options)


def analyse_quantification_run(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(
        __doc__, plot_formats=po.PLOT_FORMATS)
    options = docopt.docopt(
        docstring, argv=args,
        version="analyse_quantification_run v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Write statistics and graphs for the quantification run
    _analyse_run(logger, options)
    _analyse_resource_usage(logger, options)
