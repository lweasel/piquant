#!/usr/bin/env python

"""Usage:
    analyse_quantification_run [{log_option_spec} --scatter-max=<scatter-max-val> --log10-scatter-min=<log10-scatter-min-val> --log10-scatter-max=<log10-scatter-max-val> --plot-format=<plot-format> --grouped-threshold=<threshold>] --quant-method=<quant-method> --read-length=<read-length> --read-depth=<read-depth> --paired-end=<paired-end> --error=<errors> --bias=<bias> <tpm-file> <out-file>

-h --help                                    Show this message.
-v --version                                 Show version.
{log_option_description}
--log-level=<log-level>                      Set logging level (one of {log_level_vals}) [default: info].
--scatter-max=<scatter-max-val>              Maximum x and y values for scatter plot; a value of 0 means do not impose a maximum [default: 0].
--log10-scatter-min=<log10-scatter-min-val>  Minimum x and y values for log10 scatter plot; a value of 0 means do not impose a minimum [default: 0].
--log10-scatter-max=<log10-scatter-max-val>  Maximum x and y values for log10 scatter plot; a value of 0 means do not impose a maximum [default: 0].
--plot-format=<plot-format>                  Output format for graphs (one of {plot_formats}) [default: pdf].
--grouped-threshold=<threshold>              Minimum number of data points required for a group of transcripts to be shown on a plot [default: 300].
--quant-method=<quant-method>                Method used to quantify transcript abundances.
--read-length=<read-length>                  The length of sequence reads.
--read-depth=<read-depth>                    The depth of reads sequenced across the transcriptome.
--paired-end=<paired-end>                    Whether paired-end sequence reads were used.
--error=<errors>                             Whether the reads contain sequencing errors.
--bias=<bias>                                Whether the reads contain sequence bias.
<tpm-file>                                   File containing real and calculated TPMs.
<out-file>                                   Basename for output graph and data files.
"""

import classifiers
import collections
import docopt
import itertools
import options as opt
import pandas as pd
import parameters
import statistics
import tpms as t
import plot
import schema

TRANSCRIPT_COUNT_LABEL = "No. transcripts per gene"
TRUE_POSITIVES_LABEL = "true positive TPMs"
TPM_FILE = "<tpm-file>"
OUT_FILE_BASENAME = "<out-file>"
SCATTER_MAX = "--scatter-max"
LOG10_SCATTER_MIN = "--log10-scatter-min"
LOG10_SCATTER_MAX = "--log10-scatter-max"
PLOT_FORMAT = "--plot-format"
GROUPED_THRESHOLD = "--grouped-threshold"

TpmInfo = collections.namedtuple("TpmInfo", ["tpms", "label"])


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_list_option(
            options[PLOT_FORMAT], plot.PLOT_FORMATS, "Invalid plot format")

        opt.validate_file_option(options[TPM_FILE], "Could not open TPM file")

        options[SCATTER_MAX] = opt.validate_int_option(
            options[SCATTER_MAX],
            "Invalid maximum value for scatter plot axes")
        options[LOG10_SCATTER_MIN] = opt.validate_int_option(
            options[LOG10_SCATTER_MIN],
            "Invalid minimum value for log10 scatter plot axes")
        options[LOG10_SCATTER_MAX] = opt.validate_int_option(
            options[LOG10_SCATTER_MAX],
            "Invalid maximum value for log10 scatter plot axes")
        options[GROUPED_THRESHOLD] = opt.validate_int_option(
            options[GROUPED_THRESHOLD],
            "Invalid minimum value for number of data points for boxplots")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_non_zero_tpms(tpms):
    return tpms[(tpms[t.REAL_TPM] > 0) & (tpms[t.CALCULATED_TPM] > 0)]


def _get_tpm_infos(non_zero, tp_tpms):
    return [TpmInfo(non_zero, "non-zero real TPMs"),
            TpmInfo(tp_tpms, TRUE_POSITIVES_LABEL)]


def _add_parameter_values_to_stats(stats):
    for param in parameters.get_run_parameters():
        stats[param.name] = options[param.option_name]


def _write_overall_stats(tpms, tp_tpms, options):
    stats = t.get_stats(tpms, tp_tpms, statistics.get_statistics())
    _add_parameter_values_to_stats(stats)

    stats_file_name = statistics.get_stats_file(
        ".", options[OUT_FILE_BASENAME])
    statistics.write_stats_data(stats_file_name, stats, index=False)


def _write_stratified_stats(tpms, tp_tpms, non_zero, options):
    clsfr_stats = {}

    for classifier in classifiers.get_classifiers():
        if classifier.produces_grouped_stats():
            column_name = classifier.get_column_name()
            stats = t.get_grouped_stats(
                tpms, tp_tpms, column_name, statistics.get_statistics())
            _add_parameter_values_to_stats(stats)
            clsfr_stats[classifier] = stats

            stats_file_name = statistics.get_stats_file(
                ".", options[OUT_FILE_BASENAME], classifier)
            statistics.write_stats_data(stats_file_name, stats)

        elif classifier.produces_distribution_plots():
            for ascending in [True, False]:
                stats = t.get_distribution_stats(
                    non_zero, tp_tpms, classifier, ascending)
                _add_parameter_values_to_stats(stats)

                stats_file_name = statistics.get_stats_file(
                    ".", options[OUT_FILE_BASENAME], classifier, ascending)
                statistics.write_stats_data(
                    stats_file_name, stats, index=False)

    return clsfr_stats


def _draw_tpm_scatter_plot(tp_tpms, options):
    plot.log_tpm_scatter_plot(
        options[PLOT_FORMAT], tp_tpms,
        options[OUT_FILE_BASENAME], TRUE_POSITIVES_LABEL)


def _draw_stratified_log_ratio_boxplots(non_zero, tp_tpms, options):
    tpm_infos = _get_tpm_infos(non_zero, tp_tpms)
    clsfrs = [c for c in classifiers.get_classifiers()
              if c.produces_grouped_stats()]

    for c, ti in itertools.product(clsfrs, tpm_infos):
        plot.log_ratio_boxplot(
            options[PLOT_FORMAT], ti.tpms,
            options[OUT_FILE_BASENAME], ti.label, c,
            options[GROUPED_THRESHOLD])


def _draw_stats_vs_transcript_classifier_graphs(clsfr_stats, options):
    for classifier, stats in clsfr_stats.items():
        stats.reset_index(level=0, inplace=True)
        for statistic in statistics.get_graphable_by_classifier_statistics():
            plot.plot_statistic_vs_transcript_classifier(
                options[PLOT_FORMAT], stats,
                options[OUT_FILE_BASENAME], statistic, classifier,
                options[GROUPED_THRESHOLD])


def _draw_cumulative_transcript_distribution_graphs(
        tp_tpms, non_zero, options):

    tpm_infos = _get_tpm_infos(non_zero, tp_tpms)
    clsfrs = [c for c in classifiers.get_classifiers()
              if c.produces_distribution_plots()]
    ascending = [True, False]

    for c, asc, ti in itertools.product(clsfrs, ascending, tpm_infos):
        plot.plot_cumulative_transcript_distribution(
            options[PLOT_FORMAT], ti.tpms,
            options[OUT_FILE_BASENAME], ti.label, c, asc)


def _prepare_data(tpms):
    # Determine whether each TPM measurement is a true/false positive/negative.
    # For our purposes, marking an TPM as positive or negative is determined by
    # whether it is greater or less than an "isoform not present" cutoff value.
    logger.info("Marking positives and negatives...")
    t.mark_positives_and_negatives(tpms)

    # Calculate the percent error (positive or negative) of the calculated TPM
    # values from the real values
    logger.info("Calculating TPM percentage errors...")
    t.calculate_percent_error(tpms)

    # Calculate log ratio of calculated and real TPMs
    logger.info("Calculating TPM log ratios...")
    t.calculate_log_ratios(tpms)

    # Apply various classification measures to the TPM data
    logger.info("Applying classifiers...")
    t.apply_classifiers(tpms, classifiers.get_classifiers())


def _write_statistics(options, tpms, tp_tpms, non_zero):
    # Write statistics pertaining to the set of quantified transcripts as a
    # whole.
    logger.info("Writing overall statistics...")
    _write_overall_stats(tpms, tp_tpms, options)

    # Write statistics for TPMS stratified by various classification measures
    logger.info("Writing statistics for stratified TPMs")
    return _write_stratified_stats(tpms, tp_tpms, non_zero, options)


def _draw_graphs(options, tp_tpms, non_zero, clsfr_stats):
    # Make a scatter plot of log transformed calculated vs real TPMs
    _draw_tpm_scatter_plot(tp_tpms, options)

    # Make boxplots of log ratios stratified by various classification measures
    # (e.g. the number of transcripts per-originating gene of each transcript)
    _draw_stratified_log_ratio_boxplots(non_zero, tp_tpms, options)

    # Make plots of statistics calculated on groups of transcripts stratified
    # by classification measures
    _draw_stats_vs_transcript_classifier_graphs(clsfr_stats, options)

    # Make plots showing the percentage of isoforms above or below threshold
    # values according to various classification measures
    _draw_cumulative_transcript_distribution_graphs(tp_tpms, non_zero, options)


def _analyse_run(logger, options):
    # Read TPMs into a data frame
    logger.info("Reading TPMs...")
    tpms = pd.read_csv(options[TPM_FILE])

    _prepare_data(tpms)

    # Get data frames containing only true positive TPMs and TPMs with non-zero
    # real and estimated abundances
    logger.info("Getting filtered TPMs...")
    tp_tpms = t.get_true_positives(tpms)
    non_zero = _get_non_zero_tpms(tpms)

    clsfr_stats = _write_statistics(options, tpms, tp_tpms, non_zero)

    # Draw graphs
    logger.info("Plotting graphs...")
    _draw_graphs(options, tp_tpms, non_zero, clsfr_stats)


if __name__ == "__main__":
    # Read in command-line options
    __doc__ = opt.substitute_into_usage(__doc__).format(
        plot_formats=plot.PLOT_FORMATS)
    options = docopt.docopt(
        __doc__, version="analyse_quantification_run v0.1")

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Write statistics and graphs for the quantification run
    _analyse_run(logger, options)
