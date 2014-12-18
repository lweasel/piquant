import contextlib
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import seaborn as sb
import sys

from . import classifiers
from . import piquant_options as po
from . import resource_usage as ru
from . import statistics
from . import tpms as t

# Don't embed characters as paths when outputting SVG - assume fonts are
# installed on machine where SVG will be viewed (see
# http://matplotlib.org/users/customizing.html)
plt.rcParams['svg.fonttype'] = 'none'

# matplotlib parameters appropriate for poster output
# plt.rcParams['font.size'] = 16.0
# plt.rcParams['axes.labelsize'] = 'medium'
# plt.rcParams['xtick.labelsize'] = 'x-small'
# plt.rcParams['ytick.labelsize'] = 'x-small'
# plt.rcParams['legend.fontsize'] = 'small'


class _GroupedPlotInfo(object):
    def __init__(self, group_mqr_option, fixed_mqr_option_info):
        self.group_mqr_option = group_mqr_option
        self.fixed_mqr_option_info = fixed_mqr_option_info
        self.title = None

    def get_filename_parts(
            self, base_name, plotted, versus=None, ascending=None):

        name_elements = [base_name, plotted]
        if versus:
            name_elements += ["vs", versus]
        name_elements += ["per", self.group_mqr_option.title.lower()]
        if ascending is not None:
            name_elements.append("asc" if ascending else "desc")
        name_elements += self.fixed_mqr_option_info
        if ascending is not None:
            name_elements.append("distribution")
        return name_elements

    def set_plot_title(self, plotted, versus=None):
        title_elements = [plotted]
        if versus:
            title_elements += ["vs", _decapitalized(versus)]
        title_elements += ["per", self.group_mqr_option.title.lower()]

        title = " ".join(title_elements)
        if len(self.fixed_mqr_option_info) > 0:
            title += ": " + ", ".join(self.fixed_mqr_option_info)

        self.title = title


@contextlib.contextmanager
def _saving_new_plot(fformat, file_name_elements):
    plt.figure()
    try:
        yield
    finally:
        file_name = "_".join([str(el) for el in file_name_elements])
        file_name = file_name.replace(' ', '_')
        plt.savefig(file_name + "." + fformat, format=fformat)
        plt.close()


def _capitalized(text):
    return text[:1].upper() + text[1:]


def _decapitalized(text):
    return text[:1].lower() + text[1:]


def _get_distribution_plot_ylabel(ascending):
    return "Percentage of isoforms " + \
        ("less" if ascending else "greater") + " than threshold"


def _set_distribution_plot_bounds(xmin, xmax, ymin=None, ymax=None):
    # unused arguments
    del ymin
    del ymax

    xmargin = (xmax - xmin) / 40.0
    plt.xlim(xmin=xmin - xmargin, xmax=xmax + xmargin)
    plt.ylim(ymin=-2.5, ymax=102.5)


def _get_plot_bounds_setter(statistic):
    def _set_statistic_plot_bounds(xmin, xmax, ymin, ymax):
        xmargin = 2 * (xmax - xmin) / 100.0
        plt.xlim(xmin=xmin - xmargin, xmax=xmax + xmargin)

        stat_range = statistic.stat_range((ymin, ymax))
        if stat_range is not None:
            min_val = stat_range[0]
            max_val = stat_range[1]
            if min_val is not None:
                plt.ylim(ymin=min_val)
            if max_val is not None:
                plt.ylim(ymax=max_val)

    return _set_statistic_plot_bounds


def _set_ticks_for_classifier_plot(locations, classifier):
    plt.xticks(locations, classifier.get_value_labels(len(locations)))


def _get_group_mqr_option_values(stats_df, group_mqr_option):
    group_mqr_option_vals = \
        stats_df[group_mqr_option.name].value_counts().index.tolist()
    group_mqr_option_vals.sort()
    return group_mqr_option_vals


def _plot_grouped_statistic(
        stats_df, plot_info, xcol, ycol, xlabel, ylabel,
        plot_bounds_setter):

    group_mqr_option_vals = _get_group_mqr_option_values(
        stats_df, plot_info.group_mqr_option)

    xmin = ymin = sys.maxsize
    xmax = ymax = -sys.maxsize - 1

    for group_mqr_option_value in group_mqr_option_vals:
        group_stats = stats_df[
            stats_df[plot_info.group_mqr_option.name] == group_mqr_option_value]
        group_stats.sort(columns=xcol, axis=0, inplace=True)
        xvals = group_stats[xcol]
        yvals = group_stats[ycol]
        plt.plot(xvals, yvals, '-o',
                 label=plot_info.group_mqr_option.get_value_name(
                     group_mqr_option_value))

        group_ymin = yvals.min()
        if group_ymin < ymin:
            ymin = group_ymin

        group_ymax = yvals.max()
        if group_ymax > ymax:
            ymax = group_ymax

        group_xmin = xvals.min()
        if group_xmin < xmin:
            xmin = group_xmin

        group_xmax = xvals.max()
        if group_xmax > xmax:
            xmax = group_xmax

    plot_bounds_setter(xmin, xmax, ymin, ymax)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title=plot_info.group_mqr_option.title, loc=4)
    plt.suptitle(plot_info.title)

    return (ymin, ymax)


def _plot_grouped_stat_vs_mqr_opt(
        fformat, stats, base_name, statistic, group_mqr_option,
        varying_mqr_option, fixed_mqr_option_values):

    fixed_mqr_option_info = po.get_value_names(fixed_mqr_option_values)

    plot_info = _GroupedPlotInfo(group_mqr_option, fixed_mqr_option_info)
    name_elements = plot_info.get_filename_parts(
        base_name, statistic.name, versus=varying_mqr_option.name)

    with _saving_new_plot(fformat, name_elements):
        plot_info.set_plot_title(
            statistic.title, versus=varying_mqr_option.title)
        _plot_grouped_statistic(
            stats, plot_info, varying_mqr_option.name, statistic.name,
            varying_mqr_option.title, statistic.title,
            _get_plot_bounds_setter(statistic))


def _plot_grouped_stat_vs_clsfr(
        fformat, stats, base_name, statistic, group_mqr_option,
        classifier, fixed_mqr_option_values):

    clsfr_col = classifier.get_column_name()
    fixed_mqr_option_info = po.get_value_names(fixed_mqr_option_values)

    plot_info = _GroupedPlotInfo(group_mqr_option, fixed_mqr_option_info)
    name_elements = plot_info.get_filename_parts(
        base_name, statistic.name, versus=clsfr_col)

    with _saving_new_plot(fformat, name_elements):
        xlabel = _capitalized(classifier.get_plot_title())
        plot_info.set_plot_title(statistic.title, versus=xlabel)

        _plot_grouped_statistic(
            stats, plot_info, clsfr_col,
            statistic.name, xlabel, statistic.title,
            _get_plot_bounds_setter(statistic))

        min_xval = stats[clsfr_col].min()
        max_xval = stats[clsfr_col].max()
        _set_ticks_for_classifier_plot(
            np.arange(min_xval, max_xval + 1), classifier)


def _plot_grouped_cumulative_dist(
        fformat, stats, base_name, group_mqr_option,
        classifier, ascending, fixed_mqr_option_values):

    clsfr_col = classifier.get_column_name()
    fixed_mqr_option_info = po.get_value_names(fixed_mqr_option_values)

    plot_info = _GroupedPlotInfo(group_mqr_option, fixed_mqr_option_info)
    name_elements = plot_info.get_filename_parts(
        base_name, clsfr_col, ascending=ascending)

    with _saving_new_plot(fformat, name_elements):
        plot_info.set_plot_title(
            _capitalized(clsfr_col) + " threshold")
        _plot_grouped_statistic(
            stats, plot_info, clsfr_col, t.TRUE_POSITIVE_PERCENTAGE,
            _capitalized(clsfr_col),
            _get_distribution_plot_ylabel(ascending),
            _set_distribution_plot_bounds)


def log_tpm_scatter_plot(
        fformat, tpms, base_name, tpm_label, not_present_cutoff):

    with _saving_new_plot(fformat, [base_name, tpm_label, "log10 scatter"]):
        plt.scatter(tpms[t.LOG10_REAL_TPM].values,
                    tpms[t.LOG10_CALCULATED_TPM].values,
                    c="lightblue", alpha=0.4)

        plt.suptitle("Scatter plot of log calculated vs real TPMs: " +
                     tpm_label)
        plt.xlabel("Log10 real TPM")
        plt.ylabel("Log10 calculated TPM")

        min_val = np.log10(not_present_cutoff) - 0.2
        plt.xlim(xmin=min_val)
        plt.ylim(ymin=min_val)


def log_ratio_boxplot(
        fformat, tpms, base_name, tpm_label, classifier, threshold):

    grouping_column = classifier.get_column_name()
    grouped_tpms = tpms.groupby(grouping_column)
    tpms = grouped_tpms.filter(
        lambda x: len(x[t.REAL_TPM]) > threshold)

    with _saving_new_plot(
            fformat, [base_name, grouping_column, tpm_label, "boxplot"]):
        sb.boxplot(tpms[t.LOG10_RATIO], groupby=tpms[grouping_column],
                   sym='', color='lightblue')

        plt.suptitle("Log ratios of calculated to real TPMs: " + tpm_label)
        plt.xlabel(_capitalized(grouping_column))
        plt.ylabel("Log ratio (calculated/real TPM)")

        _set_ticks_for_classifier_plot(plt.xticks()[0], classifier)


def plot_statistic_vs_classifier(
        fformat, stats, base_name, statistic, classifier, threshold):

    stats = stats[stats[statistics.TP_NUM_TPMS] > threshold]
    clsfr_col = classifier.get_column_name()

    with _saving_new_plot(fformat, [base_name, statistic.name, "vs", clsfr_col]):
        xvals = stats[clsfr_col]
        min_xval = xvals.min()
        max_xval = xvals.max()
        yvals = stats[statistic.name]

        plt.plot(xvals, yvals, '-o')

        _get_plot_bounds_setter(statistic)(
            min_xval, max_xval, yvals.min(), yvals.max())

        plt.xlabel(_capitalized(classifier.get_plot_title()))
        plt.ylabel(statistic.title)
        plt.suptitle(statistic.title + " vs " + _decapitalized(clsfr_col))

        _set_ticks_for_classifier_plot(
            np.arange(min_xval, max_xval + 1), classifier)


def plot_transcript_cumul_dist(
        fformat, tpms, base_name, tpm_label, classifier, ascending):

    clsfr_col = classifier.get_column_name()

    with _saving_new_plot(
            fformat, [base_name, clsfr_col, tpm_label,
                      ("asc" if ascending else "desc"), "distribution"]):

        xvals, yvals = t.get_distribution(tpms, classifier, ascending)
        plt.plot(xvals, yvals, '-o')

        _set_distribution_plot_bounds(xvals[0], xvals[-1])

        plt.xlabel(_capitalized(clsfr_col))
        plt.ylabel(_get_distribution_plot_ylabel(ascending))
        plt.suptitle(_capitalized(clsfr_col) + " threshold: " + tpm_label)


def draw_prequant_usage_barplot(
        fformat, stats_dir, usage_prequant):
    pass


# Utility functions for manipulating sets of quantification run options


def _degenerate_mqr_option(mqr_option, mqr_option_values):
    return len(mqr_option_values[mqr_option]) <= 1


def _get_non_degenerate_mqr_options(mqr_options, mqr_option_values):
    return [o for o in mqr_options
            if not _degenerate_mqr_option(o, mqr_option_values)]


def _remove_from(mqr_options, to_remove):
    get_pset = lambda x: x if isinstance(x, set) \
        else (set(x) if isinstance(x, list) else set([x]))
    return get_pset(mqr_options) - get_pset(to_remove)


def _get_fixed_mqr_opts(all_mqr_options, non_fixed, mqr_option_values):
    fixed_mqr_options = [o for o in _remove_from(all_mqr_options, non_fixed)
                         if not _degenerate_mqr_option(o, mqr_option_values)]
    value_sets = [v for v in itertools.product(
                  *[mqr_option_values[o] for o in fixed_mqr_options])]
    return fixed_mqr_options, value_sets


def _get_data_for_fixed_mqr_opts(df, fixed_mqr_options, fo_values_set):
    fixed_mqr_option_values = {}
    for i, fixed_o in enumerate(fixed_mqr_options):
        fo_value = fo_values_set[i]
        df = df[df[fixed_o.name] == fo_value]
        fixed_mqr_option_values[fixed_o] = fo_value
    return df, fixed_mqr_option_values


# Making plots over multiple sets of sequencing and quantification run options


def _get_plot_subdirectory(parent_dir, sub_dir_name):
    sub_dir = os.path.join(parent_dir, sub_dir_name)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)
    return sub_dir


def draw_overall_stats_graphs(
        fformat, stats_dir, overall_stats, mqr_option_values, tpm_level):

    # Draw graphs derived from statistics calculated for the whole set of TPMs.
    # e.g. the Spearman correlation of calculated and real TPMs graphed as
    # read-depth varies, for each quantification method, in the case of
    # paired-end reads with errors and bias.
    overall_stats_dir = _get_plot_subdirectory(
        stats_dir, "overall_{l}_stats_graphs".format(l=tpm_level))

    numerical_mqr_options = \
        [o for o in po.get_multiple_quant_run_options() if o.is_numeric]

    for mqr_option in _get_non_degenerate_mqr_options(
            po.get_multiple_quant_run_options(), mqr_option_values):

        non_deg_numerical_mqr_opts = _get_non_degenerate_mqr_options(
            _remove_from(numerical_mqr_options, mqr_option), mqr_option_values)
        if len(non_deg_numerical_mqr_opts) == 0:
            continue

        mqr_option_stats_dir = _get_plot_subdirectory(
            overall_stats_dir, "per_" + mqr_option.name)

        for num_p in non_deg_numerical_mqr_opts:
            num_mqr_option_stats_dir = _get_plot_subdirectory(
                mqr_option_stats_dir, "by_" + num_p.name)

            fixed_mqr_options, fp_values_sets = \
                _get_fixed_mqr_opts(po.get_multiple_quant_run_options(),
                                    [mqr_option, num_p], mqr_option_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_mqr_option_values = \
                    _get_data_for_fixed_mqr_opts(
                        overall_stats, fixed_mqr_options, fp_values_set)

                for stat in statistics.get_graphable_statistics():
                    statistic_dir = _get_plot_subdirectory(
                        num_mqr_option_stats_dir, stat.name)

                    graph_file_basename = os.path.join(
                        statistic_dir, statistics.OVERALL_STATS_PREFIX)
                    _plot_grouped_stat_vs_mqr_opt(
                        fformat, stats_df, graph_file_basename,
                        stat, mqr_option, num_p, fixed_mqr_option_values)


def draw_grouped_stats_graphs(fformat, stats_dir, mqr_option_values, threshold):
    # Draw graphs derived from statistics calculated on groups of TPMs that
    # have been stratified into sets based on some classifier of transcripts.
    # e.g. the median percentage error of calculated vs real TPMs graphed as
    # the percentage of unique sequence per-transcript varies, for single and
    # paired-end reads, in the case of reads with errors and bias, and a
    # particular quantification method.
    grouped_stats_dir = _get_plot_subdirectory(
        stats_dir, "grouped_stats_graphs")

    num_tpms_filter = lambda x: x[statistics.TP_NUM_TPMS] > threshold

    clsfrs = classifiers.get_classifiers()
    grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]

    for clsfr in grp_clsfrs:
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, t.TRANSCRIPT, clsfr)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdirectory(
            grouped_stats_dir,
            "grouped_by_" + clsfr.get_column_name().replace(' ', '_'))

        for mqr_option in _get_non_degenerate_mqr_options(
                po.get_multiple_quant_run_options(), mqr_option_values):
            mqr_option_stats_dir = _get_plot_subdirectory(
                clsfr_dir, "per_" + mqr_option.name)

            fixed_mqr_options, fp_values_sets = \
                _get_fixed_mqr_opts(po.get_multiple_quant_run_options(),
                                    mqr_option, mqr_option_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_mqr_option_values = \
                    _get_data_for_fixed_mqr_opts(
                        clsfr_stats, fixed_mqr_options, fp_values_set)

                for stat in statistics.get_graphable_statistics():
                    statistic_dir = _get_plot_subdirectory(
                        mqr_option_stats_dir, stat.name)
                    graph_file_basename = os.path.join(
                        statistic_dir, "grouped")

                    filtered_stats_df = stats_df[num_tpms_filter(stats_df)]
                    _plot_grouped_stat_vs_clsfr(
                        fformat, filtered_stats_df, graph_file_basename,
                        stat, mqr_option, clsfr, fixed_mqr_option_values)


def draw_distribution_graphs(fformat, stats_dir, mqr_option_values):
    # Draw distributions illustrating the percentage of TPMs above or below
    # some threshold as that threshold changes. e.g. the percentage of TPMs
    # whose absolute percentage error in calculated TPM, as compared to real
    # TPM, is below a particular threshold.
    distribution_stats_dir = _get_plot_subdirectory(
        stats_dir, "distribution_stats_graphs")

    clsfrs = classifiers.get_classifiers()
    dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]
    for clsfr, asc in itertools.product(dist_clsfrs, [True, False]):
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX,
            t.TRANSCRIPT, clsfr, asc)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdirectory(
            distribution_stats_dir,
            clsfr.get_column_name().replace(' ', '_') + "_distribution")

        for mqr_option in _get_non_degenerate_mqr_options(
                po.get_multiple_quant_run_options(), mqr_option_values):

            mqr_option_stats_dir = _get_plot_subdirectory(
                clsfr_dir, "per_" + mqr_option.name)
            graph_file_basename = os.path.join(
                mqr_option_stats_dir, "distribution")

            fixed_mqr_options, fp_values_sets = \
                _get_fixed_mqr_opts(po.get_multiple_quant_run_options(),
                                    mqr_option, mqr_option_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_mqr_option_values = \
                    _get_data_for_fixed_mqr_opts(
                        clsfr_stats, fixed_mqr_options, fp_values_set)

                _plot_grouped_cumulative_dist(
                    fformat, stats_df, graph_file_basename, mqr_option,
                    clsfr, asc, fixed_mqr_option_values)


def draw_quant_res_usage_graphs(
        fformat, stats_dir, usage_data, mqr_option_values):

    usage_graphs_dir = _get_plot_subdirectory(
        stats_dir, "resource_usage_graphs")

    numerical_mqr_options = \
        [o for o in po.get_multiple_quant_run_options() if o.is_numeric]

    for mqr_option in _get_non_degenerate_mqr_options(
            po.get_multiple_quant_run_options(), mqr_option_values):

        non_deg_numerical_mqr_opts = _get_non_degenerate_mqr_options(
            _remove_from(numerical_mqr_options, mqr_option), mqr_option_values)
        if len(non_deg_numerical_mqr_opts) == 0:
            continue

        mqr_option_usage_dir = _get_plot_subdirectory(
            usage_graphs_dir, "per_" + mqr_option.name)

        for num_p in non_deg_numerical_mqr_opts:
            num_mqr_option_usage_dir = _get_plot_subdirectory(
                mqr_option_usage_dir, "by_" + num_p.name)

            fixed_mqr_options, fp_values_sets = \
                _get_fixed_mqr_opts(po.get_multiple_quant_run_options(),
                                    [mqr_option, num_p], mqr_option_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_mqr_option_values = \
                    _get_data_for_fixed_mqr_opts(
                        usage_data, fixed_mqr_options, fp_values_set)

                for usage_stat in ru.get_resource_usage_statistics():
                    usage_dir = _get_plot_subdirectory(
                        num_mqr_option_usage_dir, usage_stat.name)

                    graph_file_basename = os.path.join(usage_dir, "usage")
                    _plot_grouped_stat_vs_mqr_opt(
                        fformat, stats_df, graph_file_basename,
                        usage_stat, mqr_option, num_p, fixed_mqr_option_values)
