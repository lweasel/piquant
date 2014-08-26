import classifiers
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import parameters
import seaborn as sb
import statistics
import sys
import tpms as t

NO_FILTER_LABEL = "no filter"
GROUPED_STATS_NUM_TPMS_FILTER = 3000
ORDER_VALUES = [True, False]
PLOT_FORMATS = ["pdf", "svg", "png"]

# Don't embed characters as paths when outputting SVG - assume fonts are
# installed on machine where SVG will be viewed (see
# http://matplotlib.org/users/customizing.html)
plt.rcParams['svg.fonttype'] = 'none'

# matplotlib parameters appropriate for poster output
#plt.rcParams['font.size'] = 16.0
#plt.rcParams['axes.labelsize'] = 'medium'
#plt.rcParams['xtick.labelsize'] = 'x-small'
#plt.rcParams['ytick.labelsize'] = 'x-small'
#plt.rcParams['legend.fontsize'] = 'small'


class _NewPlot:
    def __init__(self, fformat, *file_name_elements):
        self.fformat = fformat
        self.file_name = "_".join([str(el) for el in file_name_elements])
        self.file_name = self.file_name.replace(' ', '_')

    def __enter__(self):
        plt.figure()

    def __exit__(self, type, value, traceback):
        plt.savefig(self.file_name + "." + self.fformat, format=self.fformat)
        plt.close()


def _capitalized(text):
    return text[:1].upper() + text[1:]


def _decapitalized(text):
    return text[:1].lower() + text[1:]


def _get_distribution_plot_ylabel(ascending):
    return "Percentage of isoforms " + \
        ("less" if ascending else "greater") + " than threshold"


def _set_distribution_plot_bounds(xmin, xmax, ymin=None, ymax=None):
    xmargin = (xmax - xmin) / 40.0
    plt.xlim(xmin=xmin-xmargin, xmax=xmax+xmargin)
    plt.ylim(ymin=-2.5, ymax=102.5)


def _get_statistic_plot_bounds_setter(statistic):
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


def _set_ticks_for_transcript_classifier_plot(locations, classifier):
    plt.xticks(locations, classifier.get_value_labels(len(locations)))


def _get_group_param_values(stats_df, group_param):
    group_param_vals = stats_df[group_param.name].value_counts().index.tolist()
    group_param_vals.sort()
    return group_param_vals


def _get_fixed_param_info(fixed_param_values):
    return [k.get_value_name(v) for k, v in fixed_param_values.items()]


def _get_grouped_by_param_stats_plot_file_name_elements(
        base_name, plotted, group_param, fixed_param_info,
        versus=None, ascending=None):
    name_elements = [base_name, plotted]
    if versus:
        name_elements += ["vs", versus]
    name_elements += ["per", group_param.title.lower()]
    if ascending is not None:
        name_elements.append("asc" if ascending else "desc")
    name_elements += fixed_param_info
    if ascending is not None:
        name_elements.append("distribution")
    return name_elements


def _get_stats_plot_file_name_elements(
        base_name, plotted, versus=None, ascending=None):
    name_elements = [base_name, plotted]
    if versus:
        name_elements += ["vs", versus]
    if ascending is not None:
        name_elements += ["asc" if ascending else "desc", "distribution"]
    return name_elements


def _get_grouped_by_param_stats_plot_title(
        plotted, group_param, fixed_param_info, versus=None):

    title_elements = [plotted]
    if versus:
        title_elements += ["vs", _decapitalized(versus)]
    title_elements += ["per", group_param.title.lower()]

    title = " ".join(title_elements)
    if len(fixed_param_info) > 0:
        title += ": " + ", ".join(fixed_param_info)

    return title


def _plot_statistic_grouped_by_parameter(
        stats_df, group_param, xcol, ycol, xlabel, ylabel,
        plot_bounds_setter, title):

    group_param_vals = _get_group_param_values(stats_df, group_param)

    xmin = ymin = sys.maxsize
    xmax = ymax = -sys.maxsize - 1

    for group_param_value in group_param_vals:
        group_stats = stats_df[stats_df[group_param.name] == group_param_value]
        group_stats.sort(columns=xcol, axis=0, inplace=True)
        xvals = group_stats[xcol]
        yvals = group_stats[ycol]
        plt.plot(xvals, yvals, '-o',
                 label=group_param.get_value_name(group_param_value))

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
    plt.legend(title=group_param.title, loc=4)
    plt.suptitle(title)

    return (ymin, ymax)


def _plot_statistic_vs_varying_param_grouped_by_param(
        fformat, stats, base_name, statistic,
        group_param, varying_param, fixed_param_values):

    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_grouped_by_param_stats_plot_file_name_elements(
        base_name, statistic.name,
        group_param, fixed_param_info, versus=varying_param.name)

    with _NewPlot(fformat, *name_elements):
        title = _get_grouped_by_param_stats_plot_title(
            statistic.title, group_param, fixed_param_info,
            versus=varying_param.title)
        _plot_statistic_grouped_by_parameter(
            stats, group_param, varying_param.name, statistic.name,
            varying_param.title, statistic.title,
            _get_statistic_plot_bounds_setter(statistic), title)


def _plot_statistic_vs_transcript_classifier_grouped_by_param(
        fformat, stats, base_name, statistic, group_param,
        classifier, fixed_param_values):

    clsfr_col = classifier.get_column_name()
    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_grouped_by_param_stats_plot_file_name_elements(
        base_name, statistic.name, group_param, fixed_param_info,
        versus=clsfr_col)

    with _NewPlot(fformat, *name_elements):
        xlabel = _capitalized(classifier.get_plot_title())
        title = _get_grouped_by_param_stats_plot_title(
            statistic.title, group_param, fixed_param_info, versus=xlabel)

        _plot_statistic_grouped_by_parameter(
            stats, group_param, clsfr_col, statistic.name,
            xlabel, statistic.title,
            _get_statistic_plot_bounds_setter(statistic), title)

        min_xval = stats[clsfr_col].min()
        max_xval = stats[clsfr_col].max()
        _set_ticks_for_transcript_classifier_plot(
            np.arange(min_xval, max_xval + 1), classifier)


def _plot_cumulative_transcript_distribution_grouped_by_param(
        fformat, stats, base_name, group_param,
        classifier, ascending, fixed_param_values):

    clsfr_col = classifier.get_column_name()
    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_grouped_by_param_stats_plot_file_name_elements(
        base_name, clsfr_col, group_param, fixed_param_info,
        ascending=ascending)

    with _NewPlot(fformat, *name_elements):
        title = _get_grouped_by_param_stats_plot_title(
            _capitalized(clsfr_col) + " threshold", group_param,
            fixed_param_info)
        _plot_statistic_grouped_by_parameter(
            stats, group_param, clsfr_col, t.TRUE_POSITIVE_PERCENTAGE,
            _capitalized(clsfr_col),
            _get_distribution_plot_ylabel(ascending),
            _set_distribution_plot_bounds, title)


def log_tpm_scatter_plot(fformat, tpms, base_name, quant_method, tpm_label):
    with _NewPlot(fformat, base_name, tpm_label, "log10 scatter"):
        plt.scatter(tpms[t.LOG10_REAL_TPM].values,
                    tpms[t.LOG10_CALCULATED_TPM].values,
                    c="lightblue", alpha=0.4)

        plt.suptitle("Scatter plot of log calculated vs real TPMs: " +
                     quant_method + ", " + tpm_label)
        plt.xlabel("Log10 real TPM")
        plt.ylabel("Log10 calculated TPM")

        min_val = np.log10(t.NOT_PRESENT_CUTOFF) - 0.2
        plt.xlim(xmin=min_val)
        plt.ylim(ymin=min_val)


def log_ratio_boxplot(
        fformat, tpms, base_name, quant_method, tpm_label, classifier,
        filter=None, save_to_file=True):

    grouping_column = classifier.get_column_name()
    name_elements = [tpm_label]

    if filter:
        grouped_tpms = tpms.groupby(grouping_column)
        tpms = grouped_tpms.filter(filter)
    else:
        name_elements.append(NO_FILTER_LABEL)

    title_elements = [base_name, grouping_column] + name_elements + ["boxplot"]
    with _NewPlot(fformat, *title_elements):
        sb.boxplot(tpms[t.LOG10_RATIO], groupby=tpms[grouping_column],
                   sym='', color='lightblue')

        plt.suptitle("Log ratios of calculated to real TPMs: " +
                     ", ".join([quant_method] + name_elements))

        plt.xlabel(_capitalized(grouping_column))
        plt.ylabel("Log ratio (calculated/real TPM)")

        _set_ticks_for_transcript_classifier_plot(plt.xticks()[0], classifier)


def plot_statistic_vs_transcript_classifier(
        fformat, stats, base_name, statistic, classifier):

    stats = stats[stats[statistics.NUM_TPMS] > GROUPED_STATS_NUM_TPMS_FILTER]

    clsfr_col = classifier.get_column_name()

    with _NewPlot(fformat, base_name, statistic.name, "vs", clsfr_col):
        xvals = stats[clsfr_col]
        min_xval = xvals.min()
        max_xval = xvals.max()
        yvals = stats[statistic.name]

        plt.plot(xvals, yvals, '-o')

        _get_statistic_plot_bounds_setter(statistic)(
            min_xval, max_xval, yvals.min(), yvals.max())

        plt.xlabel(_capitalized(classifier.get_plot_title()))
        plt.ylabel(statistic.title)
        plt.suptitle(statistic.title + " vs " + _decapitalized(clsfr_col))

        _set_ticks_for_transcript_classifier_plot(
            np.arange(min_xval, max_xval + 1), classifier)


def plot_cumulative_transcript_distribution(
        fformat, tpms, base_name, quant_method,
        tpm_label, classifier, ascending):

    clsfr_col = classifier.get_column_name()

    with _NewPlot(fformat, base_name, clsfr_col, tpm_label,
                 ("asc" if ascending else "desc"), "distribution"):
        xvals, yvals = t.get_distribution(tpms, classifier, ascending)
        plt.plot(xvals, yvals, '-o')

        _set_distribution_plot_bounds(xvals[0], xvals[-1])

        plt.xlabel(_capitalized(clsfr_col))
        plt.ylabel(_get_distribution_plot_ylabel(ascending))

        plt.suptitle(_capitalized(clsfr_col) + " threshold: " + quant_method +
                     ", " + tpm_label)


# Utility functions for manipulating sets of parameters


def _degenerate_param(param, param_values):
    return len(param_values[param]) <= 1


def _get_non_degenerate_params(params, param_values):
    return [p for p in params if not _degenerate_param(p, param_values)]


def _remove_from(params, to_remove):
    get_pset = lambda x: x if isinstance(x, set) \
        else (set(x) if isinstance(x, list) else set([x]))
    return get_pset(params) - get_pset(to_remove)


def _get_fixed_params(all_params, non_fixed, param_values):
    fixed_params = [p for p in _remove_from(all_params, non_fixed)
                    if not _degenerate_param(p, param_values)]
    value_sets = [v for v in itertools.product(
                  *[param_values[p] for p in fixed_params])]
    return fixed_params, value_sets


def _get_stats_for_fixed_params(stats_df, fixed_params, fp_values_set):
    fixed_param_values = {}
    for i, fp in enumerate(fixed_params):
        fp_value = fp_values_set[i]
        stats_df = stats_df[stats_df[fp.name] == fp_value]
        fixed_param_values[fp] = fp_value
    return stats_df, fixed_param_values


# Making plots over multiple sets of sequencing and quantification parameters


def _get_plot_subdirectory(parent_dir, sub_dir_name):
    sub_dir = os.path.join(parent_dir, sub_dir_name)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)
    return sub_dir


def draw_overall_stats_graphs(fformat, stats_dir, overall_stats, param_values):
    # Draw graphs derived from statistics calculated for the whole set of TPMs.
    # e.g. the Spearman correlation of calculated and real TPMs graphed as
    # read-depth varies, for each quantification method, in the case of
    # paired-end reads with errors and bias.
    overall_stats_dir = _get_plot_subdirectory(
        stats_dir, "overall_stats_graphs")

    numerical_params = \
        [p for p in parameters.get_run_parameters() if p.is_numeric]

    for param in _get_non_degenerate_params(
            parameters.get_run_parameters(), param_values):

        non_degenerate_numerical_params = _get_non_degenerate_params(
            _remove_from(numerical_params, param), param_values)
        if len(non_degenerate_numerical_params) == 0:
            continue

        param_stats_dir = _get_plot_subdirectory(
            overall_stats_dir, "per_" + param.name)

        for num_p in non_degenerate_numerical_params:
            num_param_stats_dir = _get_plot_subdirectory(
                param_stats_dir, "by_" + num_p.name)

            fixed_params, fp_values_sets = \
                _get_fixed_params(parameters.get_run_parameters(),
                                  [param, num_p], param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = _get_stats_for_fixed_params(
                    overall_stats, fixed_params, fp_values_set)

                for stat in statistics.get_graphable_statistics():
                    statistic_dir = _get_plot_subdirectory(
                        num_param_stats_dir, stat.name)

                    graph_file_basename = os.path.join(
                        statistic_dir, statistics.OVERALL_STATS_PREFIX)
                    _plot_statistic_vs_varying_param_grouped_by_param(
                        fformat, stats_df, graph_file_basename,
                        stat, param, num_p, fixed_param_values)


def draw_grouped_stats_graphs(fformat, stats_dir, param_values):
    # Draw graphs derived from statistics calculated on groups of TPMs that
    # have been stratified into sets based on some classifier of transcripts.
    # e.g. the median percentage error of calculated vs real TPMs graphed as
    # the percentage of unique sequence per-transcript varies, for single and
    # paired-end reads, in the case of reads with errors and bias, and a
    # particular quantification method.
    grouped_stats_dir = _get_plot_subdirectory(
        stats_dir, "grouped_stats_graphs")

    num_tpms_filter = \
        lambda x: x[statistics.NUM_TPMS] > GROUPED_STATS_NUM_TPMS_FILTER

    clsfrs = classifiers.get_classifiers()
    grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]

    for clsfr in grp_clsfrs:
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, clsfr)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdirectory(
            grouped_stats_dir,
            "grouped_by_" + clsfr.get_column_name().replace(' ', '_'))

        for param in _get_non_degenerate_params(
                parameters.get_run_parameters(), param_values):
            param_stats_dir = _get_plot_subdirectory(
                clsfr_dir, "per_" + param.name)

            fixed_params, fp_values_sets = \
                _get_fixed_params(parameters.get_run_parameters(),
                                  param, param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = _get_stats_for_fixed_params(
                    clsfr_stats, fixed_params, fp_values_set)

                for stat in \
                        statistics.get_graphable_by_classifier_statistics():
                    statistic_dir = _get_plot_subdirectory(
                        param_stats_dir, stat.name)
                    graph_file_basename = os.path.join(
                        statistic_dir, "grouped")

                    filtered_stats_df = stats_df[num_tpms_filter(stats_df)]
                    _plot_statistic_vs_transcript_classifier_grouped_by_param(
                        fformat, filtered_stats_df, graph_file_basename,
                        stat, param, clsfr, fixed_param_values)


def draw_distribution_graphs(fformat, stats_dir, param_values):
    # Draw distributions illustrating the percentage of TPMs above or below
    # some threshold as that threshold changes. e.g. the percentage of TPMs
    # whose absolute percentage error in calculated TPM, as compared to real
    # TPM, is below a particular threshold.
    distribution_stats_dir = _get_plot_subdirectory(
        stats_dir, "distribution_stats_graphs")

    clsfrs = classifiers.get_classifiers()
    dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]
    for clsfr, asc in itertools.product(dist_clsfrs, ORDER_VALUES):
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, clsfr, asc)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdirectory(
            distribution_stats_dir,
            clsfr.get_column_name().replace(' ', '_') + "_distribution")

        for param in _get_non_degenerate_params(
                parameters.get_run_parameters(), param_values):

            param_stats_dir = _get_plot_subdirectory(
                clsfr_dir, "per_" + param.name)
            graph_file_basename = os.path.join(
                param_stats_dir, "distribution")

            fixed_params, fp_values_sets = \
                _get_fixed_params(parameters.get_run_parameters(),
                                  param, param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = _get_stats_for_fixed_params(
                    clsfr_stats, fixed_params, fp_values_set)

                _plot_cumulative_transcript_distribution_grouped_by_param(
                    fformat, stats_df, graph_file_basename, param,
                    clsfr, asc, fixed_param_values)
