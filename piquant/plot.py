import classifiers
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import parameters
import seaborn as sb
import statistics
import tpms as t

NO_FILTER_LABEL = "no filter"
ORDER_VALUES = [True, False]


class _NewPlot:
    def __init__(self, *file_name_elements):
        self.file_name = "_".join([str(el) for el in file_name_elements])
        self.file_name = self.file_name.replace(' ', '_')

    def __enter__(self):
        plt.figure()

    def __exit__(self, type, value, traceback):
        plt.savefig(self.file_name + ".pdf", format="pdf")
        plt.close()


def _capitalized(text):
    return text[:1].upper() + text[1:]


def _get_group_param_values(stats_df, group_param):
    group_param_vals = stats_df[group_param.name].value_counts().index.tolist()
    group_param_vals.sort()
    return group_param_vals


def _get_fixed_param_info(fixed_param_values):
    return [k.get_value_name(v) for k, v in fixed_param_values.items()]


def _get_stats_plot_name_elements(
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


def _get_grouped_stats_plot_title(
        plotted, group_param, fixed_param_info, versus=None):

    title_elements = [plotted]
    if versus:
        title_elements += ["vs", versus.lower()]
    title_elements += ["per", group_param.title.lower()]

    return " ".join(title_elements) + ": " + ", ".join(fixed_param_info)


def _plot_statistic_for_grouped_param_values(
        stats_df, ycol, group_param, xcol):

    group_param_vals = _get_group_param_values(stats_df, group_param)

    for group_param_value in group_param_vals:
        group_stats = stats_df[stats_df[group_param.name] == group_param_value]
        group_stats.sort(columns=xcol, axis=0, inplace=True)
        xvals = group_stats[xcol]
        yvals = group_stats[ycol]
        plt.plot(xvals, yvals, '-o',
                 label=group_param.get_value_name(group_param_value))


def _plot_statistic(
        stats_df, statistic, group_param, xcol, xlabel, fixed_param_info):

    _plot_statistic_for_grouped_param_values(
        stats_df, statistic.name, group_param, xcol)

    plt.xlabel(xlabel)
    plt.ylabel(statistic.title)
    plt.legend(title=group_param.title, loc=4)

    stat_range = statistic.stat_range
    if stat_range is not None:
        min_val = stat_range[0]
        max_val = stat_range[1]
        if min_val is not None:
            plt.ylim(ymin=min_val)
        if max_val is not None:
            plt.ylim(ymax=max_val)

    title = _get_grouped_stats_plot_title(
        statistic.title, group_param, fixed_param_info, versus=xlabel)
    plt.suptitle(title)


def _set_distribution_plot_bounds(xmin, xmax):
    xmargin = (xmax - xmin) / 40.0
    plt.xlim(xmin=xmin-xmargin, xmax=xmax+xmargin)
    plt.ylim(ymin=-2.5, ymax=102.5)


def _set_distribution_plot_labels(clsfr_col, ascending):
    plt.xlabel(_capitalized(clsfr_col))
    plt.ylabel("Percentage of isoforms " +
               ("less" if ascending else "greater") + " than threshold")


def log_tpm_scatter_plot(tpms, base_name, quant_method, tpm_label):
    with _NewPlot(base_name, tpm_label, "log10 scatter"):
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


def log_ratio_boxplot(tpms, base_name, quant_method, tpm_label, classifier,
                      filter=None, save_to_file=True):

    grouping_column = classifier.get_column_name()
    name_elements = [tpm_label]

    if filter:
        grouped_tpms = tpms.groupby(grouping_column)
        tpms = grouped_tpms.filter(filter)
    else:
        name_elements.append(NO_FILTER_LABEL)

    title_elements = [base_name, grouping_column] + name_elements + ["boxplot"]
    with _NewPlot(*title_elements):
        sb.boxplot(tpms[t.LOG10_RATIO], groupby=tpms[grouping_column],
                   sym='', color='lightblue')

        plt.suptitle("Log ratios of calculated to real TPMs: " +
                     ", ".join([quant_method] + name_elements))

        plt.xlabel(_capitalized(grouping_column))
        plt.ylabel("Log ratio (calculated/real TPM)")

        locs, labels = plt.xticks()
        plt.xticks(locs, classifier.get_value_labels(len(labels)))


def plot_statistic_by_parameter_values(
        stats, base_name, statistic, group_param, varying_param,
        fixed_param_values):

    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_stats_plot_name_elements(
        base_name, statistic.name,
        group_param, fixed_param_info, versus=varying_param.name)

    with _NewPlot(*name_elements):
        _plot_statistic(
            stats, statistic, group_param, varying_param.name,
            varying_param.title, fixed_param_info)


def plot_statistic_by_transcript_classifier_values(
        stats, base_name, statistic, group_param,
        classifier, fixed_param_values):

    clsfr_col = classifier.get_column_name()
    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_stats_plot_name_elements(
        base_name, statistic.name, group_param, fixed_param_info,
        versus=clsfr_col)

    with _NewPlot(*name_elements):
        xlabel = _capitalized(clsfr_col)
        _plot_statistic(
            stats, statistic, group_param, clsfr_col, xlabel, fixed_param_info)

        min_xval = stats[clsfr_col].min()
        max_xval = stats[clsfr_col].max()
        plt.xlim(xmin=min_xval, xmax=max_xval)

        tick_range = np.arange(min_xval, max_xval + 1)
        plt.xticks(np.arange(min_xval, max_xval + 1),
                   classifier.get_value_labels(len(tick_range)))


def plot_cumulative_transcript_distribution(
        tpms, base_name, quant_method, tpm_label, classifier, ascending):

    clsfr_col = classifier.get_column_name()

    with _NewPlot(base_name, clsfr_col, tpm_label,
                 ("asc" if ascending else "desc"), "distribution"):
        xvals, yvals = t.get_distribution(tpms, classifier, ascending)
        plt.plot(xvals, yvals, '-o')

        _set_distribution_plot_bounds(xvals[0], xvals[-1])
        _set_distribution_plot_labels(clsfr_col, ascending)

        plt.suptitle(_capitalized(clsfr_col) + " threshold: " + quant_method +
                     ", " + tpm_label)


def plot_cumulative_transcript_distribution_grouped(
        stats, base_name, group_param,
        classifier, ascending, fixed_param_values):

    clsfr_col = classifier.get_column_name()
    fixed_param_info = _get_fixed_param_info(fixed_param_values)

    name_elements = _get_stats_plot_name_elements(
        base_name, clsfr_col, group_param, fixed_param_info,
        ascending=ascending)

    with _NewPlot(*name_elements):
        _plot_statistic_for_grouped_param_values(
            stats, t.TRUE_POSITIVE_PERCENTAGE, group_param, clsfr_col)

        _set_distribution_plot_bounds(
            stats[clsfr_col].min(), stats[clsfr_col].max())
        _set_distribution_plot_labels(clsfr_col, ascending)

        plt.legend(title=group_param.title, loc=4)

        title = _get_grouped_stats_plot_title(
            _capitalized(clsfr_col) + " threshold", group_param,
            fixed_param_info)
        plt.suptitle(title)


# Utility functions for manipulating sets of parameters


def degenerate_param(param, param_values):
    return len(param_values[param]) <= 1


def get_non_degenerate_params(params, param_values):
    return [p for p in params if not degenerate_param(p, param_values)]


def remove_from(params, to_remove):
    get_pset = lambda x: x if isinstance(x, set) \
        else (set(x) if isinstance(x, list) else set([x]))
    return get_pset(params) - get_pset(to_remove)


def get_fixed_params(all_params, non_fixed, param_values):
    fixed_params = [p for p in remove_from(all_params, non_fixed)
                    if not degenerate_param(p, param_values)]
    value_sets = [v for v in itertools.product(
                  *[param_values[p] for p in fixed_params])]
    return fixed_params, value_sets


def get_stats_for_fixed_params(stats_df, fixed_params, fp_values_set):
    fixed_param_values = {}
    for i, fp in enumerate(fixed_params):
        fp_value = fp_values_set[i]
        stats_df = stats_df[stats_df[fp.name] == fp_value]
        fixed_param_values[fp] = fp_value
    return stats_df, fixed_param_values


def draw_overall_stats_graphs(stats_dir, overall_stats, param_values):
    # Draw graphs derived from statistics calculated for the whole set of TPMs.
    # e.g. the Spearman correlation of calculated and real TPMs graphed as
    # read-depth varies, for each quantification method, in the case of
    # paired-end reads with errors and bias.
    graph_file_basename = os.path.join(
        stats_dir, statistics.OVERALL_STATS_PREFIX)

    numerical_params = \
        [p for p in parameters.get_run_parameters() if p.is_numeric]

    for param in get_non_degenerate_params(
            parameters.get_run_parameters(), param_values):
        for num_p in get_non_degenerate_params(
                remove_from(numerical_params, param), param_values):
            fixed_params, fp_values_sets = \
                get_fixed_params(parameters.get_run_parameters(),
                                 [param, num_p], param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = get_stats_for_fixed_params(
                    overall_stats, fixed_params, fp_values_set)

                for stat in statistics.get_graphable_statistics():
                    plot_statistic_by_parameter_values(
                        stats_df, graph_file_basename,
                        stat, param, num_p, fixed_param_values)


def draw_grouped_stats_graphs(stats_dir, param_values):
    # Draw graphs derived from statistics calculated on groups of TPMs that
    # have been stratified into sets based on some classifier of transcripts.
    # e.g. the median percentage error of calculated vs real TPMs graphed as
    # the percentage of unique sequence per-transcript varies, for single and
    # paired-end reads, in the case of reads with errors and bias, and a
    # particular quantification method.
    graph_file_basename = os.path.join(
        stats_dir, statistics.OVERALL_STATS_PREFIX)

    more_than_100_filter = lambda x: x[statistics.NUM_TPMS] > 100

    clsfrs = classifiers.get_classifiers()
    grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]

    for clsfr in grp_clsfrs:
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, clsfr)
        clsfr_stats = pd.read_csv(stats_file)

        for param in get_non_degenerate_params(
                parameters.get_run_parameters(), param_values):
            fixed_params, fp_values_sets = \
                get_fixed_params(parameters.get_run_parameters(),
                                 param, param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = get_stats_for_fixed_params(
                    clsfr_stats, fixed_params, fp_values_set)

                for stat in \
                        statistics.get_graphable_by_classifier_statistics():
                    filtered_stats_df = \
                        stats_df[more_than_100_filter(stats_df)]
                    plot_statistic_by_transcript_classifier_values(
                        filtered_stats_df, graph_file_basename, stat, param,
                        clsfr, fixed_param_values)


def draw_distribution_graphs(stats_dir, param_values):
    # Draw distributions illustrating the percentage of TPMs above or below
    # some threshold as that threshold changes. e.g. the percentage of TPMs
    # whose absolute percentage error in calculated TPM, as compared to real
    # TPM, is below a particular threshold.
    graph_file_basename = os.path.join(
        stats_dir, statistics.OVERALL_STATS_PREFIX)

    clsfrs = classifiers.get_classifiers()
    dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]
    for clsfr, asc in itertools.product(dist_clsfrs, ORDER_VALUES):
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, clsfr, asc)
        clsfr_stats = pd.read_csv(stats_file)

        for param in get_non_degenerate_params(
                parameters.get_run_parameters(), param_values):
            fixed_params, fp_values_sets = \
                get_fixed_params(parameters.get_run_parameters(),
                                 param, param_values)

            for fp_values_set in fp_values_sets:
                stats_df, fixed_param_values = get_stats_for_fixed_params(
                    clsfr_stats, fixed_params, fp_values_set)

                plot_cumulative_transcript_distribution_grouped(
                    stats_df, graph_file_basename, param,
                    clsfr, asc, fixed_param_values)
