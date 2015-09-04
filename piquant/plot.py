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

RESOURCE_USAGE_DIR = "resource_usage_graphs"

# Don't embed characters as paths when outputting SVG - assume fonts are
# installed on machine where SVG will be viewed (see
# http://matplotlib.org/users/customizing.html)
plt.rcParams['svg.fonttype'] = 'none'

# This is the only way I can get fliers to appear in Seaborn boxplots
plt.rcParams['lines.markeredgewidth'] = 0.5

# matplotlib parameters appropriate for poster output
# plt.rcParams['font.size'] = 16.0
# plt.rcParams['axes.labelsize'] = 'medium'
# plt.rcParams['xtick.labelsize'] = 'x-small'
# plt.rcParams['ytick.labelsize'] = 'x-small'
# plt.rcParams['legend.fontsize'] = 'small'

# matplotlib parameters appropriate for paper
#plt.rcParams['font.size'] = 16.0
#plt.rcParams['axes.labelsize'] = 'large'
#plt.rcParams['xtick.labelsize'] = 'medium'
#plt.rcParams['ytick.labelsize'] = 'medium'
#plt.rcParams['legend.fontsize'] = 'medium'


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


def _set_plot_ybounds(statistic, ymin, ymax):
    stat_range = statistic.stat_range((ymin, ymax))
    if stat_range is not None:
        min_val = stat_range[0]
        max_val = stat_range[1]
        if min_val is not None:
            plt.ylim(ymin=min_val)
        if max_val is not None:
            plt.ylim(ymax=max_val)


def _get_plot_bounds_setter(statistic):
    def _set_statistic_plot_bounds(xmin, xmax, ymin, ymax):
        xmargin = 2 * (xmax - xmin) / 100.0
        plt.xlim(xmin=xmin - xmargin, xmax=xmax + xmargin)

        _set_plot_ybounds(statistic, ymin, ymax)

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

    plt.xlabel(_capitalized(xlabel))
    plt.ylabel(_capitalized(ylabel))
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
            varying_mqr_option.get_axis_label(), statistic.get_axis_label(),
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
        plot_info.set_plot_title(
            statistic.title, versus=classifier.get_plot_title())

        _plot_grouped_statistic(
            stats, plot_info, clsfr_col, statistic.name,
            classifier.get_axis_label(), statistic.get_axis_label(),
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
            stats, plot_info, clsfr_col, t.NON_ZERO_PERCENTAGE,
            clsfr_col, _get_distribution_plot_ylabel(ascending),
            _set_distribution_plot_bounds)


def _draw_prequant_time_usage_graph(fformat, graph_file_basename, usage_data):
    with _saving_new_plot(fformat, [graph_file_basename, "time_usage"]):
        n_groups = len(usage_data.index)
        index = np.arange(n_groups)

        time_usage_stats = ru.get_time_usage_statistics()
        gap_width = 0.1
        bar_width = (1 - gap_width) / len(time_usage_stats)

        dummy, axes = plt.subplots()
        color_cycle = axes._get_lines.color_cycle

        for i, usage_stat in enumerate(time_usage_stats):
            plt.bar(index + i * bar_width,
                    usage_data[usage_stat.name].values,
                    bar_width, color=color_cycle.next(),
                    label=_capitalized(usage_stat.name.replace('-', ' ')))

        plt.xlabel('Quantification method')
        plt.ylabel('Log10 total time (s)')
        plt.title('Time taken for prequantification')
        plt.xticks(index + ((1 - gap_width) / 2),
                   usage_data["quant_method"].values)

        box = axes.get_position()
        axes.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        axes.legend(loc=6, bbox_to_anchor=(1, 0.5))


def _draw_prequant_mem_usage_graph(fformat, graph_file_basename, usage_data):
    with _saving_new_plot(fformat, [graph_file_basename, "memory_usage"]):
        n_groups = len(usage_data.index)
        index = np.arange(n_groups)

        mem_usage_stat = ru.get_memory_usage_statistics()[0]
        bar_width = 0.9

        dummy, axes = plt.subplots()
        color_cycle = axes._get_lines.color_cycle

        plt.bar(index,
                usage_data[mem_usage_stat.name].values,
                bar_width, color=color_cycle.next())

        plt.xlabel('Quantification method')
        plt.ylabel('Resident memory (Gb)')
        plt.title('Maximum resident memory during prequantification')
        plt.xticks(index + bar_width / 2,
                   usage_data["quant_method"].values)


def log_tpm_scatter_plot(
        fformat, tpms, base_name, tpm_label, not_present_cutoff):

    with _saving_new_plot(fformat, [base_name, tpm_label, "log10 scatter"]):
        plt.scatter(tpms[t.LOG10_REAL_TPM].values,
                    tpms[t.LOG10_CALCULATED_TPM].values,
                    c="lightblue", alpha=0.4)

        plt.suptitle("Scatter plot of log calculated vs real TPMs: " +
                     tpm_label)
        plt.xlabel("$log_{10}$ real TPM")
        plt.ylabel("$log_{10}$ calculated TPM")

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
        sb.boxplot(x=tpms[grouping_column],
                   y=tpms[t.LOG10_RATIO],
                   color='lightblue', sym='')

        plt.suptitle("Log ratios of calculated to real TPMs: " + tpm_label)
        plt.xlabel(_capitalized(grouping_column))
        plt.ylabel("Log ratio (calculated/real TPM)")

        _set_ticks_for_classifier_plot(plt.xticks()[0], classifier)


def plot_statistic_vs_classifier(
        fformat, stats, base_name, statistic, classifier, threshold):

    stats = stats[stats[statistics.EXPRESSED_TPMS] > threshold]
    clsfr_col = classifier.get_column_name()

    with _saving_new_plot(fformat, [base_name, statistic.name, "vs", clsfr_col]):
        xvals = stats[clsfr_col]
        min_xval = xvals.min()
        max_xval = xvals.max()
        yvals = stats[statistic.name]

        plt.plot(xvals, yvals, '-o')

        _get_plot_bounds_setter(statistic)(
            min_xval, max_xval, yvals.min(), yvals.max())

        plt.xlabel(_capitalized(classifier.get_axis_label()))
        plt.ylabel(statistic.title)
        plt.suptitle(statistic.title + " vs " + classifier.get_plot_title())

        _set_ticks_for_classifier_plot(
            np.arange(min_xval, max_xval + 1), classifier)


def _plot_statistic_boxplot(opt_dir, plot_file_prefix, fformat,
                            data, statistic, main_opt, sub_opt=None):

    output_dir = _get_plot_subdir(opt_dir, statistic.name) \
        if sub_opt else opt_dir
    graph_file_basename = os.path.join(output_dir, plot_file_prefix)

    options = [main_opt]
    if sub_opt:
        options = [sub_opt] + options

    flatten = lambda x: [item for sublist in x for item in sublist]

    name_elements = [graph_file_basename, statistic.name] + \
        flatten([["per", opt.name] for opt in options])

    with _saving_new_plot(fformat, name_elements):
        data = data.sort([opt.name for opt in options])

        for opt in options:
            data[opt.name] = data[opt.name].apply(opt.get_value_name)

        sb.boxplot(x=data[main_opt.name],
                   y=data[statistic.name],
                   hue=data[sub_opt.name] if sub_opt else None,
                   showmeans=True,
                   meanprops=dict(marker='D', markerfacecolor='none'))

        ymin = data[statistic.name].min()
        ymax = data[statistic.name].max()

        _set_plot_ybounds(statistic, ymin, ymax)

        plt.xlabel(main_opt.get_axis_label())
        plt.ylabel(statistic.get_axis_label())

        if sub_opt:
            plt.legend(title=sub_opt.title, loc=4)

        title_elements = [statistic.title, "distribution"] + \
            flatten([["per", _decapitalized(opt.title)] for opt in options])
        title = " ".join(title_elements)
        plt.suptitle(title)


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


# Making plots over multiple sets of sequencing and quantification run options


def _get_plot_subdir(parent_dir, *sub_dir_name_elems):
    sub_dir_name = "_".join(sub_dir_name_elems).replace(' ', '_')
    sub_dir = os.path.join(parent_dir, sub_dir_name)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)
    return sub_dir


def stats_graphs_vs_num_opt_drawer(
        plot_dir, fformat, grp_option, num_option, stats, plot_file_prefix):

    num_opt_dir = _get_plot_subdir(plot_dir, "by", num_option.name)

    def drawer(df, fixed_option_values):
        for stat in stats:
            stat_dir = _get_plot_subdir(num_opt_dir, stat.name)
            graph_file_basename = os.path.join(stat_dir, plot_file_prefix)

            _plot_grouped_stat_vs_mqr_opt(
                fformat, df, graph_file_basename, stat, grp_option,
                num_option, fixed_option_values)

    return drawer


def _draw_stats_graphs(
        fformat, stats_dir, sub_dir, data_frame, opt_vals_set,
        stats, plot_file_prefix):

    main_plot_dir = _get_plot_subdir(stats_dir, sub_dir)

    for option in opt_vals_set.get_non_degenerate_options():
        nondeg_numerical_opts = opt_vals_set.get_non_degenerate_options(
            numeric_only=True, opts_to_remove=option)
        if len(nondeg_numerical_opts) == 0:
            continue

        option_plot_dir = _get_plot_subdir(main_plot_dir, "per", option.name)

        for num_opt in nondeg_numerical_opts:
            opt_vals_set.exec_for_fixed_option_values_sets(
                stats_graphs_vs_num_opt_drawer(
                    option_plot_dir, fformat, option, num_opt,
                    stats, plot_file_prefix),
                [option, num_opt], data_frame)

        for stat in stats:
            _plot_statistic_boxplot(option_plot_dir, plot_file_prefix, fformat,
                                    data_frame, stat, option)

        all_nondeg_opts = opt_vals_set.get_non_degenerate_options(
            opts_to_remove=option)

        for opt in all_nondeg_opts:
            opt_dir = _get_plot_subdir(option_plot_dir, "by", opt.name)

            for stat in stats:
                _plot_statistic_boxplot(
                    opt_dir, plot_file_prefix, fformat,
                    data_frame, stat, option, opt)


def draw_quant_res_usage_graphs(
        fformat, stats_dir, usage_data, opt_vals_set):

    _draw_stats_graphs(
        fformat, stats_dir, RESOURCE_USAGE_DIR, usage_data, opt_vals_set,
        ru.get_resource_usage_statistics(), "usage")


def draw_prequant_res_usage_graphs(fformat, stats_dir, usage_data):
    plot_dir = _get_plot_subdir(stats_dir, RESOURCE_USAGE_DIR)
    graph_file_basename = os.path.join(plot_dir, "prequant")

    _draw_prequant_time_usage_graph(fformat, graph_file_basename, usage_data)
    _draw_prequant_mem_usage_graph(fformat, graph_file_basename, usage_data)


def draw_overall_stats_graphs(
        fformat, stats_dir, overall_stats, opt_vals_set, tpm_level):

    # Draw graphs derived from statistics calculated for the whole set of TPMs.
    # e.g. the Spearman correlation of calculated and real TPMs graphed as
    # read-depth varies, for each quantification method, in the case of
    # paired-end reads with errors and bias.
    sub_dir = "overall_{l}_stats_graphs".format(l=tpm_level)
    _draw_stats_graphs(
        fformat, stats_dir, sub_dir, overall_stats, opt_vals_set,
        statistics.get_graphable_statistics(), statistics.OVERALL_STATS_PREFIX)


def grouped_stats_graph_drawer(
        plot_dir, fformat, grp_option, clsfr, num_tpms_filter):

    option_stats_dir = _get_plot_subdir(plot_dir, "per", grp_option.name)

    def drawer(df, fixed_option_values):
        for stat in statistics.get_graphable_statistics():
            statistic_dir = _get_plot_subdir(option_stats_dir, stat.name)
            graph_file_basename = os.path.join(statistic_dir, "grouped")

            filtered_stats_df = df[num_tpms_filter(df)]
            _plot_grouped_stat_vs_clsfr(
                fformat, filtered_stats_df, graph_file_basename,
                stat, grp_option, clsfr, fixed_option_values)

    return drawer


def draw_grouped_stats_graphs(fformat, stats_dir, opt_vals_set, threshold):
    # Draw graphs derived from statistics calculated on groups of TPMs that
    # have been stratified into sets based on some classifier of transcripts.
    # e.g. the median percentage error of calculated vs real TPMs graphed as
    # the percentage of unique sequence per-transcript varies, for single and
    # paired-end reads, in the case of reads with errors and bias, and a
    # particular quantification method.
    grouped_stats_dir = _get_plot_subdir(stats_dir, "grouped_stats_graphs")

    num_tpms_filter = lambda x: x[statistics.EXPRESSED_TPMS] > threshold

    clsfrs = classifiers.get_classifiers()
    grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]

    for clsfr in grp_clsfrs:
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX, t.TRANSCRIPT, clsfr)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdir(
            grouped_stats_dir, "grouped_by", clsfr.get_column_name())

        for option in opt_vals_set.get_non_degenerate_options():
            opt_vals_set.exec_for_fixed_option_values_sets(
                grouped_stats_graph_drawer(
                    clsfr_dir, fformat, option, clsfr, num_tpms_filter),
                option, clsfr_stats)


def distribution_stats_graph_drawer(
        plot_dir, fformat, grp_option, clsfr, asc):

    option_stats_dir = _get_plot_subdir(plot_dir, "per", grp_option.name)
    graph_file_basename = os.path.join(option_stats_dir, "distribution")

    def drawer(df, fixed_option_values):
        _plot_grouped_cumulative_dist(
            fformat, df, graph_file_basename, grp_option,
            clsfr, asc, fixed_option_values)

    return drawer


def draw_distribution_graphs(fformat, stats_dir, opt_vals_set):
    # Draw distributions illustrating the percentage of TPMs above or below
    # some threshold as that threshold changes. e.g. the percentage of TPMs
    # whose absolute percentage error in calculated TPM, as compared to real
    # TPM, is below a particular threshold.
    distribution_stats_dir = _get_plot_subdir(
        stats_dir, "distribution_stats_graphs")

    clsfrs = classifiers.get_classifiers()
    dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]

    for clsfr, asc in itertools.product(dist_clsfrs, [True, False]):
        stats_file = statistics.get_stats_file(
            stats_dir, statistics.OVERALL_STATS_PREFIX,
            t.TRANSCRIPT, clsfr, asc)
        clsfr_stats = pd.read_csv(stats_file)

        clsfr_dir = _get_plot_subdir(
            distribution_stats_dir, clsfr.get_column_name(), "distribution")

        for option in opt_vals_set.get_non_degenerate_options():
            opt_vals_set.exec_for_fixed_option_values_sets(
                distribution_stats_graph_drawer(
                    clsfr_dir, fformat, option, clsfr, asc),
                option, clsfr_stats)
