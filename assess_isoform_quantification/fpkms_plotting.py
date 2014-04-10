import fpkms as f
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb

from collections import namedtuple

NO_FILTER_LABEL = "no filter"
CUMULATIVE_DISTRIBUTION_POINTS = 20


class NewPlot:
    def __enter__(self):
        plt.figure()

    def __exit__(self, type, value, traceback):
        plt.close()


PlotOptions = namedtuple(
    'PlotOptions', ['quant_method', 'label', 'out_file_base'])


def _savefig(base_name, name_elements, suffix=None):
    elements = [base_name] + name_elements
    if suffix:
        elements += [suffix]
    name = "_".join(elements) + ".pdf"
    name = name.replace(' ', '_')
    plt.savefig(name, format="pdf")


def log_fpkm_scatter_plot(fpkms, plot_options):
    with NewPlot():
        plt.scatter(fpkms[f.LOG10_REAL_FPKM].values,
                    fpkms[f.LOG10_CALCULATED_FPKM].values,
                    c="lightblue", alpha=0.4)

        plt.suptitle("Scatter plot of log calculated vs real FPKMs: " +
                     ", ".join([plot_options.quant_method,
                               plot_options.label]))
        plt.xlabel("Log10 real FPKM")
        plt.ylabel("Log10 calculated FPKM")

        min_val = np.log10(f.NOT_PRESENT_CUTOFF) - 0.2
        plt.xlim(xmin=min_val)
        plt.ylim(ymin=min_val)

        _savefig(plot_options.out_file_base,
                 [plot_options.label],
                 suffix="log10 scatter")


def log_ratio_boxplot(fpkms, plot_options, classifier, filter=None):

    grouping_column = classifier.get_column_name()
    name_elements = [plot_options.label]
    if filter:
        grouped_fpkms = fpkms.groupby(grouping_column)
        fpkms = grouped_fpkms.filter(filter)
    else:
        name_elements.append(NO_FILTER_LABEL)

    with NewPlot():
        sb.boxplot(fpkms[f.LOG10_RATIO], groupby=fpkms[grouping_column],
                   sym='', color='lightblue')

        plt.suptitle("Log ratios of calculated to real FPKMs: " +
                     ", ".join([plot_options.quant_method] + name_elements))

        plt.xlabel(grouping_column[:1].upper() + grouping_column[1:])
        plt.ylabel("Log ratio (calculated/real FPKM)")

        locs, labels = plt.xticks()
        plt.xticks(locs, classifier.get_value_labels(len(labels)))

        _savefig(plot_options.out_file_base,
                 [grouping_column] + name_elements,
                 suffix="boxplot")


def plot_cumulative_transcript_distribution(
        fpkms, plot_options, classifier, ascending):

    values = fpkms.apply(classifier.get_value, axis=1)
    values.sort(ascending=ascending)

    xbounds = classifier.get_distribution_plot_range()
    if xbounds is None:
        xbounds = (values.min(), values.max())

    xvals = np.linspace(xbounds[0], xbounds[1], CUMULATIVE_DISTRIBUTION_POINTS)

    size = float(len(values))
    yvals = [100 * len(values[values < x if ascending else values > x]) / size
             for x in xvals]

    with NewPlot():
        plt.plot(xvals, yvals, '-o')

        plt.ylim(ymin=-2.5, ymax=102.5)

        xmargin = (xbounds[1] - xbounds[0]) / 40.0
        plt.xlim(xmin=xbounds[0]-xmargin, xmax=xbounds[1]+xmargin)

        grouping_column = classifier.get_column_name()
        capitalized = grouping_column[:1].upper() + grouping_column[1:]
        plt.xlabel(capitalized)
        plt.ylabel("Percentage of isoforms " +
                   ("less" if ascending else "greater") + " than threshold")

        plt.suptitle(capitalized + " threshold: " +
                     ", ".join([plot_options.quant_method,
                               plot_options.label]))

        _savefig(plot_options.out_file_base,
                 [grouping_column, plot_options.label,
                  ("asc" if ascending else "desc")],
                 suffix="distribution")


def plot_statistic(stats, plot_options, statistic,
                   group_param, varying_param, fixed_param_values):

    group_param_values = stats[group_param.name].value_counts().index.tolist()
    group_param_values.sort()

    with NewPlot():

        for group_param_value in group_param_values:
            group_stats = stats[stats[group_param.name] == group_param_value]
            group_stats.sort(columns=varying_param.name, axis=0, inplace=True)
            xvals = group_stats[varying_param.name]
            yvals = group_stats[statistic.name]
            plt.plot(xvals, yvals, '-o',
                     label=group_param.get_value_name(group_param_value))

        plt.xlabel(varying_param.title)
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

        fixed_param_info = [k.get_value_name(v) for k, v
                            in fixed_param_values.items()]

        title = " ".join([statistic.title, "vs",
                          varying_param.title.lower(), "per",
                          group_param.title.lower()]) \
                + ": " + ", ".join(fixed_param_info)

        plt.suptitle(title)

        _savefig(plot_options.out_file_base,
                 [statistic.name, "vs",
                  varying_param.name, "per",
                  group_param.title.lower()]
                 + fixed_param_info)
