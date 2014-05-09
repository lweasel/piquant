import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import tpms as t

NO_FILTER_LABEL = "no filter"


class NewPlot:
    def __init__(self, *file_name_elements):
        self.file_name = "_".join([str(el) for el in file_name_elements])
        self.file_name = self.file_name.replace(' ', '_')

    def __enter__(self):
        plt.figure()

    def __exit__(self, type, value, traceback):
        plt.savefig(self.file_name, format="pdf")
        plt.close()


def log_tpm_scatter_plot(tpms, base_name, quant_method, tpm_label):
    with NewPlot(base_name, tpm_label, "log10 scatter"):
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

    title_elements = [base_name] + name_elements + ["boxplot"]
    with NewPlot(*title_elements):
        sb.boxplot(tpms[t.LOG10_RATIO], groupby=tpms[grouping_column],
                   sym='', color='lightblue')

        plt.suptitle("Log ratios of calculated to real TPMs: " +
                     ", ".join([quant_method] + name_elements))

        plt.xlabel(grouping_column[:1].upper() + grouping_column[1:])
        plt.ylabel("Log ratio (calculated/real TPM)")

        locs, labels = plt.xticks()
        plt.xticks(locs, classifier.get_value_labels(len(labels)))


def plot_cumulative_transcript_distribution(
        tpms, base_name, quant_method, tpm_label, classifier, ascending):

    xvals, yvals = t.get_distribution(tpms, classifier, ascending)
    grouping_column = classifier.get_column_name()

    with NewPlot(base_name, grouping_column, tpm_label,
                 ("asc" if ascending else "desc"), "distribution"):
        plt.plot(xvals, yvals, '-o')

        plt.ylim(ymin=-2.5, ymax=102.5)

        xmin = xvals[0]
        xmax = xvals[-1]
        xmargin = (xmax - xmin) / 40.0
        plt.xlim(xmin=xmin-xmargin, xmax=xmax+xmargin)

        capitalized = grouping_column[:1].upper() + grouping_column[1:]
        plt.xlabel(capitalized)
        plt.ylabel("Percentage of isoforms " +
                   ("less" if ascending else "greater") + " than threshold")

        plt.suptitle(capitalized + " threshold: " + quant_method +
                     ", " + tpm_label)


def plot_statistic(stats, base_name, statistic,
                   group_param, varying_param, fixed_param_values):

    group_param_values = stats[group_param.name].value_counts().index.tolist()
    group_param_values.sort()

    fixed_param_info = [k.get_value_name(v) for k, v
                        in fixed_param_values.items()]

    name_elements = [base_name, statistic.name, "vs", varying_param.name,
                     "per", group_param.title.lower()]
    name_elements += fixed_param_info

    with NewPlot(*name_elements):

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

        title = " ".join([statistic.title, "vs",
                          varying_param.title.lower(), "per",
                          group_param.title.lower()]) \
                + ": " + ", ".join(fixed_param_info)

        plt.suptitle(title)


def plot_statistic_by_classifier(
        stats, base_name, statistic, group_param,
        classifier, filter, fixed_param_values):

    clsfr_col = classifier.get_column_name()

    group_param_values = stats[group_param.name].value_counts().index.tolist()
    group_param_values.sort()

    if filter:
        stats = stats[filter(stats)]

    fixed_param_info = [k.get_value_name(v) for k, v
                        in fixed_param_values.items()]

    name_elements = [base_name, statistic.name, "vs", clsfr_col, "per",
                     group_param.title.lower()]
    name_elements += fixed_param_info

    with NewPlot(*name_elements):
        for group_param_value in group_param_values:
            group_stats = stats[stats[group_param.name] == group_param_value]
            group_stats.sort(columns=clsfr_col, axis=0, inplace=True)
            xvals = group_stats[clsfr_col]
            yvals = group_stats[statistic.name]
            plt.plot(xvals, yvals, '-o',
                     label=group_param.get_value_name(group_param_value))

        plt.xlabel(clsfr_col[:1].upper() + clsfr_col[1:])
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

        min_xval = group_stats[clsfr_col].min()
        max_xval = group_stats[clsfr_col].max()
        plt.xlim(xmin=min_xval, xmax=max_xval)

        tick_range = np.arange(min_xval, max_xval + 1)
        plt.xticks(np.arange(min_xval, max_xval + 1),
                   classifier.get_value_labels(len(tick_range)))

        title = " ".join([statistic.title, "vs",
                          clsfr_col.lower(), "per",
                          group_param.title.lower()]) \
                + ": " + ", ".join(fixed_param_info)

        plt.suptitle(title)


def plot_cumulative_transcript_distribution_grouped(
        stats, base_name, group_param,
        classifier, ascending, fixed_param_values):

    clsfr_col = classifier.get_column_name()

    group_param_values = stats[group_param.name].value_counts().index.tolist()
    group_param_values.sort()

    fixed_param_info = [k.get_value_name(v) for k, v
                        in fixed_param_values.items()]

    name_elements = [base_name, clsfr_col, "per",
                     group_param.title.lower(),
                     ("asc" if ascending else "desc")]
    name_elements += fixed_param_info

    with NewPlot(*name_elements):
        for group_param_value in group_param_values:
            group_stats = stats[stats[group_param.name] == group_param_value]
            group_stats.sort(columns=clsfr_col, axis=0, inplace=True)
            xvals = group_stats[clsfr_col]
            yvals = group_stats[t.TRUE_POSITIVE_PERCENTAGE]

            plt.plot(xvals, yvals, '-o',
                     label=group_param.get_value_name(group_param_value))

        plt.ylim(ymin=-2.5, ymax=102.5)

        xmin = xvals.values[0]
        xmax = xvals.values[-1]
        xmargin = (xmax - xmin) / 40.0
        plt.xlim(xmin=xmin-xmargin, xmax=xmax+xmargin)

        capitalized = clsfr_col[:1].upper() + clsfr_col[1:]
        plt.xlabel(capitalized)
        plt.ylabel("Percentage of isoforms " +
                   ("less" if ascending else "greater") + " than threshold")
        plt.legend(title=group_param.title, loc=4)

        title = " ".join([capitalized,
                          "threshold per",
                          group_param.title.lower()]) + \
            ": " + ", ".join(fixed_param_info)

        plt.suptitle(title)
