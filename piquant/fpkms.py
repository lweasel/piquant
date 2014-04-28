import pandas as pd
import numpy as np
import statistics as stats

TRANSCRIPT_COUNT = "num-transcripts"
LENGTH = "length"
UNIQUE_SEQ_LENGTH = "unique-length"
REAL_FPKM = "real-fpkm"
CALCULATED_FPKM = "calc-fpkm"
PERCENT_ERROR = "percent-error"
LOG10_REAL_FPKM = "log10-real-fpkm"
LOG10_CALCULATED_FPKM = "log10-calc-fpkm"
LOG10_RATIO = "log-ratio"
NON_ZERO_PERCENTAGE = "non-zero-perc"
TRUE_POSITIVE_PERCENTAGE = "tp-perc"

FALSE_POSITIVE = "false-pos"
FALSE_NEGATIVE = "false-neg"
TRUE_POSITIVE = "true-pos"
TRUE_NEGATIVE = "true-neg"

NOT_PRESENT_CUTOFF = 0.1
CUMULATIVE_DISTRIBUTION_POINTS = 20


def mark_positives_and_negatives(fpkms):
    fpkms[FALSE_NEGATIVE] = \
        (fpkms[REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
        (fpkms[CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
    fpkms[FALSE_POSITIVE] = \
        (fpkms[CALCULATED_FPKM] > NOT_PRESENT_CUTOFF) & \
        (fpkms[REAL_FPKM] <= NOT_PRESENT_CUTOFF)
    fpkms[TRUE_NEGATIVE] = \
        (fpkms[REAL_FPKM] <= NOT_PRESENT_CUTOFF) & \
        (fpkms[CALCULATED_FPKM] <= NOT_PRESENT_CUTOFF)
    fpkms[TRUE_POSITIVE] = \
        (fpkms[REAL_FPKM] > NOT_PRESENT_CUTOFF) & \
        (fpkms[CALCULATED_FPKM] > NOT_PRESENT_CUTOFF)


def get_true_positives(fpkms):
    return fpkms[fpkms[TRUE_POSITIVE]]


def calculate_percent_error(fpkms):
    fpkms[PERCENT_ERROR] = \
        100 * (fpkms[CALCULATED_FPKM] - fpkms[REAL_FPKM]) / fpkms[REAL_FPKM]


def calculate_log_ratios(fpkms):
    fpkms[LOG10_REAL_FPKM] = np.log10(fpkms[REAL_FPKM])
    fpkms[LOG10_CALCULATED_FPKM] = np.log10(fpkms[CALCULATED_FPKM])
    fpkms[LOG10_RATIO] = \
        fpkms[LOG10_CALCULATED_FPKM] - fpkms[LOG10_REAL_FPKM]


def apply_classifiers(fpkms, clsfrs):
    for classifier in clsfrs:
        column_name = classifier.get_column_name()
        fpkms[column_name] = fpkms.apply(
            classifier.get_classification_value, axis=1)


def get_stats(fpkms, tp_fpkms):
    stats_dict = {stat.name: stat.calculate(fpkms, tp_fpkms)
                  for stat in stats.get_statistics()}
    return pd.DataFrame([stats_dict])


def get_grouped_stats(fpkms, tp_fpkms, column_name):
    grouped = fpkms.groupby(column_name)
    tp_grouped = tp_fpkms.groupby(column_name)
    summary = grouped.describe()
    tp_summary = tp_grouped.describe()

    stats_dict = {stat.name: stat.calculate_grouped(
        grouped, summary, tp_grouped, tp_summary)
        for stat in stats.get_statistics()}
    return pd.DataFrame.from_dict(stats_dict)


def get_distribution(fpkms, classifier, ascending):
    values = fpkms.apply(classifier.get_value, axis=1)
    values.sort(ascending=ascending)

    xbounds = classifier.get_distribution_plot_range()
    if xbounds is None:
        xbounds = (values.min(), values.max())

    xvals = np.linspace(xbounds[0], xbounds[1], CUMULATIVE_DISTRIBUTION_POINTS)

    size = float(len(values))
    yvals = [100 * len(values[values < x if ascending else values > x]) / size
             for x in xvals]

    return xvals, yvals


def get_distribution_stats(non_zero_fpkms, tp_fpkms, classifier, ascending):
    xvals, nz_yvals = get_distribution(non_zero_fpkms, classifier, ascending)
    xvals, tp_yvals = get_distribution(tp_fpkms, classifier, ascending)

    stats_dict = {}
    stats_dict[classifier.get_column_name()] = xvals
    stats_dict[NON_ZERO_PERCENTAGE] = nz_yvals
    stats_dict[TRUE_POSITIVE_PERCENTAGE] = tp_yvals

    return pd.DataFrame.from_dict(stats_dict)
