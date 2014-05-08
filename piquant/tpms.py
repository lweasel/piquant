import pandas as pd
import numpy as np
import statistics as stats

TRANSCRIPT_COUNT = "num-transcripts"
LENGTH = "length"
UNIQUE_SEQ_LENGTH = "unique-length"
REAL_TPM = "real-tpm"
CALCULATED_TPM = "calc-tpm"
PERCENT_ERROR = "percent-error"
LOG10_REAL_TPM = "log10-real-tpm"
LOG10_CALCULATED_TPM = "log10-calc-tpm"
LOG10_RATIO = "log-ratio"
NON_ZERO_PERCENTAGE = "non-zero-perc"
TRUE_POSITIVE_PERCENTAGE = "tp-perc"

FALSE_POSITIVE = "false-pos"
FALSE_NEGATIVE = "false-neg"
TRUE_POSITIVE = "true-pos"
TRUE_NEGATIVE = "true-neg"

NOT_PRESENT_CUTOFF = 0.1
CUMULATIVE_DISTRIBUTION_POINTS = 20


def mark_positives_and_negatives(tpms):
    tpms[FALSE_NEGATIVE] = \
        (tpms[REAL_TPM] > NOT_PRESENT_CUTOFF) & \
        (tpms[CALCULATED_TPM] <= NOT_PRESENT_CUTOFF)
    tpms[FALSE_POSITIVE] = \
        (tpms[CALCULATED_TPM] > NOT_PRESENT_CUTOFF) & \
        (tpms[REAL_TPM] <= NOT_PRESENT_CUTOFF)
    tpms[TRUE_NEGATIVE] = \
        (tpms[REAL_TPM] <= NOT_PRESENT_CUTOFF) & \
        (tpms[CALCULATED_TPM] <= NOT_PRESENT_CUTOFF)
    tpms[TRUE_POSITIVE] = \
        (tpms[REAL_TPM] > NOT_PRESENT_CUTOFF) & \
        (tpms[CALCULATED_TPM] > NOT_PRESENT_CUTOFF)


def get_true_positives(tpms):
    return tpms[tpms[TRUE_POSITIVE]]


def calculate_percent_error(tpms):
    tpms[PERCENT_ERROR] = \
        100 * (tpms[CALCULATED_TPM] - tpms[REAL_TPM]) / tpms[REAL_TPM]


def calculate_log_ratios(tpms):
    tpms[LOG10_REAL_TPM] = np.log10(tpms[REAL_TPM])
    tpms[LOG10_CALCULATED_TPM] = np.log10(tpms[CALCULATED_TPM])
    tpms[LOG10_RATIO] = \
        tpms[LOG10_CALCULATED_TPM] - tpms[LOG10_REAL_TPM]


def apply_classifiers(tpms, clsfrs):
    for classifier in clsfrs:
        column_name = classifier.get_column_name()
        tpms[column_name] = tpms.apply(
            classifier.get_classification_value, axis=1)


def get_stats(tpms, tp_tpms):
    stats_dict = {stat.name: stat.calculate(tpms, tp_tpms)
                  for stat in stats.get_statistics()}
    return pd.DataFrame([stats_dict])


def get_grouped_stats(tpms, tp_tpms, column_name):
    grouped = tpms.groupby(column_name)
    tp_grouped = tp_tpms.groupby(column_name)
    summary = grouped.describe()
    tp_summary = tp_grouped.describe()

    stats_dict = {stat.name: stat.calculate_grouped(
        grouped, summary, tp_grouped, tp_summary)
        for stat in stats.get_statistics()}
    return pd.DataFrame.from_dict(stats_dict)


def get_distribution(tpms, classifier, ascending):
    values = tpms.apply(classifier.get_value, axis=1)
    values.sort(ascending=ascending)

    xbounds = classifier.get_distribution_plot_range()
    if xbounds is None:
        xbounds = (values.min(), values.max())

    xvals = np.linspace(xbounds[0], xbounds[1], CUMULATIVE_DISTRIBUTION_POINTS)

    size = float(len(values))
    yvals = [100 * len(values[values < x if ascending else values > x]) / size
             for x in xvals]

    return xvals, yvals


def get_distribution_stats(non_zero_tpms, tp_tpms, classifier, ascending):
    xvals, nz_yvals = get_distribution(non_zero_tpms, classifier, ascending)
    xvals, tp_yvals = get_distribution(tp_tpms, classifier, ascending)

    stats_dict = {}
    stats_dict[classifier.get_column_name()] = xvals
    stats_dict[NON_ZERO_PERCENTAGE] = nz_yvals
    stats_dict[TRUE_POSITIVE_PERCENTAGE] = tp_yvals

    return pd.DataFrame.from_dict(stats_dict)
