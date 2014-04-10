import pandas as pd
import numpy as np
import statistics as stats

NUM_FPKMS = "num-fpkms"
TP_NUM_FPKMS = "tp-num-fpkms"
TRANSCRIPT_COUNT = "num-transcripts"
LENGTH = "length"
UNIQUE_SEQ_LENGTH = "unique-length"
REAL_FPKM = "real-fpkm"
CALCULATED_FPKM = "calc-fpkm"
PERCENT_ERROR = "percent-error"
TP_MEDIAN_PERCENT_ERROR = "tp-median-percent-error"
TP_ERROR_FRACTION = "tp-error-frac"
LOG10_REAL_FPKM = "log10-real-fpkm"
LOG10_CALCULATED_FPKM = "log10-calc-fpkm"
LOG10_RATIO = "log-ratio"
TP_LOG10_RATIO_MEAN = "tp-log-ratio-mean"
TP_LOG10_RATIO_STD = "tp-log-ratio-std"
TP_LOG10_RATIO_MEDIAN = "tp-log-ratio-med"
TP_LOG10_FPKM_RHO = "tp-log-fpkm-rho"

FALSE_POSITIVE = "false-pos"
FALSE_NEGATIVE = "false-neg"
TRUE_POSITIVE = "true-pos"
TRUE_NEGATIVE = "true-neg"
TP_COUNT = "tp-count"

SENSITIVITY = "sensitivity"
SPECIFICITY = "specificity"

NOT_PRESENT_CUTOFF = 0.1


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


def get_sensitivity(fpkms):
    num_tp = len(fpkms[fpkms[TRUE_POSITIVE]])
    num_fn = len(fpkms[fpkms[FALSE_NEGATIVE]])
    if num_tp + num_fn == 0:
        return 1
    return float(num_tp) / (num_tp + num_fn)


def get_specificity(fpkms):
    num_fp = len(fpkms[fpkms[FALSE_POSITIVE]])
    num_tn = len(fpkms[fpkms[TRUE_NEGATIVE]])
    if num_fp + num_tn == 0:
        return 1
    return float(num_tn) / (num_tn + num_fp)


def get_error_fraction(fpkms, error_percent):
    num_errors = len(fpkms[abs(fpkms[PERCENT_ERROR]) > error_percent])
    return float(num_errors) / len(fpkms)


def get_stats(fpkms, tp_fpkms):
    stats_dict = {}

    for stat in stats.get_statistics():
        stats_dict[stat.name] = stat.calculate(fpkms, tp_fpkms)

    # Spearman correlation coefficient between real and calculated FPKMs.
    rho = tp_fpkms[LOG10_CALCULATED_FPKM].corr(
        tp_fpkms[LOG10_REAL_FPKM], method='spearman')

    # The median percent error - i.e. the median of the percent errors of
    # the calculated values from the real ones
    tp_mpe = tp_fpkms[PERCENT_ERROR].median()

    stats_dict[TP_LOG10_FPKM_RHO] = rho
    stats_dict[TP_MEDIAN_PERCENT_ERROR] = tp_mpe
    stats_dict[SENSITIVITY] = get_sensitivity(fpkms)
    stats_dict[SPECIFICITY] = get_specificity(fpkms)
    stats_dict[TP_ERROR_FRACTION] = get_error_fraction(tp_fpkms, 10)

    return pd.DataFrame([stats_dict])


def get_grouped_stats(fpkms, tp_fpkms, column_name):
        grouped = fpkms.groupby(column_name)
        tp_grouped = tp_fpkms.groupby(column_name)

        summary = grouped.describe()
        main_stats = summary[REAL_FPKM].unstack()
        main_stats = main_stats.drop(
            ["mean", "std", "min", "25%", "50%", "75%", "max"], 1)

        main_stats[SENSITIVITY] = grouped.apply(get_sensitivity)
        main_stats[SPECIFICITY] = grouped.apply(get_specificity)

        tp_summary = tp_grouped.describe()
        tp_stats = tp_summary[LOG10_RATIO].unstack()

        stats_dict = {}
        for stat in stats.get_statistics():
            stats_dict[stat.name] = stat.calculate_grouped(grouped, main_stats, tp_grouped, tp_stats)
        tp_stats = tp_stats.drop(["min", "25%", "75%", "max"], 1)

        tp_stats = tp_stats.rename(columns={
            "count": TP_COUNT,
            "mean": TP_LOG10_RATIO_MEAN,
            "std": TP_LOG10_RATIO_STD,
            "50%": TP_LOG10_RATIO_MEDIAN})

        tp_stats[TP_LOG10_FPKM_RHO] = tp_grouped.apply(
            lambda x: x[LOG10_CALCULATED_FPKM].corr(
                x[LOG10_REAL_FPKM], method="spearman"))

        pe_stats = tp_summary[PERCENT_ERROR].unstack()
        tp_stats[TP_MEDIAN_PERCENT_ERROR] = pe_stats["50%"]

        tp_stats[TP_ERROR_FRACTION] = tp_grouped.apply(get_error_fraction, 10)

        stats_df = pd.DataFrame.from_dict(stats_dict)

        return pd.concat([main_stats, tp_stats, stats_df], axis=1)
