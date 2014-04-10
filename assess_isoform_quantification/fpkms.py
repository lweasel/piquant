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

FALSE_POSITIVE = "false-pos"
FALSE_NEGATIVE = "false-neg"
TRUE_POSITIVE = "true-pos"
TRUE_NEGATIVE = "true-neg"

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
