"""
Functions and classes for calculating statistics from the results of a
transcript quantification run. Exports:

get_statistics: Return all statistic instances.
get_graphable_statistics: Return statistic instances suitable for graphing.
"""

import classifiers
import itertools
import math
import os.path

from . import piquant_options as po
from . import tpms as t

TP_NUM_TPMS = "tp-num-tpms"
OVERALL_STATS_PREFIX = "overall"

_SUMMARY_COUNT = "count"
_SUMMARY_MEDIAN = "50%"
_ZERO_TO_ONE_STAT_RANGE = (-0.025, 1.025)

_STATISTICS = []


def get_statistics():
    """Return a set of all statistics instances.

    Return a set of objects each of which can calculate a certain statistic
    from the results of a transcript quantification run. Objects for all
    statistic types are returned.
    """
    return set(_STATISTICS)


def get_graphable_statistics():
    """Return a set of statistic instances suitable for graphing.

    Return a set of objects each of which can calculate a certain statistic
    from the results of a transcript quantification run, and is interesting to
    plot across as a certain parameter (e.g. read length) is varied across
    quantification runs.
    """
    return set([s for s in get_statistics() if s.graphable])


def get_stratified_stats_types():
    clsfrs = classifiers.get_classifiers()
    grp_clsfrs = [c for c in clsfrs if c.produces_grouped_stats()]
    dist_clsfrs = [c for c in clsfrs if c.produces_distribution_plots()]

    return [{}] + \
        [{"classifier": c} for c in grp_clsfrs] + \
        [{"classifier": c, "ascending": asc}
            for c, asc in itertools.product(dist_clsfrs, [True, False])]


def get_stats_file(directory, prefix, tpm_level, classifier=None, ascending=False):
    return os.path.join(directory, "_".join([prefix, tpm_level])) + \
        (classifier.get_stats_file_suffix(ascending=ascending)
            if classifier else "_stats") + ".csv"


def write_stats_data(filename, data_frame, **kwargs):
    with open(filename, "w") as out_file:
        data_frame.to_csv(out_file, float_format="%.5f", **kwargs)


def _Statistic(cls):
    # Mark a class as capable of calculate a statistic for the results of a
    # quantification run.
    _STATISTICS.append(cls())
    return cls


class _BaseStatistic():
    # Base for classes capable of calculating a statistic

    @classmethod
    def set_options(cls, options):
        pass

    def __init__(self, name, title, graphable=True):
        self.name = name
        self.title = title
        self.graphable = graphable

    def calculate(self, tpms, tp_tpms):
        """Calculate the statistic for a set of TPMs.

        Calculate a single statistic value for the results of a quantification
        run.
        tpms: A pandas DataFrame describing the result of a quantification
        run.
        tp_tpms: A pandas DataFrame describing those results of a
        quantification run for which both real and calculated TPMs were above
        a threshold value indicating "presence" of the transcript.
        """
        raise NotImplementedError

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        """Calculate the statistic for a set of TPMs grouped by a classifier.

        Calculate a set of statistic values for the results of a quantification
        run which have been grouped according to a certain method of
        classifying transcripts. Should return a pandas Series instance.
        grouped: A pandas GroupBy instance describing the results of a
        quantification run grouped by a certain classifier of transcripts.
        grp_summary: A pandas DataFrame containing basic summary statistics
        calculated for 'grouped'.
        tp_grouped: A pandas GroupBy instance describing those results of a
        quantification run for which both real and calculated TPMs were above
        a threshold value indicating "presence" of the transcript, grouped
        by a certain classifier of transcripts.
        tp_grp_summary: A pandas DataFrame containing basic sumnmary statistics
        calculated for 'tp_grouped'.
        """
        raise NotImplementedError


@_Statistic
class _NumberOfTPMs(_BaseStatistic):
    # Calculates the total number of transcript TPMs in the results.
    def __init__(self):
        _BaseStatistic.__init__(self, "num-tpms", "No. TPMs", graphable=False)

    def calculate(self, tpms, tp_tpms):
        return len(tpms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = grp_summary[t.REAL_TPM].unstack()
        return stats[_SUMMARY_COUNT]

    def stat_range(self, vals_range):
        return (0, None)


@_Statistic
class _NumberOfTruePositiveTPMs(_BaseStatistic):
    # Calculates the total number of transcript TPMs in the results for which
    # both real and calculated TPMs are above a threshold value indicating
    # 'presence' of the transcript.
    def __init__(self):
        _BaseStatistic.__init__(self, TP_NUM_TPMS, "No. true positive TPMs")

    def calculate(self, tpms, tp_tpms):
        return len(tp_tpms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[t.REAL_TPM].unstack()
        return stats[_SUMMARY_COUNT]

    def stat_range(self, vals_range):
        return (0, None)


@_Statistic
class _SpearmanCorrelation(_BaseStatistic):
    # Calculates the Spearman rank correlation coefficient between calculated
    # and real TPMs for 'true positive' transcript TPMs (those for which both
    # real and calculated TPM were above a threshold value indicating
    # 'presence' of the transcript).
    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-log-tpm-rho", "Spearman's rho")

    @staticmethod
    def _calculate(tpms):
        return tpms[t.LOG10_CALCULATED_TPM].corr(
            tpms[t.LOG10_REAL_TPM], method='spearman')

    def calculate(self, tpms, tp_tpms):
        return _SpearmanCorrelation._calculate(tp_tpms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(_SpearmanCorrelation._calculate)

    def stat_range(self, vals_range):
        min_val = math.floor(vals_range[0] * 5) / 5.0
        return (min_val - 0.01, 1.01)


@_Statistic
class _TruePositiveErrorFraction(_BaseStatistic):
    # Calculates the percentage of 'true positive' transcript TPMs (those for
    # which both real and calculated TPMs were above a threshold value
    # indicating 'presence' of the transcript) for which the calculated TPM
    # was greater than a certain percentage above or below the real TPM.

    ERROR_FRACTION_THRESHOLD = 10

    @classmethod
    def set_options(cls, options):
        cls.ERROR_FRACTION_THRESHOLD = options[po.ERROR_FRACTION_THRESHOLD]

    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-error-frac", "True positive error fraction")

    @staticmethod
    def _calculate(tpms, error_percent):
        num_errors = len(tpms[abs(tpms[t.PERCENT_ERROR]) > error_percent])
        return float(num_errors) / len(tpms)

    def calculate(self, tpms, tp_tpms):
        return _TruePositiveErrorFraction._calculate(
            tp_tpms, _TruePositiveErrorFraction.ERROR_FRACTION_THRESHOLD)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(
            _TruePositiveErrorFraction._calculate,
            _TruePositiveErrorFraction.ERROR_FRACTION_THRESHOLD)

    def stat_range(self, vals_range):
        return _ZERO_TO_ONE_STAT_RANGE


@_Statistic
class _MedianPercentError(_BaseStatistic):
    # Calculates the median of the percent errors of the calculated compared to
    # real TPMs for 'true positive' transcript TPMs (those for which both
    # real and calculated TPMs were above a threshold value indicating
    # 'presence' of the transcript).
    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-median-percent-error", "True positive median % error")

    def calculate(self, tpms, tp_tpms):
        return tp_tpms[t.PERCENT_ERROR].median()

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[t.PERCENT_ERROR].unstack()
        return stats[_SUMMARY_MEDIAN]

    def stat_range(self, vals_range):
        division = 5.0
        closest_div = lambda x: math.floor(x / division) * division
        ymin, ymax = vals_range
        if ymin > 0 and ymax > 0:
            return (0, closest_div(ymax) + division)
        elif ymin < 0 and ymax < 0:
            return (closest_div(ymin), 0)
        else:
            return (closest_div(ymin), closest_div(ymax) + division)


@_Statistic
class _Sensitivity(_BaseStatistic):
    # Calculates the "sensitivity" of the transcript quantification method -
    # that is, the fraction of all transcripts considered to be 'present'
    # (their real TPM above a threshold value - that is, both true positives
    # and false negatives), which were correctly identified as being present
    # (just the true positives).
    def __init__(self):
        _BaseStatistic.__init__(self, "sensitivity", "Sensitivity")

    @staticmethod
    def _calculate(tpms):
        num_tp = len(tpms[tpms[t.TRUE_POSITIVE]])
        num_fn = len(tpms[tpms[t.FALSE_NEGATIVE]])
        if num_tp + num_fn == 0:
            return 1
        return float(num_tp) / (num_tp + num_fn)

    def calculate(self, tpms, tp_tpms):
        return _Sensitivity._calculate(tpms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(_Sensitivity._calculate)

    def stat_range(self, vals_range):
        min_val = math.floor(vals_range[0] * 5) / 5.0
        return (min_val - 0.01, 1.01)


@_Statistic
class _Specificity(_BaseStatistic):
    # Calculates the "specificity" of the transcript quantification method -
    # that is, the fraction of all transcripts considered to be 'not present'
    # (their real TPM below a threshold value - that is, both true negatives
    # and false positives), which were correctly identified as being present
    # (just the true negatives).
    def __init__(self):
        _BaseStatistic.__init__(self, "specificity", "Specificity")

    @staticmethod
    def _calculate(tpms):
        num_fp = len(tpms[tpms[t.FALSE_POSITIVE]])
        num_tn = len(tpms[tpms[t.TRUE_NEGATIVE]])
        if num_fp + num_tn == 0:
            return 1
        return float(num_tn) / (num_tn + num_fp)

    def calculate(self, tpms, tp_tpms):
        return _Specificity._calculate(tpms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(_Specificity._calculate)

    def stat_range(self, vals_range):
        min_val = math.floor(vals_range[0] * 5) / 5.0
        return (min_val - 0.01, 1.01)
