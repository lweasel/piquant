"""
Functions and classes for calculating statistics from the results of a
transcript quantification run. Exports:

get_statistics: Return all statistic instances.
get_graphable_statistics: Return statistic instances suitable for graphing.
"""

import itertools
import math
import os.path

from . import classifiers
from . import piquant_options as po
from . import tpms as t

EXPRESSED_TPMS = "expressed-tpms"
OVERALL_STATS_PREFIX = "overall"
GROUPED_STATS_PREFIX = "grouped"

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

    return [(None, False)] + \
        [(c, False) for c in grp_clsfrs] + \
        [(c, asc) for c, asc in itertools.product(dist_clsfrs, [True, False])]


def get_stats_file(directory, prefix, tpm_level,
                   classifier=None, ascending=False):

    return os.path.join(directory, "_".join([prefix, tpm_level])) + \
        (classifier.get_stats_file_suffix(ascending=ascending)
            if classifier else "_stats") + ".csv"


def write_stats_data(filename, data_frame, **kwargs):
    with open(filename, "w") as out_file:
        data_frame.to_csv(out_file, float_format="%.5f", **kwargs)


def _rounded_stat_range(vals_range):
    granularity = 0.2
    border = 0.01
    min_val = math.floor(vals_range[0] / granularity) * granularity
    return (min_val - border, 1 + border)


def _statistic(cls):
    # Mark a class as capable of calculate a statistic for the results of a
    # quantification run.
    _STATISTICS.append(cls())
    return cls


class _BaseStatistic(object):
    # Base for classes capable of calculating a statistic

    @classmethod
    def set_options(cls, options):
        pass

    def __init__(self, name, title, graphable=True):
        self.name = name
        self.title = title
        self.graphable = graphable

    def get_axis_label(self):
        return self.title

    def calculate(self, tpms, expressed_tpms):
        """Calculate the statistic for a set of TPMs.

        Calculate a single statistic value for the results of a quantification
        run.
        tpms: A pandas DataFrame describing the result of a quantification
        run.
        expressed_tpms: A pandas DataFrame describing the results of a
        quantification run for which the real TPMs were above a threshold value
        indicating 'presence' of the transcript.
        """
        raise NotImplementedError

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        """Calculate the statistic for a set of TPMs grouped by a classifier.

        Calculate a set of statistic values for the results of a quantification
        run which have been grouped according to a certain method of
        classifying transcripts. Should return a pandas Series instance.
        grouped: A pandas GroupBy instance describing the results of a
        quantification run grouped by a certain classifier of transcripts.
        grp_summary: A pandas DataFrame containing basic summary statistics
        calculated for 'grouped'.
        expressed_grouped: A pandas GroupBy instance describing those results
        of a quantification run for which the real TPMs were above a threshold
        value indicating 'presence' of the transcript, grouped by a certain
        classifier of transcripts.
        expressed_grp_summary: A pandas DataFrame containing basic summary
        statistics calculated for 'expressed_grouped'.
        """
        raise NotImplementedError


@_statistic
class _NumberOfExpressedTPMs(_BaseStatistic):
    # Calculates the total number of transcripts with non-zero real TPMs in the results.
    def __init__(self):
        _BaseStatistic.__init__(
            self, EXPRESSED_TPMS, "No. expressed TPMs", graphable=False)

    def calculate(self, tpms, expressed_tpms):
        return len(expressed_tpms)

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):

        stats = expressed_grp_summary[t.REAL_TPM].unstack()
        return stats[_SUMMARY_COUNT]

    def stat_range(self, vals_range):
        del vals_range
        return (0, None)


@_statistic
class _SpearmanCorrelation(_BaseStatistic):
    # Calculates the Spearman rank correlation coefficient between calculated
    # and real TPMs, for those transcripts with non-zero real TPMs.
    def __init__(self):
        _BaseStatistic.__init__(
            self, "log-tpm-rho", r"Spearman's $\rho$")

    @staticmethod
    def _calculate(tpms):
        return tpms[t.CALCULATED_TPM].corr(
            tpms[t.REAL_TPM], method='spearman')

    def calculate(self, tpms, expressed_tpms):
        return _SpearmanCorrelation._calculate(expressed_tpms)

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        return expressed_grouped.apply(_SpearmanCorrelation._calculate)

    def stat_range(self, vals_range):
        return _rounded_stat_range(vals_range)


@_statistic
class _ErrorFraction(_BaseStatistic):
    # Calculates the percentage of calculated TPMs which were greater than a
    # certain percentage above or below the real TPM, for those transcripts
    # with non-zero real TPMs.

    ERROR_FRACTION_THRESHOLD = 10

    @classmethod
    def set_options(cls, options):
        cls.ERROR_FRACTION_THRESHOLD = options[po.ERROR_FRACTION_THRESHOLD.name]

    def __init__(self):
        _BaseStatistic.__init__(
            self, "error-frac", "Error fraction")

    @staticmethod
    def _calculate(tpms, error_percent):
        num_errors = len(tpms[abs(tpms[t.PERCENT_ERROR]) > error_percent])
        return float(num_errors) / len(tpms)

    def calculate(self, tpms, expressed_tpms):
        return _ErrorFraction._calculate(
            expressed_tpms, _ErrorFraction.ERROR_FRACTION_THRESHOLD)

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        return expressed_grouped.apply(
            _ErrorFraction._calculate,
            _ErrorFraction.ERROR_FRACTION_THRESHOLD)

    def stat_range(self, vals_range):
        return _rounded_stat_range(vals_range)


@_statistic
class _MedianPercentError(_BaseStatistic):
    # Calculates the median of the percent errors of the calculated compared to
    # real TPMs, for those transcripts with non-zero real TPMs
    def __init__(self):
        _BaseStatistic.__init__(
            self, "median-percent-error", "Median % error")

    @staticmethod
    def _calculate(tpms):
        tpms = tpms[tpms[t.REAL_TPM] > 0]
        return tpms[t.PERCENT_ERROR].median()

    def calculate(self, tpms, expressed_tpms):
        return expressed_tpms[t.PERCENT_ERROR].median()

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        stats = expressed_grp_summary[t.PERCENT_ERROR].unstack()
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


@_statistic
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
        num_tp = len(
            tpms[(tpms[t.REAL_TPM] > 0) & (tpms[t.CALCULATED_TPM] > 0)])
        num_fn = len(
            tpms[(tpms[t.REAL_TPM] > 0) & (tpms[t.CALCULATED_TPM] == 0)])
        if num_tp + num_fn == 0:
            return 1
        return float(num_tp) / (num_tp + num_fn)

    def calculate(self, tpms, expressed_tpms):
        return _Sensitivity._calculate(tpms)

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        return grouped.apply(_Sensitivity._calculate)

    def stat_range(self, vals_range):
        return _rounded_stat_range(vals_range)


@_statistic
class _Specificity(_BaseStatistic):
    # Calculates the "specificity" of the transcript quantification method -
    # that is, the fraction of all transcripts considered to be 'not present'
    # (their real TPM below a threshold value - that is, both true negatives
    # and false positives), which were correctly identified as being not
    # present (just the true negatives).
    def __init__(self):
        _BaseStatistic.__init__(self, "specificity", "Specificity")

    @staticmethod
    def _calculate(tpms):
        num_fp = len(
            tpms[(tpms[t.REAL_TPM] == 0) & (tpms[t.CALCULATED_TPM] > 0)])
        num_tn = len(
            tpms[(tpms[t.REAL_TPM] == 0) & (tpms[t.CALCULATED_TPM] == 0)])
        if num_fp + num_tn == 0:
            return 1
        return float(num_tn) / (num_tn + num_fp)

    def calculate(self, tpms, expressed_tpms):
        return _Specificity._calculate(tpms)

    def calculate_grouped(self, grouped, grp_summary,
                          expressed_grouped, expressed_grp_summary):
        return grouped.apply(_Specificity._calculate)

    def stat_range(self, vals_range):
        return _rounded_stat_range(vals_range)
