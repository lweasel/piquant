"""
Functions and classes for calculating statistics from the results of a
transcript quantification run. Exports:

get_statistics: Return all statistic instances.
get_graphable_statistics: Return statistic instances suitable for graphing.
get_graphable_by_classifier_statistics: Return statistic instances suitable for
graphing by a classifier.
"""

import fpkms as f

NUM_FPKMS = "num-fpkms"

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


def get_graphable_by_classifier_statistics():
    """Return a set of statistics suitable for graphing by a classifier.

    Return a set of objects each of which can calculate a certain statistic
    from the results of a transcript quantification run, and is interesting to
    plot when FPKMs are grouped by a classifier (e.g. by the number of isoforms
    for each transcript's originating gene).
    """
    return set([s for s in get_statistics() if s.graphable_by_classifier])


def _Statistic(cls):
    # Mark a class as capable of calculate a statistic for the results of a
    # quantification run.
    _STATISTICS.append(cls())
    return cls


class _BaseStatistic():
    # Base for classes capable of calculating a statistic
    def __init__(self, name, title, graphable=True,
                 graphable_by_classifier=True, stat_range=None):
        self.name = name
        self.title = title
        self.graphable = graphable
        self.graphable_by_classifier = graphable_by_classifier
        self.stat_range = stat_range

    def calculate(self, fpkms, tp_fpkms):
        """Calculate the statistic for a set of FPKMs.

        Calculate a single statistic value for the results of a quantification
        run.
        fpkms: A pandas DataFrame describing the result of a quantification
        run.
        tp_fpkms: A pandas DataFrame describing those results of a
        quantification run for which both real and calculated FPKMs were above
        a threshold value indicating "presence" of the transcript.
        """
        raise NotImplementedError

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        """Calculate the statistic for a set of FPKMs grouped by a classifier.

        Calculate a set of statistic values for the results of a quantification
        run which have been grouped according to a certain method of
        classifying transcripts. Should return a pandas Series instance.
        grouped: A pandas GroupBy instance describing the results of a
        quantification run grouped by a certain classifier of transcripts.
        grp_summary: A pandas DataFrame containing basic summary statistics
        calculated for 'grouped'.
        tp_grouped: A pandas GroupBy instance describing those results of a
        quantification run for which both real and calculated FPKMs were above
        a threshold value indicating "presence" of the transcript, grouped
        by a certain classifier of transcripts.
        tp_grp_summary: A pandas DataFrame containing basic sumnmary statistics
        calculated for 'tp_grouped'.
        """
        raise NotImplementedError


@_Statistic
class _NumberOfFPKMs(_BaseStatistic):
    # Calculates the total number of transcript FPKMs in the results.
    def __init__(self):
        _BaseStatistic.__init__(
            self, NUM_FPKMS, "No. FPKMs",
            graphable=False, graphable_by_classifier=False,
            stat_range=(0, None))

    def calculate(self, fpkms, tp_fpkms):
        return len(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = grp_summary[f.REAL_FPKM].unstack()
        return stats[_SUMMARY_COUNT]


@_Statistic
class _NumberOfTruePositiveFPKMs(_BaseStatistic):
    # Calculates the total number of transcript FPKMs in the results for which
    # both real and calculated FPKMs are above a threshold value indicating
    # 'presence' of the transcript.
    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-num-fpkms", "No. true pos. FPKMs",
            graphable_by_classifier=False, stat_range=(0, None))

    def calculate(self, fpkms, tp_fpkms):
        return len(tp_fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[f.REAL_FPKM].unstack()
        return stats[_SUMMARY_COUNT]


@_Statistic
class _SpearmanCorrelation(_BaseStatistic):
    # Calculates the Spearman rank correlation coefficient between calculated
    # and real FPKMs for 'true positive' transcript FPKMs (those for which both
    # real and calculated FPKM were above a threshold value indicating
    # 'presence' of the transcript).
    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-log-fpkm-rho", "Spearman's rho",
            stat_range=_ZERO_TO_ONE_STAT_RANGE)

    @staticmethod
    def _calculate(fpkms):
        return fpkms[f.LOG10_CALCULATED_FPKM].corr(
            fpkms[f.LOG10_REAL_FPKM], method='spearman')

    def calculate(self, fpkms, tp_fpkms):
        return _SpearmanCorrelation._calculate(tp_fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(_SpearmanCorrelation._calculate)


@_Statistic
class _TruePositiveErrorFraction(_BaseStatistic):
    # Calculates the percentage of 'true positive' transcript FPKMs (those for
    # which both real and calculated FPKMs were above a threshold value
    # indicating 'presence' of the transcript) for which the calculated FPKM
    # was greater than a certain percentage above or below the real FPKM.

    ERROR_PERCENTAGE_THRESHOLD = 10

    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-error-frac", "True pos. error fraction",
            stat_range=_ZERO_TO_ONE_STAT_RANGE)

    @staticmethod
    def _calculate(fpkms, error_percent):
        num_errors = len(fpkms[abs(fpkms[f.PERCENT_ERROR]) > error_percent])
        return float(num_errors) / len(fpkms)

    def calculate(self, fpkms, tp_fpkms):
        return _TruePositiveErrorFraction._calculate(
            tp_fpkms, _TruePositiveErrorFraction.ERROR_PERCENTAGE_THRESHOLD)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(
            _TruePositiveErrorFraction._calculate,
            _TruePositiveErrorFraction.ERROR_PERCENTAGE_THRESHOLD)


@_Statistic
class _MedianPercentError(_BaseStatistic):
    # Calculates the median of the percent errors of the calculated compared to
    # real FPKMs for 'true positive' transcript FPKMs (those for which both
    # real and calculated FPKMs were above a threshold value indicating
    # 'presence' of the transcript).
    def __init__(self):
        _BaseStatistic.__init__(
            self, "tp-median-percent-error", "True pos. median % error")

    def calculate(self, fpkms, tp_fpkms):
        return tp_fpkms[f.PERCENT_ERROR].median()

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[f.PERCENT_ERROR].unstack()
        return stats[_SUMMARY_MEDIAN]


@_Statistic
class _Sensitivity(_BaseStatistic):
    # Calculates the "sensitivity" of the transcript quantification method -
    # that is, the fraction of all transcripts marked as 'present' (their
    # calculated FPKM above a threshold value) which truly were 'present'
    # (their real FPKM above a threshold value).
    def __init__(self):
        _BaseStatistic.__init__(self, "sensitivity", "Sensitivity",
                                stat_range=_ZERO_TO_ONE_STAT_RANGE)

    @staticmethod
    def _calculate(fpkms):
        num_tp = len(fpkms[fpkms[f.TRUE_POSITIVE]])
        num_fn = len(fpkms[fpkms[f.FALSE_NEGATIVE]])
        if num_tp + num_fn == 0:
            return 1
        return float(num_tp) / (num_tp + num_fn)

    def calculate(self, fpkms, tp_fpkms):
        return _Sensitivity._calculate(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(_Sensitivity._calculate)


@_Statistic
class _Specificity(_BaseStatistic):
    # Calculates the "specificity" of the transcript quantification method -
    # that is, the fraction of all transcripts marked as 'not present' (their
    # calculated FPKM below a threshold value) which truly were 'not present'
    # (their real FPKM below a threshold value).
    def __init__(self):
        _BaseStatistic.__init__(self, "specificity", "Specificity",
                                stat_range=_ZERO_TO_ONE_STAT_RANGE)

    @staticmethod
    def _calculate(fpkms):
        num_fp = len(fpkms[fpkms[f.FALSE_POSITIVE]])
        num_tn = len(fpkms[fpkms[f.TRUE_NEGATIVE]])
        if num_fp + num_tn == 0:
            return 1
        return float(num_tn) / (num_tn + num_fp)

    def calculate(self, fpkms, tp_fpkms):
        return _Specificity._calculate(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(_Specificity._calculate)
