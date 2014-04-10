import fpkms as f

_STATISTICS = []


def get_statistics():
    return [s() for s in _STATISTICS]


def Statistic(cls):
    _STATISTICS.append(cls)
    return cls


class BaseStatistic():
    def __init__(self, name, title, graphable=True, stat_range=None):
        self.name = name
        self.title = title
        self.graphable = graphable
        self.stat_range = stat_range


@Statistic
class NumberOfFPKMs(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(
            self, "num-fpkms", "No. FPKMs",
            graphable=False, stat_range=(0, None))

    def calculate(self, fpkms, tp_fpkms):
        return len(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = grp_summary[f.REAL_FPKM].unstack()
        return stats["count"]


@Statistic
class NumberOfTruePositiveFPKMs(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(
            self, "tp-num-fpkms", "No. true pos. FPKMs",
            stat_range=(0, None))

    def calculate(self, fpkms, tp_fpkms):
        return len(tp_fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[f.LOG10_RATIO].unstack()
        return stats["count"]


@Statistic
class SpearmanCorrelation(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(
            self, "tp-log-fpkm-rho", "Spearman's rho", stat_range=(0, 1))

    @staticmethod
    def _calculate(fpkms):
        return fpkms[f.LOG10_CALCULATED_FPKM].corr(
            fpkms[f.LOG10_REAL_FPKM], method='spearman')

    def calculate(self, fpkms, tp_fpkms):
        return SpearmanCorrelation._calculate(tp_fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(SpearmanCorrelation._calculate)


@Statistic
class TruePositiveErrorFraction(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(
            self, "tp-error-frac", "True pos. error fraction",
            stat_range=(0, 1))

    @staticmethod
    def _calculate(fpkms, error_percent):
        num_errors = len(fpkms[abs(fpkms[f.PERCENT_ERROR]) > error_percent])
        return float(num_errors) / len(fpkms)

    def calculate(self, fpkms, tp_fpkms):
        return TruePositiveErrorFraction._calculate(tp_fpkms, 10)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return tp_grouped.apply(TruePositiveErrorFraction._calculate, 10)


@Statistic
class MedianPercentError(BaseStatistic):
    # The median of the percent errors of the calculated FPKMS from the real
    # ones
    def __init__(self):
        BaseStatistic.__init__(
            self, "tp-median-percent-error", "True pos. median % error")

    def calculate(self, fpkms, tp_fpkms):
        return tp_fpkms[f.PERCENT_ERROR].median()

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        stats = tp_grp_summary[f.PERCENT_ERROR].unstack()
        return stats["50%"]


@Statistic
class Sensitivity(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(self, "sensitivity", "Sensitivity",
                               stat_range=(0, 1))

    @staticmethod
    def _calculate(fpkms):
        num_tp = len(fpkms[fpkms[f.TRUE_POSITIVE]])
        num_fn = len(fpkms[fpkms[f.FALSE_NEGATIVE]])
        if num_tp + num_fn == 0:
            return 1
        return float(num_tp) / (num_tp + num_fn)

    def calculate(self, fpkms, tp_fpkms):
        return Sensitivity._calculate(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(Sensitivity._calculate)


@Statistic
class Specificity(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(self, "specificity", "Specificity",
                               stat_range=(0, 1))

    @staticmethod
    def _calculate(fpkms):
        num_fp = len(fpkms[fpkms[f.FALSE_POSITIVE]])
        num_tn = len(fpkms[fpkms[f.TRUE_NEGATIVE]])
        if num_fp + num_tn == 0:
            return 1
        return float(num_tn) / (num_tn + num_fp)

    def calculate(self, fpkms, tp_fpkms):
        return Specificity._calculate(fpkms)

    def calculate_grouped(
            self, grouped, grp_summary, tp_grouped, tp_grp_summary):
        return grouped.apply(Specificity._calculate)
