import fpkms as f

_STATISTICS = []


def get_statistics():
    return [s() for s in _STATISTICS]


def Statistic(cls):
    _STATISTICS.append(cls)
    return cls


class BaseStatistic():
    def __init__(self, name):
        self.name = name


@Statistic
class NumberOfFPKMs(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(self, "num-fpkms")

    def calculate(self, fpkms, tp_fpkms):
        return len(fpkms)

    def calculate_grouped(self, grouped, grp_stats, tp_grouped, tp_grp_stats):
        return grp_stats["count"]


@Statistic
class NumberOfTruePositiveFPKMs(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(self, "tp-num-fpkms")

    def calculate(self, fpkms, tp_fpkms):
        return len(tp_fpkms)

    def calculate_grouped(self, grouped, grp_stats, tp_grouped, tp_grp_stats):
        return tp_grp_stats["count"]


@Statistic
class SpearmanCorrelation(BaseStatistic):
    def __init__(self):
        BaseStatistic.__init__(self, "tp-log-fpkm-rho")

    @staticmethod
    def _calculate(fpkms):
        return fpkms[f.LOG10_CALCULATED_FPKM].corr(
            fpkms[f.LOG10_REAL_FPKM], method='spearman')

    def calculate(self, fpkms, tp_fpkms):
        return SpearmanCorrelation._calculate(tp_fpkms)

    def calculate_grouped(self, grouped, grp_stats, tp_grouped, tp_grp_stats):
        return tp_grouped.apply(SpearmanCorrelation._calculate)
