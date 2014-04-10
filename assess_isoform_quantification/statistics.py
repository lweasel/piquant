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
