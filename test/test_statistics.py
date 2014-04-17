import assess_isoform_quantification.statistics as statistics
import assess_isoform_quantification.fpkms as f
import pandas as pd

REAL_FPKMS_VALS = [0.5, 1, 2, 10, 30]
CALC_FPKMS_VALS = [0.3, 3, 0.1, 5, 20]
GROUPS = [0, 0, 1, 0, 1]

GROUP_TEST_COL = "group_test"


def _get_test_fpkms():
    fpkms = pd.DataFrame.from_dict({
        f.REAL_FPKM: REAL_FPKMS_VALS,
        f.CALCULATED_FPKM: CALC_FPKMS_VALS,
        GROUP_TEST_COL: GROUPS
    })

    fpkms[f.TRUE_POSITIVE] = \
        (fpkms[f.REAL_FPKM] > f.NOT_PRESENT_CUTOFF) & \
        (fpkms[f.CALCULATED_FPKM] > f.NOT_PRESENT_CUTOFF)
    tp_fpkms = fpkms[fpkms[f.TRUE_POSITIVE]]

    return fpkms, tp_fpkms


def __get_test_grouped_fpkms():
    fpkms, tp_fpkms = _get_test_fpkms()

    grouped = fpkms.groupby(GROUP_TEST_COL)
    tp_grouped = tp_fpkms.groupby(GROUP_TEST_COL)
    summary = grouped.describe()
    tp_summary = tp_grouped.describe()

    return grouped, summary, tp_grouped, tp_summary


def test_get_statistics_returns_statistics_instances():
    stats = statistics.get_statistics()
    assert all([isinstance(s, statistics.BaseStatistic) for s in stats])


def test_get_graphable_statistics_returns_subset_of_statistics():
    stats = statistics.get_statistics()
    g_stats = statistics.get_graphable_statistics()
    assert g_stats <= stats


def test_get_graphable_statistics_returns_graphable_instances():
    g_stats = statistics.get_graphable_statistics()
    assert all([s.graphable for s in g_stats])


def test_get_graphable_by_classifier_statistics_returns_subset_of_statistics():
    stats = statistics.get_statistics()
    g_stats = statistics.get_graphable_by_classifier_statistics()
    assert g_stats <= stats


def test_get_graphable_by_classifier_statistics_return_graphable_by_classifier_instances():
    g_stats = statistics.get_graphable_by_classifier_statistics()
    assert all([s.graphable_by_classifier for s in g_stats])


def test_number_of_fpkms_statistic_calculates_correct_value():
    fpkms, tp_fpkms = _get_test_fpkms()
    stat = statistics._NumberOfFPKMs()
    assert stat.calculate(fpkms, tp_fpkms) == len(REAL_FPKMS_VALS)


def test_number_of_fpkms_statistic_calculates_correct_grouped_value():
    grouped, summary, tp_grouped, tp_summary = __get_test_grouped_fpkms()
    stat = statistics._NumberOfFPKMs()
    grouped_stats = stat.calculate_grouped(
        grouped, summary, tp_grouped, tp_summary)

    group_count = lambda x: len([v for v in GROUPS if v == x])
    group_count_test = lambda x: grouped_stats.ix[x] == group_count(x)
    assert all([group_count_test(gv) for gv in set(GROUPS)])
