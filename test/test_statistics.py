import assess_isoform_quantification.statistics as statistics
import assess_isoform_quantification.fpkms as f
import pandas as pd

REAL_FPKMS_VALS = [0.05, 0.02, 1, 2, 10, 30]
CALC_FPKMS_VALS = [0.03, 20, 3, 0.01, 5, 20]
GROUPS = [0, 1, 0, 1, 0, 1]

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


def _check_statistic_value(stat_class, correct_value):
    fpkms, tp_fpkms = _get_test_fpkms()
    stat = stat_class()
    assert stat.calculate(fpkms, tp_fpkms) == correct_value


def _check_grouped_statistic_values(stat_class, correct_value_calculator):
    g, s, tp_g, tp_s = __get_test_grouped_fpkms()
    stat = stat_class()
    grouped_stats = stat.calculate_grouped(g, s, tp_g, tp_s)
    group_count_test = \
        lambda x: grouped_stats.ix[x] == correct_value_calculator(x)
    assert all([group_count_test(gv) for gv in set(GROUPS)])


def test_get_statistics_returns_statistics_instances():
    stats = statistics.get_statistics()
    assert all([isinstance(s, statistics._BaseStatistic) for s in stats])


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
    correct_value = len(fpkms)
    _check_statistic_value(statistics._NumberOfFPKMs, correct_value)


def test_number_of_fpkms_statistic_calculates_correct_grouped_values():
    fpkms, tp_fpkms = _get_test_fpkms()
    cvc = lambda x: len(fpkms[fpkms[GROUP_TEST_COL] == x])
    _check_grouped_statistic_values(statistics._NumberOfFPKMs, cvc)


def test_number_of_true_positive_fpkms_statistic_calculates_correct_value():
    fpkms, tp_fpkms = _get_test_fpkms()
    correct_value = len(tp_fpkms)
    _check_statistic_value(
        statistics._NumberOfTruePositiveFPKMs, correct_value)


def test_number_of_true_positive_fpkms_statistic_calculates_correct_grouped_values():
    fpkms, tp_fpkms = _get_test_fpkms()
    cvc = lambda x: len(tp_fpkms[tp_fpkms[GROUP_TEST_COL] == x])
    _check_grouped_statistic_values(statistics._NumberOfTruePositiveFPKMs, cvc)
