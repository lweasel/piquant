import assess_isoform_quantification.statistics as statistics
import assess_isoform_quantification.fpkms as f
import numpy as np
import pandas as pd
import scipy.stats as scistats

REAL_FPKMS_VALS = [0.05, 0.02, 15, 2, 10, 30, 11]
CALC_FPKMS_VALS = [0.03, 20, 3, 0.01, 5, 20, 10]
GROUPS = [0, 1, 0, 1, 0, 1, 1]

GROUP_TEST_COL = "group_test"


def _get_test_fpkms():
    fpkms = pd.DataFrame.from_dict({
        f.REAL_FPKM: REAL_FPKMS_VALS,
        f.CALCULATED_FPKM: CALC_FPKMS_VALS,
        GROUP_TEST_COL: GROUPS
    })

    f.calculate_log_ratios(fpkms)
    f.calculate_percent_error(fpkms)
    f.mark_positives_and_negatives(fpkms)

    return fpkms, f.get_true_positives(fpkms)


def __get_test_grouped_fpkms():
    fpkms, tp_fpkms = _get_test_fpkms()

    grouped = fpkms.groupby(GROUP_TEST_COL)
    tp_grouped = tp_fpkms.groupby(GROUP_TEST_COL)
    summary = grouped.describe()
    tp_summary = tp_grouped.describe()

    return grouped, summary, tp_grouped, tp_summary


def _true_positive(real_fpkm, calculated_fpkm):
    return real_fpkm > f.NOT_PRESENT_CUTOFF and \
        calculated_fpkm > f.NOT_PRESENT_CUTOFF


def _true_negative(real_fpkm, calculated_fpkm):
    return real_fpkm < f.NOT_PRESENT_CUTOFF and \
        calculated_fpkm < f.NOT_PRESENT_CUTOFF


def _false_negative(real_fpkm, calculated_fpkm):
    return real_fpkm > f.NOT_PRESENT_CUTOFF and \
        calculated_fpkm < f.NOT_PRESENT_CUTOFF


def _false_positive(real_fpkm, calculated_fpkm):
    return real_fpkm < f.NOT_PRESENT_CUTOFF and \
        calculated_fpkm > f.NOT_PRESENT_CUTOFF


def _fpkm_pairs(filter=lambda r, c: True):
    return [(r, c) for r, c in zip(REAL_FPKMS_VALS, CALC_FPKMS_VALS) if
            filter(r, c)]


def _tp_fpkm_pairs():
    return _fpkm_pairs(lambda r, c: _true_positive(r, c))


def _group_fpkm_pairs(group_val, filter=lambda r, c: True):
    return [(r, c) for r, c, gv in
            zip(REAL_FPKMS_VALS, CALC_FPKMS_VALS, GROUPS) if
            (gv == group_val and filter(r, c))]


def _group_tp_fpkm_pairs(group_val):
    return _group_fpkm_pairs(group_val, lambda r, c: _true_positive(r, c))


def _check_statistic_value(stat_class, calculator, pair_func):
    fpkms, tp_fpkms = _get_test_fpkms()
    stat = stat_class()
    correct_value = calculator(pair_func())
    assert stat.calculate(fpkms, tp_fpkms) == correct_value


def _check_grouped_statistic_values(stat_class, calculator, grouped_pair_func):
    g, s, tp_g, tp_s = __get_test_grouped_fpkms()
    stat = stat_class()
    grouped_stats = stat.calculate_grouped(g, s, tp_g, tp_s)
    correct_value_calculator = lambda x: calculator(grouped_pair_func(x))
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


def _number_of_fpkms(fpkm_pairs):
    return len(fpkm_pairs)


def test_number_of_fpkms_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._NumberOfFPKMs, _number_of_fpkms, _fpkm_pairs)


def test_number_of_fpkms_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._NumberOfFPKMs, _number_of_fpkms, _group_fpkm_pairs)


def test_number_of_true_positive_fpkms_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._NumberOfTruePositiveFPKMs,
        _number_of_fpkms, _tp_fpkm_pairs)


def test_number_of_true_positive_fpkms_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._NumberOfTruePositiveFPKMs,
        _number_of_fpkms, _group_tp_fpkm_pairs)


def _spearman(fpkm_pairs):
    rs, cs = zip(*fpkm_pairs)
    return scistats.spearmanr(np.array(rs), np.array(cs))[0]


def test_spearman_correlation_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._SpearmanCorrelation, _spearman, _tp_fpkm_pairs)


def test_spearman_correlation_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._SpearmanCorrelation, _spearman, _group_tp_fpkm_pairs)


def _error_fraction(fpkm_pairs):
    error_percent = lambda r, c: abs(100 * (c - r) / float(r))
    above_threshold = \
        [r for r, c in fpkm_pairs if
         error_percent(r, c) >
            statistics._TruePositiveErrorFraction.ERROR_PERCENTAGE_THRESHOLD]
    return len(above_threshold) / float(len(fpkm_pairs))


def test_true_positive_error_fraction_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._TruePositiveErrorFraction,
        _error_fraction, _tp_fpkm_pairs)


def test_true_positive_error_fraction_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._TruePositiveErrorFraction,
        _error_fraction, _group_tp_fpkm_pairs)


def _median_percent_error(fpkm_pairs):
    error_percent = lambda r, c: 100 * (c - r) / float(r)
    percent_errors = [error_percent(r, c) for r, c in fpkm_pairs]
    return np.median(percent_errors)


def test_median_percent_error_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._MedianPercentError,
        _median_percent_error, _tp_fpkm_pairs)


def test_median_percent_error_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._MedianPercentError,
        _median_percent_error, _group_tp_fpkm_pairs)


def _sensitivity(fpkm_pairs):
    num_tp = sum([_true_positive(r, c) for r, c in fpkm_pairs])
    num_fn = sum([_false_negative(r, c) for r, c in fpkm_pairs])
    return float(num_tp) / (num_tp + num_fn)


def test_sensitivity_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._Sensitivity, _sensitivity, _fpkm_pairs)


def test_sensitivity_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._Sensitivity, _sensitivity, _group_fpkm_pairs)


def _specificity(fpkm_pairs):
    num_fp = sum([_false_positive(r, c) for r, c in fpkm_pairs])
    num_tn = sum([_true_negative(r, c) for r, c in fpkm_pairs])
    return float(num_tn) / (num_tn + num_fp)


def test_specificity_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._Specificity, _specificity, _fpkm_pairs)


def test_specificity_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._Specificity, _specificity, _group_fpkm_pairs)
