import piquant.statistics as statistics
import piquant.tpms as t
import numpy as np
import pandas as pd
import scipy.stats as scistats
import test_tpms


def _get_test_tpms():
    tpms = pd.DataFrame.from_dict({
        t.REAL_TPM: test_tpms.REAL_TPMS_VALS,
        t.CALCULATED_TPM: test_tpms.CALC_TPMS_VALS,
        test_tpms.GROUP_TEST_COL: test_tpms.GROUPS
    })

    t.calculate_log_ratios(tpms)
    t.calculate_percent_error(tpms)
    t.mark_positives_and_negatives(tpms)

    return tpms, t.get_true_positives(tpms)


def __get_test_grouped_tpms():
    tpms, tp_tpms = _get_test_tpms()

    grouped = tpms.groupby(test_tpms.GROUP_TEST_COL)
    tp_grouped = tp_tpms.groupby(test_tpms.GROUP_TEST_COL)
    summary = grouped.describe()
    tp_summary = tp_grouped.describe()

    return grouped, summary, tp_grouped, tp_summary


def _tpm_pairs(filter=lambda r, c: True):
    return [(r, c) for r, c in zip(test_tpms.REAL_TPMS_VALS,
                                   test_tpms.CALC_TPMS_VALS)
            if filter(r, c)]


def _tp_tpm_pairs():
    return _tpm_pairs(lambda r, c: test_tpms.true_positive(r, c))


def _group_tpm_pairs(group_val, filter=lambda r, c: True):
    return [(r, c) for r, c, gv in
            zip(test_tpms.REAL_TPMS_VALS,
                test_tpms.CALC_TPMS_VALS,
                test_tpms.GROUPS) if
            (gv == group_val and filter(r, c))]


def _group_tp_tpm_pairs(group_val):
    return _group_tpm_pairs(
        group_val, lambda r, c: test_tpms.true_positive(r, c))


def _check_statistic_value(stat_class, calculator, pair_func):
    tpms, tp_tpms = _get_test_tpms()
    stat = stat_class()
    correct_value = calculator(pair_func())
    assert stat.calculate(tpms, tp_tpms) == correct_value


def _check_grouped_statistic_values(stat_class, calculator, grouped_pair_func):
    g, s, tp_g, tp_s = __get_test_grouped_tpms()
    stat = stat_class()
    grouped_stats = stat.calculate_grouped(g, s, tp_g, tp_s)
    correct_value_calculator = lambda x: calculator(grouped_pair_func(x))
    group_count_test = \
        lambda x: grouped_stats.ix[x] == correct_value_calculator(x)
    assert all([group_count_test(gv) for gv in set(test_tpms.GROUPS)])


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


def _number_of_tpms(tpm_pairs):
    return len(tpm_pairs)


def test_number_of_tpms_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._NumberOfTPMs, _number_of_tpms, _tpm_pairs)


def test_number_of_tpms_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._NumberOfTPMs, _number_of_tpms, _group_tpm_pairs)


def test_number_of_true_positive_tpms_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._NumberOfTruePositiveTPMs,
        _number_of_tpms, _tp_tpm_pairs)


def test_number_of_true_positive_tpms_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._NumberOfTruePositiveTPMs,
        _number_of_tpms, _group_tp_tpm_pairs)


def _spearman(tpm_pairs):
    rs, cs = zip(*tpm_pairs)
    return scistats.spearmanr(np.array(rs), np.array(cs))[0]


def test_spearman_correlation_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._SpearmanCorrelation, _spearman, _tp_tpm_pairs)


def test_spearman_correlation_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._SpearmanCorrelation, _spearman, _group_tp_tpm_pairs)


def _error_fraction(tpm_pairs):
    error_percent = lambda r, c: abs(100 * (c - r) / float(r))
    above_threshold = \
        [r for r, c in tpm_pairs if
         error_percent(r, c) >
            statistics._TruePositiveErrorFraction.ERROR_PERCENTAGE_THRESHOLD]
    return len(above_threshold) / float(len(tpm_pairs))


def test_true_positive_error_fraction_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._TruePositiveErrorFraction,
        _error_fraction, _tp_tpm_pairs)


def test_true_positive_error_fraction_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._TruePositiveErrorFraction,
        _error_fraction, _group_tp_tpm_pairs)


def _median_percent_error(tpm_pairs):
    error_percent = lambda r, c: 100 * (c - r) / float(r)
    percent_errors = [error_percent(r, c) for r, c in tpm_pairs]
    return np.median(percent_errors)


def test_median_percent_error_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._MedianPercentError,
        _median_percent_error, _tp_tpm_pairs)


def test_median_percent_error_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._MedianPercentError,
        _median_percent_error, _group_tp_tpm_pairs)


def _sensitivity(tpm_pairs):
    num_tp = sum([test_tpms.true_positive(r, c) for r, c in tpm_pairs])
    num_fn = sum([test_tpms.false_negative(r, c) for r, c in tpm_pairs])
    return float(num_tp) / (num_tp + num_fn)


def test_sensitivity_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._Sensitivity, _sensitivity, _tpm_pairs)


def test_sensitivity_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._Sensitivity, _sensitivity, _group_tpm_pairs)


def _specificity(tpm_pairs):
    num_fp = sum([test_tpms.false_positive(r, c) for r, c in tpm_pairs])
    num_tn = sum([test_tpms.true_negative(r, c) for r, c in tpm_pairs])
    return float(num_tn) / (num_tn + num_fp)


def test_specificity_statistic_calculates_correct_value():
    _check_statistic_value(
        statistics._Specificity, _specificity, _tpm_pairs)


def test_specificity_statistic_calculates_correct_grouped_values():
    _check_grouped_statistic_values(
        statistics._Specificity, _specificity, _group_tpm_pairs)
