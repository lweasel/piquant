import numpy as np
import numpy.testing as npt
import pandas as pd
import piquant.tpms as t

REAL_TPMS_VALS = [0.05, 0.02, 15, 2, 10, 30, 11]
CALC_TPMS_VALS = [0.03, 20, 3, 0.01, 5, 20, 10]
GROUPS = [0, 1, 0, 1, 0, 1, 1]

GROUP_TEST_COL = "group_test"
NOT_PRESENT_CUTOFF = 0.1


def _get_test_tpms():
    tpms = pd.DataFrame.from_dict({
        t.REAL_TPM: REAL_TPMS_VALS,
        t.CALCULATED_TPM: CALC_TPMS_VALS,
        GROUP_TEST_COL: GROUPS
    })
    return tpms


def _get_test_tp_tpms():
    tpms = _get_test_tpms()
    return tpms[(tpms[t.REAL_TPM] > NOT_PRESENT_CUTOFF) &
                (tpms[t.CALCULATED_TPM] > NOT_PRESENT_CUTOFF)]


def _true_positive(real_tpm, calculated_tpm):
    return real_tpm > NOT_PRESENT_CUTOFF and \
        calculated_tpm > NOT_PRESENT_CUTOFF


def _true_negative(real_tpm, calculated_tpm):
    return real_tpm < NOT_PRESENT_CUTOFF and \
        calculated_tpm < NOT_PRESENT_CUTOFF


def _false_negative(real_tpm, calculated_tpm):
    return real_tpm > NOT_PRESENT_CUTOFF and \
        calculated_tpm < NOT_PRESENT_CUTOFF


def _false_positive(real_tpm, calculated_tpm):
    return real_tpm < NOT_PRESENT_CUTOFF and \
        calculated_tpm > NOT_PRESENT_CUTOFF


class _DummyStatistic:
    def __init__(self, name, true_positives):
        self.name = name
        self.true_positives = true_positives

    def calculate(self, tpms, tp_tpms):
        df = tp_tpms if self.true_positives else tpms
        return len(df)

    def calculate_grouped(self, grouped, summary, tp_grouped, tp_summary):
        df = tp_summary if self.true_positives else summary
        return df[t.REAL_TPM].unstack()["count"]


class _DummyClassifier:
    def __init__(self, name, value_func=lambda x: x[t.REAL_TPM]):
        self.name = name
        self.value_func = value_func

    def get_column_name(self):
        return self.name

    def get_classification_value(self, x):
        return self.value_func(x)


def test_mark_positives_negatives_marks_correct_entries_as_true_positive():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(NOT_PRESENT_CUTOFF, tpms)
    for index, row in tpms.iterrows():
        if _true_positive(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.TRUE_POSITIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.TRUE_POSITIVE]


def test_mark_positives_negatives_marks_correct_entries_as_false_positive():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(NOT_PRESENT_CUTOFF, tpms)
    for index, row in tpms.iterrows():
        if _false_positive(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.FALSE_POSITIVE]


def test_mark_positives_negatives_marks_correct_entries_as_true_negative():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(NOT_PRESENT_CUTOFF, tpms)
    for index, row in tpms.iterrows():
        if _true_negative(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_POSITIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.TRUE_NEGATIVE]


def test_mark_positives_negatives_marks_correct_entries_as_false_negative():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(NOT_PRESENT_CUTOFF, tpms)
    for index, row in tpms.iterrows():
        if _false_negative(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.FALSE_NEGATIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.TRUE_POSITIVE]
        else:
            assert not row[t.FALSE_NEGATIVE]


def test_get_true_positives_returns_correct_number_of_entries():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(NOT_PRESENT_CUTOFF, tpms)
    tp_tpms = t.get_true_positives(tpms)

    assert len(tp_tpms) == \
        len([x for x, y in zip(CALC_TPMS_VALS, REAL_TPMS_VALS)
            if x > NOT_PRESENT_CUTOFF and y > NOT_PRESENT_CUTOFF])


def test_calculate_percent_error_calculates_correct_values():
    tpms = _get_test_tpms()
    t.calculate_percent_error(tpms)

    for index, row in tpms.iterrows():
        val = 100 * ((CALC_TPMS_VALS[index] - REAL_TPMS_VALS[index])
                     / float(REAL_TPMS_VALS[index]))
        npt.assert_approx_equal(row[t.PERCENT_ERROR], val)


def test_calculate_log_ratios_calculates_correct_values():
    tpms = _get_test_tpms()
    t.calculate_log_ratios(tpms)

    for index, row in tpms.iterrows():
        val = np.log10(CALC_TPMS_VALS[index] / float(REAL_TPMS_VALS[index]))
        npt.assert_approx_equal(row[t.LOG10_RATIO], val)


def test_apply_classifiers_adds_correct_columns():
    classifier_names = ["a", "b", "c"]
    classifiers = [_DummyClassifier(x) for x in classifier_names]

    tpms = _get_test_tpms()
    t.apply_classifiers(tpms, classifiers)

    for name in classifier_names:
        assert tpms[name] is not None


def test_apply_classifiers_calculates_correct_values():
    name = "dummy"
    val = 5
    classifiers = [_DummyClassifier(name, lambda x: x[t.CALCULATED_TPM] + val)]

    tpms = _get_test_tpms()
    t.apply_classifiers(tpms, classifiers)

    for index, row in tpms.iterrows():
        assert row[name] == CALC_TPMS_VALS[index] + val


def test_get_stats_returns_correct_number_of_statistics():
    num_statistics = 5
    statistics = [_DummyStatistic(str(i), False)
                  for i in range(num_statistics)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_stats(tpms, tp_tpms, statistics)
    assert len(stats.columns) == num_statistics


def test_get_stats_returns_correct_column_names():
    name1 = "dummy1"
    name2 = "dummy2"
    statistics = [_DummyStatistic(name1, False), _DummyStatistic(name2, False)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_stats(tpms, tp_tpms, statistics)
    assert name1 in stats.columns
    assert name2 in stats.columns


def test_get_stats_calculates_correct_values():
    name1 = "dummy1"
    name2 = "dummy2"
    statistics = [_DummyStatistic(name1, False), _DummyStatistic(name2, True)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_stats(tpms, tp_tpms, statistics)
    assert stats[name1].ix[0] == len(tpms)
    assert stats[name2].ix[0] == len(tp_tpms)


def test_get_grouped_stats_returns_correct_number_of_statistics():
    num_statistics = 5
    statistics = [_DummyStatistic("c" + str(i), False)
                  for i in range(num_statistics)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_grouped_stats(tpms, tp_tpms, GROUP_TEST_COL, statistics)
    assert len(stats.columns) == num_statistics


def test_get_grouped_stats_returns_correct_column_names():
    name1 = "dummy1"
    name2 = "dummy2"
    statistics = [_DummyStatistic(name1, False), _DummyStatistic(name2, False)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_grouped_stats(tpms, tp_tpms, GROUP_TEST_COL, statistics)
    assert name1 in stats.columns
    assert name2 in stats.columns


def test_get_grouped_stats_calculates_correct_values():
    name1 = "dummy1"
    name2 = "dummy2"
    statistics = [_DummyStatistic(name1, False), _DummyStatistic(name2, True)]

    tpms = _get_test_tpms()
    tp_tpms = _get_test_tp_tpms()

    stats = t.get_grouped_stats(tpms, tp_tpms, GROUP_TEST_COL, statistics)
    for group in set(GROUPS):
        assert stats[name1].ix[group] == len(tpms[tpms[GROUP_TEST_COL] == group])
        assert stats[name2].ix[group] == len(tp_tpms[tp_tpms[GROUP_TEST_COL] == group])
