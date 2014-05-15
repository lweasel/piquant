import pandas as pd
import piquant.tpms as t

REAL_TPMS_VALS = [0.05, 0.02, 15, 2, 10, 30, 11]
CALC_TPMS_VALS = [0.03, 20, 3, 0.01, 5, 20, 10]
GROUPS = [0, 1, 0, 1, 0, 1, 1]

GROUP_TEST_COL = "group_test"


def _get_test_tpms():
    tpms = pd.DataFrame.from_dict({
        t.REAL_TPM: REAL_TPMS_VALS,
        t.CALCULATED_TPM: CALC_TPMS_VALS,
        GROUP_TEST_COL: GROUPS
    })
    return tpms


def true_positive(real_tpm, calculated_tpm):
    return real_tpm > t.NOT_PRESENT_CUTOFF and \
        calculated_tpm > t.NOT_PRESENT_CUTOFF


def true_negative(real_tpm, calculated_tpm):
    return real_tpm < t.NOT_PRESENT_CUTOFF and \
        calculated_tpm < t.NOT_PRESENT_CUTOFF


def false_negative(real_tpm, calculated_tpm):
    return real_tpm > t.NOT_PRESENT_CUTOFF and \
        calculated_tpm < t.NOT_PRESENT_CUTOFF


def false_positive(real_tpm, calculated_tpm):
    return real_tpm < t.NOT_PRESENT_CUTOFF and \
        calculated_tpm > t.NOT_PRESENT_CUTOFF


def test_mark_positives_negatives_marks_correct_entries_as_true_positive():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(tpms)
    for index, row in tpms.iterrows():
        if true_positive(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.TRUE_POSITIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.TRUE_POSITIVE]


def test_mark_positives_negatives_marks_correct_entries_as_false_positive():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(tpms)
    for index, row in tpms.iterrows():
        if false_positive(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.FALSE_POSITIVE]


def test_mark_positives_negatives_marks_correct_entries_as_true_negative():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(tpms)
    for index, row in tpms.iterrows():
        if true_negative(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.TRUE_NEGATIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_POSITIVE]
            assert not row[t.FALSE_NEGATIVE]
        else:
            assert not row[t.TRUE_NEGATIVE]


def test_mark_positives_negatives_marks_correct_entries_as_false_negative():
    tpms = _get_test_tpms()
    t.mark_positives_and_negatives(tpms)
    for index, row in tpms.iterrows():
        if false_negative(row[t.REAL_TPM], row[t.CALCULATED_TPM]):
            assert row[t.FALSE_NEGATIVE]
            assert not row[t.FALSE_POSITIVE]
            assert not row[t.TRUE_NEGATIVE]
            assert not row[t.TRUE_POSITIVE]
        else:
            assert not row[t.FALSE_NEGATIVE]
