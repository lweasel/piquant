import pandas as pd
import piquant.classifiers as classifiers


def _get_test_classifier(
        column_name="dummy", value_extractor=lambda x: x,
        grouped_stats=True, distribution_plot_range=None):
    return classifiers._Classifier(
        column_name, value_extractor, grouped_stats, distribution_plot_range)


def _get_test_levels_classifier(
        column_name="dummy", value_extractor=lambda x: x,
        levels=[10, 20, 30], closed=True):
    return classifiers._LevelsClassifier(
        column_name, value_extractor, levels, closed)


def test_get_classifiers_returns_classifiers_instances():
    clsfrs = classifiers.get_classifiers()
    assert all([isinstance(c, classifiers._Classifier) for c in clsfrs])


def test_classifier_get_column_name_returns_correct_name():
    name = "column name"
    c = _get_test_classifier(column_name=name)
    assert c.name == name


def test_classifier_get_value_return_correct_value():
    col_name = "column name"
    col_value = "column value"
    df = pd.DataFrame.from_dict([{col_name: col_value, "dummy": "dummy"}])
    c = _get_test_classifier(value_extractor=lambda x: x[col_name])
    assert c.get_value(df.ix[0]) == col_value


def test_classifier_get_classification_value_returns_correct_value():
    col_name = "column name"
    col_value = "column value"
    df = pd.DataFrame.from_dict([{col_name: col_value, "dummy": "dummy"}])
    c = _get_test_classifier(value_extractor=lambda x: x[col_name])
    assert c.get_classification_value(df.ix[0]) == col_value


def test_classifier_produces_grouped_stats_returns_correct_value():
    gs = False
    c = _get_test_classifier(grouped_stats=gs)
    assert c.produces_grouped_stats() == gs


def test_classifier_produces_distribution_plots_returns_correct_value():
    gs = False
    c = _get_test_classifier(grouped_stats=gs)
    assert c.produces_distribution_plots() != gs


def test_classifier_get_distribution_plot_range_returns_correct_value():
    dpr = (10, 30)
    c = _get_test_classifier(distribution_plot_range=dpr)
    assert c.get_distribution_plot_range() == dpr


def test_classifier_get_value_labels_returns_correct_labels():
    num_labels = 5
    c = _get_test_classifier()
    assert c.get_value_labels(num_labels) == range(1, num_labels + 1)


def test_classifier_get_stats_file_suffix_returns_correct_suffix_for_grouped_stats():
    c = _get_test_classifier(column_name="column name", grouped_stats=True)
    assert c.get_stats_file_suffix() == "stats_by_column_name"


def test_classifier_get_stats_file_suffix_returns_correct_suffix_for_non_grouped_stats_and_ascending_order():
    c = _get_test_classifier(column_name="column name", grouped_stats=False)
    assert c.get_stats_file_suffix(True) == \
        "distribution_stats_asc_by_column_name"


def test_classifier_get_stats_file_suffix_returns_correct_suffix_for_non_grouped_stats_and_descending_order():
    c = _get_test_classifier(column_name="column name", grouped_stats=False)
    assert c.get_stats_file_suffix(False) == \
        "distribution_stats_desc_by_column_name"


def test_levels_classifier_get_classification_value_returns_correct_value():
    col_name = "column name"
    df = pd.DataFrame.from_dict([{col_name: 25, "dummy": "dummy"}])
    c = _get_test_levels_classifier(value_extractor=lambda x: x[col_name])
    assert c.get_classification_value(df.ix[0]) == 2


def test_levels_classifier_get_value_labels_returns_correct_labels_for_closed_classifier():
    levels = [10, 20, 30, 40]
    c = _get_test_levels_classifier(levels=levels)
    assert c.get_value_labels(len(levels)) == \
        ["<= 10", "<= 20", "<= 30", "<= 40"]
    assert c.get_value_labels(len(levels) - 1) == \
        ["<= 10", "<= 20", "<= 30"]


def test_levels_classifier_get_value_labels_returns_correct_labels_for_open_classifier():
    levels = [10, 20, 30, 40]
    c = _get_test_levels_classifier(levels=levels, closed=False)
    assert c.get_value_labels(len(levels) + 1) == \
        ["<= 10", "<= 20", "<= 30", "<= 40", "> 40"]
    assert c.get_value_labels(len(levels)) == \
        ["<= 10", "<= 20", "<= 30", "<= 40"]
    assert c.get_value_labels(len(levels) - 1) == \
        ["<= 10", "<= 20", "<= 30"]
