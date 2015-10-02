import numpy as np

from . import tpms as t


class _Classifier(object):
    @classmethod
    def capitalized(cls, text):
        return text[:1].upper() + text[1:]

    def __init__(self, name, value_extractor,
                 grouped_stats=True, distribution_plot_range=None,
                 title=None, units=None):
        self.name = name
        self.value_extractor = value_extractor
        self.grouped_stats = grouped_stats
        self.distribution_plot_range = distribution_plot_range
        self.title = title if title else _Classifier.capitalized(name)
        self.units = units

    def get_value(self, row):
        return self.value_extractor(row)

    def get_classification_value(self, row):
        return self.get_value(row)

    def produces_grouped_stats(self):
        return self.grouped_stats

    def produces_distribution_plots(self):
        return not self.grouped_stats

    def get_distribution_plot_range(self):
        return self.distribution_plot_range

    def get_value_labels(self, num_labels):
        return range(1, num_labels + 1)

    def get_value_name(self, value):
        return str(value)

    def get_stats_file_suffix(self, ascending=True):
        suffix = "_stats" if self.produces_grouped_stats() else \
            "_distribution_stats_" + ("asc" if ascending else "desc")
        suffix += "_by_" + self.name.replace(' ', '_')
        return suffix

    def get_axis_label(self):
        label = self.title
        if self.units:
            label += " ({u})".format(u=self.units)
        return label


class _LevelsClassifier(_Classifier):
    def __init__(self, name, value_extractor, levels,
                 closed=False, title=None, units=None):

        _Classifier.__init__(self, name, value_extractor,
                             title=title, units=units)

        self.levels = levels
        self.closed = closed

        self.level_names = ["<= " + str(l) for l in levels]
        if not self.closed:
            self.level_names += ["> " + str(levels[-1])]

    def get_classification_value(self, row):
        row_value = _Classifier.get_classification_value(self, row)
        for i, level in enumerate(self.levels):
            if row_value <= level:
                return i
        return len(self.levels)

    def get_value_labels(self, num_labels):
        return self.level_names[:num_labels]

    def get_value_name(self, value):
        return self.level_names[value]


_CLASSIFIERS = []

_CLASSIFIERS.append(_LevelsClassifier(
    "gene transcript number", lambda x: x[t.TRANSCRIPT_COUNT],
    [5, 10, 15, 20]))

_CLASSIFIERS.append(_Classifier(
    "absolute percent error", lambda x: abs(x[t.PERCENT_ERROR]),
    grouped_stats=False, distribution_plot_range=(0, 100)))

_CLASSIFIERS.append(_LevelsClassifier(
    "log10 real TPM", lambda x: x[t.LOG10_REAL_TPM],
    [0, 0.5, 1, 1.5],
    title=r"$log_{10}$ Real TPM"))

_CLASSIFIERS.append(_LevelsClassifier(
    "transcript length", lambda x: x[t.LENGTH],
    [500, 1000, 3000], units="bp"))

_CLASSIFIERS.append(_LevelsClassifier(
    "unique sequence percentage",
    lambda x: 100 * float(x[t.UNIQUE_SEQ_LENGTH]) / x[t.LENGTH]
    if x[t.LENGTH] > 0 else 0.0,
    [20, 40, 60, 80, 100],
    closed=True))

_CLASSIFIERS.append(_LevelsClassifier(
    "unique sequence length",
    lambda x: x[t.UNIQUE_SEQ_LENGTH],
    [0, 100, 300, 1000], units="bp"))


def get_classifiers():
    return set(_CLASSIFIERS)
