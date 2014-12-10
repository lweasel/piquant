from . import tpms as t


class _Classifier(object):
    def __init__(self, column_name, value_extractor,
                 grouped_stats=True, distribution_plot_range=None,
                 plot_title=None):
        self.column_name = column_name
        self.value_extractor = value_extractor
        self.grouped_stats = grouped_stats
        self.distribution_plot_range = distribution_plot_range
        self.plot_title = plot_title if plot_title else column_name

    def get_column_name(self):
        return self.column_name

    def get_plot_title(self):
        return self.plot_title

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

    def get_stats_file_suffix(self, ascending=True):
        suffix = "_stats" if self.produces_grouped_stats() else \
            "_distribution_stats_" + ("asc" if ascending else "desc")
        suffix += "_by_" + self.get_column_name().replace(' ', '_')
        return suffix


class _LevelsClassifier(_Classifier):
    def __init__(self, column_name, value_extractor, levels,
                 closed=False, plot_title=None):

        _Classifier.__init__(self, column_name, value_extractor,
                             plot_title=plot_title)

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


_CLASSIFIERS = []

_CLASSIFIERS.append(_Classifier(
    "gene transcript number", lambda x: x[t.TRANSCRIPT_COUNT]))

_CLASSIFIERS.append(_Classifier(
    "absolute percent error", lambda x: abs(x[t.PERCENT_ERROR]),
    grouped_stats=False, distribution_plot_range=(0, 100)))

_CLASSIFIERS.append(_LevelsClassifier(
    "log10 real TPM", lambda x: x[t.LOG10_REAL_TPM],
    [0, 0.5, 1, 1.5],
    plot_title=r"$log_{10}$ real TPM"))

_CLASSIFIERS.append(_LevelsClassifier(
    "transcript length", lambda x: x[t.LENGTH],
    [1000, 3162]))

_CLASSIFIERS.append(_LevelsClassifier(
    "unique sequence percentage",
    lambda x: 100 * float(x[t.UNIQUE_SEQ_LENGTH]) / x[t.LENGTH] if x[t.LENGTH] > 0 else 0.0,
    [20, 40, 60, 80, 100],
    closed=True))


def get_classifiers():
    return set(_CLASSIFIERS)
