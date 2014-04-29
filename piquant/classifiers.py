import fpkms as f
import utils


class _Classifier():
    def __init__(self, column_name, value_extractor,
                 grouped_stats=True, distribution_plot_range=None):
        self.column_name = column_name
        self.value_extractor = value_extractor
        self.grouped_stats = grouped_stats
        self.distribution_plot_range = distribution_plot_range

    def get_column_name(self):
        return self.column_name

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
            "_distribution_stats_" + utils.get_order_string(ascending)
        suffix += "_by_" + utils.spaces_to_underscore(self.get_column_name())
        return suffix


class _LevelsClassifier(_Classifier):
    def __init__(self, column_name, value_extractor, levels, closed=False):
        _Classifier.__init__(self, column_name, value_extractor)

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
    "gene transcript number", lambda x: x[f.TRANSCRIPT_COUNT]))

_CLASSIFIERS.append(_Classifier(
    "absolute percent error", lambda x: abs(x[f.PERCENT_ERROR]),
    grouped_stats=False, distribution_plot_range=(0, 100)))

_CLASSIFIERS.append(_LevelsClassifier(
    "log10 real FPKM", lambda x: x[f.LOG10_REAL_FPKM],
    [0, 0.5, 1, 1.5]))

_CLASSIFIERS.append(_LevelsClassifier(
    "transcript length", lambda x: x[f.LENGTH],
    [1000, 3162]))

_CLASSIFIERS.append(_LevelsClassifier(
    "unique sequence percentage",
    lambda x: 100 * float(x[f.UNIQUE_SEQ_LENGTH]) / x[f.LENGTH],
    [20, 40, 60, 80, 100],
    closed=True))


def get_classifiers():
    return set(_CLASSIFIERS)
