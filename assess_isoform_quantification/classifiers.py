import fpkms as f

CLASSIFIERS = []


def Classifier(cls):
    CLASSIFIERS.append(cls)
    return cls


class BaseClassifier():
    def __init__(self, value_extractor):
        self.value_extractor = value_extractor
        self.grouped_stats = True
        self.distribution_plots = True

    def get_value(self, row):
        return self.value_extractor(row)

    def get_classification_value(self, row):
        return self.get_value(row)

    def produces_grouped_stats(self):
        return self.grouped_stats

    def produces_distribution_plots(self):
        return self.distribution_plots

    def get_distribution_plot_range(self):
        return None

    def get_value_labels(self, num_labels):
        return range(1, num_labels + 1)


@Classifier
class GeneTranscriptNumberClassifier(BaseClassifier):
    def __init__(self):
        BaseClassifier.__init__(self, lambda x: x[f.TRANSCRIPT_COUNT])
        self.distribution_plots = False

    def get_column_name(self):
        return "gene transcript number"


@Classifier
class PercentErrorClassifier(BaseClassifier):
    def __init__(self):
        BaseClassifier.__init__(self, lambda x: abs(x[f.PERCENT_ERROR]))
        self.grouped_stats = False

    def get_column_name(self):
        return "absolute percent error"

    def get_distribution_plot_range(self):
        return (0, 100)


class LevelsClassifier(BaseClassifier):
    def __init__(self, levels, value_extractor, closed):
        BaseClassifier.__init__(self, value_extractor)

        self.levels = levels
        self.closed = closed
        self.distribution_plots = False

        if self.closed:
            self.level_names = ["<= " + str(l) for l in levels]
        else:
            self.level_names = ["<= " + str(l) for l in levels] \
                + ["> " + str(levels[-1])]

    def get_classification_value(self, row):
        row_value = BaseClassifier.get_classification_value(self, row)
        for i, level in enumerate(self.levels):
            if row_value <= level:
                return i
        return len(self.levels)

    def get_value_labels(self, num_labels):
        return self.level_names[:num_labels]


@Classifier
class RealAbundanceClassifier(LevelsClassifier):
    LEVELS = [0, 0.5, 1, 1.5]

    def __init__(self):
        LevelsClassifier.__init__(
            self, RealAbundanceClassifier.LEVELS,
            lambda x: x[f.LOG10_REAL_FPKM], False)

    def get_column_name(self):
        return "log10 real FPKM"


@Classifier
class TranscriptLengthClassifier(LevelsClassifier):
    LEVELS = [1000, 3162]

    def __init__(self):
        LevelsClassifier.__init__(
            self, TranscriptLengthClassifier.LEVELS,
            lambda x: x[f.LENGTH], False)

    def get_column_name(self):
        return "transcript length"


@Classifier
class UniqueSequencePercentageClassifier(LevelsClassifier):
    LEVELS = [20, 40, 60, 80, 100]

    def __init__(self):
        LevelsClassifier.__init__(
            self, UniqueSequencePercentageClassifier.LEVELS,
            lambda x: 100 * float(x[f.UNIQUE_SEQ_LENGTH]) / x[f.LENGTH],
            True)

    def get_column_name(self):
        return "unique sequence percentage"
