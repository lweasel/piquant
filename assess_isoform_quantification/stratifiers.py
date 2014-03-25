import fpkms as f

STRATIFIERS = []


def Stratifier(cls):
    STRATIFIERS.append(cls)
    return cls


@Stratifier
class GeneTranscriptNumberStratifier:
    def get_column_name(self):
        return "gene transcript number"

    def get_stratification_value(self, row):
        return row[f.TRANSCRIPT_COUNT]

    def get_value_labels(self, num_labels):
        return range(1, num_labels + 1)


class LevelsStratifier:
    def __init__(self, levels, value_extractor):
        self.levels = levels
        self.level_names = ["< " + str(l) for l in levels] \
            + ["> " + str(levels[-1])]
        self.value_extractor = value_extractor

    def get_stratification_value(self, row):
        row_value = self.value_extractor(row)
        for i, level in enumerate(self.levels):
            if row_value < level:
                return i
        return len(self.levels)

    def get_value_labels(self, num_labels):
        return self.level_names[:num_labels]


@Stratifier
class RealAbundanceStratifier(LevelsStratifier):
    LEVELS = [0, 0.5, 1, 1.5]

    def __init__(self):
        LevelsStratifier.__init__(
            self, RealAbundanceStratifier.LEVELS,
            lambda x: x[f.LOG10_REAL_FPKM])

    def get_column_name(self):
        return "log10 real FPKM"


@Stratifier
class TranscriptLengthStratifier(LevelsStratifier):
    LEVELS = [1000, 3162]

    def __init__(self):
        LevelsStratifier.__init__(
            self, TranscriptLengthStratifier.LEVELS, lambda x: x[f.LENGTH])

    def get_column_name(self):
        return "transcript length"
