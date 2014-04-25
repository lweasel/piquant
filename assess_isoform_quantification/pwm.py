import operator

NEUTRAL_SCORE = 0.25


class PWM:
    def __init__(self, filename):
        base_weights = []
        with open(filename, 'r') as f:
            base_weights = [line.strip().split(",") for line in f]

        get_base_scores = lambda x: zip(['a', 'c', 'g', 't', 'n'],
                                        list(x) + [NEUTRAL_SCORE])
        self.pos = [{base: float(freq) for base, freq in get_base_scores(l)}
                    for l in zip(*base_weights)]
        self.length = len(self.pos)

    def score(self, sequence):
        sequence = sequence[0: self.length].lower()
        scores = [self.pos[i][base] for i, base in enumerate(sequence)]

        score_length = len(scores)
        if score_length < self.length:
            scores += [NEUTRAL_SCORE] * (self.length - score_length)

        return reduce(operator.mul, scores)
