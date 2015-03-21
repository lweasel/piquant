# pylint: disable=W0142
# pylint: disable=R0903

import functools
import operator

NEUTRAL_SCORE = 0.25


class PWM(object):
    def __init__(self, filename):
        base_weights = []
        with open(filename, 'r') as pwm_file:
            base_weights = [line.strip().split(",") for line in pwm_file]

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

        return functools.reduce(operator.mul, scores)
