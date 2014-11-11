#!/usr/bin/env python

"""
Usage:
    simulate_read_bias [{log_option_spec} --out-prefix=<out-prefix>]
        [--paired-end] --num-reads=<num-reads>
        <pwm-file> <reads_file>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
-n --num-reads=<num-reads>
    Number of reads to output.
--out-prefix=<out-prefix>
    String to be prepended to input file names for output [default: bias]
--paired-end
    Indicates the reads file contains paired-end reads.
<pwm-file>
    PWM file with positional base weights used to bias reads.
<reads_file>
    FASTA/Q file containing single or paired end reads.
"""

import collections
import docopt
import options as opt
import os.path
import pwm
import random
import schema
import sys

from __init__ import __version__

NUM_READS = "--num-reads"
OUT_PREFIX = "--out-prefix"
PAIRED_END = "--paired-end"
PWM_FILE = "<pwm-file>"
READS_FILE = "<reads_file>"


ReadScore = collections.namedtuple("ReadScore", ["read_number", "score"])


class SequenceLinePicker:
    def __init__(self, lines_per_fragment):
        self.lines_per_fragment = lines_per_fragment

    def __call__(self, line_no, line):
        return (line_no - 1) % self.lines_per_fragment == 0


class OutputPicker:
    def __init__(self, scores, lines_per_fragment):
        self.scores = scores
        self.lines_per_fragment = lines_per_fragment
        self.index = 0

    def __call__(self, line_no, line):
        if self.index >= len(self.scores):
            return False

        read_number = line_no / self.lines_per_fragment
        if self.scores[self.index].read_number < read_number:
            self.index += 1
            if self.index >= len(self.scores):
                return False

        return self.scores[self.index].read_number == read_number


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        options[NUM_READS] = opt.validate_int_option(
            options[NUM_READS],
            "Number of reads must be non-negative", nonneg=True)
        opt.validate_file_option(
            options[PWM_FILE], "PWM file should exist")
        opt.validate_file_option(
            options[READS_FILE], "Reads file should exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_fragment_counts(reads_file, num_reads, paired_end):
    with_errors = reads_file.endswith("fastq")

    num_fragments = num_reads
    lines_per_fragment = 2
    if with_errors:
        lines_per_fragment *= 2
    if paired_end:
        num_fragments /= 2
        lines_per_fragment *= 2

    return num_fragments, lines_per_fragment


def _yield_elements(enumerable, element_picker):
    return (elem for elem_no, elem in enumerate(enumerable)
            if element_picker(elem_no, elem))


def _score_fragments(reads_file, bias_pwm, num_fragments, lines_per_fragment):
    scores = []
    with open(reads_file, 'r') as rfile:
        for i, line in enumerate(
                _yield_elements(rfile, SequenceLinePicker(lines_per_fragment))):
            score = ReadScore(
                i, random.random() * bias_pwm.score(line.strip()))
            scores.append(score)

    if num_fragments > len(scores):
        sys.exit("Input file(s) did not contain enough fragments " +
                 "({ni} found, {no} required)".
                 format(ni=len(scores), no=num_fragments))

    return scores


def _select_scores(scores, num_fragments):
    selected_scores = scores[:]
    selected_scores.sort(key=lambda x: x.score, reverse=True)
    selected_scores = selected_scores[1: num_fragments]
    selected_scores.sort(key=lambda x: x.read_number)
    return selected_scores


def _write_output_file(input_file, out_prefix, scores, lines_per_fragment):
    dirname = os.path.dirname(os.path.abspath(input_file))
    basename = os.path.basename(input_file)
    output_file = os.path.join(dirname, out_prefix + "." + basename)

    with open(input_file, 'r') as in_f, open(output_file, 'w') as out_f:
        for i, line in enumerate(_yield_elements(
                in_f, OutputPicker(scores, lines_per_fragment))):
            out_f.write(line)


def _simulate_bias(logger, options):
    # Read PWM file
    logger.info("Reading PWM file " + options[PWM_FILE])
    bias_pwm = pwm.PWM(options[PWM_FILE])

    # Iterate through fragments, storing positions and scores
    logger.info("Scoring fragments according to PWM")
    num_fragments, lines_per_fragment = _get_fragment_counts(
        options[READS_FILE], options[NUM_READS], options[PAIRED_END])
    scores = _score_fragments(
        options[READS_FILE], bias_pwm, num_fragments, lines_per_fragment)
    logger.info("...scored {n} fragments.".format(n=len(scores)))

    # Sort fragments by score and select the required number of highest-scoring
    # fragments
    logger.info("Sorting and selecting {n} highest scoring fragments ".
                format(n=num_fragments))
    selected_scores = _select_scores(scores, num_fragments)

    # Write selected fragments to output file(s)
    logger.info("Writing selected fragments to output files")
    _write_output_file(options[READS_FILE], options[OUT_PREFIX],
                       selected_scores, lines_per_fragment)


def _main(docstring):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(docstring)
    options = docopt.docopt(
        docstring, version="simulate_read_bias v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Simulate bias by preferentially selecting reads according to a position
    # weight matrix
    _simulate_bias(logger, options)


if __name__ == "__main__":
    _main(__doc__)
