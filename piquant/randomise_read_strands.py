"""
Usage:
    randomise_read_strands
        [{log_option_spec} --out-prefix=<out-prefix>]
        <reads-file>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
--out-prefix=<out-prefix>
    String to be prepended to input file name for output [default: unstranded]
<reads-file>
    FASTA/Q file containing paired-end reads with the first read corresponding
    to the sense strand.

Randomly (with probability 1/2) reassign pairs of paired-end reads such that
the first read corresponds to the antisense strand.
"""

import docopt
import itertools
import os.path
import random
import schema

from . import options as opt
from .__init__ import __version__

OUT_PREFIX = "--out-prefix"
READS_FILE = "<reads-file>"


def _validate_command_line_options(options):
    try:
        opt.validate_file_option(
            options[READS_FILE], "Reads file should exist.")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_lists_per_read(reads_file):
    return 4 if reads_file.endswith("fastq") else 2


def _write_read(out_file, read_lines):
    for line in read_lines:
        out_file.write(line)


def _randomise_read_strands(reads_file, out_prefix):
    dirname = os.path.dirname(os.path.abspath(reads_file))
    basename = os.path.basename(reads_file)
    output_file = os.path.join(dirname, out_prefix + "." + basename)

    lines_per_read = _get_lists_per_read(reads_file)

    with open(reads_file, 'r') as in_f, open(output_file, 'w') as out_f:
        while True:
            read_pair_lines = list(itertools.islice(in_f, lines_per_read * 2))
            if not read_pair_lines:
                # reached end of file
                break

            read_1_lines = read_pair_lines[:lines_per_read]
            read_2_lines = read_pair_lines[lines_per_read:]

            if random.random() > 0.5:
                read_1_lines[0] = read_1_lines[0][:-2] + "2\n"
                read_2_lines[0] = read_2_lines[0][:-2] + "1\n"
                _write_read(out_f, read_2_lines)
                _write_read(out_f, read_1_lines)
            else:
                _write_read(out_f, read_1_lines)
                _write_read(out_f, read_2_lines)


def randomise_read_strands(args):
    # Read in and validate command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        docstring, argv=args,
        version="randomise_read_strands v" + __version__)

    _validate_command_line_options(options)

    # Randomly reassign pairs of reads such that the first read corresponds to
    # the antisense strand
    _randomise_read_strands(options[READS_FILE], options[OUT_PREFIX])
