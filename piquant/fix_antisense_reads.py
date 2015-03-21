"""
Usage:
    fix_antisense_reads
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
    String to be prepended to input file name for output [default: sense]
<reads-file>
    FASTA/Q file containing single-end reads.

Reverse complement any reads in the input FASTA or FASTQ file that originate from
the antisense strand.
"""

import docopt
import os.path
import schema

from . import options as opt
from .__init__ import __version__

OUT_PREFIX = "--out-prefix"
READS_FILE = "<reads-file>"

COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N"
}


def _validate_command_line_options(options):
    try:
        opt.validate_file_option(
            options[READS_FILE], "Reads file should exist.")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_lines_per_read(reads_file):
    return 4 if reads_file.endswith("fastq") else 2


def _reverse_complement(sequence):
    return ("".join(COMPLEMENT[base] for base in sequence))[::-1]


def _fix_antisense_reads(reads_file, out_prefix):
    dirname = os.path.dirname(os.path.abspath(reads_file))
    basename = os.path.basename(reads_file)
    output_file = os.path.join(dirname, out_prefix + "." + basename)

    lines_per_read = _get_lines_per_read(reads_file)
    sense = True

    with open(reads_file, 'r') as in_f, open(output_file, 'w') as out_f:
        for line_no, line in enumerate(in_f):
            if line_no % lines_per_read == 0:
                sense = line[-2:-1] == "S"
                out_f.write(line[:-2] + "S\n")
            elif line_no % lines_per_read == 1:
                out_f.write(line if sense else
                            (_reverse_complement(line.upper()[:-1]) + "\n"))
            else:
                out_f.write(line)


def fix_antisense_reads(args):
    # Read in and validate command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        docstring, argv=args,
        version="fix_antisense_reads v" + __version__)

    _validate_command_line_options(options)

    # Set up logger
    #logger = opt.get_logger_for_options(options)

    # Transform input reads file by reverse complementing all antisense reads
    _fix_antisense_reads(options[READS_FILE], options[OUT_PREFIX])
