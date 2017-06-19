"""
Usage:
    fasta_to_fastq [{log_option_spec}] <fasta-file> <fastq-file>

Options:
    {help_option_spec}
        {help_option_description}
    {ver_option_spec}
        {ver_option_description}
    {log_option_spec}
        {log_option_description}
    <fasta-file>
        input FASTA file containing reads
    <fastq-file>
        output converted FASTQ file


Convert the given fasta file to the fastq file.
"""

import docopt
import schema

import options as opt
from __init__ import __version__


FASTA_FILE = "<fasta-file>"
FASTQ_FILE = "<fastq-file>"

def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)
        opt.validate_file_option(options[FASTA_FILE],"Could not open FASTA file")
        opt.validate_file_option(options[FASTQ_FILE],"Could not create FASTQ file",False)
    except schema.SchemaError as exc:
        exit("Exiting. " + exc.code)


def _write_fastq(outfile, header, sequence, quality):
    outfile.write(header + "\n"
            + sequence + "\n"
            + r"+" + "\n"
            + quality + "\n")


def _fasta_to_fastq(logger, options):
    logger.info(
            "Converting reads from '{fasta_file}' to '{fastq_file}'...".format(fasta_file=options[FASTA_FILE],fastq_file=options[FASTQ_FILE]))
    fasta_file = open(options[FASTA_FILE])
    fastq_file = open(options[FASTQ_FILE],"a")
    header = ""
    sequence = ""
    quality = ""
    for line in fasta_file:
        line = line.strip("\n")
        if line.startswith(r">"):
            if header != "":
                _write_fastq(fastq_file, header, sequence, quality)
            header = r"@" + line[1:]
            sequence = ""
            quality = ""
        else:
            sequence += line
            length = len(sequence)
            quality = length*"I"
    fasta_file.close()
    _write_fastq(fastq_file, header, sequence, quality)
    fastq_file.close()
    logger.info("Converting complete.")


def fasta_to_fastq(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
            docstring,argv=args,
            version="fasta_to_fastq v" + __version__)

    # Validate and process command-line options
    _validate_command_line_options(options)

    # set up logger
    logger = opt.get_logger_for_options(options)

    # Convert fasta file to fastq file
    _fasta_to_fastq(logger, options)
