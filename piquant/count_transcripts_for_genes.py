"""
Usage:
    count_transcripts_for_genes [{log_option_spec}] <gtf-file>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<gtf-file>
    GTF file containing genes and transcripts.

Calculate the number of transcripts shared by the gene of origin for each
transcript in the specified GTF file.
"""

from __future__ import print_function

import collections
import docopt
import schema

from . import gtf
from . import options as opt
from . import tpms
from .__init__ import __version__

GTF_FILE = "<gtf-file>"


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)
        opt.validate_file_option(options[GTF_FILE], "Could not open GTF file")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_transcript_to_gene_map(gtf_info):
    transcript_to_gene_map = {}

    for dummy, row in gtf_info.iterrows():
        attributes_dict = gtf.get_attributes_dict(row[gtf.ATTRIBUTES_COL])
        transcript = attributes_dict[gtf.TRANSCRIPT_ID_ATTRIBUTE]
        if transcript not in transcript_to_gene_map:
            transcript_to_gene_map[transcript] = \
                attributes_dict[gtf.GENE_ID_ATTRIBUTE]

    return transcript_to_gene_map


def _get_gene_transcript_counts(transcript_to_gene_map):
    transcript_counts = collections.Counter()
    for gene in transcript_to_gene_map.values():
        transcript_counts[gene] += 1

    return transcript_counts


def _output_transcript_counts_for_genes(
        transcript_to_gene_map, transcript_counts):

    print(",".join([tpms.TRANSCRIPT, tpms.GENE, tpms.TRANSCRIPT_COUNT]))
    for transcript, gene in transcript_to_gene_map.items():
        print("{t},{g},{c}".format(
            t=transcript, g=gene, c=transcript_counts[gene]))


def _count_transcripts_for_genes(logger, options):
    logger.info("Reading GTF file {f}...".format(f=options[GTF_FILE]))
    gtf_info = gtf.read_gtf_file(options[GTF_FILE])

    logger.info("Extracting transcript to gene mappings...")
    transcript_to_gene_map = _get_transcript_to_gene_map(gtf_info)

    logger.info("Calculating transcript counts for genes...")
    transcript_counts = _get_gene_transcript_counts(
        transcript_to_gene_map)

    logger.info("Printing transcript counts for genes...")
    _output_transcript_counts_for_genes(
        transcript_to_gene_map, transcript_counts)


def count_transcripts_for_genes(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        docstring, argv=args,
        version="count_transcripts_for_genes.py v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Calculate and print per-gene transcript counts to standard out
    _count_transcripts_for_genes(logger, options)
