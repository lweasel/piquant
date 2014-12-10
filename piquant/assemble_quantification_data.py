#!/usr/bin/env python

# pylint: disable=E1103

"""
Usage:
    assemble_quantification_data [{log_option_spec}] --method=<quantification-method> --out=<output-file> <pro-file> <transcript-count-file> <unique-sequence-file>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
-m --method=<quant-method>
    Method used to quantify transcript abundances.
-o <output-file> --out=<output-file>
    Output file for real and calculated TPMs.
<pro-file>
    Flux Simulator gene expression profile file.
<transcript-count-file>
    File containing per-gene transcript counts.
<unique-sequence-file>
    File containing unique sequence lengths per-transcript.

assemble_quantification_data assembles data required to assess the accuracy of
transcript abundance estimates produced in a single quantification run, then
writes these data to an output CSV file.
"""

import pandas as pd

from . import flux_simulator as fs
from . import options as opt
from . import quantifiers as qs
from . import tpms
from .__init__ import __version__

from docopt import docopt
from schema import SchemaError

QUANT_METHOD = "--method"
OUT_FILE = "--out"
PRO_FILE = "<pro-file>"
COUNT_FILE = "<transcript-count-file>"
UNIQUE_SEQ_FILE = "<unique-sequence-file>"

SORTED_PREFIX = "sorted"
SORTED_BAM_FILE = SORTED_PREFIX + ".bam"


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_file_option(
            options[PRO_FILE], "Could not open expression profile file")
        opt.validate_file_option(
            options[COUNT_FILE], "Could not open transcript count file")
        opt.validate_file_option(
            options[UNIQUE_SEQ_FILE],
            "Could not open unique sequence lengths file")
        options[QUANT_METHOD] = opt.validate_dict_option(
            options[QUANT_METHOD], qs.get_quantification_methods(),
            "Unknown quantification method")
    except SchemaError as exc:
        exit(exc.code)


def _read_expression_profiles(pro_file):
    profiles = fs.read_expression_profiles(pro_file)
    profiles[tpms.REAL_TPM] = \
        profiles[fs.PRO_FILE_FRAC_COL].map(lambda x: 1000000 * x)
    return profiles


def _read_transcript_abundances(quantifier, profiles):
    profiles[tpms.CALCULATED_TPM] = profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].\
        map(quantifier.get_transcript_abundance)


def _read_transcript_counts(count_file, profiles):
    transcript_counts = pd.read_csv(count_file, index_col=tpms.TRANSCRIPT)

    def set_transcript_count_and_gene(t_id):
        tc_row = transcript_counts.ix[t_id]  # pylint: disable=E1103
        return pd.Series({
            tpms.GENE: tc_row[tpms.GENE],
            tpms.TRANSCRIPT_COUNT: tc_row[tpms.TRANSCRIPT_COUNT]
        })

    return profiles.merge(profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].apply(
        set_transcript_count_and_gene),
        left_index=True, right_index=True)


def _read_unique_sequence_lengths(unique_seq_file, profiles):
    unique_seqs = pd.read_csv(unique_seq_file, index_col=tpms.TRANSCRIPT)

    set_unique_length = lambda t_id: \
        unique_seqs.ix[t_id][tpms.UNIQUE_SEQ_LENGTH] \
        if t_id in unique_seqs.index else 0  # pylint: disable=E1103

    profiles[tpms.UNIQUE_SEQ_LENGTH] = \
        profiles[fs.PRO_FILE_TRANSCRIPT_ID_COL].map(set_unique_length)


def _write_quantification_data(out_file, profiles):
    profiles.rename(
        columns={
            fs.PRO_FILE_TRANSCRIPT_ID_COL: tpms.TRANSCRIPT,
            fs.PRO_FILE_LENGTH_COL: tpms.LENGTH
        },
        inplace=True)

    profiles.to_csv(
        out_file, index=False,
        cols=[tpms.TRANSCRIPT, tpms.GENE, tpms.LENGTH, tpms.UNIQUE_SEQ_LENGTH,
              tpms.TRANSCRIPT_COUNT, tpms.REAL_TPM, tpms.CALCULATED_TPM])


def _assemble_and_write_quantification_data(logger, options):
    # Read in the expression profile file, and calculate the true TPM
    # for each transcript
    logger.info("Reading expression profiles...")
    profiles = _read_expression_profiles(options[PRO_FILE])

    # Read calculated TPM values for each transcript produced by a particular
    # quantification method
    logger.info("Reading calculated TPMs...")
    _read_transcript_abundances(options[QUANT_METHOD], profiles)

    # Read per-gene transcript counts
    logger.info("Reading per-gene transcript counts...")
    profiles = _read_transcript_counts(options[COUNT_FILE], profiles)

    # Read unique sequence lengths per-transcript
    logger.info("Reading unique sequence lengths per-transcript")
    _read_unique_sequence_lengths(options[UNIQUE_SEQ_FILE], profiles)

    # Write TPMs and other relevant data to output file
    logger.info("Writing TPMs to file {out}".format(out=options[OUT_FILE]))
    _write_quantification_data(options[OUT_FILE], profiles)


def assemble_quantification_data(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt(
        docstring, argv=args,
        version="assemble_quantification_data v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Assemble and write quantification data
    _assemble_and_write_quantification_data(logger, options)
