#!/usr/bin/env python

"""Usage:
    calculate_unique_transcript_sequence [{log_option_spec}] <gtf-file>

{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<gtf-file>
    GTF file containing genes and transcripts.

Calculate the length of sequence in base pairs that is unique to each transcript
in the specified GTF file.
"""

from __future__ import print_function

import docopt
import schema

from . import gtf
from . import options as opt
from . import tpms
from .__init__ import __version__

from collections import defaultdict, namedtuple

GTF_FILE = "<gtf-file>"
NO_BASES = set()

Exon = namedtuple('Exon', ['sequence', 'start', 'end', 'strand', 'bases'])

ExonAndTranscript = namedtuple('ExonAndTranscript', ['exon', 'transcript'])


def _replace_bases(exon, new_bases):
    return Exon(exon.sequence,
                exon.start,
                exon.end,
                exon.strand,
                new_bases)


def _replace_exon(e_and_t, new_exon):
    return ExonAndTranscript(new_exon, e_and_t.transcript)


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)
        opt.validate_file_option(options[GTF_FILE], "Could not open GTF file")
    except schema.SchemaError as exc:
        exit(exc.code)


def _read_exon_info(gtf_file):
    gtf_info = gtf.read_gtf_file(gtf_file)
    exon_info = gtf_info[gtf_info[gtf.FEATURE_COL] == gtf.EXON_FEATURE]

    exon_info[gtf.TRANSCRIPT_ID_ATTRIBUTE] = exon_info[gtf.ATTRIBUTES_COL].map(
        lambda x: gtf.get_attributes_dict(x)[gtf.TRANSCRIPT_ID_ATTRIBUTE])

    return exon_info


def _get_exon_transcript_pairs(exon_info):
    column_dict = exon_info.to_dict(outtype="list")
    exon_info_list = zip(
        column_dict[gtf.SEQUENCE_COL], column_dict[gtf.START_COL],
        column_dict[gtf.END_COL], column_dict[gtf.STRAND_COL],
        column_dict[gtf.TRANSCRIPT_ID_ATTRIBUTE])

    return [ExonAndTranscript(
        Exon(str(ei[0]), int(ei[1]), int(ei[2]), str(ei[3]), None), ei[-1])
        for ei in exon_info_list]


def _get_unique_exon_transcript_pairs(exon_transcript_pairs):
    exon_transcript_map = defaultdict(list)
    for e_and_t in exon_transcript_pairs:
        exon_transcript_map[e_and_t.exon].append(e_and_t.transcript)

    unique_exon_transcript_map = \
        {k: v for k, v in exon_transcript_map.items() if len(v) == 1}

    unique_exon_transcript_list = \
        [ExonAndTranscript(exon, transcripts[0])
            for exon, transcripts in unique_exon_transcript_map.items()]

    return unique_exon_transcript_list


def _get_unique_exon_per_chromosome_map(unique_exon_transcript_pairs):
    seq_to_unique_exon_transcripts = defaultdict(list)

    for e_and_t in unique_exon_transcript_pairs:
        seq_to_unique_exon_transcripts[e_and_t.exon.sequence].\
            append(e_and_t)

    return seq_to_unique_exon_transcripts


def _get_unique_transcript_lengths(
        seq_to_unique_exon_transcripts, logger):

    transcript_lengths = defaultdict(int)

    for seq, e_and_t_list in seq_to_unique_exon_transcripts.items():
        logger.info("...processing {exons} exons for chromosome '{seq}'".
                    format(exons=len(e_and_t_list), seq=seq))

        indices = range(len(e_and_t_list))

        e_and_t_list = [
            _replace_exon(e_and_t, _replace_bases(
                e_and_t.exon,
                set(range(e_and_t.exon.start, e_and_t.exon.end + 1))))
            for e_and_t in e_and_t_list]

        for i in indices:
            exon1 = e_and_t_list[i].exon
            exon_bases = exon1.bases

            for j in indices:
                if i == j:
                    continue

                exon2 = e_and_t_list[j].exon
                if exon2.end < exon1.start or exon2.start > exon1.end:
                    # No overlap between exons
                    continue

                if exon2.start <= exon1.start and exon2.end >= exon2.end:
                    # Exon 2 completely covers exon 1
                    exon_bases = NO_BASES
                    break

                exon_bases = exon_bases - exon2.bases
                if len(exon_bases) == 0:
                    break

            transcript_lengths[e_and_t_list[i].transcript] += len(exon_bases)

    return transcript_lengths


def _output_unique_transcript_lengths(transcript_lengths):
    print(",".join([tpms.TRANSCRIPT, tpms.UNIQUE_SEQ_LENGTH]))
    for transcript, length in transcript_lengths.items():
        print("{t},{l}".format(t=transcript, l=length))


def _calculate_unique_transcript_sequence(logger, options):
    # Read exon lines information from GTF file and extract transcript ID from
    # GTF attributes.
    logger.info("Reading GTF file {f}".format(f=options[GTF_FILE]))
    exon_info = _read_exon_info(options[GTF_FILE])

    # Extract pairs of exons and transcripts IDs from GTF exon lines
    exon_transcript_pairs = _get_exon_transcript_pairs(exon_info)
    logger.info("Read {c} exon + transcripts pairs.".
                format(c=len(exon_transcript_pairs)))

    # Discard any exons which are shared between transcripts, leaving only
    # pairs of unique exons and transcripts
    unique_exon_transcript_pairs = \
        _get_unique_exon_transcript_pairs(exon_transcript_pairs)
    logger.info("Retained {c} exons unique to one transcript.".
                format(c=len(unique_exon_transcript_pairs)))

    # For those exons which originate from only one transcript, split into
    # groups per originating chromosome
    seq_to_unique_exon_transcripts = \
        _get_unique_exon_per_chromosome_map(unique_exon_transcript_pairs)

    # For the remaining exons, some will overlap - we want now calculate how
    # much of each exon is unique, and thus sum, per-transcript, the number of
    # bases unique to that transcript.
    logger.info("Removing overlaps between exons...")
    transcript_lengths = _get_unique_transcript_lengths(
        seq_to_unique_exon_transcripts, logger)

    # Write the unique number of bases per-transcript to the specified output
    # file.
    logger.info("Writing unique lengths for {n} transcripts.".
                format(n=len(transcript_lengths)))
    _output_unique_transcript_lengths(transcript_lengths)


def calculate_unique_transcript_sequence(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        docstring, argv=args,
        version="calculate_unique_transcript_sequence v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Calculate and output number of unique bases per transcript
    _calculate_unique_transcript_sequence(logger, options)
