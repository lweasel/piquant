#!/usr/bin/env python

"""Usage:
    calculate_unique_transcript_sequence [{log_option_spec}] <gtf-file>

{help_option_spec}                 {help_option_description}
{ver_option_spec}              {ver_option_description}
{log_option_spec}   {log_option_description}
<gtf-file>                GTF file containing genes and transcripts.
"""

import docopt
import gtf
import options as opt
import schema

from collections import defaultdict, namedtuple
from __init__ import __version__

GTF_FILE = "<gtf-file>"


Exon = namedtuple('Exon', ['sequence', 'start', 'end', 'strand'])

ExonAndTranscript = namedtuple('ExonAndTranscript', ['exon', 'transcript'])


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
        Exon(str(ei[0]), int(ei[1]), int(ei[2]), str(ei[3])), ei[-1])
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
    seq_to_unique_exon_transcript_map = defaultdict(list)

    for e_and_t in unique_exon_transcript_pairs:
            seq_to_unique_exon_transcript_map[e_and_t.exon.sequence].\
                append(e_and_t)

    return seq_to_unique_exon_transcript_map


def _get_exon_bases_map(unique_exon_transcript_pairs):
    return {e_and_t.exon:
            set(range(e_and_t.exon.start, e_and_t.exon.end + 1))
            for e_and_t in unique_exon_transcript_pairs}


def _get_exon_bases(index, exon_and_transcript_list, exon_bases_map):
    exon = exon_and_transcript_list[index].exon
    exon_bases = exon_bases_map[exon]
    return exon, exon_bases


def _calculate_exon_difference(exon1, exon1_bases, exon2, exon2_bases):
    if exon2.end < exon1.start or exon2.start > exon1.end:
        # No overlap between exons
        return exon1_bases

    if exon2.start <= exon1.start and exon2.end >= exon2.end:
        # Exon 2 completely covers exon 1
        return set()

    return exon1_bases - exon2_bases


def _get_unique_transcript_lengths(
        unique_exon_transcript_pairs,
        seq_to_unique_exon_transcript_map, logger):

    exon_bases_map = _get_exon_bases_map(unique_exon_transcript_pairs)

    transcript_lengths = defaultdict(int)

    for seq, e_and_t_list in seq_to_unique_exon_transcript_map.items():
        logger.info("...processing chromosome '{seq}'".format(seq=seq))
        indices = range(len(e_and_t_list))

        for i in indices:
            e1, e1_bases = _get_exon_bases(i, e_and_t_list, exon_bases_map)

            for j in indices:
                if i == j:
                    continue

                e2, e2_bases = _get_exon_bases(j, e_and_t_list, exon_bases_map)

                e1_bases = _calculate_exon_difference(
                    e1, e1_bases, e2, e2_bases)
                if len(e1_bases) == 0:
                    break

            transcript = e_and_t_list[i].transcript
            transcript_lengths[transcript] += len(e1_bases)

    return transcript_lengths


def _output_unique_transcript_lengths(transcript_lengths):
    print("transcript,unique-length")
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
    seq_to_unique_exon_transcript_map = \
        _get_unique_exon_per_chromosome_map(unique_exon_transcript_pairs)

    # For the remaining exons, some will overlap - we want now calculate how
    # much of each exon is unique, and thus sum, per-transcript, the number of
    # bases unique to that transcript.
    logger.info("Removing overlaps between exons...")
    transcript_lengths = _get_unique_transcript_lengths(
        unique_exon_transcript_pairs,
        seq_to_unique_exon_transcript_map, logger)

    # Write the unique number of bases per-transcript to the specified output
    # file.
    logger.info("Writing unique lengths for {n} transcripts.".
                format(n=len(transcript_lengths)))
    _output_unique_transcript_lengths(transcript_lengths)


if __name__ == "__main__":
    # Read in command-line options
    __doc__ = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(
        __doc__, version="calculate_unique_transcript_sequence v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    # Calculate and output number of unique bases per transcript
    _calculate_unique_transcript_sequence(logger, options)
