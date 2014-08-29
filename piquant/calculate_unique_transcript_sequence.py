#!/usr/bin/env python

"""Usage:
    calculate_unique_transcript_sequence [{log_option_spec}] <gtf-file>

-h --help                 Show this message.
-v --version              Show version.
{log_option_spec}   {log_option_description}
<gtf-file>                GTF file containing genes and transcripts.
"""

import docopt
import gtf
import options as opt
import schema

from collections import defaultdict, namedtuple

GTF_FILE = "<gtf-file>"

# Read in command-line options
__doc__ = opt.substitute_into_usage(__doc__)
options = docopt.docopt(
    __doc__, version="calculate_unique_transcript_sequence v0.1")

# Validate command-line options
try:
    opt.validate_log_level(options)
    opt.validate_file_option(options[GTF_FILE], "Could not open GTF file")
except schema.SchemaError as exc:
    exit(exc.code)

# Set up logger
logger = opt.get_logger_for_options(options)

# Read exon lines information from GTF file and extract transcript ID from GTF
# attributes.
logger.info("Reading GTF file {f}".format(f=options[GTF_FILE]))

gtf_info = gtf.read_gtf_file(options[GTF_FILE])
exon_info = gtf_info[gtf_info[gtf.FEATURE_COL] == gtf.EXON_FEATURE]

exon_info[gtf.TRANSCRIPT_ID_ATTRIBUTE] = exon_info[gtf.ATTRIBUTES_COL].map(
    lambda x: gtf.get_attributes_dict(x)[gtf.TRANSCRIPT_ID_ATTRIBUTE])

# Extract pairs of exons and transcripts IDs from GTF exon lines
column_dict = exon_info.to_dict(outtype="list")
exon_info_list = zip(
    column_dict[gtf.SEQUENCE_COL], column_dict[gtf.START_COL],
    column_dict[gtf.END_COL], column_dict[gtf.STRAND_COL],
    column_dict[gtf.TRANSCRIPT_ID_ATTRIBUTE])

Exon = namedtuple('Exon', ['sequence', 'start', 'end', 'strand'])
ExonAndTranscript = namedtuple('ExonAndTranscript', ['exon', 'transcript'])
exon_transcript_pairs = \
    [ExonAndTranscript(
        Exon(str(ei[0]), int(ei[1]), int(ei[2]), str(ei[3])), ei[-1])
        for ei in exon_info_list]

logger.info("Read {c} exon + transcripts pairs.".
            format(c=len(exon_transcript_pairs)))

# Construct a map from exons to the transcripts that contain that exon. Then
# discard any exons which are shared between transcripts.
exon_transcript_map = defaultdict(list)
for e_and_t in exon_transcript_pairs:
    exon_transcript_map[e_and_t.exon].append(e_and_t.transcript)

unique_exon_transcript_map = \
    {k: v for k, v in exon_transcript_map.items() if len(v) == 1}

logger.info("Retained {c} exons unique to one transcript.".
            format(c=len(unique_exon_transcript_map)))

# For those exons which originate from only one transcript, split into groups
# per originating chromosome
unique_exon_transcript_list = \
    [ExonAndTranscript(exon, transcripts[0])
        for exon, transcripts in unique_exon_transcript_map.items()]

seq_to_unique_exon_transcript_map = defaultdict(list)

for e_and_t in unique_exon_transcript_list:
        seq_to_unique_exon_transcript_map[e_and_t.exon.sequence].\
            append(e_and_t)

# For the remaining exons, some will overlap - we want now calculate how much
# of each exon is unique, and thus sum, per-transcript, the number of bases
# unique to that transcript.
logger.info("Removing overlaps between exons...")

exon_bases = {e_and_t.exon:
              set(range(e_and_t.exon.start, e_and_t.exon.end + 1))
              for e_and_t in unique_exon_transcript_list}

transcript_lengths = defaultdict(int)

for seq, e_and_t_list in seq_to_unique_exon_transcript_map.items():
    logger.info("...processing chromosome '{seq}'".format(seq=seq))
    indices = range(len(e_and_t_list))
    for i in indices:
        exon1 = e_and_t_list[i].exon
        exon1_bases = exon_bases[exon1]

        for j in indices:
            if i == j:
                continue

            exon2 = e_and_t_list[j].exon
            if exon2.end < exon1.start or exon2.start > exon1.end:
                # No overlap between exons
                continue

            if exon2.start <= exon1.start and exon2.end >= exon2.end:
                # Exon 2 completely covers exon 1
                exon1_bases = None
                break

            exon1_bases = exon1_bases - exon_bases[exon2]
            if len(exon1_bases) == 0:
                break

        transcript = e_and_t_list[i].transcript
        transcript_lengths[transcript] += \
            0 if exon1_bases is None else len(exon1_bases)

# Write the unique number of bases per-transcript to the specified output file.
logger.info("Writing unique lengths for {n} transcripts.".
            format(n=len(transcript_lengths)))

print("transcript,unique-length")
for transcript, length in transcript_lengths.items():
    print("{t},{l}".format(t=transcript, l=length))
