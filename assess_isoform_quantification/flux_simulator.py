import pandas as pd

# TODO: numeric literals in get_fragment_bounds

# n.b. a Flux Simulator read identifier takes a form like:
#
# 17:23917458-23919441W:ENSMUST00000115428:3404:610:576:768:S
#
# Here "17" is the region, "23917458" and 23919441 are the start and end
# positions of the originating transcript (both in 1-based co-ordinates),
# "610" is the length of the processed transcript, "576" is the start
# position of the fragment relative to the transcript and "768" is the
# end position of the fragment relative to the transcript.

PRO_FILE_LOCUS_COL = 0
PRO_FILE_TRANSCRIPT_ID_COL = 1
PRO_FILE_CODING_COL = 2
PRO_FILE_LENGTH_COL = 3
PRO_FILE_FRAC_COL = 4
PRO_FILE_NUM_COL = 5


def read_expression_profiles(pro_file):
    profiles = pd.read_csv(pro_file, delim_whitespace=True, header=None)
    #profiles.set_index(PRO_FILE_TRANSCRIPT_ID_COL, inplace=True)
    return profiles


def _get_read_identifier_elems(read_identifier):
    return read_identifier.split(":")


def get_transcript_id(read_identifier):
    """Return originating transcript ID for a Flux Simulator read."""
    rid_elems = _get_read_identifier_elems(read_identifier)
    return rid_elems[2]


def get_transcript_length(read_identifier):
    """Return the processed length of the originating transcript of a read"""
    rid_elems = _get_read_identifier_elems(read_identifier)
    return int(rid_elems[4])


def get_transcript_bounds(read_identifier):
    """Return transcript bounds for a Flux Simulator read.

    Along with the region of the originating transcript, this function
    will return (in 0-based, half-open co-ordinates) the start and end
    positions of the read's originating transcript."""
    rid_elems = _get_read_identifier_elems(read_identifier)
    region = rid_elems[0]
    start_str, end_str = rid_elems[1][:-1].split("-")
    return region, int(start_str) - 1, int(end_str)


def get_fragment_bounds(read_identifier):
    """Return fragment bounds for a Flux Simulator read.

    This function behaves the same as get_transcript_bounds above, unless:
    i)  The fragment start position relative to the transcript is negative,
        in which case the bounds will be extended at the start to accommodate
        this, or
    ii) The fragment end position relative to the transcript is beyond the
        end of the transcript (i.e. end pos >= len), in which case the bounds
        will be extended at the end to accommodate this.
    """
    rid_elems = _get_read_identifier_elems(read_identifier)
    region = rid_elems[0]

    t_start_str, t_end_str = rid_elems[1][:-1].split("-")
    t_start = int(t_start_str) - 1
    t_end = int(t_end_str)
    t_len = int(rid_elems[4])

    f_start = int(rid_elems[5])
    if f_start < 0:
        t_start += f_start

    f_end = int(rid_elems[6])
    if f_end >= t_len:
        t_end += f_end - t_len + 1

    return region, t_start, t_end
