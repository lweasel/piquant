def get_fragment_bounds(read_identifier):
    """Return region bounds for a Flux Simulator read identifier.

    A Flux Simulator read identifier takes a form like:

    17:23917458-23919441W:ENSMUST00000115428:3404:610:576:768:S

    Here "17" is the region, "23917458" and 23919441 are the start and end
    positions of the originating transcript (both in 1-based co-ordinates),
    "610" is the length of the processed transcript, "576" is the start
    position of the fragment relative to the transcript and "768" is the
    end position of the fragment relative to the transcript.

    Along with the region of the originating transcript, this function
    will return (in 0-based, half-open co-ordinates) the start and end
    positions of the originating transcript, unless:
    i)  The fragment start position relative to the transcript is negative,
        in which case the bounds will be extended at the start to accommodate
        this, or
    ii) The fragment end position relative to the transcript is beyond the
        end of the transcript (i.e. end pos >= len), in which case the bounds
        will be extended at the end to accommodate this.
    """
    rid_elems = read_identifier.split(":")
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
