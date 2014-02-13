from assess_isoform_quantification.flux_simulator import get_fragment_bounds


def _get_read_identifier(t_region, t_start, t_end, t_len, f_start, f_end):
    return "{r}:{ts}-{te}W:ENSMUST00000115428:1:{tl}:{fs}:{fe}:S".format(
        r=t_region, ts=t_start, te=t_end, tl=t_len, fs=f_start, fe=f_end)


def test_get_fragment_bounds_returns_correct_region():
    region = 17
    rid = _get_read_identifier(region, 1000, 2000, 500, 0, 100)
    b_region, b_start, b_end = get_fragment_bounds(rid)
    assert b_region == str(region)


def test_get_fragment_bounds_returns_correct_start():
    start = 1000
    rid = _get_read_identifier(17, start, 2000, 500, 0, 100)
    b_region, b_start, b_end = get_fragment_bounds(rid)
    assert b_start == start - 1


def test_get_fragment_bounds_returns_correct_end():
    end = 2000
    rid = _get_read_identifier(17, 1000, end, 500, 0, 100)
    b_region, b_start, b_end = get_fragment_bounds(rid)
    assert b_end == end


def test_get_fragment_bounds_returns_correct_start_when_fragment_starts_before_transcript():
    t_start = 1000
    f_start = -100
    rid = _get_read_identifier(17, t_start, 2000, 500, f_start, 100)
    b_region, b_start, b_end = get_fragment_bounds(rid)
    assert b_start == t_start - 1 + f_start


def test_get_fragment_bounds_returns_correct_start_when_fragment_ends_after_transcript():
    t_end = 2000
    t_len = 500
    f_end = 550
    rid = _get_read_identifier(17, 1000, t_end, t_len, 400, f_end)
    b_region, b_start, b_end = get_fragment_bounds(rid)
    assert b_end == t_end + (f_end - t_len + 1)
