import assess_isoform_quantification.flux_simulator as fs


def _get_read_identifier(t_region=1, t_start=1000, t_end=2000,
                         t_len=500, f_start=0, f_end=100,
                         t_id="ENSMUST00000115428"):
    return "{r}:{ts}-{te}W:{tid}:1:{tl}:{fs}:{fe}:S".format(
        r=t_region, ts=t_start, te=t_end,
        tl=t_len, fs=f_start, fe=f_end,
        tid=t_id)


def test_get_transcript_id_returns_correct_identifier():
    t_id = "abc"
    rid = _get_read_identifier(t_id=t_id)
    identifier = fs.get_transcript_id(rid)
    assert identifier == t_id


def test_get_transcript_length_returns_correct_value():
    t_len = 123
    rid = _get_read_identifier(t_len=t_len)
    length = fs.get_transcript_length(rid)
    assert length == t_len


def test_get_transcript_bounds_returns_correct_region():
    region = "X"
    rid = _get_read_identifier(t_region=region)
    t_region, t_start, t_end = fs.get_transcript_bounds(rid)
    assert t_region == region


def test_get_transcript_bounds_returns_correct_start():
    start = 1000
    rid = _get_read_identifier(t_start=start)
    t_region, t_start, t_end = fs.get_transcript_bounds(rid)
    assert t_start == start - 1


def test_get_transcript_bounds_returns_correct_end():
    end = 2000
    rid = _get_read_identifier(t_end=end)
    t_region, t_start, t_end = fs.get_transcript_bounds(rid)
    assert t_end == end


def test_get_fragment_bounds_returns_correct_region():
    region = 17
    rid = _get_read_identifier(t_region=region)
    b_region, b_start, b_end = fs.get_fragment_bounds(rid)
    assert b_region == str(region)


def test_get_fragment_bounds_returns_correct_start():
    start = 1000
    rid = _get_read_identifier(t_start=start)
    b_region, b_start, b_end = fs.get_fragment_bounds(rid)
    assert b_start == start - 1


def test_get_fragment_bounds_returns_correct_end():
    end = 2000
    rid = _get_read_identifier(t_end=end)
    b_region, b_start, b_end = fs.get_fragment_bounds(rid)
    assert b_end == end


def test_get_fragment_bounds_returns_correct_start_when_fragment_starts_before_transcript():
    t_start = 1000
    f_start = -100
    rid = _get_read_identifier(t_start=t_start, f_start=f_start)
    b_region, b_start, b_end = fs.get_fragment_bounds(rid)
    assert b_start == t_start - 1 + f_start


def test_get_fragment_bounds_returns_correct_start_when_fragment_ends_after_transcript():
    t_end = 2000
    t_len = 500
    f_end = 550
    rid = _get_read_identifier(t_end=t_end, t_len=t_len, f_end=f_end)
    b_region, b_start, b_end = fs.get_fragment_bounds(rid)
    assert b_end == t_end + (f_end - t_len + 1)
