from parse_uniprot_relnotes import parse_relnotes


def test_parses_real_relnotes_format():
    text = open("tests/data/accumulation/relnotes_excerpt.txt").read()
    info = parse_relnotes(text)
    assert info["release"] == "2026_03"
    assert info["release_date"] == "2026-09-02"


def test_raises_on_unparseable_text():
    import pytest
    with pytest.raises(ValueError):
        parse_relnotes("nothing recognizable here")
