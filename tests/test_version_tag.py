import json
from version_tag import build_version_tag, manifest_hash

def test_manifest_hash_excludes_version_tag_field():
    manifest_without_tag = {"date": "20260909", "n": 100}
    manifest_with_tag = dict(manifest_without_tag, version_tag="whatever")
    assert manifest_hash(manifest_without_tag) == manifest_hash(manifest_with_tag)

def test_manifest_hash_is_7_hex_chars():
    h = manifest_hash({"a": 1})
    assert len(h) == 7
    assert all(c in "0123456789abcdef" for c in h)

def test_manifest_hash_changes_with_content():
    assert manifest_hash({"a": 1}) != manifest_hash({"a": 2})

def test_build_version_tag_format():
    fields = {"date": "20260909", "n": "11024", "params": "mmseqs95c90"}
    tag = build_version_tag(fields, prefix="bfd", manifest={"date": "20260909", "n": 11024})
    parts = tag.split("-")
    assert parts[0] == "bfd"
    assert parts[1] == "v20260909"
    assert parts[2] == "n11024"
    assert parts[3] == "mmseqs95c90"
    assert len(parts[4]) == 7
