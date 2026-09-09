#!/usr/bin/env python3
"""Short hash-based version-tag builder shared by UniProt reference and BFD
dataset provenance manifests (design doc: Data provenance & versioning).

The hash is always computed over the manifest with any existing
`version_tag` key removed first -- hashing a manifest that already embeds
its own tag would make the tag circular (its value would depend on itself).
"""
import hashlib
import json


def manifest_hash(manifest: dict) -> str:
    clean = {k: v for k, v in manifest.items() if k != "version_tag"}
    serialized = json.dumps(clean, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:7]


def build_version_tag(fields: dict, prefix: str, manifest: dict) -> str:
    """fields: ordered {label: value} pairs rendered as `label+value`
    segments (e.g. {"date": "20260909"} -> "v20260909"). The first field is
    rendered with a leading 'v', the rest bare -- matching the spec's
    example `bfd-v20260909-n11024-mmseqs95c90-a1b2c3d`."""
    segments = [prefix]
    for i, (label, value) in enumerate(fields.items()):
        if i == 0:
            segments.append(f"v{value}")
        elif len(label) == 1:
            segments.append(f"{label}{value}")
        else:
            segments.append(f"{value}")
    segments.append(manifest_hash(manifest))
    return "-".join(segments)
