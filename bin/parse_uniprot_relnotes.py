#!/usr/bin/env python3
"""Parse UniProt's relnotes.txt for the current release version and date.

Real format confirmed 2026-09-09 against
https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt:
    Header of the UniProt Knowledgebase Release 2026_03 (02-Sept-2026)
"""
import re
from datetime import datetime

_RELEASE_RE = re.compile(
    r"UniProt Knowledgebase Release (\d{4}_\d{2}) \((\d{1,2}-[A-Za-z]+-\d{4})\)"
)
_MONTHS = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
     "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"], start=1)}


def parse_relnotes(text: str) -> dict:
    match = _RELEASE_RE.search(text)
    if not match:
        raise ValueError("could not find a release line in relnotes.txt")
    release, date_str = match.groups()
    day, month_str, year = date_str.split("-")
    release_date = datetime(int(year), _MONTHS[month_str], int(day)).strftime("%Y-%m-%d")
    return {"release": release, "release_date": release_date}
