"""
herd/validation.py — sequence cleaning and validation, all on CPU.

Runs BEFORE anything touches the GPU. Cleaning lives here rather than inside
`embed_sequences` so that a bad sequence can be reported on its own instead of
killing the whole batch, and so `esm.embed_sequences` only ever receives
sequences that are already clean.
"""

import re
from dataclasses import dataclass, field
from typing import List, Tuple

from .config import VALID_AA, UNKNOWN_AA, SEQ_LENGTH

# Compiled once: anything outside the 20 canonical residues.
_NON_CANONICAL = re.compile(f"[^{VALID_AA}]")


def clean_sequence(seq: str) -> Tuple[str, List[str]]:
    """
    Non-canonical residues -> 'X', then truncate to SEQ_LENGTH. Exactly the
    recipe used to build the training embeddings; do not change it.
    Returns (clean_sequence, warnings).
    """
    warnings: List[str] = []
    up = seq.upper()

    n_bad = len(_NON_CANONICAL.findall(up))
    if n_bad:
        warnings.append(f"{n_bad} non-canonical residue(s) -> '{UNKNOWN_AA}'")
    clean = _NON_CANONICAL.sub(UNKNOWN_AA, up)

    if len(clean) > SEQ_LENGTH:
        warnings.append(f"truncated from {len(clean)} to {SEQ_LENGTH} residues")
        clean = clean[:SEQ_LENGTH]

    return clean, warnings


@dataclass
class ValidatedRecord:
    id: str
    seq: str                       # cleaned, ready for ESM
    orig_length: int
    warnings: List[str] = field(default_factory=list)


@dataclass
class InvalidRecord:
    id: str
    reason: str


def validate_records(
    records: List[Tuple[str, str]]
) -> Tuple[List[ValidatedRecord], List[InvalidRecord]]:
    """
    Validates and cleans a list of (id, seq) -> (valid, invalid). Invalid
    records are returned separately instead of aborting the whole run.

    No minimum length is enforced: the 40 aa floor was a dataset filter, not
    an inference one.
    """
    valid: List[ValidatedRecord] = []
    invalid: List[InvalidRecord] = []

    for sid, seq in records:
        raw = (seq or "").strip()
        if not raw:
            invalid.append(InvalidRecord(sid, "empty sequence"))
            continue
        clean, warns = clean_sequence(raw)
        if not clean:
            invalid.append(InvalidRecord(sid, "no valid residues"))
            continue
        valid.append(ValidatedRecord(id=sid, seq=clean, orig_length=len(raw), warnings=warns))

    return valid, invalid