"""
Minimal tests against the INSTALLED package.

With the src/ layout these cannot import the repository folder by accident:
they always import the installed `herd`, which is the point.

    pip install -e .
    pytest
"""

from herd.validation import clean_sequence, validate_records
from herd.fasta import parse_fasta_text
from herd.config import SEQ_LENGTH, CORE_HOSTS, SPECIES_ORDER
from herd.models import resolve_host


def test_non_canonical_residues_become_x():
    clean, warnings = clean_sequence("MKVBBBZZZ")
    assert "B" not in clean and "Z" not in clean
    assert clean.count("X") == 6
    assert "non-canonical" in warnings[0]


def test_truncation_to_1022():
    clean, warnings = clean_sequence("M" * 1315)
    assert len(clean) == SEQ_LENGTH
    assert warnings == [f"truncated from 1315 to {SEQ_LENGTH} residues"]


def test_empty_sequence_is_invalid():
    valid, invalid = validate_records([("a", "")])
    assert valid == []
    assert invalid[0].reason == "empty sequence"


def test_fasta_text_multiline_is_joined():
    ids, seqs = parse_fasta_text(">p1\nMKV\nLAA\n>p2\nMSE")
    assert ids == ["p1", "p2"]
    assert seqs == ["MKVLAA", "MSE"]


def test_bare_sequence_without_header():
    ids, seqs = parse_fasta_text("MKVLAA")
    assert ids == ["sequence_1"]
    assert seqs == ["MKVLAA"]


def test_host_aliases_resolve():
    assert resolve_host("Canis_familiaris") == "Canis_sp"
    assert resolve_host("Homo") == "Homo_sapiens"
    assert resolve_host("Tyrannosaurus") is None


def test_core_hosts_are_a_subset_of_the_heads():
    assert set(CORE_HOSTS) <= set(SPECIES_ORDER)
    assert len(SPECIES_ORDER) == 15
