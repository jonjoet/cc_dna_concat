"""Tests for dna_concat.io."""

import tempfile
from pathlib import Path

from dna_concat.io import (
    parse_inline,
    read_csv,
    read_fasta,
    write_csv,
    write_fasta,
)
from dna_concat.models import Construct, NamedSequence

FIXTURES = Path(__file__).parent / "fixtures"


def test_read_fasta():
    seqs = read_fasta(FIXTURES / "orfs.fasta")
    assert len(seqs) == 3
    assert seqs[0].name == "ORF1"
    assert seqs[0].sequence == "ATGAAAGGG"


def test_read_csv():
    seqs = read_csv(
        FIXTURES / "barcodes.csv",
        name_column="barcode_id",
        sequence_column="barcode_seq",
    )
    assert len(seqs) == 3
    assert seqs[0].name == "BC01"
    assert seqs[0].sequence == "AACCGGTT"


def test_parse_inline():
    entries = [
        {"name": "seq1", "sequence": "AAAA"},
        {"name": "seq2", "sequence": "TTTT"},
    ]
    seqs = parse_inline(entries)
    assert len(seqs) == 2
    assert seqs[1].name == "seq2"


def test_write_and_read_fasta():
    constructs = [
        Construct(
            assignments={"a": NamedSequence("x", "AAA"), "b": NamedSequence("y", "TTT")},
            slot_order=["a", "b"],
        ),
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "out.fasta"
        write_fasta(constructs, path, "{a}_{b}")
        seqs = read_fasta(path)
        assert len(seqs) == 1
        assert seqs[0].name == "x_y"
        assert seqs[0].sequence == "AAATTT"


def test_write_csv():
    constructs = [
        Construct(
            assignments={"orf": NamedSequence("GFP", "ATG"), "bc": NamedSequence("BC01", "AA")},
            slot_order=["orf", "bc"],
        ),
        Construct(
            assignments={"orf": NamedSequence("RFP", "ATG"), "bc": NamedSequence("BC02", "CC")},
            slot_order=["orf", "bc"],
        ),
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "out.csv"
        write_csv(constructs, path, "{orf}_{bc}")
        seqs = read_csv(path, name_column="name", sequence_column="full_sequence")
        assert len(seqs) == 2
        assert seqs[0].name == "GFP_BC01"
        assert seqs[0].sequence == "ATGAA"
