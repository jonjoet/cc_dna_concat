"""Tests for dna_concat.assembly."""

import pytest

from dna_concat.assembly import assemble
from dna_concat.models import AssemblySpec, NamedSequence, Slot, SlotBehavior


def test_all_fixed():
    spec = AssemblySpec(
        slots=[
            Slot(name="a", behavior=SlotBehavior.FIXED, sequences=[NamedSequence("x", "AAA")]),
            Slot(name="b", behavior=SlotBehavior.FIXED, sequences=[NamedSequence("y", "TTT")]),
        ]
    )
    constructs = assemble(spec)
    assert len(constructs) == 1
    assert constructs[0].full_sequence == "AAATTT"


def test_fixed_plus_zip():
    spec = AssemblySpec(
        slots=[
            Slot(name="flank", behavior=SlotBehavior.FIXED, sequences=[NamedSequence("f", "GC")]),
            Slot(
                name="orf",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("A", "AAA"), NamedSequence("B", "BBB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("X", "XXX"), NamedSequence("Y", "YYY")],
            ),
        ]
    )
    constructs = assemble(spec)
    assert len(constructs) == 2
    assert constructs[0].full_sequence == "GCAAAXXX"
    assert constructs[1].full_sequence == "GCBBBYYY"


def test_product_only():
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="orf",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("G", "GGG"), NamedSequence("R", "RRR"), NamedSequence("B", "BBB")],
            ),
        ]
    )
    constructs = assemble(spec)
    assert len(constructs) == 6  # 2 * 3
    seqs = [c.full_sequence for c in constructs]
    assert "AAGGG" in seqs
    assert "AARRR" in seqs
    assert "AABBB" in seqs
    assert "BBGGG" in seqs
    assert "BBRRR" in seqs
    assert "BBBBB" in seqs


def test_product_plus_zip():
    # 2 product slots * (2, 3) = 6 combos; zip must have 6 entries
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="orf",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("G", "GGG"), NamedSequence("R", "RRR"), NamedSequence("B", "BBB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[
                    NamedSequence(f"BC{i}", f"BC{i}") for i in range(6)
                ],
            ),
        ]
    )
    constructs = assemble(spec)
    assert len(constructs) == 6


def test_product_zip_mismatch():
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("BC1", "X"), NamedSequence("BC2", "Y"), NamedSequence("BC3", "Z")],
            ),
        ]
    )
    with pytest.raises(ValueError, match="Zip count.*must equal product count"):
        assemble(spec)


def test_construct_names():
    spec = AssemblySpec(
        slots=[
            Slot(
                name="orf",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("GFP", "ATG"), NamedSequence("RFP", "ATG")],
            ),
        ],
        name_template="{orf}_construct",
    )
    constructs = assemble(spec)
    names = [c.format_name(spec.name_template) for c in constructs]
    assert names == ["GFP_construct", "RFP_construct"]


# --- Trim behavior tests ---


def test_zip_only_trim_to_shortest():
    """Zip-only: slots with lengths 3 and 5, trim to 3."""
    spec = AssemblySpec(
        slots=[
            Slot(
                name="a",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence(f"a{i}", f"A{i}") for i in range(3)],
            ),
            Slot(
                name="b",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence(f"b{i}", f"B{i}") for i in range(5)],
            ),
        ],
        allow_zip_trim=True,
    )
    constructs = assemble(spec)
    assert len(constructs) == 3


def test_product_plus_zip_trim():
    """Product count 6, zip has 8 entries, trim=True → 6 constructs."""
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="orf",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("G", "GGG"), NamedSequence("R", "RRR"), NamedSequence("B", "BBB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence(f"BC{i}", f"BC{i}") for i in range(8)],
            ),
        ],
        allow_zip_trim=True,
    )
    constructs = assemble(spec)
    assert len(constructs) == 6


def test_product_plus_zip_trim_false_still_errors():
    """Product count 6, zip has 8 entries, trim=False → ValueError."""
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="orf",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("G", "GGG"), NamedSequence("R", "RRR"), NamedSequence("B", "BBB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence(f"BC{i}", f"BC{i}") for i in range(8)],
            ),
        ],
        allow_zip_trim=False,
    )
    with pytest.raises(ValueError, match="Zip count.*must equal product count"):
        assemble(spec)


def test_zip_trim_too_few_still_errors():
    """Zip count < product count even with trim → error."""
    spec = AssemblySpec(
        slots=[
            Slot(
                name="promoter",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("pA", "AA"), NamedSequence("pB", "BB")],
            ),
            Slot(
                name="orf",
                behavior=SlotBehavior.PRODUCT,
                sequences=[NamedSequence("G", "GGG"), NamedSequence("R", "RRR"), NamedSequence("B", "BBB")],
            ),
            Slot(
                name="bc",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence(f"BC{i}", f"BC{i}") for i in range(4)],
            ),
        ],
        allow_zip_trim=True,
    )
    with pytest.raises(ValueError, match="less than product count"):
        assemble(spec)
