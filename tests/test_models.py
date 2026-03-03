"""Tests for dna_concat.models."""

import pytest

from dna_concat.models import (
    AssemblySpec,
    Construct,
    NamedSequence,
    Slot,
    SlotBehavior,
)


def test_fixed_slot_requires_one_sequence():
    with pytest.raises(ValueError, match="exactly 1 sequence"):
        Slot(
            name="test",
            behavior=SlotBehavior.FIXED,
            sequences=[
                NamedSequence("a", "AAA"),
                NamedSequence("b", "BBB"),
            ],
        )


def test_fixed_slot_with_one_sequence():
    slot = Slot(
        name="test",
        behavior=SlotBehavior.FIXED,
        sequences=[NamedSequence("a", "AAA")],
    )
    assert slot.sequences[0].sequence == "AAA"


def test_slot_requires_at_least_one_sequence():
    with pytest.raises(ValueError, match="at least one sequence"):
        Slot(name="test", behavior=SlotBehavior.ZIP, sequences=[])


def test_zip_slot_mismatch():
    with pytest.raises(ValueError, match="mismatched counts"):
        AssemblySpec(
            slots=[
                Slot(
                    name="a",
                    behavior=SlotBehavior.ZIP,
                    sequences=[NamedSequence("x", "A"), NamedSequence("y", "B")],
                ),
                Slot(
                    name="b",
                    behavior=SlotBehavior.ZIP,
                    sequences=[NamedSequence("x", "A")],
                ),
            ]
        )


def test_zip_slots_same_size_ok():
    spec = AssemblySpec(
        slots=[
            Slot(
                name="a",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("x", "A"), NamedSequence("y", "B")],
            ),
            Slot(
                name="b",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("p", "C"), NamedSequence("q", "D")],
            ),
        ]
    )
    assert len(spec.slots) == 2


def test_zip_slot_mismatch_allowed_with_trim():
    """Mismatched zip slot lengths are allowed when allow_zip_trim=True."""
    spec = AssemblySpec(
        slots=[
            Slot(
                name="a",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("x", "A"), NamedSequence("y", "B"), NamedSequence("z", "C")],
            ),
            Slot(
                name="b",
                behavior=SlotBehavior.ZIP,
                sequences=[NamedSequence("p", "D"), NamedSequence("q", "E")],
            ),
        ],
        allow_zip_trim=True,
    )
    assert len(spec.slots) == 2


def test_construct_full_sequence():
    c = Construct(
        assignments={
            "a": NamedSequence("x", "AAA"),
            "b": NamedSequence("y", "TTT"),
        },
        slot_order=["a", "b"],
    )
    assert c.full_sequence == "AAATTT"


def test_construct_format_name():
    c = Construct(
        assignments={
            "orf": NamedSequence("GFP", "ATG"),
            "barcode": NamedSequence("BC01", "AACC"),
        },
        slot_order=["orf", "barcode"],
    )
    assert c.format_name("{orf}_{barcode}") == "GFP_BC01"


def test_auto_name_template():
    spec = AssemblySpec(
        slots=[
            Slot(name="flank", behavior=SlotBehavior.FIXED, sequences=[NamedSequence("f", "GC")]),
            Slot(name="orf", behavior=SlotBehavior.ZIP, sequences=[NamedSequence("x", "A")]),
            Slot(name="bc", behavior=SlotBehavior.PRODUCT, sequences=[NamedSequence("y", "B")]),
        ]
    )
    assert spec.name_template == "{orf}_{bc}"


def test_auto_name_template_all_fixed():
    spec = AssemblySpec(
        slots=[
            Slot(name="a", behavior=SlotBehavior.FIXED, sequences=[NamedSequence("x", "A")]),
        ]
    )
    assert spec.name_template == "construct"
