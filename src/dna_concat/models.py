"""Data models for dna_concat."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field


class SlotBehavior(enum.Enum):
    FIXED = "fixed"
    ZIP = "zip"
    PRODUCT = "product"
    VARIABLE_STUFFER = "variable_stuffer"


@dataclass
class NamedSequence:
    name: str
    sequence: str


@dataclass
class Slot:
    name: str
    behavior: SlotBehavior
    sequences: list[NamedSequence]
    # Variable-stuffer-only fields (must be None for every other behavior):
    counterpart: str | None = None
    truncate_side: str | None = None

    def __post_init__(self) -> None:
        if not self.sequences:
            raise ValueError(f"Slot '{self.name}' must have at least one sequence")
        if self.behavior is SlotBehavior.FIXED and len(self.sequences) != 1:
            raise ValueError(
                f"Fixed slot '{self.name}' must have exactly 1 sequence, "
                f"got {len(self.sequences)}"
            )
        if self.behavior is SlotBehavior.VARIABLE_STUFFER:
            if len(self.sequences) != 1:
                raise ValueError(
                    f"Variable stuffer slot '{self.name}' must have exactly 1 "
                    f"sequence, got {len(self.sequences)}"
                )
            if not self.counterpart:
                raise ValueError(
                    f"Variable stuffer slot '{self.name}' requires a "
                    f"'counterpart' slot name"
                )
            if self.truncate_side not in ("left", "right"):
                raise ValueError(
                    f"Variable stuffer slot '{self.name}' requires truncate_side "
                    f"'left' or 'right', got {self.truncate_side!r}"
                )
        elif self.counterpart is not None or self.truncate_side is not None:
            raise ValueError(
                f"Slot '{self.name}': 'counterpart' and 'truncate_side' are only "
                f"valid for variable_stuffer slots"
            )


@dataclass
class Construct:
    assignments: dict[str, NamedSequence]
    slot_order: list[str]

    @property
    def full_sequence(self) -> str:
        return "".join(self.assignments[name].sequence for name in self.slot_order)

    def format_name(self, template: str) -> str:
        values = {name: ns.name for name, ns in self.assignments.items()}
        return template.format(**values)


@dataclass
class AssemblySpec:
    slots: list[Slot]
    name_template: str | None = None
    allow_zip_trim: bool = False

    def __post_init__(self) -> None:
        self._validate_zip_slots()
        self._validate_stuffer_slots()
        if self.name_template is None:
            self.name_template = self._auto_name_template()

    def _validate_zip_slots(self) -> None:
        zip_slots = [s for s in self.slots if s.behavior is SlotBehavior.ZIP]
        if not zip_slots:
            return
        sizes = {len(s.sequences) for s in zip_slots}
        if self.allow_zip_trim:
            # With trimming, zip slots can differ — they'll be trimmed to the
            # minimum length (or to product count) during assembly.
            return
        if len(sizes) > 1:
            detail = ", ".join(
                f"'{s.name}'={len(s.sequences)}" for s in zip_slots
            )
            raise ValueError(
                f"Zip slots have mismatched counts: {detail}"
            )

    def _validate_stuffer_slots(self) -> None:
        stuffers = [
            s for s in self.slots if s.behavior is SlotBehavior.VARIABLE_STUFFER
        ]
        if not stuffers:
            return
        by_name = {s.name: s for s in self.slots}
        for stuffer in stuffers:
            counterpart = stuffer.counterpart
            if counterpart == stuffer.name:
                raise ValueError(
                    f"Variable stuffer slot '{stuffer.name}' cannot use itself as "
                    f"its counterpart"
                )
            if counterpart not in by_name:
                raise ValueError(
                    f"Variable stuffer slot '{stuffer.name}' references unknown "
                    f"counterpart slot '{counterpart}'"
                )
            if by_name[counterpart].behavior is SlotBehavior.VARIABLE_STUFFER:
                raise ValueError(
                    f"Variable stuffer slot '{stuffer.name}' cannot use another "
                    f"variable stuffer ('{counterpart}') as its counterpart"
                )

    def _auto_name_template(self) -> str:
        excluded = {SlotBehavior.FIXED, SlotBehavior.VARIABLE_STUFFER}
        varying = [s.name for s in self.slots if s.behavior not in excluded]
        if not varying:
            return "construct"
        return "_".join(f"{{{name}}}" for name in varying)

    @property
    def slot_order(self) -> list[str]:
        return [s.name for s in self.slots]
