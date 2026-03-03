"""Data models for dna_concat."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field


class SlotBehavior(enum.Enum):
    FIXED = "fixed"
    ZIP = "zip"
    PRODUCT = "product"


@dataclass
class NamedSequence:
    name: str
    sequence: str


@dataclass
class Slot:
    name: str
    behavior: SlotBehavior
    sequences: list[NamedSequence]

    def __post_init__(self) -> None:
        if not self.sequences:
            raise ValueError(f"Slot '{self.name}' must have at least one sequence")
        if self.behavior is SlotBehavior.FIXED and len(self.sequences) != 1:
            raise ValueError(
                f"Fixed slot '{self.name}' must have exactly 1 sequence, "
                f"got {len(self.sequences)}"
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

    def _auto_name_template(self) -> str:
        non_fixed = [
            s.name for s in self.slots if s.behavior is not SlotBehavior.FIXED
        ]
        if not non_fixed:
            return "construct"
        return "_".join(f"{{{name}}}" for name in non_fixed)

    @property
    def slot_order(self) -> list[str]:
        return [s.name for s in self.slots]
