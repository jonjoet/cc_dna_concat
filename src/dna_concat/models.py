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
    zip_group: str = "default"

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

    def __post_init__(self) -> None:
        self._validate_zip_groups()
        if self.name_template is None:
            self.name_template = self._auto_name_template()

    def _validate_zip_groups(self) -> None:
        zip_groups: dict[str, list[Slot]] = {}
        for slot in self.slots:
            if slot.behavior is SlotBehavior.ZIP:
                zip_groups.setdefault(slot.zip_group, []).append(slot)
        for group_name, group_slots in zip_groups.items():
            sizes = {len(s.sequences) for s in group_slots}
            if len(sizes) > 1:
                detail = ", ".join(
                    f"'{s.name}'={len(s.sequences)}" for s in group_slots
                )
                raise ValueError(
                    f"Zip group '{group_name}' has mismatched counts: {detail}"
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
