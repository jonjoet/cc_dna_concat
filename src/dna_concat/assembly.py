"""Core combination algorithm for dna_concat."""

from __future__ import annotations

import itertools

from dna_concat.models import AssemblySpec, Construct, NamedSequence, Slot, SlotBehavior


def assemble(spec: AssemblySpec) -> list[Construct]:
    """Generate all constructs from an AssemblySpec."""
    fixed_slots = [s for s in spec.slots if s.behavior is SlotBehavior.FIXED]
    product_slots = [s for s in spec.slots if s.behavior is SlotBehavior.PRODUCT]
    zip_slots = [s for s in spec.slots if s.behavior is SlotBehavior.ZIP]
    stuffer_slots = [
        s for s in spec.slots if s.behavior is SlotBehavior.VARIABLE_STUFFER
    ]

    # Build fixed assignments (same for every construct)
    fixed_assignments = {s.name: s.sequences[0] for s in fixed_slots}

    # Phase 1: Cartesian product of product slots
    if product_slots:
        product_combos = list(
            itertools.product(*(s.sequences for s in product_slots))
        )
    else:
        product_combos = [()]  # single empty combo

    product_count = len(product_combos)

    # Phase 2: Determine zip count and validate/trim
    if zip_slots:
        zip_count = min(len(s.sequences) for s in zip_slots)
        if product_slots:
            if spec.allow_zip_trim:
                if zip_count < product_count:
                    raise ValueError(
                        f"Zip count ({zip_count}) is less than product count "
                        f"({product_count}); cannot trim to fill missing entries"
                    )
                zip_count = product_count
            else:
                if zip_count != product_count:
                    raise ValueError(
                        f"Zip count ({zip_count}) must equal product count "
                        f"({product_count}) when both are present"
                    )
        total = zip_count
    else:
        total = product_count

    # Phase 3: Build constructs
    constructs = []
    for i in range(total):
        assignments = dict(fixed_assignments)

        # Add product assignments
        if product_slots:
            for slot, seq in zip(product_slots, product_combos[i]):
                assignments[slot.name] = seq

        # Add zip assignments
        for slot in zip_slots:
            assignments[slot.name] = slot.sequences[i]

        # Resolve variable stuffers last: each depends on its counterpart's
        # length in THIS construct (counterparts are never stuffers, so they
        # are already assigned above).
        for slot in stuffer_slots:
            counterpart_len = len(assignments[slot.counterpart].sequence)
            assignments[slot.name] = _truncate_stuffer(slot, counterpart_len)

        constructs.append(
            Construct(assignments=assignments, slot_order=spec.slot_order)
        )

    return constructs


def _truncate_stuffer(slot: Slot, counterpart_len: int) -> NamedSequence:
    """Truncate a variable-stuffer slot so that it plus its counterpart equal
    the untruncated stuffer's length.

    ``truncate_side`` selects which end is trimmed, independently of where the
    counterpart sits.
    """
    full = slot.sequences[0]
    keep = len(full.sequence) - counterpart_len
    if keep < 0:
        raise ValueError(
            f"Variable stuffer slot '{slot.name}': counterpart length "
            f"({counterpart_len}) exceeds stuffer length "
            f"({len(full.sequence)}); cannot hold the combined length constant"
        )
    if slot.truncate_side == "right":
        truncated = full.sequence[:keep]  # trim the right end, keep the left
    else:  # "left"
        # trim the left end, keep the right; index by length so keep == 0 -> ""
        truncated = full.sequence[len(full.sequence) - keep:]
    return NamedSequence(name=full.name, sequence=truncated)
