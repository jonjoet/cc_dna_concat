"""Core combination algorithm for dna_concat."""

from __future__ import annotations

import itertools

from dna_concat.models import AssemblySpec, Construct, NamedSequence, Slot, SlotBehavior


def assemble(spec: AssemblySpec) -> list[Construct]:
    """Generate all constructs from an AssemblySpec."""
    fixed_slots = [s for s in spec.slots if s.behavior is SlotBehavior.FIXED]
    product_slots = [s for s in spec.slots if s.behavior is SlotBehavior.PRODUCT]
    zip_slots = [s for s in spec.slots if s.behavior is SlotBehavior.ZIP]

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

        constructs.append(
            Construct(assignments=assignments, slot_order=spec.slot_order)
        )

    return constructs
