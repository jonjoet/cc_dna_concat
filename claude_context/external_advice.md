# initial concerns from external source

- Zip groups can silently diverge in length: _validate_zip_groups in src/dna_concat/models.py only checks uniformity within each group, but assemble() in src/dna_concat/assembly.py indexes every zip slot with the same i. If two groups have lengths 3 and 4 (and no product slots), iteration runs for total=len(zip_slots[0]) and raises an IndexError once the shorter group is exhausted. Add a cross-group length check (or iterate per group) before assembly to fail fast with a clear message.
- assemble() materializes the full Cartesian product (product_combos = list(itertools.product(...))) in src/dna_concat/assembly.py, which explodes memory for realistic combinatorial libraries (e.g., 6 slots × 20 parts each → 64M tuples). Consider streaming via nested loops/generator, or at least add a guardrail/estimate to warn users before allocating the entire product.

# Summary of conversation with external source


- **ZIP validation:** `_validate_zip_groups()` enforces equal lengths within each zip group, but `assemble()` currently assumes a single global length and forces it to match product-count when both behaviors coexist.
- **Group semantics clarification:** Intended design is per-group zipping (independent lengths) combined multiplicatively via product; current implementation doesn’t honor group distinctions.
- **Observed behavior:** All-ZIP assemblies succeed, but mixing product and ZIP slots fails unless counts match, because `assemble()` inspects only the first ZIP slot and applies that length globally.
- **Next steps:** Refactor `assemble()` to bucket slots by `zip_group`, drop the global length tie-in, and iterate per group to align runtime behavior with intended design.