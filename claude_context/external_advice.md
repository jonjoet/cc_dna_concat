# Design Review Notes

Issues identified during code review, retained here for context on design decisions.

## Zip length validation

`AssemblySpec._validate_zip_slots()` enforces equal lengths across zip slots unless `allow_zip_trim` is enabled. `assemble()` pairs zip entries with product combinations 1:1 when both behaviors coexist. Without trimming, their counts must match. With trimming, each zip slot must have at least the product count; extra entries are ignored. If there are no product slots, trimming uses the shortest zip slot's length. These are intentional validation constraints.

## Memory considerations

`assemble()` materializes the full Cartesian product (`product_combos = list(itertools.product(...))`) which could be expensive for very large combinatorial libraries (e.g., 6 slots x 20 parts = 64M tuples). A streaming/generator approach would address this, but it hasn't been needed for realistic use cases so far.

## Zip group semantics

The original design considered per-group zipping (independent lengths combined multiplicatively), but the current implementation uses a simpler model: all zip slots share one global length. This covers the primary use cases without the complexity of group-level semantics. See `hierarchical_assembly.md` for the future design that would generalize this.
