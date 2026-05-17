# Design Review Notes

Issues identified during code review, retained here for context on design decisions.

## Zip length validation

`_validate_zip_groups()` enforces equal lengths within zip slots, but `assemble()` assumes a single global length and forces it to match product-count when both behaviors coexist. When mixing product and zip slots, assembly fails unless counts match — this is intentional (it's a validation constraint, not a bug).

## Memory considerations

`assemble()` materializes the full Cartesian product (`product_combos = list(itertools.product(...))`) which could be expensive for very large combinatorial libraries (e.g., 6 slots x 20 parts = 64M tuples). A streaming/generator approach would address this, but it hasn't been needed for realistic use cases so far.

## Zip group semantics

The original design considered per-group zipping (independent lengths combined multiplicatively), but the current implementation uses a simpler model: all zip slots share one global length. This covers the primary use cases without the complexity of group-level semantics. See `hierarchical_assembly.md` for the future design that would generalize this.
