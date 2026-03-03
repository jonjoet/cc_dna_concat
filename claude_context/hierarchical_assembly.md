# Hierarchical Assembly: Future Design

## Motivation

The current flat model assigns a behavior (`fixed`, `zip`, `product`) to each slot independently. This works well for common cases — combinatorial libraries with optional barcoding — but can't express scenarios where **subgroups of slots combine differently and then compose at a higher level**.

Example: "product of (promoter × ORF), product of (terminator × RBS), zip those two groups with a barcode." The flat model can't represent this because `product` and `zip` are global: all product slots form one Cartesian product, and all zip slots lock-step together.

## Why we're deferring

The flat model covers our primary use cases:

1. **Full combinatorial** — all product slots → Cartesian product
2. **Paired iteration** — all zip slots lock-step together
3. **Combinatorial + barcoding** — product expansion paired 1:1 with zip entries
4. **Trimming** — `allow_zip_trim` handles the case where zip slots have more entries than needed

Adding hierarchy before we have concrete user demand would be premature complexity. The flat model is easy to understand, validate, and debug.

## Proposed tree-based design

### Core idea

Replace flat per-slot behaviors with a **tree of combination nodes**. Each node is either a leaf (a slot) or an internal node that specifies how to combine its children (`product` or `zip`).

### YAML syntax

```yaml
# Simple case (backwards-compatible): flat list with behaviors
slots:
  - name: promoter
    behavior: product
    source: ...
  - name: orf
    behavior: product
    source: ...
  - name: barcode
    behavior: zip
    source: ...

# Hierarchical case: explicit grouping
assembly:
  strategy: zip
  allow_trim: true
  children:
    - strategy: product
      children:
        - slot: promoter
        - slot: orf
    - strategy: product
      children:
        - slot: terminator
        - slot: rbs
    - slot: barcode

# Slot definitions stay flat (just source data, no behavior)
slots:
  - name: promoter
    source: ...
  - name: orf
    source: ...
  - name: terminator
    source: ...
  - name: rbs
    source: ...
  - name: barcode
    source: ...
```

### Algorithm sketch

The assembly algorithm becomes a recursive evaluator over the tree:

```python
def evaluate(node: AssemblyNode) -> list[dict[str, NamedSequence]]:
    """Recursively evaluate an assembly tree node.

    Returns a list of partial assignment dicts. Each dict maps
    slot names to their chosen NamedSequence for that combination.
    """
    if isinstance(node, SlotLeaf):
        # Leaf: return one partial assignment per sequence
        return [{node.slot_name: seq} for seq in node.sequences]

    child_results = [evaluate(child) for child in node.children]

    if node.strategy == "product":
        # Cartesian product of all children's assignment lists
        result = []
        for combo in itertools.product(*child_results):
            merged = {}
            for partial in combo:
                merged.update(partial)
            result.append(merged)
        return result

    elif node.strategy == "zip":
        # Lockstep across children's assignment lists
        lengths = [len(r) for r in child_results]
        if node.allow_trim:
            n = min(lengths)
        else:
            if len(set(lengths)) > 1:
                raise ValueError(f"Zip children have mismatched counts: {lengths}")
            n = lengths[0]
        result = []
        for i in range(n):
            merged = {}
            for child_result in child_results:
                merged.update(child_result[i])
            result.append(merged)
        return result
```

### Flat-to-tree desugaring

For backwards compatibility, a flat config with per-slot behaviors desugars into a tree:

1. All `product` slots become children of an implicit `product` root
2. All `zip` slots become children of an implicit `zip` wrapper
3. If both exist, the root is a `zip` node containing:
   - A `product` node wrapping all product slots
   - Each zip slot as a direct child
4. `fixed` slots are handled separately (injected into every assignment)

This means existing configs work unchanged.

### Implementation considerations

1. **YAML ergonomics** — deeply nested YAML is hard to read. The flat syntax should remain the default for simple cases, with the tree syntax opt-in for complex scenarios.

2. **Validation and error messages** — errors like "group A produced 12 tuples but group B produced 8" need to trace back to specific slots and groups in the config. The tree structure should carry source location info.

3. **Name templates** — auto-generating names from a tree requires deciding which slots to include. A reasonable default: all leaf slots that aren't fixed.

4. **Memory** — the recursive evaluator materializes intermediate results. For very large trees, a streaming/generator approach would be needed, but this is an optimization concern, not a design concern.

5. **Slot ordering** — constructs need a deterministic slot order for sequence concatenation. The tree's depth-first traversal order is the natural choice, matching left-to-right reading of the YAML.

### Migration path

1. Keep the flat model as the primary interface
2. Add tree-based assembly as an alternative config format
3. Internally, convert flat configs to trees (desugaring)
4. Eventually, the tree evaluator handles all assembly, with flat configs as syntactic sugar
