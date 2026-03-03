# dna-concat

Flexible DNA sequence concatenation with slot-based combination logic.

Define multi-part DNA constructs in YAML, and dna-concat generates all combinations as FASTA or CSV — handling fixed regions, paired iteration (zip), and full combinatorial (product) assembly in a single config.

## Quick example

```yaml
slots:
  - name: promoter
    behavior: product
    source:
      inline:
        - name: pCMV
          sequence: GTTGACATTGATTATTGACTAG
        - name: pEF1a
          sequence: CCTCGGCCTCTGAGCTATTCC

  - name: orf
    behavior: product
    source:
      inline:
        - name: GFP
          sequence: ATGGTGAGCAAGGGCGAGGAG
        - name: RFP
          sequence: ATGGTGAGCAAGGGCGAGGAC

output:
  fasta: output/constructs.fasta
```

This produces 4 constructs (2 promoters x 2 ORFs) written to a FASTA file.

## Installation

### pip

```bash
pip install -e .
```

### Docker

```bash
docker build -t dna-concat .
```

## Usage

### CLI

```bash
# Run assembly
dna-concat run config.yaml

# Dry run — validate config and print summary
dna-concat run config.yaml --dry-run

# Override output paths
dna-concat run config.yaml --output-fasta out.fasta --output-csv out.csv
```

### Python API

```python
from dna_concat import AssemblySpec, Slot, SlotBehavior, NamedSequence, assemble

spec = AssemblySpec(slots=[
    Slot(name="flank", behavior=SlotBehavior.FIXED,
         sequences=[NamedSequence("flank", "GCTAGC")]),
    Slot(name="insert", behavior=SlotBehavior.PRODUCT,
         sequences=[NamedSequence("A", "ATGAAA"), NamedSequence("B", "ATGCCC")]),
])

constructs = assemble(spec)
for c in constructs:
    print(c.format_name(spec.name_template), c.full_sequence)
```

### Docker

```bash
# Run with a config file (mount your data directory)
docker run -v $(pwd)/data:/data dna-concat -c "dna-concat run /data/config.yaml"

# Run the bundled example
docker run dna-concat -c "dna-concat run examples/paired_iteration.yaml --dry-run"
```

## Slot behaviors

| Behavior  | Description | Construct count |
|-----------|-------------|-----------------|
| `fixed`   | Same single sequence in every construct | 1 (constant) |
| `zip`     | Sequences paired by index across all zip slots | N (zip length) |
| `product` | Full Cartesian product across all product slots | N1 x N2 x ... |

When both `product` and `zip` slots are present, the zip count must equal the product count — each product combination is paired 1:1 with a zip entry. Set `allow_zip_trim: true` to silently truncate zip slots that have more entries than needed (zip slots with too few entries still produce an error).

## YAML config reference

```yaml
# Optional: name template using slot names as placeholders
name_template: "{promoter}_{orf}"

# Optional: allow trimming zip slots to match product count (default: false)
allow_zip_trim: false

# Optional: output paths (relative to config file location)
output:
  fasta: output/constructs.fasta
  csv: output/constructs.csv

# Required: ordered list of slots
slots:
  - name: slot_name          # unique slot identifier
    behavior: fixed|zip|product
    source:
      # One of:
      inline:
        - name: seq_name
          sequence: ATGCCC
      fasta: path/to/seqs.fasta
      csv:
        path: path/to/seqs.csv
        name_column: name        # default: "name"
        sequence_column: sequence # default: "sequence"
        delimiter: ","            # default: ","
```

## Docker usage

The Dockerfile uses a multi-stage build. Stage 1 (`test`) runs pytest; stage 2 (`production`) creates a clean image with only runtime dependencies.

```bash
# Run tests only
docker build --target test -t dna-concat:test .

# Build production image (tests run as intermediate stage)
docker build -t dna-concat .

# Run CLI via bash entrypoint
docker run dna-concat -c "dna-concat run examples/paired_iteration.yaml --dry-run"
```

The production image includes `procps` for Nextflow compatibility and uses `ENTRYPOINT ["/bin/bash"]` so Nextflow can pass commands via CMD.

## Examples

- `examples/paired_iteration.yaml` — 5 slots with fixed flanks and zipped ORF/barcode pairs (3 constructs)
- `examples/full_combinatorial.yaml` — 2 product slots for full Cartesian product (6 constructs)
- `examples/combinatorial_barcodes.yaml` — product + zip: 6 promoter/ORF combos each paired with a unique barcode
