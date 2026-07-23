# dna-concat

Flexible DNA sequence concatenation with slot-based combination logic.

Define multi-part DNA constructs in YAML, and dna-concat generates all combinations as FASTA or CSV — handling fixed regions, paired iteration (zip), and full combinatorial (product) assembly in a single config. Includes a CLI, a Python API, and a Streamlit web interface (work in progress).

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

### pip with Streamlit GUI

```bash
pip install -e ".[streamlit]"
```

### Docker (CLI)

```bash
docker build -t dna-concat .
```

### Docker (Streamlit GUI)

```bash
docker build -f Dockerfile.streamlit -t dna-concat-streamlit:latest .
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

### Docker (CLI)

```bash
# Run with a config file (mount your data directory)
docker run -v $(pwd)/data:/data dna-concat -c "dna-concat run /data/config.yaml"

# Run the bundled example
docker run dna-concat -c "dna-concat run examples/paired_iteration.yaml --dry-run"
```

### Streamlit GUI

Run locally:

```bash
pip install -e ".[streamlit]"
streamlit run src/dna_concat/gui.py
# Access at http://localhost:8501
```

Or via Docker:

```bash
docker build -f Dockerfile.streamlit -t dna-concat-streamlit:latest .
docker run -p 8501:8501 dna-concat-streamlit:latest
# Access at http://localhost:8501
```

The GUI provides an interactive interface for building slot configurations, running assemblies, and downloading results as FASTA or CSV. You can configure slots manually or upload a YAML config file via the sidebar. Note that YAML configs referencing external fasta/csv files will need those sequences re-added via the file upload UI, since file paths from your local machine are not accessible inside the browser.

## Slot behaviors

| Behavior  | Description | Construct count |
|-----------|-------------|-----------------|
| `fixed`   | Same single sequence in every construct | 1 (constant) |
| `zip`     | Sequences paired by index across all zip slots | N (zip length) |
| `product` | Full Cartesian product across all product slots | N1 x N2 x ... |
| `variable_stuffer` | Single sequence truncated so it plus a named counterpart slot always equal the untruncated stuffer length | 1 (constant) |

When both `product` and `zip` slots are present, the zip count must equal the product count — each product combination is paired 1:1 with a zip entry. Set `allow_zip_trim: true` to silently truncate zip slots that have more entries than needed (zip slots with too few entries still produce an error).

### Variable stuffer

A `variable_stuffer` slot holds a single sequence and truncates it per construct so that the stuffer plus its `counterpart` slot always sum to the untruncated stuffer's length — only the stuffer is trimmed. This keeps a construct region at a fixed total length as the counterpart (e.g. an ORF) varies.

- `counterpart` names the slot whose length is subtracted. It may be any non-stuffer slot, in any position — other slots are allowed to sit between the stuffer and its counterpart (they add length on top of the fixed stuffer+counterpart total).
- `truncate_side` (`left` or `right`) selects which end of the stuffer is trimmed, independently of where the counterpart sits.

If the counterpart is ever longer than the untruncated stuffer, assembly raises an error (the combined length can't be held constant).

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
    behavior: fixed|zip|product|variable_stuffer
    # For variable_stuffer slots only:
    counterpart: other_slot   # name of the slot whose length is subtracted
    truncate_side: left|right # which end of the stuffer to trim
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
- `examples/variable_stuffer.yaml` — a variable stuffer that shrinks as a zipped ORF grows, holding ORF + stuffer at a constant length
