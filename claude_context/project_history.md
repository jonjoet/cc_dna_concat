# dna-concat: Project History and Context

## Motivation

Existing bioinformatics tools for DNA construct assembly either assume simple concatenation or require full-blown cloning simulation. dna-concat fills the gap: a lightweight tool for generating all combinations of multi-part DNA constructs from a YAML config, with explicit control over how slots combine.

The core use case is designing synthetic biology libraries — promoter/ORF/barcode combinations where some parts are fixed, some iterate in lockstep (zip), and others form Cartesian products.

## Core abstraction: Slots and behaviors

The central design idea is the **slot**: a named position in a construct that holds one or more candidate sequences. Each slot has a **behavior** that determines how its sequences combine with other slots:

- **fixed** — exactly one sequence, present in every construct (e.g., flanking regions, spacers)
- **zip** — sequences paired by index; all zip slots must have the same count (unless `allow_zip_trim` is set)
- **product** — full Cartesian product across all product slots

When both product and zip behaviors are present, the zip count must equal the product count, creating a 1:1 pairing between combinatorial results and zip entries. This enables patterns like "6 promoter/ORF combos, each assigned a unique barcode." With `allow_zip_trim: true`, zip slots with more entries than the product count are silently truncated; zip slots with too few entries still produce an error.

## Architecture and modules

```
src/dna_concat/
  __init__.py    — public API exports (AssemblySpec, Slot, Construct, assemble, etc.)
  models.py      — dataclasses: SlotBehavior, NamedSequence, Slot, Construct, AssemblySpec
  assembly.py    — core algorithm: fixed/product/zip combination logic
  config.py      — YAML parsing, source resolution (inline, FASTA, CSV)
  io.py          — file I/O: read/write FASTA and CSV, parse inline sequences
  cli.py         — Click CLI: `dna-concat run` with --dry-run, --output-fasta, --output-csv
```

Key design decisions:
- **Dataclasses over Pydantic** — minimal dependencies, sufficient for config validation
- **BioPython for FASTA** — standard library for sequence I/O
- **Click for CLI** — clean subcommand structure (`dna-concat run`)
- **Path resolution relative to config** — all file paths in YAML resolve from the config file's directory
- **Auto-generated name templates** — if no `name_template` is specified, one is built from non-fixed slot names

## Current state

- All 6 modules implemented
- 26 tests passing across 5 test files (test_models, test_assembly, test_config, test_io, test_cli)
- Multi-stage Dockerfile: test stage runs pytest, production stage is Nextflow-compatible
- 3 example YAML configs demonstrating paired iteration, full combinatorial, and mixed product+zip patterns
- CLI supports dry-run validation, output path overrides, and both FASTA and CSV output

## File layout

```
pyproject.toml
Dockerfile
README.md
src/dna_concat/          — package source (6 modules)
tests/                   — pytest suite (26 tests)
  fixtures/              — test FASTA/CSV files
examples/                — 3 YAML example configs
claude_context/          — this file
```
