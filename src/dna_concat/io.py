"""File I/O for dna_concat: FASTA, CSV, and inline sequence reading/writing."""

from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dna_concat.models import Construct, NamedSequence


def read_fasta(path: str | Path) -> list[NamedSequence]:
    """Read sequences from a FASTA file."""
    path = Path(path)
    results = []
    for record in SeqIO.parse(str(path), "fasta"):
        results.append(NamedSequence(name=record.id, sequence=str(record.seq)))
    return results


def read_csv(
    path: str | Path,
    name_column: str = "name",
    sequence_column: str = "sequence",
    delimiter: str = ",",
) -> list[NamedSequence]:
    """Read sequences from a CSV file."""
    path = Path(path)
    results = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            results.append(
                NamedSequence(
                    name=row[name_column], sequence=row[sequence_column]
                )
            )
    return results


def parse_inline(entries: list[dict]) -> list[NamedSequence]:
    """Parse inline sequence definitions from config dicts."""
    return [NamedSequence(name=e["name"], sequence=e["sequence"]) for e in entries]


def write_fasta(
    constructs: list[Construct], path: str | Path, name_template: str
) -> None:
    """Write constructs to a FASTA file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for construct in constructs:
        name = construct.format_name(name_template)
        record = SeqRecord(Seq(construct.full_sequence), id=name, description="")
        records.append(record)
    SeqIO.write(records, str(path), "fasta")


def write_csv(
    constructs: list[Construct], path: str | Path, name_template: str
) -> None:
    """Write constructs to a CSV file with per-slot breakdown columns."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not constructs:
        return

    slot_names = constructs[0].slot_order
    fieldnames = ["name", "full_sequence"]
    for sn in slot_names:
        fieldnames.extend([f"{sn}_name", f"{sn}_seq"])

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for construct in constructs:
            row: dict[str, str] = {
                "name": construct.format_name(name_template),
                "full_sequence": construct.full_sequence,
            }
            for sn in slot_names:
                ns = construct.assignments[sn]
                row[f"{sn}_name"] = ns.name
                row[f"{sn}_seq"] = ns.sequence
            writer.writerow(row)
