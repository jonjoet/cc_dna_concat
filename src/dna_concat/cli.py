"""CLI for dna_concat."""

from __future__ import annotations

import click

from dna_concat.assembly import assemble
from dna_concat.config import load_config
from dna_concat.io import write_csv, write_fasta


@click.group()
def main() -> None:
    """dna-concat: Flexible DNA sequence concatenation."""


@main.command()
@click.argument("config", type=click.Path(exists=True))
@click.option("--output-fasta", type=click.Path(), default=None, help="Override FASTA output path")
@click.option("--output-csv", type=click.Path(), default=None, help="Override CSV output path")
@click.option("--dry-run", is_flag=True, help="Validate config and print summary without writing")
def run(config: str, output_fasta: str | None, output_csv: str | None, dry_run: bool) -> None:
    """Generate concatenated constructs from a YAML config."""
    spec, output_config = load_config(config)
    constructs = assemble(spec)

    if dry_run:
        click.echo(f"Slots: {len(spec.slots)}")
        for slot in spec.slots:
            click.echo(f"  {slot.name}: {slot.behavior.value} ({len(slot.sequences)} sequences)")
        click.echo(f"Constructs: {len(constructs)}")
        click.echo(f"Name template: {spec.name_template}")
        if constructs:
            click.echo(f"First: {constructs[0].format_name(spec.name_template)}")
            click.echo(f"  Length: {len(constructs[0].full_sequence)} bp")
        return

    fasta_path = output_fasta or output_config.get("fasta")
    csv_path = output_csv or output_config.get("csv")

    if fasta_path:
        write_fasta(constructs, fasta_path, spec.name_template)
        click.echo(f"Wrote {len(constructs)} constructs to {fasta_path}")

    if csv_path:
        write_csv(constructs, csv_path, spec.name_template)
        click.echo(f"Wrote {len(constructs)} constructs to {csv_path}")

    if not fasta_path and not csv_path:
        click.echo("No output paths specified. Use --output-fasta or --output-csv, or set in config.")
