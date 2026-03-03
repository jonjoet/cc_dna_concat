"""YAML configuration parsing for dna_concat."""

from __future__ import annotations

from pathlib import Path

import yaml

from dna_concat.io import parse_inline, read_csv, read_fasta
from dna_concat.models import AssemblySpec, Slot, SlotBehavior


def load_config(config_path: str | Path) -> tuple[AssemblySpec, dict]:
    """Load a YAML config file and return (AssemblySpec, output_config).

    File paths in the config are resolved relative to the config file's directory.
    """
    config_path = Path(config_path)
    config_dir = config_path.parent

    with open(config_path) as f:
        raw = yaml.safe_load(f)

    slots = []
    for slot_def in raw["slots"]:
        behavior = SlotBehavior(slot_def["behavior"])
        sequences = _resolve_source(slot_def["source"], config_dir)
        slots.append(
            Slot(
                name=slot_def["name"],
                behavior=behavior,
                sequences=sequences,
            )
        )

    spec = AssemblySpec(
        slots=slots,
        name_template=raw.get("name_template"),
        allow_zip_trim=raw.get("allow_zip_trim", False),
    )

    output_config = raw.get("output", {})
    # Resolve output paths relative to config dir
    for key in ("fasta", "csv"):
        if key in output_config:
            output_config[key] = str(config_dir / output_config[key])

    return spec, output_config


def _resolve_source(source: dict, config_dir: Path) -> list:
    """Resolve a slot source block to a list of NamedSequence."""
    if "inline" in source:
        return parse_inline(source["inline"])
    elif "fasta" in source:
        return read_fasta(config_dir / source["fasta"])
    elif "csv" in source:
        csv_conf = source["csv"]
        return read_csv(
            path=config_dir / csv_conf["path"],
            name_column=csv_conf.get("name_column", "name"),
            sequence_column=csv_conf.get("sequence_column", "sequence"),
            delimiter=csv_conf.get("delimiter", ","),
        )
    else:
        raise ValueError(f"Unknown source type: {list(source.keys())}")
