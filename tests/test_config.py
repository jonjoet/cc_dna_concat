"""Tests for dna_concat.config."""

from pathlib import Path

import yaml

from dna_concat.config import load_config
from dna_concat.models import SlotBehavior

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_zip_config():
    spec, output_config = load_config(FIXTURES / "simple_zip.yaml")
    assert len(spec.slots) == 4
    assert spec.slots[0].behavior is SlotBehavior.FIXED
    assert spec.slots[1].behavior is SlotBehavior.ZIP
    assert len(spec.slots[1].sequences) == 3  # 3 ORFs from FASTA
    assert len(spec.slots[2].sequences) == 3  # 3 barcodes from CSV
    assert spec.name_template == "{orf}_{barcode}"
    assert "fasta" in output_config


def test_load_product_config():
    spec, output_config = load_config(FIXTURES / "product.yaml")
    assert len(spec.slots) == 2
    assert spec.slots[0].behavior is SlotBehavior.PRODUCT
    assert len(spec.slots[0].sequences) == 2
    assert len(spec.slots[1].sequences) == 3
    # Auto-generated name template
    assert spec.name_template == "{promoter}_{orf}"


def test_allow_zip_trim_parsing(tmp_path):
    """allow_zip_trim in YAML config is parsed and passed to AssemblySpec."""
    config = {
        "allow_zip_trim": True,
        "slots": [
            {
                "name": "a",
                "behavior": "zip",
                "source": {"inline": [{"name": "x", "sequence": "AAA"}]},
            },
        ],
    }
    config_path = tmp_path / "trim.yaml"
    config_path.write_text(yaml.dump(config))
    spec, _ = load_config(config_path)
    assert spec.allow_zip_trim is True


def test_allow_zip_trim_defaults_false(tmp_path):
    """allow_zip_trim defaults to False when not specified."""
    config = {
        "slots": [
            {
                "name": "a",
                "behavior": "zip",
                "source": {"inline": [{"name": "x", "sequence": "AAA"}]},
            },
        ],
    }
    config_path = tmp_path / "no_trim.yaml"
    config_path.write_text(yaml.dump(config))
    spec, _ = load_config(config_path)
    assert spec.allow_zip_trim is False


def test_load_variable_stuffer_config(tmp_path):
    """A variable_stuffer slot parses into a Slot with its counterpart/side."""
    config = {
        "slots": [
            {
                "name": "orf",
                "behavior": "zip",
                "source": {"inline": [{"name": "x", "sequence": "ATG"}]},
            },
            {
                "name": "stuffer",
                "behavior": "variable_stuffer",
                "counterpart": "orf",
                "truncate_side": "left",
                "source": {"inline": [{"name": "pad", "sequence": "AAAAAA"}]},
            },
        ],
    }
    config_path = tmp_path / "stuffer.yaml"
    config_path.write_text(yaml.dump(config))
    spec, _ = load_config(config_path)
    assert spec.slots[1].behavior is SlotBehavior.VARIABLE_STUFFER
    assert spec.slots[1].counterpart == "orf"
    assert spec.slots[1].truncate_side == "left"
