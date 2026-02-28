"""Tests for dna_concat.config."""

from pathlib import Path

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
