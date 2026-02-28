"""Tests for dna_concat.cli."""

from pathlib import Path

from click.testing import CliRunner

from dna_concat.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


def test_dry_run_zip():
    runner = CliRunner()
    result = runner.invoke(main, ["run", str(FIXTURES / "simple_zip.yaml"), "--dry-run"])
    assert result.exit_code == 0
    assert "Constructs: 3" in result.output
    assert "Slots: 4" in result.output


def test_dry_run_product():
    runner = CliRunner()
    result = runner.invoke(main, ["run", str(FIXTURES / "product.yaml"), "--dry-run"])
    assert result.exit_code == 0
    assert "Constructs: 6" in result.output


def test_output_fasta(tmp_path):
    out = tmp_path / "out.fasta"
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["run", str(FIXTURES / "simple_zip.yaml"), "--output-fasta", str(out)],
    )
    assert result.exit_code == 0
    assert out.exists()
    content = out.read_text()
    assert ">ORF1_BC01" in content


def test_output_csv(tmp_path):
    out = tmp_path / "out.csv"
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["run", str(FIXTURES / "product.yaml"), "--output-csv", str(out)],
    )
    assert result.exit_code == 0
    assert out.exists()
    assert "full_sequence" in out.read_text()
