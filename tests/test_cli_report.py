from __future__ import annotations

from pathlib import Path

from paraxial_optics_analyzer.cli import main
from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.report import write_report


def test_cli_numeric_output(capsys):
    rc = main(["examples/singlet_bk7.yaml", "--no-report", "--rings", "2"])

    out = capsys.readouterr().out
    assert rc == 0
    assert "Plano-convex singlet" in out
    assert "EFL:" in out
    assert "best-focus RMS spot:" in out


def test_write_png_report(tmp_path: Path):
    prescription = load_prescription("examples/singlet_bk7.yaml")
    out = tmp_path / "report.png"

    result = write_report(prescription, out, n_rings=2)

    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_cli_missing_file_returns_error(tmp_path: Path, capsys):
    rc = main([str(tmp_path / "no_such_file.yaml"), "--no-report"])
    err = capsys.readouterr().err
    assert rc != 0
    assert "error" in err.lower()


def test_cli_field_angle_too_large_returns_error(capsys):
    # sin(90°) = 1 -> launch direction has dz = 0; CLI should report it cleanly.
    rc = main([
        "examples/singlet_bk7.yaml",
        "--no-report",
        "--field-angle-deg", "90",
        "--rings", "2",
    ])
    err = capsys.readouterr().err
    assert rc != 0
    assert "error" in err.lower()
    assert "field angle too large" in err
