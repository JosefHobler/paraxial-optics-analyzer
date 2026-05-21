from __future__ import annotations

from pathlib import Path

from paraxial_optics_analyzer.cli import main
from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.report import write_report


def test_cli_info_numeric_output(capsys):
    rc = main(["info", "examples/singlet_bk7.yaml", "--rings", "2"])

    out = capsys.readouterr().out
    assert rc == 0
    assert "Plano-convex singlet" in out
    assert "EFL:" in out
    assert "best-focus RMS spot:" in out


def test_cli_report_writes_pdf(tmp_path: Path, capsys):
    out_path = tmp_path / "report.pdf"
    rc = main([
        "report", "examples/singlet_bk7.yaml",
        "-o", str(out_path),
        "--rings", "2",
    ])
    stdout = capsys.readouterr().out
    assert rc == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    assert "report:" in stdout


def test_write_png_report_direct(tmp_path: Path):
    """write_report is the library-facing entry — keep covering it directly."""
    prescription = load_prescription("examples/singlet_bk7.yaml")
    out = tmp_path / "report.png"

    result = write_report(prescription, out, n_rings=2)

    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_cli_missing_file_returns_error(tmp_path: Path, capsys):
    rc = main(["info", str(tmp_path / "no_such_file.yaml")])
    err = capsys.readouterr().err
    assert rc != 0
    assert "error" in err.lower()


def test_cli_field_angle_too_large_returns_error(capsys):
    # sin(90°) = 1 -> launch direction has dz = 0; CLI should report it cleanly.
    rc = main([
        "info", "examples/singlet_bk7.yaml",
        "--field-angle-deg", "90",
        "--rings", "2",
    ])
    err = capsys.readouterr().err
    assert rc != 0
    assert "error" in err.lower()
    assert "field angle too large" in err
