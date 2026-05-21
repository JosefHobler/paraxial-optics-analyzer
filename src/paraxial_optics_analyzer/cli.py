"""CLI entry point — subcommand-based interface (argparse)."""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from paraxial_optics_analyzer import __version__
from paraxial_optics_analyzer.analysis import find_best_focus, spot_diagram
from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.paraxial import trace_paraxial
from paraxial_optics_analyzer.prescription import Prescription, PrescriptionError
from paraxial_optics_analyzer.raytrace import TraceError
from paraxial_optics_analyzer.report import write_report
from paraxial_optics_analyzer.validate import run_all as run_validation

_DESCRIPTION = (
    "Paraxial optics analyzer — sequential ray tracing for centered lens systems."
)

_EPILOG = """\
examples:
  analyze info examples/singlet_bk7.yaml
  analyze info examples/cooke_triplet.yaml --field-angle-deg 5
  analyze report examples/cooke_triplet.yaml -o cooke.pdf
  analyze validate

run `analyze <command> --help` for command-specific options.
"""

_TRACE_HINT = "lower --field-angle-deg or widen the aperture (semi_diameter)"


# --- parser ---


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="analyze",
        description=_DESCRIPTION,
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )

    sub = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        required=True,
        title="commands",
    )

    info = sub.add_parser(
        "info",
        help="print first-order properties and spot statistics",
        description="Print paraxial EFL/BFL/f-number, spot RMS, and best-focus offset.",
    )
    _add_prescription_arg(info)
    _add_analysis_args(info)
    info.set_defaults(func=_cmd_info)

    rep = sub.add_parser(
        "report",
        help="write a 4-panel PDF/PNG report",
        description="Run the full analysis and write a 4-panel summary report.",
    )
    _add_prescription_arg(rep)
    rep.add_argument(
        "-o", "--output", metavar="PATH", default=None,
        help="output path (.pdf or .png). default: <stem>_report.pdf",
    )
    _add_analysis_args(rep)
    rep.set_defaults(func=_cmd_report)

    val = sub.add_parser(
        "validate",
        help="run built-in physics self-checks",
        description=(
            "Run built-in physics self-checks: lensmaker equation, paraxial limit "
            "of the real trace, and a Cooke-triplet EFL cross-check."
        ),
    )
    val.set_defaults(func=_cmd_validate)

    return parser


def _add_prescription_arg(p: argparse.ArgumentParser) -> None:
    p.add_argument("prescription", help="path to a lens prescription YAML file")


def _add_analysis_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--field-angle-deg", metavar="DEG", type=float, default=0.0,
        help="field angle in degrees (default: %(default)s)",
    )
    p.add_argument(
        "--rings", metavar="N", type=int, default=6,
        help="hexapolar pupil rings (default: %(default)s)",
    )


# --- helpers ---


def _fail(msg: str, *, hint: str | None = None, code: int = 1) -> int:
    print(f"error: {msg}", file=sys.stderr)
    if hint:
        print(f"  hint: {hint}", file=sys.stderr)
    return code


def _load(path_str: str) -> Prescription | int:
    try:
        return load_prescription(Path(path_str))
    except PrescriptionError as e:
        return _fail(str(e), code=2)


# --- commands ---


def _cmd_info(args: argparse.Namespace) -> int:
    pre = _load(args.prescription)
    if isinstance(pre, int):
        return pre

    fa = math.radians(args.field_angle_deg)
    try:
        para = trace_paraxial(pre)
        spot = spot_diagram(pre, field_angle_rad=fa, n_rings=args.rings)
        if spot.n_failed and len(spot.points) == 0:
            raise TraceError("every sampled ray missed an aperture or hit TIR")
        best = find_best_focus(pre, field_angle_rad=fa, n_rings=args.rings)
    except TraceError as e:
        return _fail(f"ray trace failed — {e}", hint=_TRACE_HINT, code=3)
    except ValueError as e:
        return _fail(str(e), hint=_TRACE_HINT, code=2)

    u = pre.units
    print(f"{pre.name}")
    print(f"  EFL: {para.efl:.12g} {u}")
    print(f"  BFL: {para.bfl:.12g} {u}")
    print(f"  image distance: {para.image_distance:.12g} {u}")
    print(f"  f-number: f/{para.f_number:.6g}")
    print(f"  nominal RMS spot: {spot.rms:.12g} {u}")
    print(f"  best-focus offset: {best.image_plane_offset:.12g} {u}")
    print(f"  best-focus RMS spot: {best.rms_at_best:.12g} {u}")
    if spot.n_failed:
        print(f"  failed rays: {spot.n_failed}", file=sys.stderr)
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    pre = _load(args.prescription)
    if isinstance(pre, int):
        return pre

    path = Path(args.prescription)
    out = Path(args.output) if args.output else path.with_name(f"{path.stem}_report.pdf")
    try:
        report_path = write_report(
            pre, out,
            field_angle_deg=args.field_angle_deg,
            n_rings=args.rings,
        )
    except TraceError as e:
        return _fail(f"could not render report — {e}", hint=_TRACE_HINT, code=3)
    except ValueError as e:
        return _fail(str(e), hint=_TRACE_HINT, code=2)
    print(f"report: {report_path}")
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    failed = 0
    for r in run_validation():
        print(r.format_line())
        if not r.passed:
            failed += 1
    if failed:
        print(f"\n{failed} check(s) failed", file=sys.stderr)
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
