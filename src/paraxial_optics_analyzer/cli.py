"""CLI entry point."""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from paraxial_optics_analyzer.analysis import find_best_focus, spot_diagram
from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.paraxial import trace_paraxial
from paraxial_optics_analyzer.prescription import PrescriptionError
from paraxial_optics_analyzer.raytrace import TraceError
from paraxial_optics_analyzer.report import write_report


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="analyze",
        description="Analyze a sequential spherical lens prescription",
    )
    p.add_argument("prescription", help="Path to a lens prescription YAML file.")
    p.add_argument("-o", "--output", default=None,
                   help="Report path (.pdf or .png). Default: <stem>_report.pdf")
    p.add_argument("--field-angle-deg", type=float, default=0.0,
                   help="Field angle in degrees for spot/ray-fan analysis.")
    p.add_argument("--rings", type=int, default=6,
                   help="Hexapolar pupil rings for spot and focus analysis.")
    p.add_argument("--no-report", action="store_true",
                   help="Print numeric results only; do not write a matplotlib report.")
    return p


def _fail(msg: str, *, hint: str | None = None, code: int = 1) -> int:
    print(f"error: {msg}", file=sys.stderr)
    if hint:
        print(f"  hint: {hint}", file=sys.stderr)
    return code


_TRACE_HINT = "try lowering --field-angle-deg or widening the aperture (semi_diameter)"


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    path = Path(args.prescription)
    fa = math.radians(args.field_angle_deg)

    try:
        pre = load_prescription(path)
    except PrescriptionError as e:
        return _fail(str(e), code=2)

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

    if not args.no_report:
        out = Path(args.output) if args.output else path.with_name(f"{path.stem}_report.pdf")
        try:
            report_path = write_report(pre, out, field_angle_deg=args.field_angle_deg, n_rings=args.rings)
        except TraceError as e:
            return _fail(f"could not render report — {e}", hint=_TRACE_HINT, code=3)
        print(f"  report: {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
