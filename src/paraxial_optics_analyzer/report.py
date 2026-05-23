from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from paraxial_optics_analyzer.analysis import find_best_focus, ray_fan, spot_diagram
from paraxial_optics_analyzer.paraxial import trace_paraxial
from paraxial_optics_analyzer.prescription import Prescription
from paraxial_optics_analyzer.raytrace import image_plane_z


def write_report(
    pre: Prescription,
    output_path: str | Path,
    *,
    field_angle_deg: float = 0.0,
    n_rings: int = 6,
) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fa = float(np.deg2rad(field_angle_deg))
    para = trace_paraxial(pre)
    best = find_best_focus(pre, fa, n_rings=n_rings)
    spot = spot_diagram(pre, fa, n_rings=n_rings, image_plane_offset=best.image_plane_offset)
    fan_t = ray_fan(pre, fa, axis="tangential", image_plane_offset=best.image_plane_offset)
    fan_s = ray_fan(pre, fa, axis="sagittal", image_plane_offset=best.image_plane_offset)

    fig = plt.figure(figsize=(11.0, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax_sum = fig.add_subplot(gs[0, 0])
    ax_spot = fig.add_subplot(gs[0, 1])
    ax_fan = fig.add_subplot(gs[1, 0])
    ax_focus = fig.add_subplot(gs[1, 1])

    _summary_panel(ax_sum, pre, para, best)
    _spot_panel(ax_spot, spot)
    _fan_panel(ax_fan, fan_t, fan_s)
    _focus_panel(ax_focus, pre, fa, n_rings, best.image_plane_offset)

    fig.suptitle(pre.name, fontsize=15)
    fig.savefig(out)
    plt.close(fig)
    return out


def _summary_panel(ax, pre, para, best) -> None:
    ax.axis("off")
    u = pre.units
    lines = [
        "Paraxial result",
        f"EFL: {para.efl:.6g} {u}",
        f"BFL: {para.bfl:.6g} {u}",
        f"Image plane z: {image_plane_z(pre):.6g} {u}",
        f"Image distance: {para.image_distance:.6g} {u}",
        f"f-number: f/{para.f_number:.4g}",
        "",
        "Best focus",
        f"Defocus: {best.image_plane_offset:.6g} {u}",
        f"RMS at nominal: {best.rms_at_nominal:.6g} {u}",
        f"RMS at best: {best.rms_at_best:.6g} {u}",
    ]
    ax.text(0.0, 1.0, "\n".join(lines), va="top", family="monospace", fontsize=10)


def _spot_panel(ax, spot) -> None:
    pts = spot.points
    if len(pts):
        ax.scatter(pts[:, 0], pts[:, 1], s=12, color="#1f77b4", alpha=0.85)
        ax.scatter([spot.centroid[0]], [spot.centroid[1]], s=45, marker="+", color="#d62728")
    ax.set_title(f"Spot diagram, RMS={spot.rms:.4g}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="datalim")


def _fan_panel(ax, tan, sag) -> None:
    ax.plot(tan.pupil_coords, tan.transverse, label="Tangential")
    ax.plot(sag.pupil_coords, sag.transverse, label="Sagittal")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_title("Ray fan")
    ax.set_xlabel("pupil coordinate")
    ax.set_ylabel("transverse aberration")
    ax.grid(True, alpha=0.3)
    ax.legend()


def _focus_panel(ax, pre, fa, n_rings, best_offset) -> None:
    span = max(1.0, abs(best_offset) * 2.0)
    offsets = np.linspace(best_offset - span, best_offset + span, 41)
    rms = [
        spot_diagram(pre, fa, n_rings=n_rings, image_plane_offset=float(o)).rms
        for o in offsets
    ]
    ax.plot(offsets, rms, color="#2ca02c")
    ax.axvline(best_offset, color="#d62728", linewidth=1.0, label="best")
    ax.set_title("Focus search")
    ax.set_xlabel("image-plane offset")
    ax.set_ylabel("RMS spot radius")
    ax.grid(True, alpha=0.3)
    ax.legend()
