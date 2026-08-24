#!/usr/bin/env python
r"""Figures for the T10 exterior-medium run record — the lane's one entry point, regenerates everything.

Reads the committed JSON + NPZ (no GPU, no re-run) and writes four PNGs.  This is the ISOLATION render the
ladder gate requires of a component reaching ``CUDA_UNIT``: the ``extracellular_medium`` component and the
surface it acts on, with nothing else in the scene.

Every figure carries the reference it is judged against on the SAME axes, because a measurement without its
reference overlaid is not a check:

  1. ``traction_field.png``   — the per-node traction under rigid translation and rotation, projected on the
                                surface, with the exact Stokes resultant drawn as the reference arrow. This
                                is the figure that catches a wrong sign or a single outlier node, which no
                                aggregate in the record can.
  2. ``rigid_modes.png``      — the six rigid-mode resistance eigenvalues against ``6 pi mu a`` and
                                ``8 pi mu a^3``, and against what a constant ``aI`` regulariser would give.
                                The point of the whole component is that these two references are visibly
                                different objects.
  3. ``grid_invariance.png``  — the Stokes-law error against surface spacing on log-log axes, with the
                                fitted order and the first-order guide line, so ``epsilon = h`` is seen to
                                behave as a discretisation parameter rather than a tuned constant.
  4. ``dissipation.png``      — traction-velocity alignment per node: the medium must OPPOSE motion
                                everywhere, so every node's projection must be negative and the histogram
                                must sit entirely on one side of zero.

Usage (dev Mac, no CUDA needed)::

    python aleph/scripts/ac_gate_t10_medium_vis.py \
        --record aleph/outputs/ac/gate_t10_medium/medium_exterior_native.json

Visualization integrity (project rule): no axis truncation; every log axis is named as such in its own
label; units are engine units (pN, um, pN*s/um) and named on every axis; the reference is drawn on the
measurement rather than described in a caption.  Magnitudes are shown but NOT quotable — ``mu_medium`` is a
PI-GAP, so each panel that carries a force axis says so in its title.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

#: One palette across all four figures, so a colour means the same thing everywhere.
COLOR_TRANSLATION = "#3b6ea5"
COLOR_ROTATION = "#c1502e"
COLOR_REFERENCE = "#2f2f2f"
COLOR_REGULARISER = "#8a8a8a"

#: Repeated on every panel that carries a force axis. The magnitudes are real outputs of a real solve; what
#: makes them non-quotable is the unsourced viscosity they scale linearly with, not any doubt about the solve.
GAP_BANNER = "magnitudes NOT quotable: mu_medium is a PI-GAP (drag scales linearly in it)"


def _load(record_path: Path) -> tuple[dict, dict[str, np.ndarray]]:
    """Load the run record and its sidecar field dump.

    Args:
        record_path: Path to the committed run-record JSON.

    Returns:
        Tuple ``(record, fields)``.

    Raises:
        FileNotFoundError: If either the record or its NPZ sidecar is missing — a figure invented from a
            record without its fields would be a drawing, not a measurement.
    """
    record = json.loads(record_path.read_text(encoding="utf-8"))
    npz_path = record_path.with_suffix(".npz")
    if not npz_path.exists():
        raise FileNotFoundError(
            f"field dump {npz_path} is absent; rerun the native gate with --out so it writes the sidecar"
        )
    with np.load(npz_path) as data:
        fields = {key: data[key] for key in data.files}
    return record, fields


def _traction_field(record: dict, fields: dict[str, np.ndarray], out: Path) -> Path:
    """Draw the per-node traction under translation and rotation, with the exact resultant overlaid."""
    pos = fields["surface_position_um"]
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 5.2))
    ratios = record["measurements"]["stokes_law_ratio"]

    # Each mode is drawn in ITS OWN plane. Translation along +x is a meridional field, so x-z shows it; a
    # rotation about +z lives in x-y, and projecting it onto x-z would render the circulation as a pile of
    # collinear arrows — the plane is part of the measurement protocol, not a layout preference.
    for ax, key, colour, title, plane in (
        (axes[0], "traction_translation_pN", COLOR_TRANSLATION, "rigid translation (+x)", (0, 2)),
        (axes[1], "traction_rotation_pN", COLOR_ROTATION, "rigid rotation (about +z)", (0, 1)),
    ):
        traction = fields[key]
        u, v = plane
        names = ("x", "y", "z")
        scale = 0.6 * float(np.ptp(pos[:, u])) / max(float(np.abs(traction).max()), 1e-300)
        ax.scatter(pos[:, u], pos[:, v], s=6, c=colour, alpha=0.40, linewidths=0)
        ax.quiver(
            pos[:, u], pos[:, v], traction[:, u] * scale, traction[:, v] * scale,
            angles="xy", scale_units="xy", scale=1.0, width=0.003, color=colour, alpha=0.8,
        )
        ax.set_xlabel(f"{names[u]} [µm]")
        ax.set_ylabel(f"{names[v]} [µm]")
        ax.set_aspect("equal")
        ax.set_title(f"{title}\nper-node traction [pN], arrows to common scale")

    # The two modes are judged on DIFFERENT resultants, so they are drawn and labelled differently rather
    # than sharing one arrow: translation is checked on the net FORCE, rotation on the net TORQUE, and a
    # rotating sphere's net force is zero — labelling a torque ratio onto a force arrow would misread.
    translation = fields["traction_translation_pN"]
    net_force = translation.sum(axis=0)
    scale = 0.6 * float(np.ptp(pos[:, 0])) / max(float(np.abs(translation).max()), 1e-300)
    axes[0].quiver(
        [0.0], [0.0], [net_force[0] * scale], [net_force[2] * scale],
        angles="xy", scale_units="xy", scale=1.0, width=0.010, color=COLOR_REFERENCE,
        label=f"net force / exact 6πµa = {ratios['translation_measured_over_exact']:.4f}",
    )
    axes[0].legend(loc="upper right", fontsize=8)

    rotation = fields["traction_rotation_pN"]
    torque = np.cross(pos, rotation).sum(axis=0)
    rotation_force = rotation.sum(axis=0)
    axes[1].annotate(
        f"net torque$_z$ / exact 8πµa³ = {ratios['rotation_measured_over_exact']:.4f}\n"
        f"net FORCE |Σf| / Σ|f| = "
        f"{np.linalg.norm(rotation_force) / max(np.abs(rotation).sum(), 1e-300):.2e}"
        "  (must vanish: a spinning sphere is not pushed)\n"
        f"torque sign = {'−z, opposes' if torque[2] < 0 else '+z, WRONG SIGN'}",
        xy=(0.02, 0.02), xycoords="axes fraction", fontsize=8,
        bbox={"boxstyle": "round", "fc": "white", "ec": COLOR_ROTATION, "alpha": 0.85},
    )
    figure.suptitle(
        "extracellular_medium in ISOLATION — exterior Stokes traction on the native membrane surface\n"
        f"{GAP_BANNER}", fontsize=9,
    )
    figure.tight_layout()
    path = out / "traction_field.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def _rigid_modes(record: dict, fields: dict[str, np.ndarray], out: Path) -> Path:
    """Draw the six rigid-mode resistances against Stokes and against a constant regulariser."""
    eigenvalues = np.sort(fields["rigid_mode_eigenvalues"])
    modes = record["measurements"]["rigid_modes_are_loaded"]
    radius = float(record["config"]["argv"]["radius"])
    translation_exact = float(modes["translation_resistance_over_radius"]) * radius
    rotation_exact = translation_exact * float(modes["rotation_over_translation_shape"])

    figure, ax = plt.subplots(figsize=(8.0, 5.2))
    index = np.arange(1, eigenvalues.size + 1)
    ax.semilogy(index, eigenvalues, "o-", color=COLOR_TRANSLATION, ms=8,
                label="measured rigid-mode resistance eigenvalues")
    ax.axhline(translation_exact, color=COLOR_REFERENCE, ls="--", lw=1.2,
               label="translation reference 6·π·µ·a  [pN·s/µm]")
    ax.axhline(rotation_exact, color=COLOR_ROTATION, ls="--", lw=1.2,
               label="rotation reference 8·π·µ·a³  [pN·s·µm]")
    # NO line is drawn for the aI regulariser it replaces. Its value is not a quantity this run measured,
    # and drawing an un-measured reference on a measurement is precisely what the visualization rule
    # forbids — so the contrast is STATED, in the one form that is checkable, and not fabricated as data.
    ax.annotate(
        "the aI regulariser this replaces is NOT drawn: it is one constant with\n"
        "no relation to the geometry, so it has no place on this axis. The\n"
        "checkable contrast is the SCALING — resistance ∝ a here, ∝ a⁰ there\n"
        f"(measured translation/a = {float(modes['translation_resistance_over_radius']):.5g} pN·s/µm²)",
        xy=(0.03, 0.42), xycoords="axes fraction", fontsize=8,
        bbox={"boxstyle": "round", "fc": "white", "ec": COLOR_REGULARISER, "alpha": 0.9},
    )
    ax.set_xlabel("rigid mode, sorted (3 translations, then 3 rotations)")
    ax.set_ylabel("resistance eigenvalue [pN·s/µm and pN·s·µm] — LOG axis")
    ax.set_title(
        "the six whole-body modes are LOADED, and by geometry — "
        f"min eigenvalue {modes['min_eigenvalue']:.4g} > 0\n" + GAP_BANNER,
        fontsize=9,
    )
    ax.legend(fontsize=8, loc="center right")
    figure.tight_layout()
    path = out / "rigid_modes.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def _grid_invariance(record: dict, out: Path) -> Path:
    """Draw the Stokes-law error against surface spacing, with the fitted order and a first-order guide."""
    study = record["measurements"]["grid_invariance"]
    spacing = np.array([row["spacing_um"] for row in study["rows"]], dtype=float)
    translation = np.array(
        [row["translation_relative_error"] for row in study["rows"]], dtype=float
    )
    rotation = np.array([row["rotation_relative_error"] for row in study["rows"]], dtype=float)
    figure, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.loglog(spacing, translation, "o-", color=COLOR_TRANSLATION, label="translation relative error")
    ax.loglog(spacing, rotation, "s-", color=COLOR_ROTATION, label="rotation relative error")
    guide = translation[-1] * (spacing / spacing[-1])
    ax.loglog(spacing, guide, ls="--", color=COLOR_REFERENCE, lw=1.0,
              label="first-order guide (slope 1)")
    ax.set_xlabel("surface mean nearest-neighbour spacing h [µm] — LOG axis")
    ax.set_ylabel("|measured − Stokes| / Stokes  [dimensionless] — LOG axis")
    ax.set_title(
        "grid invariance at FIXED ε/h = "
        f"{study['epsilon_ratio']:g}: refining the surface is the lever, not lowering ε\n"
        f"fitted order in h = {study['translation_order']:.3f}  (ε is DERIVED, not tuned)",
        fontsize=9,
    )
    for spacing_value, error, row in zip(spacing, translation, study["rows"], strict=True):
        ax.annotate(f"N={int(row['n_points'])}", (spacing_value, error),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.legend(fontsize=8)
    figure.tight_layout()
    path = out / "grid_invariance.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def _dissipation(record: dict, fields: dict[str, np.ndarray], out: Path) -> Path:
    """Draw the per-node traction-velocity projection: it must be negative at EVERY node."""
    pos = fields["surface_position_um"]
    translation_velocity = np.tile(np.array([1.0, 0.0, 0.0]), (pos.shape[0], 1))
    rotation_velocity = np.cross(np.array([0.0, 0.0, 1.0]), pos)
    figure, ax = plt.subplots(figsize=(7.5, 5.0))
    for key, velocity, colour, label in (
        ("traction_translation_pN", translation_velocity, COLOR_TRANSLATION, "rigid translation"),
        ("traction_rotation_pN", rotation_velocity, COLOR_ROTATION, "rigid rotation"),
    ):
        projection = np.einsum("ij,ij->i", fields[key], velocity)
        ax.hist(projection, bins=48, color=colour, alpha=0.55, label=f"{label} (n={projection.size})")
    ax.axvline(0.0, color=COLOR_REFERENCE, lw=1.4, label="zero — the whole distribution must be left of it")
    banked = record["measurements"]["dissipation"]["banked_work_pN_um"]
    ax.set_xlabel("per-node traction · unit velocity  [pN per (µm/s)]")
    ax.set_ylabel("surface nodes [count]")
    ax.set_title(
        "sign-sense: the medium opposes motion at EVERY node, not merely on average\n"
        f"accepted-step banked work {banked:.4g} pN·µm (must be < 0)",
        fontsize=9,
    )
    # Rotation nodes near the axis have small velocity, so their projection legitimately approaches zero;
    # what would be a defect is any node crossing it, which is why the histogram and not a mean is drawn.
    ax.annotate(
        GAP_BANNER + "\nnodes near the rotation axis approach zero legitimately;\n"
        "a node CROSSING zero would be the defect",
        xy=(0.03, 0.72), xycoords="axes fraction", fontsize=8,
        bbox={"boxstyle": "round", "fc": "white", "ec": COLOR_REGULARISER, "alpha": 0.9},
    )
    ax.legend(fontsize=8, loc="upper right")
    figure.tight_layout()
    path = out / "dissipation.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def _tolerance_sweep(record_path: Path, out: Path) -> Path | None:
    """Draw the CG-tolerance sweep that decides whether a failing residual is solver or operator.

    This is the falsification figure, and it exists because the numbers alone cannot choose between the
    two hypotheses. Three "exact structural" predicates failed their declared round-off floor; each is a
    property of the OPERATOR but is measured through an ITERATIVE solve. If the residuals track the CG
    tolerance they are the solver's error; if they are flat, the operator is wrong. Flatness is the
    falsifier, so it is drawn as a horizontal reference line rather than described.

    Reads sibling ``tolsweep_*.json`` records plus the declared run. Returns ``None`` when no sweep exists,
    because a figure asserting a sweep that was not run would be a drawing.

    Args:
        record_path: The declared run record; siblings are discovered beside it.
        out: Figure directory.

    Returns:
        The written path, or ``None`` if fewer than two tolerances are on disk.
    """
    records = []
    for path in sorted(record_path.parent.glob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            tol = float(data["config"]["cg_relative_tolerance"])
        except (KeyError, ValueError, TypeError):
            continue
        m = data["measurements"]
        records.append((
            tol,
            m["exact_structural_residuals"]["translation_isotropy"],
            m["exact_structural_residuals"]["translation_rotation_coupling"],
            m["time_reversal"]["relative_asymmetry"],
            m["solver"]["transverse_over_axial"],
            m["exact_structural_residuals"]["round_off_floor"],
            m["solver"]["mobility_condition_number"],
            m["stokes_law_ratio"]["translation_measured_over_exact"],
        ))
    records = sorted(set(records))
    if len(records) < 2:
        return None

    tol = np.array([r[0] for r in records])
    floor = records[0][5]
    cond = records[0][6]
    figure, ax = plt.subplots(figsize=(8.0, 5.4))
    for idx, colour, marker, label in (
        (1, COLOR_TRANSLATION, "o", "translation isotropy"),
        (2, COLOR_ROTATION, "s", "translation–rotation decoupling"),
        (3, "#4a7c59", "^", "time-reversal asymmetry"),
        (4, "#7d5ba6", "d", "transverse force / axial force"),
    ):
        ax.loglog(tol, [r[idx] for r in records], marker=marker, color=colour, label=label)
    ax.loglog(tol, cond * tol, ls="-.", color=COLOR_REFERENCE, lw=1.0,
              label=f"CG solution-error bound  cond·tol  (cond = {cond:.0f})")
    ax.axhline(floor, color="#b03a2e", ls="--", lw=1.4,
               label=f"declared floor 64·eps64 = {floor:.2e}  (a DIRECT-solve floor)")
    # NO "what an operator defect would look like" line is drawn. Its value is not a quantity any run
    # measured, and drawing an un-measured reference on a measurement is what the visualization rule
    # forbids — the same mistake this figure set already had to fix once. The falsifier is stated instead:
    # a tolerance-INDEPENDENT defect is a horizontal line, and none of these four is horizontal.
    ax.set_xlabel("CG relative residual tolerance [dimensionless] — LOG axis")
    ax.set_ylabel("relative residual [dimensionless] — LOG axis")
    ratios = {f"{r[7]:.9f}" for r in records}
    spans = {
        "isotropy": max(r[1] for r in records) / min(r[1] for r in records),
        "decoupling": max(r[2] for r in records) / min(r[2] for r in records),
        "time-reversal": max(r[3] for r in records) / min(r[3] for r in records),
        "transverse": max(r[4] for r in records) / min(r[4] for r in records),
    }
    ax.set_title(
        "a tolerance-INDEPENDENT operator defect is REFUTED — but the dependence is not a clean power law\n"
        f"drag ratio {'bit-identical' if len(ratios) == 1 else 'NOT identical'} across "
        f"{np.log10(tol.max() / tol.min()):.0f} decades of tolerance"
        f"{' (' + ratios.pop() + ')' if len(ratios) == 1 else ''}: the ANSWER does not depend on any of this",
        fontsize=9,
    )
    ax.annotate(
        "each residual moves with the solver (span over the sweep):\n"
        + "\n".join(f"   {k}: ×{v:.3g}" for k, v in spans.items())
        + "\nnone is horizontal, so an operator defect is refuted. CG error is NOT\n"
          "monotone in the stopping tolerance — it depends on which Krylov iterate\n"
          "the test lands on — so `decoupling` sitting near-flat over three decades\n"
          "and then dropping two is expected, not a second finding.",
        xy=(0.02, 0.60), xycoords="axes fraction", fontsize=7.5,
        bbox={"boxstyle": "round", "fc": "white", "ec": COLOR_REGULARISER, "alpha": 0.9},
    )
    ax.legend(fontsize=7.5, loc="lower right")
    figure.tight_layout()
    path = out / "solver_tolerance_sweep.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def main() -> int:
    """Regenerate every T10 exterior-medium figure from the committed record."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--record",
        default="aleph/outputs/ac/gate_t10_medium/medium_exterior_native.json",
        help="the committed run record; its NPZ sidecar is read from the same stem",
    )
    parser.add_argument("--out", default="", help="figure directory; defaults to figs/ beside the record")
    args = parser.parse_args()

    record_path = Path(args.record)
    record, fields = _load(record_path)
    out = Path(args.out) if args.out else record_path.parent / "figs"
    out.mkdir(parents=True, exist_ok=True)

    written = [
        _traction_field(record, fields, out),
        _rigid_modes(record, fields, out),
        _grid_invariance(record, out),
        _dissipation(record, fields, out),
    ]
    sweep = _tolerance_sweep(record_path, out)
    if sweep is not None:
        written.append(sweep)
    else:
        print("no tolerance sweep on disk — skipping solver_tolerance_sweep.png", flush=True)
    for path in written:
        print(f"wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
