"""H.7 manifold FA-traction FIELD — driver + viz (PI direction (b)).

Wires the surface manifold (a shared cell-surface coordinate system, GEOMETRY
ONLY) as the spatial substrate for focal-adhesion traction: it builds an
FA-adhered MCF7 cell at the physiological operating point, attaches a
:class:`~ffn_sim.cortex.surface_manifold.SurfaceManifold` icosphere at
``R_cell``, and uses :func:`ffn_sim.bridge.manifold_traction.measure_fa_traction_field`
to turn the EXPLICIT integrin↔ligand bond forces into a per-contact-patch
**traction vector field** (the design's first-class output for an adherent cell).

The manifold supplies geometry (patches + local frames + areas); ALL forces come
from the explicit FA bonds (``F_mag = k_int_bare·clip(|Δr|−integrin_r0, 0)`` along
the integrin→ligand unit vector). The contact-conservation gate (Σ per-patch
force ≡ Σ per-bond force) is asserted; a failure exits nonzero.

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_manifold_traction.png`` — (a) basal-cap
    traction-magnitude heatmap (Lambert azimuthal projection) + (b) basal-patch
    traction vectors + (c) the φ-azimuthal traction distribution (FA polarity).
  * ``ffn_sim/outputs/h7/figs/h7_manifold_traction.npz`` — per-patch fields
    (ParaView/Blender export handle).
  * ``ffn_sim/outputs/h7/figs/h7_manifold_traction.json`` — scalar summary.

Usage (smoke, CPU dev — small cell, some integrins bind):
    python -m ffn_sim.scripts.h7_manifold_traction --n-filaments 160 \
        --n-nuc-beads 400 --device cpu --allow-cpu-dev --warmup 300 --sample 200
Usage (full production scale, gbook GPU):
    python -m ffn_sim.scripts.h7_manifold_traction --device gpu
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

from ffn_sim.bridge.fa import (  # cross-check the type literals against the SoT
    BOND_TYPE_INTEGRIN as FA_BOND_INTEGRIN,
)
from ffn_sim.bridge.manifold_traction import (
    BOND_TYPE_INTEGRIN,
    measure_fa_traction_field,
)
from ffn_sim.cell.manifest import (
    build_baseline_cell,
    load_manifest,
    resolve_baseline,
)
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.surface_manifold import (
    SurfaceManifold,
    n_tri_for_subdivisions,
)

# The driver depends on the bridge/cell type literals agreeing. Assert it once so
# a future rename in fa.py cannot silently desync the measurement module.
assert BOND_TYPE_INTEGRIN == FA_BOND_INTEGRIN, (
    "manifold_traction.BOND_TYPE_INTEGRIN must match bridge.fa.BOND_TYPE_INTEGRIN"
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
_PA_PER_PA = 1.0  # SI default; kept explicit for the axis labels


# ---------------------------------------------------------------------------
# Build + measure
# ---------------------------------------------------------------------------
def build_fa_cell(
    *,
    n_filaments: int | None,
    n_nuc_beads: int | None,
    warmup: int,
    sample: int,
    device,
    seed: int,
    softstart: int | None,
    constrained: bool,
):
    """Build + settle an FA-adhered MCF7 cell (mirrors h7_gate_b's FA-enable).

    Unconstrained-soft (``constrained=False``) is the stable default here — we
    only need the FA bond forces, not the rigid Lagrange channel. Equilibration
    (with a force-ramped softstart) settles the cell onto the substrate so some
    integrins bind before the measurement.

    Returns:
        ``(cell, p_fa)`` — the built :class:`Cell` and the resolved FA params
        (``ResolvedH4``) read from the same manifest, so ``k_int_bare`` /
        ``integrin_r0`` match the bonds the cell was wired with.
    """
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest["optional_subsystems"]["fa"]["enabled"] = True  # adhered op point
    if n_filaments is not None:
        co = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
        co["n_filaments"] = int(n_filaments)
        co["demo_mode"] = True
    if n_nuc_beads is not None:
        manifest["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)

    # Resolve the FA params from the SAME manifest the cell is built from. The
    # built ``Cell`` does not surface ``p_fa``, so we read the resolved FA params
    # here (build_baseline_cell resolves the identical baseline internally, so
    # k_int_bare / integrin_r0 match the bonds the cell was wired with bit-for-bit).
    p_fa = resolve_baseline(manifest).p_fa

    cell = build_baseline_cell(
        manifest=manifest,
        device=device,
        seed=seed,
        constrained=constrained,
        equilibrate=True,
        equilibrate_steps=warmup,
        equilibrate_softstart_steps=(
            softstart if softstart is not None else max(100, warmup // 4)
        ),
    )
    if sample > 0:
        cell.simulation.run(sample)
    return cell, p_fa


def _choose_subdivisions(R_cell: float, p_fa) -> int:
    """DERIVE the icosphere subdivision so a patch resolves an FA footprint.

    Not a magic number: we want the mean patch edge ``≈ R_cell·(edge/R for an
    icosphere)`` to be at or below the FA footprint scale ``√A_mature`` so each
    contact patch resolves (at most) ~one FA. The level-s icosphere mean edge
    scales ≈ ``1.05·R / 2ˢ`` (geodesic dome heuristic), so

        ``2ˢ ≳ 1.05·R_cell / √A_mature``  →  ``s = ceil(log2(...))``,

    clamped to the canonical {2,3,4} resolution grid (320/1280/5120 faces). The
    FA footprint comes from the FA params (geometry of the physics), the mesh
    scale from R_cell — both derived, neither tuned to a gate. Resolution
    sensitivity is reported (``subdivisions_hint``) and can be swept via
    ``--subdivisions``.
    """
    import math

    A_mature = float(getattr(p_fa, "A_mature", 0.0))
    fa_footprint = math.sqrt(A_mature) if A_mature > 0.0 else R_cell
    ratio = 1.05 * R_cell / max(fa_footprint, 1.0e-12)
    s = max(2, math.ceil(math.log2(max(ratio, 1.0))))
    return int(min(s, 4))  # cap at the canonical fine grid (5120 faces)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------
def _lambert_azimuthal_south(centroids: np.ndarray, R: float) -> np.ndarray:
    """Lambert azimuthal equal-area projection of the SOUTH cap (z<0) to 2D.

    The basal contact sits on the south pole (the cell rests on the z=0
    substrate; the cortex south pole is at z ≈ −R_cell). Equal-area so patch
    areas are not distorted in the heatmap. Returns ``(n, 2)`` projected coords
    (NaN for north-hemisphere patches, which the basal cap never uses).
    """
    x, y, z = centroids[:, 0], centroids[:, 1], centroids[:, 2]
    rho = np.linalg.norm(centroids, axis=1)
    rho = np.where(rho > 0.0, rho, 1.0)
    zc = z / rho  # cos(colatitude from +z); south pole zc = −1
    # South-cap Lambert (project from the NORTH pole): k = sqrt(2/(1 − zc)).
    denom = 1.0 - zc
    denom = np.where(denom > 1.0e-12, denom, np.nan)
    k = np.sqrt(2.0 / denom)
    X = k * (x / rho) * R
    Y = k * (y / rho) * R
    return np.stack([X, Y], axis=1)


def _make_figure(field: dict, *, R_cell: float, out_png: Path, title: str) -> None:
    """Three-panel FA-traction figure (heatmap / vectors / φ-distribution)."""
    centroids = field["patch_centroids"]
    basal_mask = field["basal_patch_mask"]
    tmag = field["patch_traction_mag_Pa"]
    tn = field["patch_normal_traction_Pa"]
    fvec = field["patch_force_vec_N"]

    fig = plt.figure(figsize=(18.0, 6.2), constrained_layout=True)
    fig.suptitle(title, fontsize=12)

    # --- (a) basal-cap traction-magnitude heatmap (Lambert equal-area) -----
    ax0 = fig.add_subplot(1, 3, 1)
    proj = _lambert_azimuthal_south(centroids, R_cell)
    bm = basal_mask & np.isfinite(proj[:, 0])
    if np.any(bm):
        sc = ax0.scatter(
            proj[bm, 0] * 1e6,
            proj[bm, 1] * 1e6,
            c=tmag[bm],
            s=40,
            cmap="viridis",
            edgecolors="k",
            linewidths=0.3,
        )
        cb = fig.colorbar(sc, ax=ax0)
        cb.set_label("|traction|  [Pa]")
    else:
        ax0.text(0.5, 0.5, "no basal contact patches", ha="center",
                 va="center", transform=ax0.transAxes)
    ax0.set_aspect("equal")
    ax0.set_xlabel("Lambert X  [µm]")
    ax0.set_ylabel("Lambert Y  [µm]")
    ax0.set_title("(a) basal-cap traction magnitude\n(equal-area south-pole projection)")

    # --- (b) basal-patch traction VECTORS (in-plane component, top view) ----
    ax1 = fig.add_subplot(1, 3, 2)
    if np.any(basal_mask):
        cx = centroids[basal_mask, 0] * 1e6
        cy = centroids[basal_mask, 1] * 1e6
        # In-plane (xy) force component → the shear the cell feels.
        fx = fvec[basal_mask, 0]
        fy = fvec[basal_mask, 1]
        col = tn[basal_mask]  # colour by signed normal traction
        q = ax1.quiver(
            cx, cy, fx, fy, col,
            cmap="coolwarm", angles="xy", scale_units="xy",
            pivot="tail", width=0.005,
        )
        cb1 = fig.colorbar(q, ax=ax1)
        cb1.set_label("normal traction  n̂·ΣF / A  [Pa]\n(+ outward / − inward)")
        # Annotate the dominant vector scale.
        fmax = float(np.max(np.hypot(fx, fy))) if fx.size else 0.0
        ax1.set_title(
            f"(b) basal-patch traction vectors (xy force)\nmax |F_xy| = {fmax:.2e} N"
        )
    else:
        ax1.text(0.5, 0.5, "no basal contact patches", ha="center",
                 va="center", transform=ax1.transAxes)
        ax1.set_title("(b) basal-patch traction vectors")
    ax1.set_aspect("equal")
    ax1.set_xlabel("x  [µm]")
    ax1.set_ylabel("y  [µm]")

    # --- (c) φ-azimuthal traction distribution (FA polarity) ---------------
    ax2 = fig.add_subplot(1, 3, 3, projection="polar")
    edges = field["phi_bin_edges_rad"]
    phi_t = field["phi_traction_mag_Pa"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    ax2.bar(centers, phi_t, width=width, bottom=0.0, color="tab:orange",
            edgecolor="k", linewidth=0.3, alpha=0.85)
    ax2.set_theta_zero_location("E")
    ax2.set_theta_direction(1)
    ax2.set_title(
        "(c) φ-azimuthal |traction|  [Pa]\n"
        f"polarity={field['polarity_magnitude']:.3f} "
        f"@ {np.degrees(field['polarity_angle_rad']):.0f}°",
        pad=18,
    )
    # Overlay the net polarity heading as a radial arrow.
    pol = field["polarity_magnitude"]
    if pol > 0.0 and np.any(phi_t > 0):
        ax2.annotate(
            "", xy=(field["polarity_angle_rad"], float(np.max(phi_t))),
            xytext=(field["polarity_angle_rad"], 0.0),
            arrowprops=dict(arrowstyle="-|>", color="navy", lw=2.0),
        )

    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _report(field: dict, *, manifold: SurfaceManifold, subdivisions: int) -> bool:
    """Print the gate result + totals; return the conservation pass flag."""
    cc = field["contact_conservation"]
    print("=" * 70, flush=True)
    print("H.7 MANIFOLD FA-TRACTION FIELD  (manifold = GEOMETRY only;", flush=True)
    print("                                 forces = explicit integrin bonds)", flush=True)
    print("-" * 70, flush=True)
    print(f"  manifold icosphere: subdivisions={subdivisions}, "
          f"n_tri={manifold.n_tri} (expected {n_tri_for_subdivisions(subdivisions)}), "
          f"R={manifold.R*1e6:.2f} µm", flush=True)
    print(f"  mean patch edge   : {manifold.mean_edge_length*1e6:.3f} µm", flush=True)
    print(f"  resolution note   : {field['subdivisions_hint']}", flush=True)
    print("-" * 70, flush=True)
    print(f"  integrin_ligand bonds      : {field['n_integrin_bonds']} "
          f"({field['n_bonds_loaded']} load-bearing, F>0)", flush=True)
    print(f"  basal contact patches      : {field['n_basal_patches']} "
          f"/ {manifold.n_tri}", flush=True)
    print(f"  contact area               : {field['contact_area_m2']*1e12:.3f} µm²", flush=True)
    print("-" * 70, flush=True)
    print("  TOTALS (cell-side traction, normal n̂·ΣF signed + / out, − / in):", flush=True)
    print(f"    total NORMAL  force    = {field['total_normal_force_N']:+.4e} N", flush=True)
    print(f"    total TANGENT force    = {field['total_tangential_force_N']:+.4e} N", flush=True)
    print(f"    total NORMAL  traction = {field['total_normal_traction_Pa']:+.4e} Pa", flush=True)
    print(f"    total TANGENT traction = {field['total_tangential_traction_Pa']:+.4e} Pa", flush=True)
    fv = field["total_force_vec_N"]
    print(f"    Σ bond force vector    = "
          f"[{fv[0]:+.3e}, {fv[1]:+.3e}, {fv[2]:+.3e}] N", flush=True)
    print("-" * 70, flush=True)
    print(f"  FA polarity (in-plane)     : |p|={field['polarity_magnitude']:.4f} "
          f"@ {np.degrees(field['polarity_angle_rad']):+.1f}°  (0=isotropic, 1=fully polar)",
          flush=True)
    print("-" * 70, flush=True)
    status = "PASS ✓" if cc["passed"] else "FAIL ✗"
    print(f"  CONTACT-CONSERVATION GATE  : {status}", flush=True)
    print(f"    Σ per-patch force = [{cc['sum_patch_force_N'][0]:+.3e}, "
          f"{cc['sum_patch_force_N'][1]:+.3e}, {cc['sum_patch_force_N'][2]:+.3e}] N", flush=True)
    print(f"    Σ per-bond  force = [{cc['sum_bond_force_N'][0]:+.3e}, "
          f"{cc['sum_bond_force_N'][1]:+.3e}, {cc['sum_bond_force_N'][2]:+.3e}] N", flush=True)
    print(f"    max |residual|    = {cc['max_abs_residual_N']:.3e} N "
          f"(atol={cc['atol']:.1e}, rtol={cc['rtol']:.1e})", flush=True)
    print("    (the manifold only RE-BINS the explicit bond forces — it must", flush=True)
    print("     neither create nor lose force.)", flush=True)
    print("=" * 70, flush=True)
    return bool(cc["passed"])


def _save_exports(field: dict, *, out_npz: Path, out_json: Path,
                  manifold: SurfaceManifold, subdivisions: int, meta: dict) -> None:
    """Write per-patch fields (.npz, ParaView/Blender handle) + scalar summary (.json)."""
    np.savez_compressed(
        out_npz,
        patch_centroids_m=field["patch_centroids"],
        patch_normals=np.asarray(manifold.tri_normals),
        patch_e1=np.asarray(manifold.tri_e1),
        patch_e2=np.asarray(manifold.tri_e2),
        patch_area_m2=field["patch_area_m2"],
        patch_force_vec_N=field["patch_force_vec_N"],
        patch_normal_force_N=field["patch_normal_force_N"],
        patch_tangential_force_N=field["patch_tangential_force_N"],
        patch_normal_traction_Pa=field["patch_normal_traction_Pa"],
        patch_tangential_traction_Pa=field["patch_tangential_traction_Pa"],
        patch_traction_mag_Pa=field["patch_traction_mag_Pa"],
        basal_patch_mask=field["basal_patch_mask"],
        basal_patch_ids=field["basal_patch_ids"],
        phi_bin_edges_rad=field["phi_bin_edges_rad"],
        phi_traction_mag_Pa=field["phi_traction_mag_Pa"],
        phi_force_mag_N=field["phi_force_mag_N"],
        tris=np.asarray(manifold.tris),
        verts_m=np.asarray(manifold.verts),
    )
    summary = {
        "subdivisions": int(subdivisions),
        "n_tri": int(manifold.n_tri),
        "R_cell_m": float(manifold.R),
        "mean_patch_edge_m": float(manifold.mean_edge_length),
        "n_integrin_bonds": int(field["n_integrin_bonds"]),
        "n_bonds_loaded": int(field["n_bonds_loaded"]),
        "n_basal_patches": int(field["n_basal_patches"]),
        "contact_area_m2": float(field["contact_area_m2"]),
        "total_normal_force_N": float(field["total_normal_force_N"]),
        "total_tangential_force_N": float(field["total_tangential_force_N"]),
        "total_normal_traction_Pa": float(field["total_normal_traction_Pa"]),
        "total_tangential_traction_Pa": float(field["total_tangential_traction_Pa"]),
        "total_force_vec_N": [float(v) for v in field["total_force_vec_N"]],
        "polarity_magnitude": float(field["polarity_magnitude"]),
        "polarity_angle_deg": float(np.degrees(field["polarity_angle_rad"])),
        "contact_conservation_passed": bool(field["contact_conservation"]["passed"]),
        "contact_conservation_max_abs_residual_N": float(
            field["contact_conservation"]["max_abs_residual_N"]
        ),
        "subdivisions_hint": field["subdivisions_hint"],
        **meta,
    }
    out_json.write_text(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=None,
                    help="cortex filament count override (omit for full x40 scale)")
    ap.add_argument("--n-nuc-beads", type=int, default=None)
    ap.add_argument("--warmup", type=int, default=300, help="equilibration baoab steps")
    ap.add_argument("--sample", type=int, default=200, help="sample steps at op point")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--softstart", type=int, default=None,
                    help="force-ramped softstart steps (default warmup//4)")
    ap.add_argument("--constrained", action="store_true",
                    help="use the rigid M-SHAKE backbone (default unconstrained-soft, "
                         "which is stable and sufficient — we only read FA bond forces)")
    ap.add_argument("--subdivisions", type=int, default=None,
                    help="icosphere subdivision level (default DERIVED from R_cell vs the "
                         "FA footprint √A_mature; sweep to check resolution sensitivity)")
    ap.add_argument("--phi-bins", type=int, default=36,
                    help="azimuthal histogram bins for the FA-polarity distribution")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    cell, p_fa = build_fa_cell(
        n_filaments=args.n_filaments, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, sample=args.sample, device=dev, seed=args.seed,
        softstart=args.softstart, constrained=args.constrained,
    )
    if p_fa is None:
        raise RuntimeError(
            "resolved p_fa is None — FA adhesion did not resolve; cannot measure "
            "the FA traction field. (Manifest fa.enabled must be True.)"
        )
    R_cell = float(cell.p_cortex.R_cell)

    subdivisions = (
        args.subdivisions if args.subdivisions is not None
        else _choose_subdivisions(R_cell, p_fa)
    )
    # Manifold = GEOMETRY ONLY (icosphere at the cell radius). The slaved-shell
    # set_verts(bead_cloud) update is intentionally NOT applied here: it needs a
    # vertex-count-matched deformed cloud (manifold has 10·4ˢ+2 verts, which does
    # not match the cortex actin bead count), and for an isotropic spread cell the
    # R_cell icosphere is the correct contact geometry. set_verts remains available
    # for a future slaved shell; it carries no mechanics either way.
    manifold = SurfaceManifold.icosphere(subdivisions=subdivisions, radius=R_cell)

    field = measure_fa_traction_field(
        cell.simulation, manifold,
        k_int_bare=float(p_fa.k_int_bare),
        integrin_r0=float(p_fa.integrin_r0),
        n_phi_bins=int(args.phi_bins),
    )

    passed = _report(field, manifold=manifold, subdivisions=subdivisions)

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = _OUT_DIR / "h7_manifold_traction.png"
    out_npz = _OUT_DIR / "h7_manifold_traction.npz"
    out_json = _OUT_DIR / "h7_manifold_traction.json"
    scale = ("FULL x40 scale" if args.n_filaments is None
             else f"smoke (n_filaments={args.n_filaments})")
    title = (
        f"H.7 manifold FA-traction field — {scale}\n"
        f"{field['n_integrin_bonds']} integrin bonds, "
        f"{field['n_basal_patches']} basal patches, "
        f"conservation={'PASS' if passed else 'FAIL'}"
    )
    _make_figure(field, R_cell=R_cell, out_png=out_png, title=title)
    meta = {
        "n_filaments": int(cell.p_cortex.n_filaments),
        "is_full_scale": args.n_filaments is None,
        "device": args.device,
        "seed": int(args.seed),
        "k_int_bare_N_per_m": float(p_fa.k_int_bare),
        "integrin_r0_m": float(p_fa.integrin_r0),
    }
    _save_exports(field, out_npz=out_npz, out_json=out_json,
                  manifold=manifold, subdivisions=subdivisions, meta=meta)
    print(f"  wrote figure : {out_png}", flush=True)
    print(f"  wrote fields : {out_npz}", flush=True)
    print(f"  wrote summary: {out_json}", flush=True)

    # CONTACT-CONSERVATION GATE: a failure means the binning lost/created force.
    # Do NOT silently accept it — exit nonzero (gate, not a fudge).
    if not passed:
        print("  CONTACT-CONSERVATION GATE FAILED — exiting nonzero.", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
