"""H.7 cell↔ECM CONTACT manifold — driver + viz (PI direction (b)).

The geometry layer beneath the FA-traction field: builds an FA-adhered MCF7 cell
at the physiological operating point, attaches a surface manifold, and measures
the cell↔ECM **contact relation** (``bridge/ecm_contact.py``) — which patches are
in adhesive contact with the substrate ligands, the cell-surface↔ligand
correspondence, contact gap + area — then measures the FA traction field
(``bridge/manifold_traction.py``) on the SAME manifold and runs the
cross-consistency GATE: every FA-engaged patch must lie inside the geometric
contact footprint (a bond can only form where a ligand is reachable). A gate
failure exits nonzero.

Geometry only: the manifold supplies frames/areas, the ECM supplies ligand
positions, the FA bonds supply forces; nothing here creates a force.

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_ecm_contact.png`` — south-cap Lambert maps:
    (a) contact patches + ligand home patches, (b) per-patch contact gap, and
    (c) the FA-engaged patches overlaid on the contact footprint (the gate).
  * ``ffn_sim/outputs/h7/figs/h7_ecm_contact.json`` — scalar summary + gate.

Usage (smoke, CPU dev):
    python -m ffn_sim.scripts.h7_ecm_contact --n-filaments 160 --n-nuc-beads 400 \
        --device cpu --allow-cpu-dev --warmup 300 --sample 200
Usage (full scale, gbook GPU):
    python -m ffn_sim.scripts.h7_ecm_contact --device gpu
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ffn_sim.archive.hoomd_legacy.bridge.ecm_contact import (  # noqa: E402
    TYPE_LIGAND,
    assert_traction_within_contact,
    measure_ecm_contact_manifold,
)
from ffn_sim.archive.hoomd_legacy.bridge.fa import TYPE_LIGAND as FA_TYPE_LIGAND  # noqa: E402
from ffn_sim.archive.hoomd_legacy.bridge.manifold_traction import (  # noqa: E402
    measure_fa_traction_field,
)
from ffn_sim.common.production_policy import (  # noqa: E402
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.common.surface_manifold import SurfaceManifold  # noqa: E402
from ffn_sim.scripts.h7_manifold_traction import (  # noqa: E402
    _choose_subdivisions,
    _lambert_azimuthal_south,
    build_fa_cell,
)

# The contact module's ligand-type literal must match the FA single-source-of-truth.
assert TYPE_LIGAND == FA_TYPE_LIGAND, (
    "ecm_contact.TYPE_LIGAND must match bridge.fa.TYPE_LIGAND"
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"


def _derive_contact_gap(p_fa, manifold: SurfaceManifold) -> float:
    """DERIVE the contact reach from FA geometry + manifold slack (no magic #).

    A ligand is "reachable" by a patch when it lies within the integrin capture
    radius + the integrin standoff above the substrate (the physics reach), plus
    the geometric slack that a ligand can sit a patch circumradius from the patch
    centroid (the broad-phase home-patch metric). All three are derived: the FA
    reach from the FA params, the slack from the manifold geometry.
    """
    R_FA = float(getattr(p_fa, "capture_radius_R_FA", 0.0))
    h_int = float(getattr(p_fa, "h_integrin_above_substrate", 0.0))
    fa_reach = R_FA + h_int
    return fa_reach + float(manifold.max_circumradius)


def _make_figure(
    contact: dict,
    traction: dict,
    *,
    manifold: SurfaceManifold,
    R_cell: float,
    out_png: Path,
) -> None:
    centroids = manifold.tri_centroids
    proj = _lambert_azimuthal_south(centroids, R_cell)
    south = centroids[:, 2] < 0.0
    ok = south & np.isfinite(proj[:, 0])

    fig, axes = plt.subplots(1, 3, figsize=(18.0, 6.2), constrained_layout=True)

    # (a) contact patches + ligand home patches.
    ax = axes[0]
    ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.85", label="patch")
    cmask = contact["contact_mask"] & ok
    ax.scatter(proj[cmask, 0], proj[cmask, 1], s=20, c="tab:green",
               label=f"contact ({contact['n_contact_patches']})")
    home = np.unique(contact["ligand_home_patch"]) if contact["n_ligands"] else np.array([], int)
    home = home[np.isin(home, np.flatnonzero(ok))]
    if home.size:
        ax.scatter(proj[home, 0], proj[home, 1], s=24, facecolors="none",
                   edgecolors="k", linewidths=0.8, label="ligand home patch")
    ax.set_title(f"(a) ECM contact footprint — area={contact['contact_area_m2']*1e12:.1f} µm²")
    ax.legend(loc="upper right", fontsize=8)

    # (b) per-patch contact gap [nm] (finite only).
    ax = axes[1]
    gap = contact["patch_min_gap_m"].copy()
    finite = np.isfinite(gap) & ok
    sc = ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.9")
    if np.any(finite):
        sc = ax.scatter(proj[finite, 0], proj[finite, 1], s=20,
                        c=gap[finite] * 1e9, cmap="viridis")
        fig.colorbar(sc, ax=ax, label="nearest-ligand gap [nm]")
    ax.set_title(f"(b) contact gap (mean={contact['mean_contact_gap_m']*1e9:.1f} nm)")

    # (c) FA-engaged patches ON the contact footprint (the gate).
    ax = axes[2]
    ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.85")
    ax.scatter(proj[cmask, 0], proj[cmask, 1], s=20, c="tab:green", label="contact")
    fa = traction["basal_patch_mask"] & ok
    ax.scatter(proj[fa, 0], proj[fa, 1], s=26, c="tab:red", marker="x",
               label=f"FA-engaged ({int(np.count_nonzero(traction['basal_patch_mask']))})")
    ax.set_title("(c) FA-engaged ⊆ contact (cross-consistency gate)")
    ax.legend(loc="upper right", fontsize=8)

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("Lambert X (south cap)")
    axes[0].set_ylabel("Lambert Y")

    fig.suptitle(
        f"H.7 cell↔ECM contact manifold  (R_cell={R_cell*1e6:.1f} µm, "
        f"{manifold.n_tri} patches, {contact['n_ligands']} ligands)",
        fontsize=11,
    )
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160)
    ap.add_argument("--n-nuc-beads", type=int, default=400)
    ap.add_argument("--warmup", type=int, default=300)
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--softstart", type=int, default=None)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--subdivisions", type=int, default=None)
    ap.add_argument("--constrained", action="store_true")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd

    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    cell, p_fa = build_fa_cell(
        n_filaments=args.n_filaments, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, sample=args.sample, device=dev, seed=args.seed,
        softstart=args.softstart, constrained=bool(args.constrained),
    )
    if p_fa is None:
        raise RuntimeError(
            "resolved p_fa is None — FA adhesion did not resolve (manifest "
            "fa.enabled must be True); cannot measure the contact manifold."
        )
    R_cell = float(cell.p_cortex.R_cell)
    subdivisions = (
        int(args.subdivisions) if args.subdivisions is not None
        else _choose_subdivisions(R_cell, p_fa)
    )
    manifold = SurfaceManifold.icosphere(subdivisions=subdivisions, radius=R_cell)

    contact_gap = _derive_contact_gap(p_fa, manifold)
    contact = measure_ecm_contact_manifold(
        cell.simulation, manifold, contact_gap=contact_gap,
    )
    traction = measure_fa_traction_field(
        cell.simulation, manifold,
        k_int_bare=p_fa.k_int_bare, integrin_r0=p_fa.integrin_r0,
    )
    gate = assert_traction_within_contact(contact, traction)

    out_png = _OUT_DIR / "h7_ecm_contact.png"
    _make_figure(contact, traction, manifold=manifold, R_cell=R_cell, out_png=out_png)

    summary = {
        "R_cell_m": R_cell,
        "subdivisions": subdivisions,
        "n_tri": manifold.n_tri,
        "contact_gap_m": contact_gap,
        "n_ligands": contact["n_ligands"],
        "n_ligands_in_contact": contact["n_ligands_in_contact"],
        "n_contact_patches": contact["n_contact_patches"],
        "contact_area_m2": contact["contact_area_m2"],
        "contact_area_fraction": contact["contact_area_fraction"],
        "mean_contact_gap_m": contact["mean_contact_gap_m"],
        "partition_passed": contact["partition"]["passed"],
        "n_fa_patches": gate["n_fa_patches"],
        "n_integrin_bonds": traction["n_integrin_bonds"],
        "cross_consistency_gate": {
            "passed": gate["passed"],
            "n_fa_outside_contact": gate["n_fa_outside_contact"],
            "fa_outside_ids": gate["fa_outside_ids"].tolist(),
        },
        "figure": str(out_png),
    }
    out_json = _OUT_DIR / "h7_ecm_contact.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    ok = bool(contact["partition"]["passed"] and gate["passed"])
    print(
        f"\npartition gate: {'PASS' if contact['partition']['passed'] else 'FAIL'}; "
        f"FA⊆contact gate: {'PASS' if gate['passed'] else 'FAIL'} "
        f"({gate['n_fa_outside_contact']} FA patches outside contact)"
    )
    print(f"figure → {out_png}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
