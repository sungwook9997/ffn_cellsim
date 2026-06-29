"""H.7 surface-manifold substrate — integrated demo on a REAL cell (capstone).

Ties together the (b) spatial-substrate increments on a real FA-adhered MCF7 cell:

  1. build a small FA-adhered cell at the physiological operating point;
  2. **fit** the surface manifold to the cell's actual cortex beads
     (``SurfaceManifold.fit_to_cloud`` — the slaved-to-beads deformable shell, so
     the substrate tracks the REAL deformed surface, not a static icosphere);
  3. carve the lamellipodium **region masks** (basal_ring / polarized_patch,
     ``cortex/manifold_regions.py``) on the fitted shell;
  4. measure the cell↔ECM **contact manifold** (``bridge/ecm_contact.py``) and the
     FA **traction field** (``bridge/manifold_traction.py``) on the SAME fitted
     shell;
  5. run the cross-consistency gate (FA-engaged patches ⊆ geometric contact).

All four layers (deformable geometry · region masks · contact · traction) live on
ONE co-registered manifold — the shared-frame contact infrastructure the H.7
surface-manifold design calls for. Geometry/measurement only; no force created.

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_manifold_substrate_demo.png`` — south-cap Lambert
    maps on the FITTED shell: regions · ECM contact · FA traction.
  * ``ffn_sim/outputs/h7/figs/h7_manifold_substrate_demo.json`` — integrated metrics.

Usage (smoke, CPU dev):
    python -m ffn_sim.scripts.h7_manifold_substrate_demo --n-filaments 160 \
        --n-nuc-beads 400 --device cpu --allow-cpu-dev --warmup 300 --sample 200
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ffn_sim.bridge.ecm_contact import (  # noqa: E402
    assert_traction_within_contact,
    measure_ecm_contact_manifold,
)
from ffn_sim.bridge.manifold_traction import measure_fa_traction_field  # noqa: E402
from ffn_sim.common.production_policy import (  # noqa: E402
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.manifold_regions import (  # noqa: E402
    basal_ring_region,
    polarized_patch_region,
)
from ffn_sim.common.surface_manifold import SurfaceManifold  # noqa: E402
from ffn_sim.scripts.h7_ecm_contact import _derive_contact_gap  # noqa: E402
from ffn_sim.scripts.h7_manifold_traction import (  # noqa: E402
    _choose_subdivisions,
    _lambert_azimuthal_south,
    build_fa_cell,
)

_FIG_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"


def _cortex_positions(cell) -> np.ndarray:
    """Tag-ordered cortex-actin bead positions [0, n_cortex_actin)."""
    nca = int(cell.p_cortex.n_filaments * cell.p_cortex.beads_per_filament)
    with cell.simulation.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv][:nca].copy()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160)
    ap.add_argument("--n-nuc-beads", type=int, default=400)
    ap.add_argument("--warmup", type=int, default=300)
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--softstart", type=int, default=None)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--subdivisions", type=int, default=None)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    cell, p_fa = build_fa_cell(
        n_filaments=args.n_filaments, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, sample=args.sample, device=dev, seed=args.seed,
        softstart=args.softstart, constrained=False,
    )
    if p_fa is None:
        raise RuntimeError("FA did not resolve (manifest fa.enabled must be True).")
    R_cell = float(cell.p_cortex.R_cell)
    subdiv = int(args.subdivisions) if args.subdivisions is not None else _choose_subdivisions(R_cell, p_fa)

    # --- fit the manifold to the REAL cortex beads (deformable shell) ---
    cortex = _cortex_positions(cell)
    manifold = SurfaceManifold.icosphere(subdivisions=subdiv, radius=R_cell)
    r_before = float(np.mean(np.linalg.norm(manifold.tri_centroids, axis=1)))
    manifold.fit_to_cloud(cortex)
    r_after = float(np.mean(np.linalg.norm(manifold.tri_centroids, axis=1)))

    # --- regions + contact + traction on the SAME fitted shell ---
    basal = basal_ring_region(manifold, R_cell=R_cell, rest_length=p_fa_rest(cell),
                              collar_half_angle=3.0 * manifold.mean_edge_length / R_cell)
    polar = polarized_patch_region(manifold, R_cell=R_cell, half_angle_azimuth=math.pi / 6.0)
    contact_gap = _derive_contact_gap(p_fa, manifold)
    contact = measure_ecm_contact_manifold(cell.simulation, manifold, contact_gap=contact_gap)
    traction = measure_fa_traction_field(
        cell.simulation, manifold, k_int_bare=p_fa.k_int_bare, integrin_r0=p_fa.integrin_r0)
    gate = assert_traction_within_contact(contact, traction)

    _make_figure(manifold, basal, polar, contact, traction, R_cell=R_cell,
                 out_png=_FIG_DIR / "h7_manifold_substrate_demo.png")

    summary = {
        "R_cell_m": R_cell, "subdivisions": subdiv, "n_tri": manifold.n_tri,
        "fit": {"mean_centroid_r_before_m": r_before, "mean_centroid_r_after_m": r_after,
                "cortex_beads": int(cortex.shape[0])},
        "regions": {"basal_ring_patches": int(basal.patch_ids.size),
                    "polarized_patch_patches": int(polar.patch_ids.size),
                    "basal_ring_area_frac": basal.area_fraction,
                    "polarized_patch_area_frac": polar.area_fraction},
        "contact": {"n_ligands": contact["n_ligands"],
                    "n_contact_patches": contact["n_contact_patches"],
                    "contact_area_m2": contact["contact_area_m2"],
                    "partition_passed": contact["partition"]["passed"]},
        "traction": {"n_integrin_bonds": traction["n_integrin_bonds"],
                     "n_basal_patches": traction["n_basal_patches"],
                     "total_tangential_traction_Pa": traction["total_tangential_traction_Pa"]},
        "cross_gate_FA_subset_contact": {"passed": gate["passed"],
                                         "n_fa_outside_contact": gate["n_fa_outside_contact"]},
        "figure": str(_FIG_DIR / "h7_manifold_substrate_demo.png"),
    }
    (_FIG_DIR / "h7_manifold_substrate_demo.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    ok = contact["partition"]["passed"] and gate["passed"]
    print(f"\nsubstrate demo: fit {r_before*1e6:.2f}→{r_after*1e6:.2f}µm; "
          f"regions br={basal.patch_ids.size}/pp={polar.patch_ids.size}; "
          f"contact {contact['n_contact_patches']} patches; FA⊆contact={'PASS' if gate['passed'] else 'FAIL'}")
    return 0 if ok else 1


def p_fa_rest(cell) -> float:
    """Cortex backbone rest length ℓ₀ [m] (region collar scale)."""
    return float(cell.p_cortex.rest_length)


def _make_figure(manifold, basal, polar, contact, traction, *, R_cell, out_png):
    centroids = manifold.tri_centroids
    proj = _lambert_azimuthal_south(centroids, R_cell)
    south = centroids[:, 2] < 0.0
    ok = south & np.isfinite(proj[:, 0])
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 6.2), constrained_layout=True)

    ax = axes[0]
    ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.85")
    ax.scatter(proj[basal.mask & ok, 0], proj[basal.mask & ok, 1], s=18,
               c="tab:blue", label=f"basal_ring ({basal.patch_ids.size})")
    ax.scatter(proj[polar.mask & ok, 0], proj[polar.mask & ok, 1], s=18,
               c="tab:red", label=f"polarized_patch ({polar.patch_ids.size})")
    ax.set_title("(a) regions on fitted shell"); ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.85")
    cm = contact["contact_mask"] & ok
    ax.scatter(proj[cm, 0], proj[cm, 1], s=18, c="tab:green",
               label=f"ECM contact ({contact['n_contact_patches']})")
    ax.set_title(f"(b) ECM contact ({contact['n_ligands']} ligands)"); ax.legend(fontsize=8)

    ax = axes[2]
    tmag = traction["patch_traction_mag_Pa"].copy()
    ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.9")
    bm = traction["basal_patch_mask"] & ok
    if np.any(bm):
        sc = ax.scatter(proj[bm, 0], proj[bm, 1], s=24, c=tmag[bm], cmap="magma")
        fig.colorbar(sc, ax=ax, label="|traction| [Pa]")
    ax.set_title(f"(c) FA traction ({traction['n_integrin_bonds']} bonds)")

    for ax in axes:
        ax.set_aspect("equal"); ax.set_xlabel("Lambert X (south cap)")
    axes[0].set_ylabel("Lambert Y")
    fig.suptitle("H.7 surface-manifold substrate on a real cell — deformable shell + "
                 "regions + ECM contact + FA traction (co-registered)", fontsize=11)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
