"""Native full-extent ECM build + GPU indentation (5R×5R, R=7.5µm → 75×75 µm) on the A5000.

The CPU validation uses converged REVs (modulus is intensive). This builds the ECM at the PI-specified
5R×5R lateral extent with a practical thickness and runs the indentation at NATIVE scale on cuda:0, where
a large bead (R ≫ mesh ξ) engages many fibers — so the fibrillar indentation E_eff converges up toward the
bulk shear modulus (resolving the local-vs-bulk softening seen at small bead/mesh ratio). Saves the
full-extent network npz for the standalone ECM viewer.

Run on gbook:  cd ~/ff_scratch && python -m aleph.scripts.ff_ecm_native --device cuda:0
Out:  aleph/outputs/ff/ecm_lib/native/{ecm_native_metrics.json, collagen_iso_native.npz, collagen_aln_native.npz}
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from aleph.laws import ecm_library as L
from aleph.laws import ecm_mechanics as M

OUT = "aleph/outputs/ff/ecm_lib/native"
R_CELL = 7.5


def save_npz(ecm, path):
    """Persist a full-extent ECM network for the viewer (positions, segments, crosslinks, pinned, box)."""
    np.savez_compressed(
        path, pos=ecm.net.pos.astype(np.float32), segments=ecm.net.segments.astype(np.int32),
        fiber_offsets=ecm.net.fiber_offsets.astype(np.int32),
        xl_i=ecm.xl_i.astype(np.int32), xl_j=ecm.xl_j.astype(np.int32),
        seg_i=ecm.seg_i.astype(np.int32), seg_j=ecm.seg_j.astype(np.int32),
        pinned=ecm.pinned, box_lo=ecm.box_lo, box_hi=ecm.box_hi,
        material=ecm.material, mesh=ecm.mesh_size_um, z=ecm.connectivity_z, S=ecm.S_measured,
        is_gel=(ecm.net.bend_triples.shape[0] == 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--conc", type=float, default=2.0, help="collagen mg/mL")
    ap.add_argument("--extent", type=float, default=5.0, help="lateral extent in cell radii (5 = 5R×5R)")
    ap.add_argument("--thick", type=float, default=24.0, help="slab thickness [µm]")
    ap.add_argument("--steps", type=int, default=5000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    # PI "세포 크기의 5배x5배" → lateral side = extent × cell diameter (5 × 15µm = 75µm)
    side = a.extent * 2.0 * R_CELL
    lo, hi = [0, 0, 0], [side, side, a.thick]
    spec = L.get_spec("collagen_I")
    res = {"device": a.device, "side_um": side, "thick_um": a.thick, "conc": a.conc, "R_cell": R_CELL}

    # ── build isotropic + aligned full-extent environments ──
    t0 = time.time()
    rng = np.random.default_rng(1)
    iso = L.build_fibrillar_ecm(spec, lo, hi, concentration=a.conc, dim=3, alignment_S=0.0,
                                pin_faces=("z_lo",), pin_margin_um=2.0, rng=rng)
    rng = np.random.default_rng(2)
    aln = L.build_fibrillar_ecm(spec, lo, hi, concentration=a.conc, dim=3, alignment_S=0.6,
                                director=(1, 0, 0), pin_faces=("z_lo",), pin_margin_um=2.0, rng=rng)
    save_npz(iso, f"{OUT}/collagen_iso_native.npz")
    save_npz(aln, f"{OUT}/collagen_aln_native.npz")
    res["build"] = {"iso_nodes": iso.meta["n_nodes"], "iso_nfib": iso.meta["n_fibers"],
                    "iso_mesh": iso.mesh_size_um, "iso_z": iso.connectivity_z, "iso_S": iso.S_measured,
                    "aln_nodes": aln.meta["n_nodes"], "aln_S": aln.S_measured, "build_s": time.time() - t0}
    print(f"built iso {iso.meta['n_nodes']} nodes / aln {aln.meta['n_nodes']} nodes  ({time.time()-t0:.0f}s)",
          flush=True)

    # ── native bulk shear modulus (converged) ──
    t0 = time.time()
    G = M.shear_modulus(iso, gamma=0.02, n_steps=a.steps, device=a.device)["G_Pa"]
    res["bulk_shear_G_Pa"] = G
    E_expect = 2 * (1 + spec.poisson) * G
    print(f"native bulk shear G={G:.1f} Pa  (Hertz-expect E≈{E_expect:.0f})  ({time.time()-t0:.0f}s)", flush=True)

    # ── indenter-size sweep: E_eff → bulk as R_ind/mesh grows (local→continuum) ──
    res["indent_sweep"] = []
    for R in (6.0, 12.0, 20.0, 30.0):
        t0 = time.time()
        ind = M.indentation_modulus(iso, indenter_R_um=R, max_depth_um=min(0.1 * R, 0.15 * a.thick),
                                    n_depths=5, k_ind=3.0e3, relax_steps=a.steps, device=a.device)
        row = {"R_ind_um": R, "R_over_mesh": R / iso.mesh_size_um, "E_eff_Pa": ind["E_eff_Pa"],
               "E_over_bulkE": ind["E_eff_Pa"] / max(E_expect, 1e-9), "contacts": ind["contacts"]}
        res["indent_sweep"].append(row)
        print(f"  indent R={R:4.0f} (R/ξ={R/iso.mesh_size_um:.0f}): E_eff={ind['E_eff_Pa']:.1f} Pa "
              f"({ind['E_eff_Pa']/E_expect:.2f}× bulk-E)  ({time.time()-t0:.0f}s)", flush=True)

    # ── CONFINED indentation: thin slab so the bead compresses the network against the substrate.
    # Sub-isostatic networks are soft under LOCAL (thick-slab) indentation but stiffen toward bulk when
    # confined — the shear/compression mode-decoupling of semiflexible networks (van Oosten 2019). Sweep
    # slab thickness to show E_eff climb from the soft local response toward the bulk modulus.
    res["confinement_sweep"] = []
    for thick in (6.0, 10.0, 16.0, 24.0):
        rng = np.random.default_rng(4)
        slab = L.build_fibrillar_ecm(spec, [0, 0, 0], [side, side, thick], concentration=a.conc, dim=3,
                                     alignment_S=0.0, pin_faces=("z_lo",), pin_margin_um=1.5, rng=rng)
        ind = M.indentation_modulus(slab, indenter_R_um=20.0, max_depth_um=0.18 * thick, n_depths=5,
                                    k_ind=3.0e3, relax_steps=a.steps, device=a.device)
        row = {"thick_um": thick, "E_eff_Pa": ind["E_eff_Pa"], "E_over_bulkE": ind["E_eff_Pa"] / max(E_expect, 1e-9)}
        res["confinement_sweep"].append(row)
        print(f"  confined thick={thick:4.0f}µm: E_eff={ind['E_eff_Pa']:.1f} Pa ({ind['E_eff_Pa']/E_expect:.2f}× bulk-E)",
              flush=True)

    with open(f"{OUT}/ecm_native_metrics.json", "w") as f:
        json.dump(res, f, indent=2, default=float)
    print(f"wrote {OUT}/ecm_native_metrics.json + 2 npz", flush=True)


if __name__ == "__main__":
    main()
