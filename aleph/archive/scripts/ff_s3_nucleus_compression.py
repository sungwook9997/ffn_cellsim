"""S3 — nucleus + cytoplasm: the layered compression response (native).

Professor's flow S3: add the NUCLEUS to the validated cortex and reveal its mechanics. FF is
cortex-native, so we honour the "isolate each compartment" intent by running the SAME plate
compression on cortex-only (S1 baseline) and on cortex+nucleus, and comparing — the difference IS
the nucleus contribution. The nucleus (R_nuc = N:C·R0, N:C = 0.68 sourced, E_nuc = 399 Pa MCF7, n_beads = 3000)
couples two ways (network_warp docstring): (1) hydrostatically at all strains — it displaces
V_nuc so the cytoplasm is smaller and ΔP rises faster; (2) by DIRECT plate contact once the
half-gap drops below R_nuc, i.e. strain > 1 − R_nuc/R_cell ≈ 0.30. RESULT (Finding 13): with
physiological MCF7 params the nucleus (E=399 Pa) is SOFTER than the bare cortex (~3624 Pa), so its
contribution is MODEST (+~11% at ε=0.45), not a stiff dominator; the contact onset is still detectable
as a force-ratio acceleration at ε≈0.30.

Gates (physically-correct observables): (G1) CONTACT ONSET — the cortex+nucleus/cortex-only
force-ratio ACCELERATES past ε ≈ 1−R_nuc/R0 (increment above ≥1.5× below); (G2) the nucleus is a NET
stiffener (ratio>1, monotonic). Magnitude is REPORTED, not thresholded. (An earlier equatorial-BULGE
gate was retired — structurally impossible for this volume-non-conserving shell nucleus; see Finding 13.)

Native: --nfil 70686 --device cuda:0. Each strain run is quasi-static (converged, 15000 steps).

Usage: python -m aleph.scripts.ff_s3_nucleus_compression --nfil 70686 --device cuda:0
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from aleph.laws.compartments import resolve_nucleus
from aleph.laws.cell_geometry import MCF7_GEOMETRY, sample_cell_geometry
from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.scripts.ff_s1_sphere import LOAD_PHYSIO
from aleph.scripts.ff_viewer_html import build_viewer

R_NUC_FRAC = MCF7_GEOMETRY.nc_ratio   # 0.68 = SOURCED MCF7 nucleus:cell RADIUS ratio (N:C), from the
#                            measured distribution (BioNumbers vol anchor + IFC n=2164 + photoacoustic;
#                            common/cell_geometry.py). Replaces the prior 0.70 whose "Moore 2016" cite is
#                            NOT in the KB (phantom). Central run uses this mean; --geom-seed samples the dist.
N_BEADS = 3000             # production nucleus beads (keeps the nucleus CFL below the cortex CFL)
E_NUC_PA = 399.0           # MCF7 in-situ (Fischer/Hayn/Mierke 2020), resolve_nucleus default
NU_NUC = 0.499             # nucleus Poisson ratio: near-incompressible stand-in for ν→½ (PI decision A,
#                            2026-07-15); volume-conserving flattening (literature). None → legacy compressible.


def _build(nfil, seed=1, R_um=None):
    params = CortexParams(R_um=R_um) if R_um is not None else CortexParams()
    cx = build_crosslinked_cortex(params, n_filaments=nfil, n_xl=nfil,
                                  n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx


def run_strain(nfil, strain, n_steps, device, with_nucleus, seed=1, geom_seed=None):
    """One plate-compression run at a fixed strain; returns (metrics, cortex_pos, nucleus_pos, R0, seg).

    ``geom_seed=None`` → CENTRAL geometry (R_cell = 7.5 µm volume anchor via CortexParams, R_nuc = N:C·R0,
    back-compat). An int → draw (R_cell, R_nuc) from the SOURCED MCF7 size distribution
    (``common/cell_geometry``) for one ensemble realization (quenched geometry disorder), building the
    cortex at the sampled R_cell and the nucleus at the co-varying sampled R_nuc."""
    if geom_seed is not None:
        R_cell_s, R_nuc_s = sample_cell_geometry(np.random.default_rng(geom_seed), MCF7_GEOMETRY)
        cx = _build(nfil, seed, R_um=R_cell_s)
    else:
        cx = _build(nfil, seed)
    R0 = cx.R0_mean
    R_nuc_um = R_nuc_s if geom_seed is not None else R_NUC_FRAC * R0
    nuc = (resolve_nucleus(R_nuc_um=R_nuc_um, n_beads=N_BEADS, E_nuc_Pa=E_NUC_PA, nu_nuc=NU_NUC)
           if with_nucleus else None)                # nu_nuc → incompressible (PI decision A)
    pos_all, m = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=float(strain), nucleus=nuc, n_steps=n_steps,
        turgor_every=(50 if device.startswith("cuda") else 20), device=device, **LOAD_PHYSIO)
    pos = np.asarray(pos_all)
    Nc, Ne = int(m["Nc"]), int(m["Ne"])
    pcx = pos[:Nc]
    pnuc = pos[Ne:] if (with_nucleus and int(m["n_nuc"]) > 0) else np.zeros((0, 3))
    seg = cx.net.segments                                    # cortex line topology (same across strains)
    return m, pcx, pnuc, R0, seg


def main() -> None:
    ap = argparse.ArgumentParser(description="S3 nucleus layered-compression response.")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--steps", type=int, default=15000)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--strains", default="0.1,0.2,0.3,0.4,0.5")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere/s3_nucleus_compression.json")
    ap.add_argument("--fig-out", default="aleph/outputs/mech_hier/figs/s3_nucleus_layered.png")
    ap.add_argument("--html-out", default="aleph/outputs/mech_hier/figs/s3_nucleus_compression_3d.html")
    ap.add_argument("--baseline", default=None, help="reuse cortex-only F from this prior sweep json "
                    "(cortex-only is nucleus-independent → skip re-running it; halves the sweep)")
    ap.add_argument("--geom-seed", type=int, default=None, help="draw (R_cell, R_nuc) from the SOURCED "
                    "MCF7 size distribution (common/cell_geometry) for one ensemble realization; "
                    "omit → central geometry (7.5 µm vol anchor, N:C=0.68). Vary across seeds for a spread.")
    a = ap.parse_args()

    strains = [float(x) for x in a.strains.split(",")]
    rows = []
    R0 = None
    if a.baseline:                                          # reuse the nucleus-independent cortex-only baseline
        base = json.load(open(a.baseline))
        for r in base["rows"]:
            if not r["with_nucleus"] and r["strain"] in strains:
                rows.append(r)
        print(f"# reusing {len(rows)} cortex-only baseline rows from {a.baseline} (skipping cortex-only runs)", flush=True)
    configs = (True,) if a.baseline else (False, True)
    cxf, nucf, seg = [], [], None                            # cortex-line + nucleus-point frames (cortex+nuc runs)
    for s in strains:
        for with_nuc in configs:                             # cortex-only baseline (unless reused), then cortex+nucleus
            m, pcx, pnuc, R0, seg_i = run_strain(a.nfil, s, a.steps, a.device, with_nuc, geom_seed=a.geom_seed)
            if with_nuc:                                     # capture geometry for the deformation viewer
                seg = seg_i if seg is None else seg
                c = pcx.mean(0)
                cxf.append(((pcx - c)[seg].reshape(-1, 3))[:, [0, 2, 1]])   # compression z → viewer-y (visible flatten)
                nucf.append((pnuc - c)[:, [0, 2, 1]])
            rows.append({"strain": s, "with_nucleus": with_nuc,
                         "F_plate_pN": float(m["F_plate_pN"]), "dP_turgor_Pa": float(m["dP_turgor_Pa"]),
                         "nucleus_contact": bool(m.get("nucleus_contact", False)),
                         "R_nuc_eq_um": float(m.get("R_nuc_eq_um", 0.0)),
                         "R_nuc_z_um": float(m.get("R_nuc_z_um", 0.0)),
                         "nuc_vol_conserved": float(m.get("nuc_vol_conserved", 1.0)),
                         "K_vol_nuc_Pa": float(m.get("K_vol_nuc_Pa", 0.0)),
                         "V_over_V0": float(m["V_over_V0"]), "n_nuc": int(m.get("n_nuc", 0))})
            tag = "cortex+nuc" if with_nuc else "cortex-only"
            print(f"strain={s:4.2f} {tag:11s}  F={m['F_plate_pN']:9.0f} pN  "
                  f"contact={m.get('nucleus_contact', False)}  R_nuc_eq={m.get('R_nuc_eq_um', 0.0):.2f}um  "
                  f"dP={m['dP_turgor_Pa']:.1f} Pa", flush=True)

    # --- gates ---
    # The nucleus's mechanical engagement is detected via the cortex+nucleus/cortex-only FORCE RATIO.
    # The signature of DIRECT plate contact (vs the always-on hydrostatic channel) is that the ratio
    # ACCELERATES once the plate passes the nucleus (ε > 1−R_nuc/R0): below contact only the hydrostatic
    # V_cyto reduction acts (weak), above it the plate directly loads the nucleus shell.
    #  G1 (contact onset at theory): the per-step ratio INCREMENT above the contact strain exceeds the
    #     increment below by ≥1.5× — a convex kink at the geometrically-predicted ε_contact.
    #  G2 (nucleus is a net stiffener): the ratio is >1 and rises monotonically with strain.
    # NOTE (2026-07-15, adversarial review wf_b7a1f531): an earlier G1 tested an equatorial BULGE
    # (R_nuc_eq > R_nuc). That is STRUCTURALLY IMPOSSIBLE for this nucleus model — the shell is
    # independent per-bead RADIAL springs with NO volume conservation, so plate contact FLATTENS the
    # poles (xy<R_nuc) but never pushes the equator past R_nuc (R_nuc_eq_max ≡ R_nuc). The correct
    # contact observable is the force-ratio acceleration (or bead z-flattening), NOT a bulge. See Finding 13.
    s_contact_theory = 1.0 - R_NUC_FRAC                       # ≈ 0.30
    R_nuc = R_NUC_FRAC * R0
    ratios, rnuc_eq = {}, {}
    for s in strains:
        Fc = next(r["F_plate_pN"] for r in rows if r["strain"] == s and not r["with_nucleus"])
        Fn = next(r["F_plate_pN"] for r in rows if r["strain"] == s and r["with_nucleus"])
        ratios[s] = Fn / max(abs(Fc), 1e-9)
        rnuc_eq[s] = next(r["R_nuc_eq_um"] for r in rows if r["strain"] == s and r["with_nucleus"])
    # per-interval ratio increments, split at the contact strain
    incs = [(0.5 * (strains[i] + strains[i + 1]), ratios[strains[i + 1]] - ratios[strains[i]])
            for i in range(len(strains) - 1)]
    inc_lo = [dr for mid, dr in incs if mid < s_contact_theory]
    inc_hi = [dr for mid, dr in incs if mid >= s_contact_theory]
    inc_lo_m = float(np.mean(inc_lo)) if inc_lo else float("nan")
    inc_hi_m = float(np.mean(inc_hi)) if inc_hi else float("nan")
    accel = (inc_hi_m / inc_lo_m) if (inc_lo and inc_hi and inc_lo_m > 0) else float("nan")
    ratio_max = ratios[strains[-1]]
    monotonic = all(ratios[strains[i + 1]] >= ratios[strains[i]] - 1e-6 for i in range(len(strains) - 1))
    g1_contact = bool(inc_lo and inc_hi) and accel >= 1.5       # contact-onset acceleration at ε_contact
    g2_active = ratio_max > 1.0 and monotonic                  # nucleus is a net (if modest) stiffener
    rnuc_eq_max = max(rnuc_eq.values(), default=0.0)
    print(f"# S3 G1 contact onset (ratio accel @ ε≈{s_contact_theory:.2f}): inc_lo={inc_lo_m:.4f} "
          f"inc_hi={inc_hi_m:.4f} accel={accel:.2f}x -> {'PASS' if g1_contact else 'FAIL'}", flush=True)
    print(f"# S3 G2 nucleus net stiffener: ratio_max={ratio_max:.3f} monotonic={monotonic} -> "
          f"{'PASS' if g2_active else 'FAIL'}", flush=True)
    _bulge = "bulges R_nuc_eq→%.2f" % rnuc_eq_max if rnuc_eq_max > R_nuc * 1.01 else \
             "R_nuc_eq=%.2f≈R_nuc (compressible, no bulge)" % rnuc_eq_max
    print(f"# S3 nucleus contribution MODEST: +{(ratio_max-1)*100:.0f}% at ε={strains[-1]:.2f} "
          f"(soft nucleus E={E_NUC_PA:.0f}<cortex → shear-limited; {_bulge}, R_nuc={R_nuc:.2f}). "
          f"See Finding 13.", flush=True)

    out = {"nfil": a.nfil, "R0_um": R0, "R_nuc_um": R_nuc, "E_nuc_Pa": E_NUC_PA,
           "n_beads": N_BEADS, "s_contact_theory": s_contact_theory,
           "force_ratio": {str(s): ratios[s] for s in strains},
           "R_nuc_eq": {str(s): rnuc_eq[s] for s in strains}, "R_nuc_eq_max": rnuc_eq_max,
           "ratio_increments": incs, "inc_lo_mean": inc_lo_m, "inc_hi_mean": inc_hi_m,
           "contact_accel_x": accel, "ratio_max": ratio_max, "nucleus_pct_at_max": (ratio_max - 1) * 100,
           "g1_contact_onset_pass": bool(g1_contact), "g2_nucleus_active_pass": bool(g2_active),
           "s3_pass": bool(g1_contact and g2_active), "rows": rows}
    with open(a.out, "w") as f:
        json.dump(out, f, indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fc = [next(r["F_plate_pN"] for r in rows if r["strain"] == s and not r["with_nucleus"]) for s in strains]
        fn = [next(r["F_plate_pN"] for r in rows if r["strain"] == s and r["with_nucleus"]) for s in strains]
        fig, ax = plt.subplots(figsize=(7.0, 4.8))
        ax.plot(strains, fc, "o-", color="#0072B2", lw=1.8, label="cortex only (S1)")
        ax.plot(strains, fn, "s-", color="#D55E00", lw=1.8, label="cortex + nucleus (S3)")
        ax.axvline(s_contact_theory, color="#999", ls="--", lw=1.0,
                   label=f"nucleus contact ε≈{s_contact_theory:.2f} (R_nuc/R={R_NUC_FRAC})")
        ax.set(xlabel="compression strain ε", ylabel="plate force F [pN]",
               title=f"S3 layered compression — nucleus stiffens past contact (NF={a.nfil})\n"
                     f"nucleus bulge R_nuc_eq_max {rnuc_eq_max:.2f}µm (R_nuc {R_nuc:.2f}) · "
                     f"force ratio lo {ratio_lo:.2f} → hi {ratio_hi:.2f}")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(a.fig_out, dpi=130)
        plt.close(fig)
        print(f"wrote {a.fig_out}", flush=True)
    except Exception as e:
        print(f"[fig skip] {e}", flush=True)

    # deformation viewer: cortex (lines) + nucleus (points), frames = strain sweep; the cut slider
    # reveals the nucleus inside being squeezed as the cell flattens (compression axis → viewer-y).
    if cxf and seg is not None:
        sc = f"cell compression + nucleus (ε {strains[0]:.2f}→{strains[-1]:.2f}, cut to reveal ▶)"
        scenes = {sc: [
            {"name": f"cortex (NF={a.nfil})", "kind": "lines", "verts": cxf[0], "frames": cxf,
             "color": "#8fbff0", "size": 1.1, "opacity": 0.35},
            {"name": f"nucleus ({N_BEADS} beads, R_nuc={R_NUC_FRAC*R0:.2f}µm)", "kind": "points",
             "verts": nucf[0], "frames": nucf, "color": "#D55E00", "size": 3.0, "opacity": 0.9},
        ]}
        try:
            build_viewer(scenes, out=a.html_out,
                         title=f"S3 nucleus layered compression (NF={a.nfil}) — cortex flattens, "
                               f"nucleus (soft E={E_NUC_PA:.0f}Pa) squeezed past ε≈{s_contact_theory:.2f}")
            print(f"wrote {a.html_out}", flush=True)
        except Exception as e:
            print(f"[html skip] {e}", flush=True)

    print(f"# S3 OVERALL: {'PASS' if (g1_contact and g2_active) else 'FAIL'}", flush=True)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
