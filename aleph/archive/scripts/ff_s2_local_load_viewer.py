"""S2 localized-load viewer + sign-reversal gate.

Applies a localized patch force (top pole) and ramps it INWARD (compression, 2A → a dimple)
and OUTWARD (tension, 2B → an outward bulge), capturing the deformed cortex + von-Mises stress
per force level. Two animated stress-coloured scenes + the sign-reversal gate (inward force →
inward deflection, outward → outward). Full-res, self-contained HTML (three.js).

Usage: python -m aleph.scripts.ff_s2_local_load_viewer --nfil 70686 --device cuda:0
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from scipy.spatial import ConvexHull

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.laws.ff_virial_stress import cortex_node_stress, cortex_node_areal_strain, face_areas
from aleph.scripts.ff_s1_sphere import LOAD_PHYSIO
from aleph.scripts.ff_resting_viewer import _turbo, _turbo_gradient
from aleph.scripts.ff_viewer_html import build_viewer


def load_ramp(nfil, fnodes, n_steps, device, cos_thresh=0.92, seed=1):
    """Ramp the patch force over `fnodes`.

    Returns (geom_frames, svm_frames, strain_frames, seg, R0, deflections). Both stress (σ_vm)
    and areal strain are captured per force level so the viewer can show the localized-load
    response as stress AND strain fields (per PI: FF-engine-style stress+strain viz).
    """
    geom, svm_frames, strain_frames, defl = [], [], [], []
    seg = faces = area0 = None
    last_pcx = rest_pcx = centre = None
    R0 = float(np.linalg.norm(
        (b := build_crosslinked_cortex(CortexParams(), n_filaments=nfil, n_xl=nfil,
                                       n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))).net.pos
        - b.net.pos.mean(0), axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    for f in fnodes:
        cx = build_crosslinked_cortex(CortexParams(), n_filaments=nfil, n_xl=nfil,
                                      n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))
        cx.R0_mean = R0
        pos_all, m = simulate_whole_cell_compression_on_device(
            cx, NMIIA_MINIFIL_STALL_PN, n_steps=n_steps,
            turgor_every=(50 if device.startswith("cuda") else 20), device=device,
            local_load={"axis": (0.0, 0.0, 1.0), "cos_thresh": cos_thresh, "f_node_pN": float(f)}, **LOAD_PHYSIO)
        Nc = cx.net.n_nodes
        pcx = np.asarray(pos_all)[:Nc]
        if seg is None:
            seg = cx.net.segments
            faces = ConvexHull(pcx).simplices.astype(np.int64)   # fixed triangulation from the near-rest (f=0) frame
            area0 = face_areas(pcx, faces)
        geom.append(pcx[seg].reshape(-1, 3))          # load axis = +z (top pole); diagonal camera shows the dimple/bulge
        svm, _ = cortex_node_stress(
            pcx, np.full(Nc, V0 / Nc),
            xl_ij=np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64), k_xl=cx.xl_k, r0_xl=cx.xl_rest,
            myo_ij=np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int64),
            f_myo=float(NMIIA_MINIFIL_STALL_PN), dP=float(m["dP_turgor_Pa"]))
        svm_frames.append(svm)
        strain_frames.append(cortex_node_areal_strain(pcx, faces, area0))
        defl.append(float(m["patch_deflection_um"]))
        if rest_pcx is None:
            rest_pcx, centre = pcx.copy(), pcx.mean(0)   # f≈0 frame = per-node rest reference
        last_pcx = pcx                                   # keep the max-|force| frame
        print(f"f_node={f:+8.0f} pN  deflection={m['patch_deflection_um']:+.4f} um  svm_p95={np.percentile(svm,95):.0f} Pa", flush=True)
    return geom, svm_frames, strain_frames, seg, R0, defl, last_pcx, rest_pcx, centre


def localization_profile(loaded_pcx, rest_pcx, centre, axis=(0.0, 0.0, 1.0), nbins=18):
    """Load-induced radial displacement vs polar angle θ from the load axis (patch-load localization).

    Differences the loaded frame against the **per-node f=0 rest frame** (not the mean radius R0),
    so the static shell roughness cancels and only the load response survives — critical because at
    any scale the resting shell radius varies node-to-node by O(the load deflection). Returns
    (theta_deg_centres, mean_radial_displacement_um, count_per_bin). A *localized* load confines the
    displacement to small θ (the load pole), decaying toward the equator/antipode; a global mode is
    ≈flat in θ. This is S2's second quantitative check beyond sign reversal.
    """
    rel_rest = np.asarray(rest_pcx) - np.asarray(centre)
    rad_rest = np.maximum(np.linalg.norm(rel_rest, axis=1), 1e-12)
    rhat = rel_rest / rad_rest[:, None]                            # per-node rest radial unit vector
    disp = np.asarray(loaded_pcx) - np.asarray(rest_pcx)
    disp = disp - disp.mean(0)                                    # remove rigid-body COM translation (drift)
    radial_disp = np.einsum("ij,ij->i", disp, rhat)               # signed radial displacement [µm]
    ax = np.asarray(axis, float); ax = ax / np.linalg.norm(ax)
    theta = np.degrees(np.arccos(np.clip(rhat @ ax, -1.0, 1.0)))
    edges = np.linspace(0.0, 180.0, nbins + 1)
    idx = np.clip(np.digitize(theta, edges) - 1, 0, nbins - 1)
    prof = np.array([radial_disp[idx == b].mean() if np.any(idx == b) else np.nan for b in range(nbins)])
    cnt = np.array([int(np.sum(idx == b)) for b in range(nbins)])
    return 0.5 * (edges[:-1] + edges[1:]), prof, cnt


def _stress_scene(geom, svm_frames, seg, hi, name):
    cf = [_turbo(s[seg].reshape(-1), 0.0, hi) for s in svm_frames]
    return [{"name": name, "kind": "lines", "verts": geom[0], "frames": geom, "color": "#8fbff0",
             "size": 1.5, "opacity": 0.7, "color_frames": cf, "cbar": "σ_vm [Pa]"}]


def _strain_scene(geom, strain_frames, seg, smag, name):
    cf = [_turbo(s[seg].reshape(-1), -smag, smag) for s in strain_frames]  # signed → diverging
    return [{"name": name, "kind": "lines", "verts": geom[0], "frames": geom, "color": "#8fbff0",
             "size": 1.5, "opacity": 0.7, "color_frames": cf, "cbar": "areal strain"}]


def main() -> None:
    ap = argparse.ArgumentParser(description="S2 localized-load sign-reversal viewer + gate.")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--fmax", type=float, default=4000.0, help="max per-node patch force [pN]")
    ap.add_argument("--nramp", type=int, default=6)
    ap.add_argument("--out", default="aleph/outputs/mech_hier/figs/s2_local_load_3d.html")
    ap.add_argument("--gate-out", default="aleph/outputs/mech_hier/s1_sphere/s2_sign_reversal.json")
    ap.add_argument("--fig-out", default="aleph/outputs/mech_hier/figs/s2_localization_decay.png")
    a = ap.parse_args()

    ramp = np.linspace(0.0, a.fmax, a.nramp)
    gin, svin, stin, seg, R0, din, pin, rin, cin = load_ramp(a.nfil, -ramp, a.steps, a.device)    # inward (2A)
    gout, svout, stout, _, _, dout, pout, rout, cout = load_ramp(a.nfil, ramp, a.steps, a.device)  # outward (2B)
    hi = float(max(np.percentile(np.concatenate(svin + svout), 98), 1e-6))
    # areal-strain colour scale: the fixed rest triangulation folds into a few DEGENERATE faces at the
    # dimple under large deformation, whose |strain| spikes to O(10) and would dominate a p98. Use a
    # robust p90 capped at a physical bound so the real strain gradient (~σ/E ~ 0.2) shows detail and the
    # degenerate outliers merely saturate. (A cleaner fix filters degenerate faces in cortex_node_areal_strain.)
    smag = float(np.clip(np.percentile(np.abs(np.concatenate(stin + stout)), 90), 1e-6, 1.0))
    s_in = "INWARD force → dimple: STRESS σ_vm (2A, ▶)"
    s_out = "OUTWARD force → bulge: STRESS σ_vm (2B, ▶)"
    s_in_e = "INWARD force → dimple: areal STRAIN (2A, ▶)"
    s_out_e = "OUTWARD force → bulge: areal STRAIN (2B, ▶)"
    scenes = {s_in: _stress_scene(gin, svin, seg, hi, f"inward load (NF={a.nfil})"),
              s_out: _stress_scene(gout, svout, seg, hi, f"outward load (NF={a.nfil})"),
              s_in_e: _strain_scene(gin, stin, seg, smag, f"inward strain (NF={a.nfil})"),
              s_out_e: _strain_scene(gout, stout, seg, smag, f"outward strain (NF={a.nfil})")}
    _sc = {"grad": _turbo_gradient()}
    cbars = {s_in: {"lo": 0.0, "hi": hi, "label": "cortex σ_vm", "unit": "Pa", **_sc},
             s_out: {"lo": 0.0, "hi": hi, "label": "cortex σ_vm", "unit": "Pa", **_sc},
             s_in_e: {"lo": -smag, "hi": smag, "label": "areal strain", "unit": "", **_sc},
             s_out_e: {"lo": -smag, "hi": smag, "label": "areal strain", "unit": "", **_sc}}
    build_viewer(scenes, out=a.out, cbars=cbars,
                 title=f"S2 localized load (NF={a.nfil}) — inward=dimple, outward=bulge; stress + strain (sign reversal)")

    # --- Gate 1: sign reversal (inward force → inward deflection, outward → outward) ---
    sign_pass = din[-1] < -1e-3 and dout[-1] > 1e-3

    # --- Gate 2: spatial localization (deflection confined to the load pole, decays with θ) ---
    th, prof_in, cnt_in = localization_profile(pin, rin, cin)         # inward, max |force| vs its rest frame
    _, prof_out, cnt_out = localization_profile(pout, rout, cout)
    pole = th <= 20.0                                                 # load-pole cap
    far = (th >= 70.0) & (th <= 110.0)                                # equatorial band
    pole_in = float(np.nanmean(prof_in[pole])); far_in = float(np.nanmean(prof_in[far]))
    loc_ratio = abs(pole_in) / max(abs(far_in), 1e-9)                 # localized ⇒ ≫ 1
    loc_pass = loc_ratio >= 3.0
    print(f"# S2 LOCALIZATION: pole|Δr|={abs(pole_in):.3f}um  equator|Δr|={abs(far_in):.3f}um  "
          f"ratio={loc_ratio:.1f} -> {'PASS' if loc_pass else 'FAIL'}", flush=True)

    gate = {"inward_deflection_um": din[-1], "outward_deflection_um": dout[-1],
            "sign_reversal_pass": bool(sign_pass), "f_max_pN": a.fmax, "nfil": a.nfil,
            "inward_ramp": din, "outward_ramp": dout,
            "localization": {"theta_deg": th.tolist(), "profile_inward_um": prof_in.tolist(),
                             "profile_outward_um": prof_out.tolist(),
                             "pole_defl_um": pole_in, "equator_defl_um": far_in,
                             "loc_ratio": loc_ratio, "localization_pass": bool(loc_pass)},
            "s2_pass": bool(sign_pass and loc_pass)}
    with open(a.gate_out, "w") as f:
        json.dump(gate, f, indent=2)

    # decay figure: signed radial deflection vs polar angle (localization is a 1-D profile → line plot)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax_ = plt.subplots(figsize=(6.4, 4.2))
        ax_.axhline(0, color="#999", lw=0.8)
        ax_.plot(th, prof_in, "o-", color="#0072B2", label=f"inward (f={a.fmax:.0f} pN/node)")
        ax_.plot(th, prof_out, "s-", color="#D55E00", label=f"outward (f={a.fmax:.0f} pN/node)")
        ax_.axvspan(0, 20, color="#0072B2", alpha=0.07, label="load pole (θ≤20°)")
        ax_.set_xlabel("polar angle θ from load axis [deg]")
        ax_.set_ylabel("signed radial deflection Δr [µm]")
        ax_.set_title(f"S2 localization — deflection decays from load pole (NF={a.nfil}, ratio={loc_ratio:.1f})")
        ax_.legend(fontsize=8); fig.tight_layout()
        fig.savefig(a.fig_out, dpi=130); plt.close(fig)
        print(f"wrote {a.fig_out}", flush=True)
    except Exception as e:                                            # figure is best-effort
        print(f"[fig skip] {e}", flush=True)

    print(f"# S2 SIGN-REVERSAL: inward={din[-1]:+.3f}um outward={dout[-1]:+.3f}um -> "
          f"{'PASS' if sign_pass else 'FAIL'}", flush=True)
    print(f"# S2 OVERALL: {'PASS' if (sign_pass and loc_pass) else 'FAIL'} "
          f"(sign={sign_pass}, localization={loc_pass})", flush=True)
    print(f"wrote {a.out} + {a.gate_out}", flush=True)


if __name__ == "__main__":
    main()
