r"""Internal-displacement DIAGNOSTIC — is the native cell's interior moving physically? (PI 2026-07-16)

The PI flagged that the cell's INTERNAL displacement "looks odd". This diagnoses it on the NATIVE full
cell (CLAUDE.md HARD: native + full + interactive 3D, never coarse/matplotlib). Build the physiological
cell (cortex NF + 3000-bead nucleus + reservoir membrane + 40-tube MT aster), relax to the pre-stressed
baseline (strain 0), then AFM-compress to `strain`, and measure the per-node displacement
Δ = pos(loaded) − pos(rest), split by compartment (cortex shell / MT aster+MTOC / nucleus). Renders a
CUT-AWAY interactive 3D cell (a hemisphere of cortex removed so the interior is visible) coloured by |Δ|,
so the interior kinematics are checkable by eye. Reports per-compartment displacement + non-physical
signatures (rigid COM drift, MTOC drift, nucleus centroid shift, interior vs shell magnitude).

Run on the gbook A5000:
  ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.ff_internal_displacement_diag \
     --nf 70686 --device cuda:0 --strain 0.10 --out ~/ff_scratch/internal_disp.html
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex, TURGOR_DP0
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.laws.microtubule import build_microtubule_aster
from aleph.laws.compartments import resolve_nucleus, resolve_membrane


def build_cell(nf, seed=1, mt_reach=False):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=max(1, nf // 10),
                                  rng=np.random.default_rng(seed))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    R0 = cx.R0_mean
    nuc = resolve_nucleus(R_nuc_um=0.68 * R0, n_beads=3000, nu_nuc=0.499)   # incompressible (clean FSI coupling)
    mem = resolve_membrane(f_excess=0.25)
    # mt_reach → MT tips REACH the cortex (compression-strut engagement fix); else the old 6µm (tips short → frozen)
    aster = build_microtubule_aster(centre=cx.net.pos.mean(axis=0), n_mt=40,
                                    L_mt_um=6.0, reach_R_um=(R0 if mt_reach else None))
    return cx, nuc, mem, aster, R0


def run_state(nf, strain, n_steps, device, seed=1, mt_reach=False, fsi=False, m_biot=300.0, if_cage=False,
              n_fil=400, mem_surface=False, mem_subdiv=4, mem_gap=0.15):
    """One relaxed state at `strain` (pre-stressed physiological turgor). Returns (pos_all, m).

    ``fsi=True`` runs the spatial-Biot two-way FSI path; ``if_cage=True`` builds the radial-spoke IF cage
    (nucleus↔cortex load path) so the diff rest→load isolates whether the nucleus now FOLLOWS the cortex.
    ``mem_surface=True`` wires the REAL plasma-membrane sheet (Helfrich bending + area tension + ERM tether),
    replacing the lumped 2γ/R tension — the defining bilayer bending mechanic the lumped model lacked."""
    cx, nuc, mem, aster, R0 = build_cell(nf, seed, mt_reach=mt_reach)
    _if = dict(n_fil=int(n_fil)) if if_cage else None
    _msurf = dict(subdiv=int(mem_subdiv), gap_um=float(mem_gap)) if mem_surface else None
    if fsi:
        biot = dict(D_um2_s=50.0, dx_um=0.6, M_biot_Pa=float(m_biot), inner_steps=30, n_phys=120)
        pos_all, m = simulate_whole_cell_compression_on_device(
            cx, NMIIA_MINIFIL_STALL_PN, strain=float(strain), nucleus=nuc, membrane=mem, microtubule=aster,
            pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=300.0, rigid_plate=True,
            v_press_um_s=1.0, biot_fsi=biot, if_cage=_if, membrane_surface=_msurf, device=device)
    else:
        pos_all, m = simulate_whole_cell_compression_on_device(
            cx, NMIIA_MINIFIL_STALL_PN, strain=float(strain), nucleus=nuc, membrane=mem, microtubule=aster,
            pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=300.0, n_steps=n_steps,
            turgor_every=50, rigid_plate=True, if_cage=_if, membrane_surface=_msurf, device=device)
    return pos_all, m, R0


def diagnose(nf, strain, n_steps, device, out_html, mt_reach=False, fsi=False, m_biot=300.0, if_cage=False,
             n_fil=400, mem_surface=False, mem_subdiv=4, mem_gap=0.15):
    t0 = time.time()
    print(f"# internal-displacement diag: NF={nf} strain={strain} device={device} mt_reach={mt_reach} fsi={fsi} if_cage={if_cage} mem_surface={mem_surface}", flush=True)
    _mk = dict(mt_reach=mt_reach, fsi=fsi, m_biot=m_biot, if_cage=if_cage, n_fil=n_fil,
               mem_surface=mem_surface, mem_subdiv=mem_subdiv, mem_gap=mem_gap)
    pos_rest, m0, R0 = run_state(nf, 0.0, n_steps, device, **_mk)
    print(f"#   rest built: Nc={m0['Nc']} Ne={m0['Ne']} n_nuc={m0['n_nuc']} R_eq={m0['R_eq_um']:.3f} ({time.time()-t0:.0f}s)", flush=True)
    if mem_surface:
        print(f"#   [membrane REST] n_mem={m0.get('n_mem',0)} R_mem={m0.get('R_mem_um',0):.3f}µm "
              f"bend={m0.get('mem_bend_over_8pikappa',0):.3f}×8πκ areal_strain={m0.get('mem_areal_strain',0):+.4f} "
              f"ERM_bleb={m0.get('erm_ruptured_bleb',0)}/{m0.get('n_mem',0)} f_rupt={m0.get('erm_f_rupt_pN',0):.2f}pN", flush=True)
    pos_load, m1, _ = run_state(nf, strain, n_steps, device, **_mk)
    Nc, Ne, n_nuc = int(m1["Nc"]), int(m1["Ne"]), int(m1["n_nuc"])
    print(f"#   loaded: R_eq={m1['R_eq_um']:.3f} F={m1['F_plate_pN']:.1f}pN ({time.time()-t0:.0f}s)", flush=True)

    # align on the plate axis: the compression is symmetric about the centroid; remove rigid COM drift so
    # the displacement is the true internal kinematics, not translation.
    c_rest = pos_rest[:Nc].mean(0); c_load = pos_load[:Nc].mean(0)
    com_drift = float(np.linalg.norm(c_load - c_rest))
    disp = (pos_load - c_load) - (pos_rest - c_rest)          # centroid-referenced per-node displacement
    dmag = np.linalg.norm(disp, axis=1)
    raw = pos_load - pos_rest                                 # RAW displacement (un-referenced): frozen ⇒ ≈ COM drift
    raw_mag = np.linalg.norm(raw, axis=1)

    def stat(sl, name):
        d = dmag[sl]; rw = raw_mag[sl]
        return {"compartment": name, "n": int(d.size), "mean_um": float(d.mean()) if d.size else 0.0,
                "p50_um": float(np.median(d)) if d.size else 0.0, "max_um": float(d.max()) if d.size else 0.0,
                "raw_mean_um": float(rw.mean()) if d.size else 0.0, "raw_std_um": float(rw.std()) if d.size else 0.0}
    cortex_sl = slice(0, Nc); mt_sl = slice(Nc, Ne); nuc_sl = slice(Ne, Ne + n_nuc)
    comps = [stat(cortex_sl, "cortex_shell"), stat(mt_sl, "MT_aster+MTOC"), stat(nuc_sl, "nucleus")]
    mtoc_idx = int(m1.get("mtoc_idx", -1))
    mtoc_disp = float(dmag[mtoc_idx]) if 0 <= mtoc_idx < dmag.size else float("nan")
    nuc_centroid_shift = float(np.linalg.norm((pos_load[nuc_sl] - c_load).mean(0) - (pos_rest[nuc_sl] - c_rest).mean(0))) if n_nuc else 0.0

    # FROZEN-INTERIOR test: does the interior DEFORM at all? (raw_std ≈ 0 ⇒ rigid; nucleus should FLATTEN under load)
    def shape(sl_pos):
        p = sl_pos - sl_pos.mean(0)
        return {"z_half_um": float(np.abs(p[:, 2]).max()), "xy_rad_um": float(np.hypot(p[:, 0], p[:, 1]).max())}
    nuc_shape = None
    if n_nuc:
        nr = shape(pos_rest[nuc_sl]); nl = shape(pos_load[nuc_sl])
        nuc_shape = {"rest": nr, "load": nl,
                     "z_flatten_frac": (nr["z_half_um"] - nl["z_half_um"]) / max(nr["z_half_um"], 1e-9),
                     "xy_bulge_frac": (nl["xy_rad_um"] - nr["xy_rad_um"]) / max(nr["xy_rad_um"], 1e-9)}
    # MT tip ↔ cortex minimum gap (are the tips even touching the cortex at this strain?)
    from scipy.spatial import cKDTree
    mt_tips_load = pos_load[mt_sl]
    tip_gap = float(cKDTree(pos_load[:Nc]).query(mt_tips_load)[0].min()) if Ne > Nc else float("nan")
    interior_frozen = bool((comps[1]["raw_std_um"] < 0.01) and (comps[2]["raw_std_um"] < 0.01))

    # axial (plate z) vs radial (xy) split — compression should push interior mostly RADIALLY outward at the
    # equator + AXIALLY inward at the poles. An interior that drifts/one-sidedly translates is the anomaly.
    zc = c_load[2]
    axial = np.abs(disp[:, 2]); radial = np.hypot(disp[:, 0], disp[:, 1])
    diag = {
        "nf": nf, "strain": strain, "R0_um": R0, "R_eq_rest_um": m0["R_eq_um"], "R_eq_load_um": m1["R_eq_um"],
        "Nc": Nc, "Ne": Ne, "n_nuc": n_nuc, "n_mt": int(m1.get("n_mt", 0)),
        "F_plate_pN": m1["F_plate_pN"], "com_drift_um": com_drift,
        "per_compartment": comps, "mtoc_disp_um": mtoc_disp, "nucleus_centroid_shift_um": nuc_centroid_shift,
        "interior_over_shell_mean_ratio": (comps[2]["mean_um"] / comps[0]["mean_um"]) if comps[0]["mean_um"] > 1e-9 else None,
        "axial_mean_um": float(axial.mean()), "radial_mean_um": float(radial.mean()),
        "nucleus_shape": nuc_shape, "mt_tip_cortex_gap_um": tip_gap, "interior_frozen": interior_frozen,
        "membrane": ({k: m1[k] for k in ("n_mem", "R_mem_um", "mem_area_um2", "mem_A0_um2", "mem_areal_strain",
                     "mem_bend_energy_pN_um", "mem_bend_over_8pikappa", "kappa_tilde_pN_um", "erm_bound",
                     "erm_ruptured_bleb", "erm_f_rupt_pN", "k_erm_pN_um") if k in m1}
                     if m1.get("membrane_surface_on") else None),
        "wall_s": time.time() - t0,
        "verdict": ("INTERIOR DECOUPLED — MT+nucleus rigid (raw_std≈0), not deforming under cortex load: "
                    "no cytoplasm stress-transmission / no IF cage / MT tip-contact not engaged"
                    if interior_frozen else "interior deforms with load"),
    }
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    (Path(out_html).with_suffix(".json")).write_text(json.dumps(diag, indent=1))
    _render(pos_load, disp, dmag, Nc, Ne, n_nuc, c_load, R0, strain, diag, out_html,
            mem_subdiv=mem_subdiv, n_mem=int(m1.get("n_mem", 0)) if mem_surface else 0)
    print("# --- diagnosis ---", flush=True)
    for c in comps:
        print(f"#   {c['compartment']:16s} n={c['n']:7d}  ref|Δ|mean={c['mean_um']*1e3:7.1f}nm  RAW mean={c['raw_mean_um']*1e3:7.1f}nm  RAW std={c['raw_std_um']*1e3:6.1f}nm", flush=True)
    print(f"#   COM drift={com_drift*1e3:.1f}nm  MTOC|Δ|={mtoc_disp*1e3:.1f}nm  MT-tip↔cortex gap={tip_gap:.3f}µm", flush=True)
    if nuc_shape:
        print(f"#   nucleus flatten: z {nuc_shape['rest']['z_half_um']:.3f}→{nuc_shape['load']['z_half_um']:.3f}µm "
              f"({nuc_shape['z_flatten_frac']*100:.1f}%)  xy bulge {nuc_shape['xy_bulge_frac']*100:.1f}%", flush=True)
    print(f"#   VERDICT: {diag['verdict']}  ({diag['wall_s']:.0f}s)", flush=True)
    return diag


def _render(pos, disp, dmag, Nc, Ne, n_nuc, c, R0, strain, diag, out_html, mem_subdiv=0, n_mem=0):
    """CUT-AWAY interactive 3D: remove the y>c_y hemisphere of the cortex so the interior (nucleus + MT) is
    visible; colour every node by |Δ| (turbo). The interior kinematics are then checkable by eye. The plasma
    membrane (if wired) is drawn as a TRANSLUCENT triangle sheet enclosing the cortex, coloured by |Δ| so a
    bleb (ERM rupture, high |Δ|) is visible on the bilayer."""
    import matplotlib.cm as cm
    from aleph.scripts.ff_viewer_html import build_viewer
    turbo = cm.get_cmap("turbo")
    hi = float(np.percentile(dmag, 99)) or 1e-9

    def cols(idx):
        return (turbo(np.clip(dmag[idx] / hi, 0, 1))[:, :3] * 255).astype(np.uint8)

    cortex = pos[:Nc]; keep = cortex[:, 1] <= c[1]              # cut the near hemisphere → interior visible
    ci = np.where(keep)[0]
    cortex_layer = {"name": f"cortex shell (cut-away hemisphere, {ci.size} pts) — |Δ| turbo", "kind": "points",
                    "verts": cortex[ci], "size": 1.4, "color_frames": [cols(ci)]}
    mt_idx = np.arange(Nc, Ne)
    mt_layer = {"name": f"MT aster + MTOC ({mt_idx.size} pts) — |Δ| turbo", "kind": "points",
                "verts": pos[mt_idx], "size": 2.6, "color_frames": [cols(mt_idx)]}
    layers = [cortex_layer, mt_layer]
    if n_nuc:
        nuc_idx = np.arange(Ne, Ne + n_nuc)
        layers.append({"name": f"nucleus ({n_nuc} beads) — |Δ| turbo", "kind": "points",
                       "verts": pos[nuc_idx], "size": 2.0, "color_frames": [cols(nuc_idx)]})
    if n_mem:                                                  # translucent plasma-membrane sheet (Helfrich bilayer)
        from aleph.laws.membrane_surface import build_membrane_mesh
        mem_faces = build_membrane_mesh(1.0, subdivisions=int(mem_subdiv)).faces   # topology only (radius irrelevant)
        mem_idx = np.arange(pos.shape[0] - n_mem, pos.shape[0])
        layers.append({"name": f"plasma membrane (Helfrich sheet, {n_mem} nodes) — translucent, |Δ| turbo",
                       "kind": "mesh", "verts": pos[mem_idx], "faces": mem_faces, "opacity": 0.28,
                       "color_frames": [cols(mem_idx)]})
    grad = ["#%02x%02x%02x" % tuple((np.asarray(turbo(x)[:3]) * 255).astype(int)) for x in np.linspace(0, 1, 16)]
    cbars = {"internal displacement |Δ| [µm]": {"lo": 0.0, "hi": hi, "label": "|Δ| node displacement [µm]", "grad": grad}}
    cm_ = diag["per_compartment"]
    title = (f"FF MCF7 cell — INTERNAL displacement under AFM strain={strain} (cut-away, |Δ| turbo) · "
             f"cortex {cm_[0]['mean_um']*1e3:.0f}nm / MT {cm_[1]['mean_um']*1e3:.0f}nm / nucleus {cm_[2]['mean_um']*1e3:.0f}nm mean · "
             f"COM drift {diag['com_drift_um']*1e3:.0f}nm · MTOC {diag['mtoc_disp_um']*1e3:.0f}nm · NATIVE full (Nc={Nc})")
    build_viewer({"internal displacement |Δ| [µm]": layers}, out_html, title=title, cbars=cbars)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nf", type=int, default=70686)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--strain", type=float, default=0.10)
    ap.add_argument("--n-steps", type=int, default=6000)
    ap.add_argument("--out", default="aleph/outputs/mech_hier/figs/internal_displacement.html")
    ap.add_argument("--mt-reach", action="store_true", help="build MT tips reaching the cortex (compression-strut fix)")
    ap.add_argument("--fsi", action="store_true", help="run the spatial-Biot two-way FSI (fluid couples to cortex+nucleus+MT)")
    ap.add_argument("--m-biot", type=float, default=300.0, help="Biot coupling modulus M_fsi [Pa]")
    ap.add_argument("--if-cage", action="store_true", help="build the radial-spoke IF cage (nucleus↔cortex load path)")
    ap.add_argument("--n-fil", type=int, default=400, help="IF cage spoke count (native ~7600, PI magic number)")
    ap.add_argument("--membrane-surface", action="store_true", help="wire the REAL plasma-membrane sheet (Helfrich bending + area tension + ERM), replacing the lumped 2γ/R")
    ap.add_argument("--mem-subdiv", type=int, default=4, help="membrane icosphere subdivisions (4→2562 nodes)")
    ap.add_argument("--mem-gap", type=float, default=0.15, help="membrane offset outside the cortex [µm] (ERM tether length)")
    a = ap.parse_args()
    diagnose(a.nf, a.strain, a.n_steps, a.device, a.out, mt_reach=a.mt_reach, fsi=a.fsi, m_biot=a.m_biot,
             if_cage=a.if_cage, n_fil=a.n_fil, mem_surface=a.membrane_surface, mem_subdiv=a.mem_subdiv, mem_gap=a.mem_gap)


if __name__ == "__main__":
    main()
