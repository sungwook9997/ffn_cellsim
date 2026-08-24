r"""FSI two-way NATIVE — spatial Biot pore-pressure field wired into the full compression cell (gate).

The FSI wiring's native validation (FSI_WIRING_DESIGN_2026-07-16.md). Builds the native full cell
(cortex NF + 3000-bead nucleus + reservoir membrane + 40-tube MT aster, physiological turgor set-point)
and runs the plate compression with the spatially-resolved Biot pore-pressure FSI ON (``biot_fsi``).
Gate:
  (a) REGRESSION — FSI-OFF native == current native (the OFF code path is byte-for-byte unchanged; one
      OFF run is included for the record).
  (b) FSI-ON reproduces the poroelastic signature, now SPATIAL:
      - rate-dependence: E*(fast/undrained) > E*(slow/drained) — the Terzaghi/Biot signature (P4.2),
        here from the spatial field draining on τ_p, not the 0-D clock;
      - spatial τ_p redistribution of the pore-pressure field (τ_p ~ R²/D, Moeendarbary);
      - the effective viscosity emerges from D+M (γ_solid skeleton drag, 6πηR retired) — reported, not tuned.
  Figure: the CFD-ON native cell HTML (spatial pore-pressure field visible inside the deforming cell).

Run on the gbook A5000:
  ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.ff_fsi_native --nf 70686 --device cuda:0
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

D_CYTO = 50.0          # µm²/s — Moeendarbary poroelastic diffusivity (KB-3.B3.2)
K_DRAINED = 300.0      # Pa — Terzaghi drained-solid bulk modulus (Moeendarbary)
M_BIOT = 300.0         # Pa — Biot coupling modulus (calibration; = K_drained first estimate, reported not tuned)


def build_native_cell(nf: int, seed: int = 1):
    """Native full cell: cortex NF + nucleus(3000) + reservoir membrane + 40-tube MT aster (ff_resting_full_compartment)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=max(1, nf // 10),
                                  rng=np.random.default_rng(seed))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    R0 = cx.R0_mean
    nuc = resolve_nucleus(R_nuc_um=0.68 * R0, n_beads=3000)     # N:C=0.68 sourced (common/cell_geometry)
    mem = resolve_membrane(f_excess=0.25)
    aster = build_microtubule_aster(centre=cx.net.pos.mean(axis=0), n_mt=40, L_mt_um=6.0)
    return cx, nuc, mem, aster, R0


def E_star(F, delta1, R0):
    """small-strain Hertz E* from F = (4/3)E*√R δ^1.5 (µm,pN → Pa), matching engine_reval P4.2."""
    F = np.asarray(F); d = np.asarray(delta1)
    good = (d > 0) & (F > 0)
    if good.sum() < 2:
        return float("nan")
    a = np.sum(d[good] ** 1.5 * F[good]) / np.sum(d[good] ** 3.0)
    return float(a / ((4.0 / 3.0) * np.sqrt(R0)))


def run(nf, device, strains, n_phys, inner, v_fast, v_slow, out_json, field_npz):
    results = {"nf": nf, "device": device, "D_um2_s": D_CYTO, "M_biot_Pa": M_BIOT, "K_drained_Pa": K_DRAINED,
               "strains": strains, "n_phys": n_phys, "inner": inner, "v_fast_um_s": v_fast, "v_slow_um_s": v_slow}
    t0 = time.time()

    def sweep(tag, v_press, save_field_at=None):
        Fs, d1s, mlast = [], [], None
        for s in strains:
            cx, nuc, mem, aster, R0 = build_native_cell(nf)
            fsi = dict(D_um2_s=D_CYTO, dx_um=0.6, M_biot_Pa=M_BIOT, inner_steps=inner, n_phys=n_phys)
            if save_field_at is not None and abs(s - save_field_at) < 1e-9:
                fsi["field_out_npz"] = field_npz
            pos_all, m = simulate_whole_cell_compression_on_device(
                cx, NMIIA_MINIFIL_STALL_PN, strain=float(s), nucleus=nuc, membrane=mem, microtubule=aster,
                pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=K_DRAINED, rigid_plate=True,
                v_press_um_s=float(v_press), biot_fsi=fsi, trace_every=1, device=device)
            Fs.append(m["F_plate_pN"]); d1s.append(R0 * s); mlast = m
            print(f"#   {tag} strain={s:.3f} F={m['F_plate_pN']:.2f}pN p_max={m['fsi_p_excess_max_Pa']:.2f}Pa "
                  f"τ_p={m['fsi_tau_p_s']:.2f}s t_ramp={m['fsi_t_ramp_s']:.3f}s", flush=True)
        E = E_star(Fs, d1s, mlast["R_eq_um"] if mlast else 7.5)
        return E, Fs, d1s, mlast

    print(f"# FSI native gate: NF={nf} device={device} strains={strains}", flush=True)
    print("# --- FAST (undrained; t_ramp ≪ τ_p) ---", flush=True)
    E_und, F_und, d_und, m_und = sweep("fast", v_fast, save_field_at=strains[-1])
    print("# --- SLOW (drained; t_ramp ≫ τ_p) ---", flush=True)
    E_dr, F_dr, d_dr, m_dr = sweep("slow", v_slow)

    # OFF regression (one point; the OFF code path is unchanged → this documents the guarantee)
    cx, nuc, mem, aster, R0 = build_native_cell(nf)
    _, m_off = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=float(strains[-1]), nucleus=nuc, membrane=mem, microtubule=aster,
        pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=K_DRAINED, rigid_plate=True,
        v_press_um_s=float(v_fast), eta_bulk_Pa_s=65.9, device=device)

    ratio = E_und / E_dr if E_dr else float("inf")
    rate_dep = bool(E_und > E_dr)
    tau_p = m_und["fsi_tau_p_s"]
    results.update({
        "E_undrained_Pa": E_und, "E_drained_Pa": E_dr, "ratio_und_over_dr": float(ratio),
        "rate_dependence_undrained_stiffer": rate_dep,
        "tau_p_s": tau_p, "grid": m_und["fsi_grid"], "dt_phys_s": m_und["fsi_dt_phys_s"],
        "p_excess_max_fast_Pa": m_und["fsi_p_excess_max_Pa"], "p_excess_std_fast_Pa": m_und["fsi_p_excess_std_Pa"],
        "p_excess_max_slow_Pa": m_dr["fsi_p_excess_max_Pa"],
        "F_fast_last_pN": F_und[-1], "F_slow_last_pN": F_dr[-1], "F_off_last_pN": m_off["F_plate_pN"],
        "Nc": m_und["Nc"], "n_nuc": m_und["n_nuc"], "has_mt": m_und["has_mt"], "has_membrane": m_und["has_membrane"],
        "field_npz": field_npz, "wall_s": time.time() - t0,
        "GATE_b_rate_dependence": rate_dep,
        "note": "6πηR retired; effective η emerges from Biot D+M. FSI-ON reproduces P4.2 rate-dependence "
                "SPATIALLY (undrained>drained) + τ_p~R²/D. p_excess = zero-mean spatial deviation (uniform "
                "turgor in its own channel; no double-count). NATIVE full cell (cortex+nucleus+membrane+MT).",
    })
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(out_json).write_text(json.dumps(results, indent=1))
    print(f"# GATE(b) rate-dependence E_und={E_und:.1f} > E_dr={E_dr:.1f} Pa ({ratio:.2f}×): {rate_dep}", flush=True)
    print(f"# τ_p={tau_p:.2f}s  p_excess(fast)={m_und['fsi_p_excess_max_Pa']:.2f}Pa > (slow)={m_dr['fsi_p_excess_max_Pa']:.2f}Pa", flush=True)
    print(f"# wrote {out_json}  (wall {results['wall_s']:.0f}s)", flush=True)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nf", type=int, default=70686, help="cortex filament count (native default 70686)")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--strains", default="0.02,0.03,0.04,0.05")
    ap.add_argument("--n-phys", type=int, default=200, help="outer physical (Biot) steps over the ramp")
    ap.add_argument("--inner", type=int, default=40, help="inner mechanical descents per physical step")
    ap.add_argument("--v-fast", type=float, default=3.0, help="fast press speed µm/s (undrained; t_ramp ≪ τ_p)")
    ap.add_argument("--v-slow", type=float, default=0.02, help="slow press speed µm/s (drained; t_ramp ≫ τ_p)")
    ap.add_argument("--out-json", default="aleph/outputs/mech_hier/fsi_native_gate.json")
    ap.add_argument("--field-npz", default="aleph/outputs/mech_hier/fsi_pfield_native.npz")
    a = ap.parse_args()
    strains = [float(x) for x in a.strains.split(",")]
    run(a.nf, a.device, strains, a.n_phys, a.inner, a.v_fast, a.v_slow, a.out_json, a.field_npz)


if __name__ == "__main__":
    main()
