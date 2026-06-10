"""H.7 active-myosin force-budget audit — split the active-γ floor into its causes.

Gate-A (REFUTE) + Gate-B (REFUTE, connected mesh) leave one frontier
(H7_GATE_B_RESULT_2026-06-09.md §5 PRIMARY): the active-myosin cortical tension is
~100-400× below the turgor-driven structural tension and the Hosseini MCF7 band
[0.18,0.40] mN/m. WHY does myosin contribute ~0 spanning tension even though
grip_walk genuinely walks (Gate-A: s_grip→0.5, real 16× contraction)?

This audit decomposes the floor into the THREE §5 candidates with ONE built-cell
measurement (no re-implementation of the binder — it reads the live updater state and
the snapshot the integrator actually sees, so γ_soft here == the gate's γ_soft):

  (1) ENGAGEMENT      — n_engaged / n_heads_total. The kinetic equilibrium bound
                        fraction is k_on/(k_on+k_off0) ≈ 0.99 (k_on=50, k_off0=0.35),
                        so any large shortfall is a GEOMETRIC/availability limit
                        (heads cannot reach actin), not kinetics.
  (2) PER-HEAD FORCE  — mean |T| per engaged attach bond, T = k_head_actin·(r−r0),
                        r0≈0 in grip_walk so T = k·r delivered straight from geometry
                        (myosin.py:1428 CH1). vs F_stall_per_head (the ceiling).
  (3) TRANSMISSION    — generated force Σ|T| over all engaged attach bonds vs the
                        force that MOP actually reads as spanning tension
                        (γ_soft·2πR). The ratio is the geometric/cancellation
                        efficiency: how much local contraction becomes hoop tension.

It also reports the COHERENT CEILING — γ if every engaged bond's in-tangent-plane
tension counted with no cancellation — and the BAND-CLOSURE budget: the engaged count
× per-head force needed to reach γ=0.27 mN/m, so the dominant lever is named with a
number, not a guess.

NO tuning: every constant is config-resolved (mcf7_baseline.yaml / phase1_h3.yaml).
The build is the Gate-B operating point (suspended/rounded, FA OFF, turgor ON,
connected mesh) — identical harness to h7_gate_b_probe._build_settled_cell.

Usage (smoke, fast):
    python -m ffn_sim.scripts.h7_active_force_budget --n-filaments 160 \
        --warmup 1500 --contract-steps 20000 --device cpu --allow-cpu-dev
Usage (full ×40, gbook GPU, contraction plateau):
    python -m ffn_sim.scripts.h7_active_force_budget --device gpu \
        --warmup 4000 --contract-steps 2000000
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension

_PN = 1.0e12   # N -> pN
_MNM = 1.0e3   # N/m -> mN/m
_UM = 1.0e6    # m -> µm
_HOSSEINI_BAND = (0.18e-3, 0.40e-3)  # N/m, MCF7 interphase IQR (contract §7)


def _build_settled_cell(*, n_filaments, n_nuc_beads, warmup, softstart, device, seed,
                        areal_density=None, cm_z_struct=None, cm_bundle_mult=None,
                        faithful=False, stepping_mode=None, xlink_k=None,
                        xl_n=None, xl_koff0=None, xl_xbeta=None):
    """Suspended/rounded MCF7 (FA OFF, turgor ON), connected mesh, grip_walk myosin.

    Mirrors h7_gate_b_probe._build_settled_cell exactly so this audit measures the
    SAME operating point the gate measures. ``areal_density`` overrides the myosin
    minifilament areal density [1/µm²] for the SENSITIVITY sweep ONLY — it is NOT a
    production config change; the production datum stays at the literature 0.6/µm²
    pending the PI datum decision (this is a mechanism-confirmation, not gate-chasing).

    ``cm_z_struct`` / ``cm_bundle_mult`` / ``faithful`` vary the cortex CONSTRUCTION
    load-path for the architecture-as-lever WALL-A sensitivity test (h7 ARCH thrust):
    they change topology only, NOT a production datum, and do NOT touch any band/gate.
    Defaults (None / False) → the production build (uniform mesh, z=3.7, bundle=2)."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest["optional_subsystems"]["fa"]["enabled"] = False
    if n_filaments is not None:
        co = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
        co["n_filaments"] = int(n_filaments)
        co["demo_mode"] = True
    if areal_density is not None:
        myo = (manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
               .setdefault("myosin", {}))
        myo["areal_density_per_um2"] = float(areal_density)
    if stepping_mode is not None:
        myo = (manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
               .setdefault("myosin", {}))
        myo["stepping_mode"] = str(stepping_mode)
    if xlink_k is not None:
        # A/B SENSITIVITY ONLY (NOT a production config change): override the dynamic
        # crosslinker harmonic-bond stiffness k_intra/k_attach to isolate the
        # transmission lever's contribution to γ. Production stays at the loop18
        # re-anchored 1.0e-3 N/m (KB-1.28); pass 1.0e-7 to reproduce the PRE-re-anchor
        # transmission-floored state for a clean A/B (PI_DECISION_BRIEF_2026-06-10 §1).
        xl = (manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
              .setdefault("dynamic_crosslinkers", {}))
        xl["k_intra"] = float(xlink_k)
        xl["k_attach"] = float(xlink_k)
    if xl_n is not None or xl_koff0 is not None or xl_xbeta is not None:
        # ab_xlink_v2 transmission redesign (Slater 2021 Soft Matter 17:10274): the
        # actual stress-transmission knobs are crosslink DENSITY + force-dependent
        # unbinding (Bell-Evans k_off0, x_beta), NOT stiffness (v1 tested stiffness →
        # γ ~5%). A/B SENSITIVITY ONLY (not a production change): production stays at
        # the phase1_h3 anchors (n_xl, alpha_k_off0=0.066 Ferrer2008, alpha_x_beta=0.4nm).
        # Sweep values must be KU-anchored ±multiples (no magic numbers); see
        # docs/AB_XLINK_TRANSMISSION_REDESIGN_2026-06-11.md.
        xl = (manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
              .setdefault("dynamic_crosslinkers", {}))
        if xl_n is not None:
            xl["n_xl"] = int(xl_n)
        if xl_koff0 is not None:
            xl["alpha_k_off0"] = float(xl_koff0)
        if xl_xbeta is not None:
            xl["alpha_x_beta"] = float(xl_xbeta)
    if n_nuc_beads is not None:
        manifest["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)
    cm_kw = {}
    if cm_z_struct is not None:
        cm_kw["cm_z_struct"] = float(cm_z_struct)
    if cm_bundle_mult is not None:
        cm_kw["cm_bundle_mult"] = int(cm_bundle_mult)
    # continuous_stroke uses the stiff cross-bridge k (=1e-3) → k_backbone≈1e-2
    # → τ_backbone≈9 ns < the legacy cortex dt 13 ns. reconcile_dt folds the myosin
    # k_backbone into the global CFL and lowers the integrator dt accordingly (NO
    # edit to the frozen integrator/ — only the dt scalar). Additive: for the soft
    # legacy modes τ_backbone≫dt so dt is unchanged (bit-identical).
    reconcile_dt = (stepping_mode == "continuous_stroke")
    cell = build_baseline_cell(
        manifest=manifest, device=device, seed=seed,
        constrained=False, with_baoab=True,
        equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=softstart,
        connected_mesh=True, faithful_connected_mesh=bool(faithful),
        reconcile_dt=reconcile_dt,
        **cm_kw,
    )
    return cell


def _engaged_attach_bonds(snap, n_cortex_actin):
    """Return (head_tags, actin_tags, bin_names) for every live myosin attach bond.

    Attach bonds are the dynamic head→cortex-actin bonds (types
    'cortex_myosin_attach_b{i}'); they exist only when a head is engaged."""
    types = list(snap.bonds.types)
    tids = np.asarray(snap.bonds.typeid)
    group = np.asarray(snap.bonds.group)
    is_attach = np.array(
        [types[i].startswith("cortex_myosin_attach_b") for i in range(len(types))]
    )
    if len(tids) == 0:
        return np.zeros((0, 2), dtype=np.int64), np.zeros(0, dtype=np.int64)
    mask = is_attach[tids]
    g = group[mask].astype(np.int64)
    # Orient each bond as (head, actin): the cortex-actin end is the bead < n_cortex_actin.
    head = np.where(g[:, 0] < n_cortex_actin, g[:, 1], g[:, 0])
    actin = np.where(g[:, 0] < n_cortex_actin, g[:, 0], g[:, 1])
    return np.stack([head, actin], axis=1), tids[mask]


def audit(*, cell, n_contract_steps, sample_every):
    """Run the contraction phase, sampling the force-budget decomposition."""
    sim = cell.simulation
    R = float(cell.p_cortex.R_cell)
    k_ha = float(cell.p_myosin.k_head_actin)         # SCALED (mesoscale) if enabled
    F_stall = float(cell.p_myosin.F_stall_per_head)  # SCALED
    n_heads_total = int(2 * cell.p_myosin.n_heads_per_side * cell.p_myosin.n_motors_per_cell)
    n_cortex_actin = int(cell.n_cortex_actin)
    meso = getattr(cell.p_myosin, "extras", None) or {}
    factor = float(meso.get("mesoscale_force_factor", 1.0))
    # §9 continuous_stroke: the delivered per-head force is NOT the attach bond's
    # k·r (the bond is k=0); it is F=min(k·s_grip, F_stall) read off the updater's
    # s_grip. Use that as the per-head tension so the force budget reports the TRUE
    # cross-bridge load, not the geometric k·r the harmonic-bond modes carried.
    stepping_mode = str(getattr(cell.p_myosin, "stepping_mode", "binned_r0"))
    _continuous = stepping_mode == "continuous_stroke"
    if _continuous:
        from ffn_sim.cortex.myosin import continuous_stroke_force

    def _sample(tick):
        snap = sim.state.get_snapshot()
        if snap.communicator.rank != 0:
            return None
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        bonds, _ = _engaged_attach_bonds(snap, n_cortex_actin)
        n_eng = int(bonds.shape[0])
        # Per-bond geometry + delivered tension T = k·(r − r0), r0≈0 (grip_walk).
        # A = 4πR² is the cortex surface area (the correct 2D-stress normalisation).
        A_surf = 4.0 * np.pi * R * R
        if n_eng:
            d = pos[bonds[:, 0]] - pos[bonds[:, 1]]      # head − actin
            r = np.linalg.norm(d, axis=1)
            if _continuous:
                # Delivered force = min(k·s_grip, F_stall) (the MyosinHeadForce
                # law) PER BOND, aligned with r/uhat so g_ik stays elementwise-
                # correct. Map each bond's head tag → head local → s_grip.
                act = cell.myosin_action
                H = int(cell.p_myosin.n_heads_per_side)
                N = int(cell.p_myosin.n_backbone)
                per_motor = int(cell.p_myosin.n_particles_per_motor)
                mt0 = int(cell.myosin_layout.motor_tag_start)
                off = bonds[:, 0].astype(np.int64) - mt0          # head tags
                head_locals = (off // per_motor) * (2 * H) + ((off % per_motor) - N)
                T = continuous_stroke_force(
                    act._head_grip_s[head_locals], k_ha, F_stall)
            else:
                T = k_ha * r                              # delivered force [N]
            uhat = d / r[:, None].clip(min=1e-30)
            # radial direction at the actin-bead end (outward normal):
            normals = pos[bonds[:, 1]] / np.linalg.norm(
                pos[bonds[:, 1]], axis=1, keepdims=True).clip(min=1e-30)
            radial_comp = np.sum(uhat * normals, axis=1)           # û·n̂
            inplane2 = np.clip(1.0 - radial_comp**2, 0.0, 1.0)     # 1−(û·n̂)²
            gen_force = float(np.sum(np.abs(T)))                    # Σ|T| generated [N]
            mean_T = float(np.mean(np.abs(T)))
            mean_r = float(np.mean(r))
            # Irving-Kirkwood 2D surface-stress estimate over the engaged attach
            # bonds (same form cortical_tension.py uses for its cross-check:
            # γ_IK = Σ T·L·(1−(r̂·û)²)/(8πR²) = Σ T·L·inplane² / (2·A)). This is the
            # PHYSICAL isotropic-network surface tension these bonds carry — the
            # legitimate cross-check on the MOP γ_soft (no hidden "ceiling" above it).
            g_ik = float(np.sum(np.abs(T) * r * inplane2)) / (2.0 * A_surf)
        else:
            gen_force = mean_T = mean_r = g_ik = 0.0
        # The gate's own γ measurement (same entry point Gate-B uses).
        ct = measure_cortical_tension(
            sim, R_cell=R, p_enclosed_volume=cell.p_enclosed_volume,
        )
        g_soft = float(ct["gamma_soft"])               # N/m
        transmitted_force = g_soft * 2.0 * np.pi * R   # Σ F_cut implied by γ_soft [N]
        # Bond-type-resolved soft γ: split the active channel into the MYOSIN bonds
        # (the direct dipole) vs the cortex-ACTIN backbone bonds (the network the
        # contraction must load to amplify beyond the direct dipole). If myosin
        # tension reaches the actin network, g_soft_actin rises with myoON; if it is
        # shunted locally (Gate-A r/r0=1), g_soft_actin is myosin-insensitive. This
        # is the empirical generation-vs-propagation test (vs assuming ℓ).
        all_types = set(sim.state.bond_types)
        myo_types = {t for t in all_types if t.startswith("cortex_myosin")}
        actin_types = {"cortex-bond"} & all_types
        g_soft_myo = float(measure_cortical_tension(
            sim, R_cell=R, cortical_bond_types=myo_types)["gamma_soft"]) if myo_types else 0.0
        g_soft_actin = float(measure_cortical_tension(
            sim, R_cell=R, cortical_bond_types=actin_types)["gamma_soft"]) if actin_types else 0.0
        return {
            "tick": int(tick),
            "n_engaged": n_eng,
            "n_heads_total": n_heads_total,
            "engaged_frac": n_eng / n_heads_total if n_heads_total else 0.0,
            "s_grip_over_l0": float(
                cell.myosin_action._head_grip_s.mean() / float(cell.p_cortex.rest_length)
            ) if cell.myosin_action is not None else 0.0,
            "mean_T_pN": mean_T * _PN,
            "mean_attach_r_nm": mean_r * 1e9,
            "F_stall_pN": F_stall * _PN,
            "gen_force_nN": gen_force * 1e9,
            "transmitted_force_nN": transmitted_force * 1e9,
            "g_soft_mN_m": g_soft * _MNM,
            "g_soft_myosin_bonds_mN_m": g_soft_myo * _MNM,
            "g_soft_actin_bonds_mN_m": g_soft_actin * _MNM,
            "g_ik_estimate_mN_m": g_ik * _MNM,
            "g_rigid_mN_m": float(ct["gamma_rigid"]) * _MNM,
            "g_passive_mN_m": float(ct["gamma_passive"]) * _MNM,
        }

    rows = []
    # initial sample (loading phase, s_grip≈0)
    s0 = _sample(0)
    if s0:
        rows.append(s0)
        print(f"  [tick 0] n_eng={s0['n_engaged']}/{s0['n_heads_total']} "
              f"({s0['engaged_frac']*100:.2f}%)  g_soft={s0['g_soft_mN_m']:.3e} mN/m",
              flush=True)
    done = 0
    while done < n_contract_steps:
        chunk = min(sample_every, n_contract_steps - done)
        sim.run(chunk)
        done += chunk
        s = _sample(done)
        if s:
            rows.append(s)
            print(f"  [tick {done}] n_eng={s['n_engaged']} "
                  f"({s['engaged_frac']*100:.2f}%)  s_grip/l0={s['s_grip_over_l0']:.3f}  "
                  f"g_soft={s['g_soft_mN_m']:.3e}  g_ik={s['g_ik_estimate_mN_m']:.3e}  "
                  f"meanT={s['mean_T_pN']:.2f}pN  gen={s['gen_force_nN']:.3f}nN", flush=True)

    # Summary: take the last quartile as the plateau.
    tail = rows[max(1, 3 * len(rows) // 4):] or rows[-1:]

    def _avg(key):
        return float(np.mean([r[key] for r in tail]))

    g_soft_plateau = _avg("g_soft_mN_m")
    band_lo, band_hi = _HOSSEINI_BAND[0] * _MNM, _HOSSEINI_BAND[1] * _MNM

    # --- Analytical parameter-implied active tension (the "what do these params
    # predict" number, independent of the sim measurement). Anchored to the
    # project-canonical active-gel relation γ_active = σ_active·h (KB-3.5), with the
    # microscopic active stress of an isotropic force-dipole population
    #   σ_active = n_3D · P,   P = f_dipole · ℓ   (dipole moment; NO 1/2 — the virial
    #   sum over distinct dipoles counts each once),  n_3D = n_2D / h.
    # Since γ = σ·h, the cortex thickness h CANCELS:
    #   γ_active ≈ n_2D · f_dipole · ℓ_minifil
    # f_dipole = the ONE-SIDED head force = n_heads_per_side·F_stall (a bipolar
    # minifilament is a force dipole; the two 28-head sets are the +/− of the SAME
    # dipole, not additive). ℓ = backbone (dipole arm). FULL engage+stall = upper
    # envelope. CONVENTION NOTE: some active-gel papers carry an extra 1/2 (dipole
    # defined as f·ℓ/2); that would HALVE this envelope and DOUBLE the density needed
    # — i.e. the gap is 18× (this form) to 36× (½ form). The σ·h anchor (σ=1.35 kPa
    # ↔ γ=0.27 mN/m, h=200 nm) is satisfied by THIS form at n_2D≈16/µm², so it is the
    # primary; the ½ ambiguity is reported as a factor-2 band, not hidden.
    n2d = float((meso.get("native_n_motors", cell.p_myosin.n_motors_per_cell))
                / (4.0 * np.pi * float(cell.p_cortex.R_cell) ** 2))  # 1/m²
    f_stall_native = F_stall / max(factor, 1.0)                       # un-scaled per-head
    f_minifil = cell.p_myosin.n_heads_per_side * f_stall_native       # one-sided dipole force
    ell_minifil = float(cell.p_myosin.backbone_length)
    g_analytic = n2d * f_minifil * ell_minifil                        # N/m (σ·h-anchored, no ½)

    # Band-closure budget at the MEASURED per-head T (independent of the analytic).
    mean_T_N = _avg("mean_T_pN") * 1e-12  # pN→N
    # γ_soft scales ~linearly with engaged-head count at fixed geometry, so the
    # head count that would lift the measured plateau into band_lo:
    need_engaged_for_band = (
        _avg("n_engaged") * band_lo / g_soft_plateau if g_soft_plateau > 0 else None)

    summary = {
        "operating_point": "suspended/rounded, FA OFF, turgor ON, connected mesh, grip_walk",
        "scale": {
            "n_filaments": int(cell.p_cortex.n_filaments),
            "n_cortex_actin": int(cell.n_cortex_actin),
            "n_motors": int(cell.p_myosin.n_motors_per_cell),
            "n_heads_per_side": int(cell.p_myosin.n_heads_per_side),
            "n_heads_total": int(2 * cell.p_myosin.n_heads_per_side
                                 * cell.p_myosin.n_motors_per_cell),
            "mesoscale_force_factor": factor,
            "native_n_motors": float(meso.get("native_n_motors", cell.p_myosin.n_motors_per_cell)),
            "native_areal_density_per_um2": n2d / 1e12,
        },
        "band_mN_m": [band_lo, band_hi],
        "phase": ("loading (s_grip≈0)" if _avg("s_grip_over_l0") < 0.02
                  else "contraction (s_grip developed)"),
        "plateau": {
            "s_grip_over_l0": _avg("s_grip_over_l0"),
            "g_soft_mN_m": g_soft_plateau,
            "g_soft_myosin_bonds_mN_m": _avg("g_soft_myosin_bonds_mN_m"),
            "g_soft_actin_bonds_mN_m": _avg("g_soft_actin_bonds_mN_m"),
            "g_ik_estimate_mN_m": _avg("g_ik_estimate_mN_m"),
            "g_rigid_mN_m": _avg("g_rigid_mN_m"),
            "g_passive_mN_m": _avg("g_passive_mN_m"),
            "engaged_frac": _avg("engaged_frac"),
            "n_engaged": _avg("n_engaged"),
            "mean_T_pN": _avg("mean_T_pN"),
            "F_stall_pN": rows[-1]["F_stall_pN"],
            "mean_attach_r_nm": _avg("mean_attach_r_nm"),
            "gen_force_nN": _avg("gen_force_nN"),
            "transmitted_force_nN": _avg("transmitted_force_nN"),
        },
        "analytic": {
            "g_active_analytic_mN_m": g_analytic * _MNM,
            "native_areal_density_per_um2": n2d / 1e12,
            "f_minifilament_pN": f_minifil * _PN,
            "ell_minifilament_nm": ell_minifil * 1e9,
            "gap_factor_analytic_to_band_lo": band_lo / (g_analytic * _MNM)
            if g_analytic > 0 else None,
            "note": "γ≈½·n2D·f_minifil·ℓ at FULL engagement+stall — param upper envelope",
        },
        "band_closure_budget": {
            "need_engaged_heads_for_band_lo": need_engaged_for_band,
            "have_engaged_heads": _avg("n_engaged"),
            "have_heads_total": int(2 * cell.p_myosin.n_heads_per_side
                                    * cell.p_myosin.n_motors_per_cell),
            "gap_factor_g_soft_to_band_lo": band_lo / g_soft_plateau if g_soft_plateau > 0 else None,
        },
        "rows": rows,
    }
    return summary


def _verdict(s):
    p = s["plateau"]; b = s["band_closure_budget"]; a = s["analytic"]
    print("=" * 76, flush=True)
    print("H.7 ACTIVE-MYOSIN FORCE-BUDGET AUDIT — plateau decomposition", flush=True)
    print("-" * 76, flush=True)
    print(f"  scale: n_filaments={s['scale']['n_filaments']}, "
          f"{s['scale']['n_cortex_actin']} cortex beads, "
          f"{s['scale']['n_heads_total']} myosin heads, "
          f"meso×{s['scale']['mesoscale_force_factor']:.2f} "
          f"(native {s['scale']['native_n_motors']:.0f} motors @ "
          f"{s['scale']['native_areal_density_per_um2']:.2f}/µm²)", flush=True)
    print(f"  band (Hosseini IQR) = [{s['band_mN_m'][0]:.3f}, {s['band_mN_m'][1]:.3f}] mN/m",
          flush=True)
    print("-" * 76, flush=True)
    print(f"  (1) ENGAGEMENT   : {p['n_engaged']:.0f}/{s['scale']['n_heads_total']} "
          f"= {p['engaged_frac']*100:.2f}%  (kinetic equil ≈99% → shortfall is geometric)",
          flush=True)
    print(f"  (2) PER-HEAD T   : {p['mean_T_pN']:.3f} pN  (F_stall={p['F_stall_pN']:.2f} pN, "
          f"mean attach r={p['mean_attach_r_nm']:.0f} nm)", flush=True)
    print(f"  (3) GENERATED Σ|T| = {p['gen_force_nN']:.3f} nN", flush=True)
    print("-" * 76, flush=True)
    print(f"  γ_soft (active, MOP) = {p['g_soft_mN_m']:.4e} mN/m   "
          f"[{b['gap_factor_g_soft_to_band_lo']:.0f}× under band_lo]", flush=True)
    print(f"    ├ myosin bonds (direct dipole) = {p['g_soft_myosin_bonds_mN_m']:.4e} mN/m",
          flush=True)
    print(f"    └ actin  bonds (network propag) = {p['g_soft_actin_bonds_mN_m']:.4e} mN/m",
          flush=True)
    print(f"  γ_IK  (active, virial cross-check) = {p['g_ik_estimate_mN_m']:.4e} mN/m", flush=True)
    print(f"  γ_rigid (turgor) = {p['g_rigid_mN_m']:.4e} mN/m   "
          f"γ_passive(YL) = {p['g_passive_mN_m']:.4e} mN/m", flush=True)
    print("-" * 76, flush=True)
    print(f"  ANALYTIC param-implied γ_active (½·n2D·f_minifil·ℓ, FULL engage+stall):",
          flush=True)
    print(f"    = {a['g_active_analytic_mN_m']:.4e} mN/m   "
          f"[{(a['gap_factor_analytic_to_band_lo'] or 0):.1f}× under band_lo]", flush=True)
    print(f"    (n2D={a['native_areal_density_per_um2']:.2f}/µm², "
          f"f_minifil={a['f_minifilament_pN']:.0f} pN, ℓ={a['ell_minifilament_nm']:.0f} nm)",
          flush=True)
    print("-" * 76, flush=True)
    print("  BAND-CLOSURE BUDGET (γ_soft ~∝ engaged count at fixed geometry):", flush=True)
    if b["need_engaged_heads_for_band_lo"] is not None:
        print(f"    need engaged heads ≈ {b['need_engaged_heads_for_band_lo']:.0f} "
              f"(have {b['have_engaged_heads']:.0f} engaged, "
              f"{b['have_heads_total']} total)", flush=True)
    print("=" * 76, flush=True)
    print("  READING (two-wall decomposition):", flush=True)
    an_gap = a["gap_factor_analytic_to_band_lo"] or 0
    g_myo = p["g_soft_myosin_bonds_mN_m"]
    g_actin = p["g_soft_actin_bonds_mN_m"]
    propag = (g_actin / g_myo) if g_myo > 0 else 0.0
    print(f"   WALL A (propagation): actin-network γ / myosin-dipole γ = {propag*100:.1f}%", flush=True)
    if propag < 0.5:
        print(f"     → myosin tension does NOT load the actin network (no prestress", flush=True)
        print(f"       amplification). Band needs the long actin load-path (ℓ_path≫ℓ_minifil);", flush=True)
        print(f"       it is absent → TRANSMISSION wall (Gate-A verdict, located).", flush=True)
    else:
        print(f"     → actin network carries myosin tension; propagation present.", flush=True)
    print(f"   WALL B (direct-dipole generation): full-engage+stall envelope "
          f"{an_gap:.1f}× under band_lo", flush=True)
    print(f"     → even the direct dipole (ℓ=minifilament) is sub-band; band needs WALL-A", flush=True)
    print(f"       amplification AND/OR a higher force budget (density/stall datum, PI-gated).", flush=True)
    print("   NOTE: if s_grip≈0 this is the LOADING phase — confirm with a contraction run", flush=True)
    print("     (does actin-network γ rise as s_grip→0.5?). Do NOT tune params (band LOCKED).",
          flush=True)
    print("=" * 76, flush=True)


def _figure(s, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = s["rows"]
    ticks = [r["tick"] for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    band_lo, band_hi = s["band_mN_m"]

    ax = axes[0, 0]
    ax.plot(ticks, [r["g_soft_mN_m"] for r in rows], "o-", color="#c0392b",
            label="γ_soft (active, MOP)")
    ax.plot(ticks, [r["g_ik_estimate_mN_m"] for r in rows], "s--", color="#e67e22",
            label="γ_IK (virial cross-check)")
    ax.axhline(s["analytic"]["g_active_analytic_mN_m"], color="#8e44ad", ls=":", lw=1.6,
               label="analytic γ (full engage+stall)")
    ax.axhspan(band_lo, band_hi, color="green", alpha=0.15, label="Hosseini band")
    ax.set_yscale("log"); ax.set_xlabel("contraction step"); ax.set_ylabel("γ (mN/m)")
    ax.set_title("active γ: MOP vs virial vs analytic vs band")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(ticks, [r["engaged_frac"] * 100 for r in rows], "o-", color="#2e86c1")
    ax.set_xlabel("contraction step"); ax.set_ylabel("engaged heads (%)")
    ax.axhline(99, color="grey", ls="--", lw=1, label="kinetic equil ≈99%")
    ax.set_title("(1) engagement fraction"); ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(ticks, [r["mean_T_pN"] for r in rows], "o-", color="#8e44ad")
    ax.axhline(s["plateau"]["F_stall_pN"], color="grey", ls="--", lw=1,
               label=f"F_stall={s['plateau']['F_stall_pN']:.1f} pN")
    ax.set_xlabel("contraction step"); ax.set_ylabel("mean |T| per engaged head (pN)")
    ax.set_title("(2) per-head delivered force"); ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(ticks, [r["gen_force_nN"] for r in rows], "o-", color="#16a085",
            label="generated Σ|T| (engaged bonds)")
    ax.set_xlabel("contraction step"); ax.set_ylabel("generated force (nN)")
    ax.set_title("(3) total generated active force"); ax.legend(fontsize=8)

    fig.suptitle("H.7 active-myosin force-budget audit "
                 f"({'smoke' if s['scale']['n_filaments'] < 900 else 'full ×40'})",
                 fontweight="bold")
    fig.savefig(out_png, dpi=130)
    print(f"  figure → {out_png}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160,
                    help="smoke scale; omit/0 for full ×40 production scale")
    ap.add_argument("--n-nuc-beads", type=int, default=400)
    ap.add_argument("--warmup", type=int, default=1500)
    ap.add_argument("--softstart", type=int, default=200)
    ap.add_argument("--contract-steps", type=int, default=20000)
    ap.add_argument("--sample-every", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--areal-density", type=float, default=None,
                    help="override myosin minifilament areal density [1/µm²] "
                         "(sensitivity ONLY — not a production config change)")
    ap.add_argument("--cm-z-struct", type=float, default=None,
                    help="cortex construction anchor density z_struct (ARCH WALL-A "
                         "sensitivity; topology only, default prod 3.7)")
    ap.add_argument("--cm-bundle-mult", type=int, default=None,
                    help="cortex construction bundle_mult (ARCH WALL-A sensitivity; "
                         "topology only, default prod 2)")
    ap.add_argument("--stepping-mode", type=str, default=None,
                    choices=["binned_r0", "grip_walk", "continuous_stroke"],
                    help="override myosin stepping_mode (§9 continuous_stroke = the "
                         "continuous per-head custom force; default = manifest grip_walk)")
    ap.add_argument("--faithful", action="store_true",
                    help="use the bimodal faithful cortex (Arp2/3 branches + bimodal "
                         "lengths) instead of the uniform production mesh (ARCH test)")
    ap.add_argument("--xlink-k", type=float, default=None,
                    help="override dynamic crosslinker k_intra/k_attach [N/m] for the "
                         "transmission-lever A/B (sensitivity ONLY — production = loop18 "
                         "re-anchored 1.0e-3; pass 1.0e-7 for the pre-re-anchor state)")
    ap.add_argument("--xl-n", type=int, default=None,
                    help="ab_xlink_v2: override crosslink COUNT n_xl (density knob; prod "
                         "1000). Slater-redesign transmission A/B — sensitivity ONLY.")
    ap.add_argument("--xl-koff0", type=float, default=None,
                    help="ab_xlink_v2: override α-actinin Bell zero-force off-rate "
                         "alpha_k_off0 [1/s] (prod 0.066, Ferrer2008). Lower = more "
                         "connected/longer transmission. KU-anchored ±multiples only.")
    ap.add_argument("--xl-xbeta", type=float, default=None,
                    help="ab_xlink_v2: override α-actinin Bell length alpha_x_beta [m] "
                         "(prod 0.4e-9). Force sensitivity of unbinding.")
    ap.add_argument("--sweep-densities", type=str, default=None,
                    help="comma-list of densities [1/µm²] to sweep (mechanism confirm: γ∝ρ); "
                         "writes a γ-vs-density curve instead of a single audit")
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_active_force_budget.json")
    ap.add_argument("--out-png", type=str,
                    default="ffn_sim/outputs/h7/figs/h7_active_force_budget.png")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments

    if args.sweep_densities:
        return _run_density_sweep(args, dev, n_fil)

    cell = _build_settled_cell(
        n_filaments=n_fil, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, softstart=args.softstart, device=dev, seed=args.seed,
        areal_density=args.areal_density,
        cm_z_struct=args.cm_z_struct, cm_bundle_mult=args.cm_bundle_mult,
        faithful=args.faithful, stepping_mode=args.stepping_mode,
        xlink_k=args.xlink_k,
        xl_n=args.xl_n, xl_koff0=args.xl_koff0, xl_xbeta=args.xl_xbeta,
    )
    s = audit(cell=cell, n_contract_steps=args.contract_steps,
              sample_every=args.sample_every)
    _verdict(s)
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(s, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    _figure(s, args.out_png)
    return 0


def _run_density_sweep(args, dev, n_fil):
    """Mechanism-confirmation: γ_active vs myosin areal density (expect ~linear).

    Quantifies the density at which the active envelope/measure reaches the active
    target — turns the '~10-20× gap' into a concrete curve for the PI datum decision.
    NOT a production change: production density stays at the literature 0.6/µm²."""
    densities = [float(x) for x in args.sweep_densities.split(",")]
    band_lo = _HOSSEINI_BAND[0] * _MNM
    # MCF7-specific ACTIVE cortical tension (Hosseini/Fischer-Friedrich 2021, Biophys J
    # 120(16):3516, PMID 34022239): γ_act = 0.39-0.41 mN/m interphase (AFM confinement, the
    # active component directly). This is the MCF7 active anchor — supersedes the 0.135 estimate
    # (≈50% of the 0.27 suspended-total). See H7_MYOSIN_DENSITY_LITERATURE_2026-06-09.md.
    active_target = 0.40
    pts = []
    for rho in densities:
        print(f"\n######## areal_density = {rho:.2f} /µm² ########", flush=True)
        cell = _build_settled_cell(
            n_filaments=n_fil, n_nuc_beads=args.n_nuc_beads,
            warmup=args.warmup, softstart=args.softstart, device=dev, seed=args.seed,
            areal_density=rho,
        )
        s = audit(cell=cell, n_contract_steps=args.contract_steps,
                  sample_every=args.sample_every)
        p = s["plateau"]; a = s["analytic"]
        pts.append({
            "areal_density_per_um2": rho,
            "mesoscale_force_factor": s["scale"]["mesoscale_force_factor"],
            "g_soft_mN_m": p["g_soft_mN_m"],
            "g_analytic_envelope_mN_m": a["g_active_analytic_mN_m"],
            "engaged_frac": p["engaged_frac"],
            "mean_T_pN": p["mean_T_pN"],
        })
        print(f"  → ρ={rho:.2f}: g_soft={p['g_soft_mN_m']:.3e}  "
              f"envelope={a['g_active_analytic_mN_m']:.3e} mN/m", flush=True)
    out = {
        "note": "myosin density sensitivity (mechanism confirm γ∝ρ); production stays 0.6/µm²",
        "band_lo_mN_m": band_lo,
        "active_target_mN_m": active_target,
        "active_target_basis": "MCF7 γ_act interphase, Hosseini/Fischer-Friedrich 2021 BiophysJ PMID 34022239",
        "phase": "loading (s_grip≈0)" if args.contract_steps < 100000 else "contraction",
        "points": pts,
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"\n  sweep json → {args.out_json}", flush=True)
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rhos = [p["areal_density_per_um2"] for p in pts]
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.plot(rhos, [p["g_analytic_envelope_mN_m"] for p in pts], "s--", color="#8e44ad",
            label="analytic envelope (full engage+stall)")
    ax.plot(rhos, [p["g_soft_mN_m"] for p in pts], "o-", color="#c0392b",
            label=f"γ_soft measured ({out['phase']})")
    ax.axhline(band_lo, color="green", ls="-", lw=1.5, label="Hosseini band_lo 0.18")
    ax.axhline(active_target, color="darkgreen", ls=":", lw=1.5,
               label="MCF7 γ_act ~0.40 (Hosseini/FF 2021)")
    ax.axvline(0.6, color="grey", ls="--", lw=1, label="literature 0.6/µm² (HeLa proxy)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("myosin minifilament areal density (1/µm²)")
    ax.set_ylabel("active cortical tension γ (mN/m)")
    ax.set_title("H.7 active-γ vs myosin density — generation lever (mechanism confirm)")
    ax.legend(fontsize=8)
    fig.savefig(args.out_png, dpi=130)
    print(f"  sweep figure → {args.out_png}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
