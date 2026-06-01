#!/usr/bin/env python
"""KU-3.5 — myosin-active cortex under constrained-BD (foundational driver).

Builds a cortex + myosin cell, warms it up with the standard BAOAB at the
CFL dt to relax construction overlaps, then hands the relaxed configuration
to the constrained integrator (rigid actin backbone, M-SHAKE + Fixman) for a
fast-dt production run where myosin actively contracts the network.

Saves the cortex-actin trajectory (for visualisation) plus per-snapshot
diagnostics: myosin engagement, total binding/stepping, mean cortex radius
(a contraction proxy), constraint drift.

NOTE — cortical tension γ: the rigorous Young-Laplace / method-of-planes
measurement with RIGID constraints (the actin "tension" is a Lagrange
multiplier, not a harmonic-bond force) is a Sanity-Gate measurement-protocol
decision pending PI ratification (like the L_p band). This driver reports a
CONTRACTION PROXY (mean radius vs time) and myosin engagement, NOT a
ratified γ — KU-3.5 PASS is not declared here.

ERM is OFF by default: at the LJ-limited production dt (~0.001·τ_bend = 0.7 μs)
the ratified k_ERM=1e-4 (τ_ERM=3.9 μs) re-violates CFL; the k_ERM for this dt
is a PI item. Myosin-active contraction is demonstrable without ERM.

Usage:
    python ffn_sim/scripts/h3_ku35_tension.py --n-fil 150 --dt-factor 0.001
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.erm import resolve_erm, attach_erm_to_simulation
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.common import checkpoint as _ckpt
from ffn_sim.common import integrity as _integrity

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def _tension_method_of_planes(sim, R_cell: float, n_planes: int = 12) -> float:
    """Mechanical cortical tension γ by method-of-planes (KU-3.5).

    For each of n_planes diametral planes through the cell center (axes
    spread on the sphere), sum the bond tensions of every harmonic bond
    that crosses the plane, projected onto the cut normal, and divide by
    the cut circumference 2πR_cell. Average across plane orientations.

    Includes ALL harmonic bonds present in the snap (myosin head springs +
    head-actin attach + xlink intra + xlink attach + any soft cortex
    component). Rigid backbone constraints (M-SHAKE Lagrange multipliers)
    are NOT included in this v1 — their shell-tension contribution is a
    separate calculation flagged for PI follow-up. Soft-bond tension is
    the dominant contractile signal once myosin engages.
    """
    btypes = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv]
        bg = np.asarray(s.bonds.group)
        bt = np.asarray(s.bonds.typeid)
    # Harmonic bond params: read from the integrator's bond force.
    bond_force = None
    for f in sim.operations.integrator.forces:
        if hasattr(f, "params") and any(t in btypes for t in getattr(f, "params", {})):
            bond_force = f; break
    if bond_force is None:
        return float("nan")
    # Per-bond k, r0 (in tag-ordered position frame).
    k_arr = np.zeros(bg.shape[0]); r0_arr = np.zeros(bg.shape[0])
    for i, tname in enumerate(btypes):
        try:
            kp = bond_force.params[tname]
            mask = (bt == i)
            k_arr[mask] = float(kp["k"]); r0_arr[mask] = float(kp["r0"])
        except Exception:
            pass
    # Bond endpoints + direction + scalar tension T = k·(|r|−r0).
    rA = pos_byTag[bg[:, 0]]; rB = pos_byTag[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1); L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)   # +T = tension (extension), −T = compression
    # Sample plane normals isotropically (Fibonacci-like).
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (i + 0.5) / n_planes
    rxy = np.sqrt(1 - z * z); theta = 2 * np.pi * i / phi
    normals = np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)
    gammas = []
    for n_hat in normals:
        a_side = rA @ n_hat; b_side = rB @ n_hat
        crossing = (a_side * b_side) < 0
        if not crossing.any():
            gammas.append(0.0); continue
        # Force along cut normal: |T| · |u · n_hat| (signed by tension dir).
        # The total force across the cut (tensile) = sum over crossing bonds.
        f_cut = (T[crossing] * (u[crossing] @ n_hat))
        # The cortex tension γ = force / circumference of the great circle.
        gammas.append(float(np.sum(f_cut)) / (2.0 * np.pi * R_cell))
    arr = np.asarray(gammas)
    # Use magnitude (sign tracks tension/compression); report mean of |γ|
    # — the tension scalar (Codex method-of-planes recommendation).
    return float(np.mean(np.abs(arr)))


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def _run_fingerprint(
    *, n_fil: int, seed: int, dt_factor: float, with_xlinks: bool,
    with_erm: bool, k_erm_fast: float, n_warmup: int, n_sample: int,
    interval: int,
) -> str:
    """Resume-gating fingerprint for this run's physics + schedule parameters.

    Covers every parameter that, if changed, must invalidate an existing
    checkpoint (filament count, seed, timestep, physics toggles, sampling
    schedule). Output-only / operational flags (``out``, checkpoint cadence,
    integrity toggle, ``device``) are deliberately excluded so the same
    physical run resumes regardless of where it is written or which machine
    continues it.

    Returns:
        Hex fingerprint from :func:`ffn_sim.common.checkpoint.compute_fingerprint`.
    """
    return _ckpt.compute_fingerprint(
        {
            "ku": "KU-3.5",
            "n_fil": int(n_fil),
            "seed": int(seed),
            "dt_factor": float(dt_factor),
            "with_xlinks": bool(with_xlinks),
            "with_erm": bool(with_erm),
            "k_erm_fast": float(k_erm_fast) if with_erm else None,
            "n_warmup": int(n_warmup),
            "n_sample": int(n_sample),
            "interval": int(interval),
        }
    )


def _restore_positions(sim, positions_by_tag: np.ndarray) -> None:
    """Write checkpointed per-tag positions back into a freshly built sim.

    Used on resume: the sim is rebuilt from scratch (same topology / forces /
    constraints), then the last-checkpointed configuration is written into the
    state via tag-ordered ``set_snapshot``. Velocities keep the rebuild
    defaults — the integrator thermostat re-randomises momenta, the documented
    (minor) statistical seam at the resume boundary.

    Args:
        sim: Active HOOMD ``Simulation`` (already built, before sampling).
        positions_by_tag: ``(N, 3)`` array of positions to restore (tag order).

    Raises:
        ValueError: If the saved particle count does not match the rebuilt sim.
    """
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        n_now = int(snap.particles.N)
        saved = np.asarray(positions_by_tag, dtype=np.float64)
        if saved.shape != (n_now, 3):
            raise ValueError(
                f"checkpoint position shape {saved.shape} != rebuilt sim "
                f"({n_now}, 3); refusing to restore"
            )
        # Mirror the warm-up restore path exactly (see below in run()): the
        # checkpointed positions come from _tagpos(), the same ordering in
        # which pos_warm is written straight into snap.particles.position, so
        # a direct assignment is correct. (The aggregate get_snapshot()
        # SnapshotParticleData exposes no per-particle .tag, hence no reindex.)
        snap.particles.position[:] = saved
    sim.state.set_snapshot(snap)


def _tension_method_of_planes_rigid(
    sim, action, R_cell: float, dt: float, n_planes: int = 12
) -> float:
    """Rigid-bond Lagrange contribution to cortical tension γ (KU-3.5).

    Companion to :func:`_tension_method_of_planes` (soft-bond M-OP), which
    misses the dominant tension source — the rigid actin backbone under
    M-SHAKE constraint. Per ``RIGID_LAGRANGE_TENSION_DESIGN.md §2`` (R1):
    the SHAKE Lagrange multiplier ``λ`` already carries the bond's
    constraint impulse; conversion via ``T_bond = λ · r₀ / Δt`` recovers
    the scalar bond tension (signed: + stretched, − compressed) which is
    then projected onto each cut normal exactly as in the soft-bond M-OP.

    Returns 0.0 if ``action.lambda_buf`` is None (``record_lambda=False``
    or no production step has run yet) — keeps the caller safe to enable
    R1 incrementally without crashing when the buffer hasn't filled.
    """
    lam = action.lambda_buf
    if lam is None or not isinstance(lam, np.ndarray):
        return 0.0
    # Host-numpy view (the Action's private buffer is device-resident cupy on a
    # GPU device after the GPU-main port; this property copies to host).
    chains = action.chains_tag_stacked
    if chains is None or chains.ndim != 2:
        return 0.0  # ragged path: not implemented here (cortex is uniform).
    r0 = float(action._chain_rest_length)

    # Bond endpoints in tag-ordered position frame (mirrors the soft M-OP).
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv]

    rA = pos_byTag[chains[:, :-1]]    # (F, m, 3) bond endpoint a
    rB = pos_byTag[chains[:, 1:]]     # (F, m, 3) bond endpoint b
    d = rB - rA                       # (F, m, 3) bond vector
    L = np.linalg.norm(d, axis=-1)    # (F, m)
    L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[..., None]         # (F, m, 3) unit bond direction
    # T_bond = λ · r₀ / Δt   (signed scalar tension, units: N)
    T = lam * (r0 / dt)               # (F, m)

    # Fibonacci-like isotropic plane normals (same convention as the
    # soft-bond M-OP so the two contributions are sampled consistently).
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (i + 0.5) / n_planes
    rxy = np.sqrt(1 - z * z); theta = 2 * np.pi * i / phi
    normals = np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)

    # Flatten (F, m) bond list for the plane-crossing reduction.
    rA_flat = rA.reshape(-1, 3)
    rB_flat = rB.reshape(-1, 3)
    u_flat = u.reshape(-1, 3)
    T_flat = T.reshape(-1)

    gammas = []
    for n_hat in normals:
        a_side = rA_flat @ n_hat
        b_side = rB_flat @ n_hat
        crossing = (a_side * b_side) < 0
        if not crossing.any():
            gammas.append(0.0); continue
        # Force along the cut normal = T · (u · n̂); cortex γ = sum / (2π R).
        f_cut = T_flat[crossing] * (u_flat[crossing] @ n_hat)
        gammas.append(float(np.sum(f_cut)) / (2.0 * np.pi * R_cell))
    arr = np.asarray(gammas)
    # Match the soft-bond M-OP convention (mean of |γ| over orientations).
    return float(np.mean(np.abs(arr)))


def run(n_fil: int, seed: int, *, dt_factor: float = 0.001, with_xlinks: bool = True,
        with_erm: bool = True, k_erm_fast: float = 5.6e-5,
        n_warmup: int = 40_000, n_sample: int = 80, interval: int = 5_000,
        device: str = "cpu", out: str | None = None,
        checkpoint_every: int = 10, skip_integrity: bool = False,
        resume: bool = True) -> dict:
    # Pre-run module-integrity guard (best-effort): a cosmetic Syncthing issue
    # must never abort an otherwise valid run, so any failure is caught and
    # downgraded to a warning here. Skipped entirely with skip_integrity=True.
    if not skip_integrity:
        try:
            _integrity.assert_clean()
            print("PROGRESS integrity check passed", flush=True)
        except _integrity.IntegrityError as exc:
            print(
                "PROGRESS WARNING integrity check flagged issues "
                f"(continuing):\n{exc}",
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001  (never crash on the guard)
            print(
                f"PROGRESS WARNING integrity check errored ({exc}); continuing",
                flush=True,
            )

    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament; F = p.n_filaments; nca = F * N
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=dtc) if with_xlinks else None
    # ERM at fast dt: re-derive k_ERM from CFL (PI 2026-05-29).
    # k_ERM = cfl_safety · γ_b / dt → k=0.1·3.9e-10/dt ≈ 5.6e-5 at dt=0.7μs.
    p_erm = None
    if with_erm:
        p_erm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell)
        from dataclasses import replace
        try: p_erm = replace(p_erm, k_ERM=k_erm_fast)
        except Exception: p_erm.k_ERM = k_erm_fast
    dev = hoomd.device.GPU() if device == "gpu" else hoomd.device.CPU(notice_level=0)

    # Resume probe (default ON): if a fingerprint-matching checkpoint exists we
    # seed the production sim directly from it and SKIP the (expensive) warm-up,
    # since the checkpointed configuration is already past equilibration. The
    # fingerprint covers the physics + schedule, so a mismatch falls back to a
    # full fresh run with no risk of splicing incompatible state.
    fingerprint = _run_fingerprint(
        n_fil=n_fil, seed=seed, dt_factor=dt_factor, with_xlinks=with_xlinks,
        with_erm=with_erm, k_erm_fast=k_erm_fast, n_warmup=n_warmup,
        n_sample=n_sample, interval=interval,
    )
    ckpt_state = _ckpt.load_checkpoint(out, fingerprint) if (resume and out) else None

    # Phase 1 — warm-up (standard BAOAB, CFL dt) to relax overlaps. Skipped on
    # a checkpoint resume (the saved configuration is already equilibrated).
    pos_warm = None
    if ckpt_state is None:
        hw = build_cortex_full_simulation(
            p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
            constrained=False, rng=np.random.default_rng(seed))
        # Attach ERM at WARM-UP dt (cfl) with config k_ERM (=1e-4); stable here.
        if with_erm:
            from dataclasses import replace as _replace
            p_erm_warm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell)
            attach_erm_to_simulation(
                hw["sim"], p_erm_warm, actin_cortex_tag_range=(0, nca),
                gamma_b=p.gamma_b, cfl_safety_factor=p.cfl_safety_factor,
                cfl_strict=True,
            )
        hw["sim"].run(0); hw["sim"].run(n_warmup)
        pos_warm = _tagpos(hw["sim"])
        del hw

    # Phase 2 — constrained production (rigid actin backbone, fast dt).
    hc = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        constrained=True, constrained_dt=dtc, rng=np.random.default_rng(seed))
    # Attach ERM at fast dt with re-derived k_ERM.
    if with_erm:
        attach_erm_to_simulation(
            hc["sim"], p_erm, actin_cortex_tag_range=(0, nca),
            gamma_b=p.gamma_b, cfl_safety_factor=p.cfl_safety_factor,
            cfl_strict=True,
        )
    sim = hc["sim"]; ma = hc["myosin_action"]; act = hc["baoab_action"]
    # R1 — RIGID_LAGRANGE_TENSION_DESIGN.md: enable Lagrange-multiplier
    # capture on the constrained Action so the rigid-bond shell-tension
    # contribution becomes available to the production-loop M-OP. No
    # behavioural change (λ is already computed inside SHAKE every step);
    # only the discarded result is now captured into act.lambda_buf.
    if hasattr(act, "record_lambda"):
        act.record_lambda = True

    # Seed Phase-2 positions: from the resumed checkpoint, else from warm-up.
    if ckpt_state is not None:
        _restore_positions(sim, ckpt_state["positions_by_tag"])
    else:
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            snap.particles.position[:] = pos_warm
        sim.state.set_snapshot(snap)
    sim.run(0)

    frames = np.empty((n_sample, F, N, 3))
    diag = []
    start_sample = 0
    # Preload accumulated history on resume so the saved JSON/npz are complete
    # (resumed work only computes the remaining samples).
    if ckpt_state is not None:
        start_sample = int(ckpt_state["sample_index"])
        diag = list(ckpt_state["diag_list"])
        for j, fr in enumerate(ckpt_state["frames_so_far"]):
            if j < n_sample:
                frames[j] = np.asarray(fr, dtype=frames.dtype).reshape(F, N, 3)
        print(f"RESUMED from sample {start_sample}/{n_sample}", flush=True)

    checkpoint_every = max(1, int(checkpoint_every))
    r_prev = _tagpos(sim)
    r0 = float(np.linalg.norm(r_prev[:nca], axis=1).mean())
    t0 = time.time()
    print(
        f"[start] device={device} n_fil={F} n_part={r_prev.shape[0]} "
        f"dt_s={dtc:.3e} dt_factor={dt_factor} n_warmup={n_warmup} "
        f"n_sample={n_sample} interval={interval} "
        f"hoomd={hoomd.version.version} gpu_build={hoomd.version.gpu_enabled}",
        flush=True,
    )
    for k in range(start_sample, n_sample):
        sim.run(interval)
        r = _tagpos(sim)
        frames[k] = r[:nca].reshape(F, N, 3)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        # LJ-CFL early-warning: max per-bead displacement over the interval
        # (canonical seed2 crash 2026-05-28 was runaway LJ overlap after motor
        # saturation → int32 image-guard overflow with no warning). Track both
        # min-image disp (physical signal for diag) and RAW disp (catches the
        # actual runaway: canonical crash had |frac coord|~1e9·box_L; legitimate
        # wraps ≤ box_L; runaway ≫ box_L).
        box_L = float(sim.state.box.L[0])
        dvec_raw = r - r_prev
        dvec = dvec_raw - box_L * np.round(dvec_raw / box_L)
        max_disp = float(np.linalg.norm(dvec, axis=1).max())
        max_disp_raw = float(np.linalg.norm(dvec_raw, axis=1).max())
        r_prev = r
        # KU-3.5 method-of-planes cortical tension.
        # γ_soft = soft-bond contribution (xlinks + myosin head-actin
        # attach + ERM + myosin internal); γ_rigid = M-SHAKE Lagrange
        # contribution (R1, RIGID_LAGRANGE_TENSION_DESIGN.md) — without
        # the latter KU-3.5 systematically under-reports tension by ~200×
        # (dominant cortex stress propagates through the actin backbone
        # which is invisible to soft-bond summation). γ_total = γ_soft + γ_rigid.
        gamma_soft = _tension_method_of_planes(sim, p.R_cell)
        gamma_rigid = _tension_method_of_planes_rigid(sim, act, p.R_cell, dtc)
        gamma_total = gamma_soft + gamma_rigid
        diag.append(dict(
            step=int(sim.timestep), r_cortex_um=rmean * 1e6,
            r_over_r0=rmean / r0, myosin_engaged=int(ma.n_engaged),
            bind_total=int(ma.n_bind_total), step_advances=int(ma.n_step_advances_total),
            tension_soft_mN_per_m=gamma_soft * 1e3,
            tension_rigid_mN_per_m=gamma_rigid * 1e3,
            tension_total_mN_per_m=gamma_total * 1e3,
            # Backward-compat: existing sweep_analysis reads tension_mN_per_m
            # as the gate quantity; promote γ_total there so analysis pipes
            # the corrected value through with no changes.
            tension_mN_per_m=gamma_total * 1e3,
            max_drift=float(act.max_constraint_drift),
            max_disp_um=max_disp * 1e6,
        ))
        wall = time.time() - t0
        eta_min = wall * (n_sample - (k + 1)) / max(k + 1, 1) / 60.0
        warn = " WARN_LJ_CFL" if max_disp_raw > 5.0 * box_L else ""
        print(
            f"PROGRESS sample={k + 1}/{n_sample} step={int(sim.timestep)} "
            f"wall={wall:.1f}s eta={eta_min:.1f}min "
            f"r/r0={rmean / r0:.3f} engaged={int(ma.n_engaged)} "
            f"steps_adv={int(ma.n_step_advances_total)} "
            f"gamma_total_mN/m={gamma_total * 1e3:.3e} gamma_soft={gamma_soft * 1e3:.3e} gamma_rigid={gamma_rigid * 1e3:.3e} drift={float(act.max_constraint_drift):.2e} "
            f"max_disp_um={max_disp * 1e6:.3f}{warn}",
            flush=True,
        )

        # Periodic position-level checkpoint (every checkpoint_every completed
        # samples, plus a final one). Best-effort: a checkpoint write failure
        # must NOT crash a valid production run. Needs an --out path to derive
        # the checkpoint location; the no-out (default-path) case skips it.
        if out and (((k + 1) % checkpoint_every == 0) or ((k + 1) == n_sample)):
            try:
                _ckpt.save_checkpoint(
                    out,
                    sample_index=k + 1,
                    positions_by_tag=r,
                    diag_list=diag,
                    frames_so_far=[frames[j] for j in range(k + 1)],
                    params_fingerprint=fingerprint,
                    rng_seed=int(seed),
                )
                print(
                    f"PROGRESS checkpoint saved at sample {k + 1}/{n_sample}",
                    flush=True,
                )
            except Exception as exc:  # noqa: BLE001  (ckpt failure non-fatal)
                print(f"WARN checkpoint save failed: {exc}", flush=True)

    outdir = PKG / "outputs" / "h3" / "production" / "ku35"
    outdir.mkdir(parents=True, exist_ok=True)
    if out is None:
        out_json = outdir / f"ku35_n{n_fil}_s{seed}.json"
        frames_npz = outdir / f"ku35_frames_n{n_fil}_s{seed}.npz"
    else:
        out_json = Path(out)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        frames_npz = out_json.with_suffix(".npz")
    np.savez_compressed(frames_npz, frames=frames)
    # Plateau (last third) tension mean — γ_total (soft + rigid).
    last_third = diag[max(0, len(diag) * 2 // 3):]
    g_plateau_total = float(np.mean([d["tension_total_mN_per_m"] for d in last_third]))
    g_plateau_soft = float(np.mean([d["tension_soft_mN_per_m"] for d in last_third]))
    g_plateau_rigid = float(np.mean([d["tension_rigid_mN_per_m"] for d in last_third]))
    result = dict(
        n_fil=int(F), seed=int(seed), dt_s=float(dtc), dt_factor=float(dt_factor),
        dt_speedup_vs_cfl=float(dtc / p.dt_cfl), r0_um=r0 * 1e6,
        with_xlinks=bool(with_xlinks), with_erm=bool(with_erm),
        k_erm_fast=float(k_erm_fast) if with_erm else None,
        r_final_over_r0=diag[-1]["r_over_r0"], myosin_engaged_final=diag[-1]["myosin_engaged"],
        step_advances_final=diag[-1]["step_advances"],
        # γ split (R1): soft + rigid + total. Backward-compat field
        # tension_plateau_mN_per_m now refers to γ_total (the gate quantity).
        tension_plateau_mN_per_m=g_plateau_total,
        tension_plateau_soft_mN_per_m=g_plateau_soft,
        tension_plateau_rigid_mN_per_m=g_plateau_rigid,
        tension_final_mN_per_m=diag[-1]["tension_total_mN_per_m"],
        tension_final_soft_mN_per_m=diag[-1]["tension_soft_mN_per_m"],
        tension_final_rigid_mN_per_m=diag[-1]["tension_rigid_mN_per_m"],
        ku35_target_mN_per_m=0.5, ku35_band_mN_per_m=[0.35, 0.65],
        n_motors=int(p_myo.n_motors_per_cell), max_drift=diag[-1]["max_drift"],
        note=("method-of-planes γ_total = γ_soft + γ_rigid (R1 in "
              "RIGID_LAGRANGE_TENSION_DESIGN.md). γ_rigid via SHAKE λ · r₀/dt; "
              "γ_soft = harmonic bond stretch sum."),
        diag=diag,
    )
    out_json.write_text(json.dumps(result, indent=2))
    print("RESULT " + json.dumps({k: v for k, v in result.items() if k != "diag"}), flush=True)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--no-xlinks", dest="with_xlinks", action="store_false", default=True)
    ap.add_argument("--no-erm", dest="with_erm", action="store_false", default=True)
    ap.add_argument("--k-erm", type=float, default=5.6e-5)
    ap.add_argument("--n-warmup", type=int, default=40_000)
    ap.add_argument("--n-sample", type=int, default=80)
    ap.add_argument("--interval", type=int, default=5_000)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument(
        "--checkpoint-every", type=int, default=10,
        help="Write a position-level checkpoint every N completed samples "
             "(requires --out to derive the checkpoint path)",
    )
    ap.add_argument(
        "--skip-integrity", action="store_true",
        help="Skip the pre-run module-integrity guard",
    )
    ap.add_argument(
        "--resume", action=argparse.BooleanOptionalAction, default=True,
        help="Resume from a fingerprint-matching checkpoint if one exists "
             "(default: on; disable with --no-resume)",
    )
    args = ap.parse_args()
    run(args.n_fil, args.seed, dt_factor=args.dt_factor,
        with_xlinks=args.with_xlinks, with_erm=args.with_erm, k_erm_fast=args.k_erm,
        n_warmup=args.n_warmup, n_sample=args.n_sample, interval=args.interval,
        device=args.device, out=args.out,
        checkpoint_every=args.checkpoint_every, skip_integrity=args.skip_integrity,
        resume=args.resume)

    # Visualize-at-closeout (CLAUDE.md hard rule + memory
    # feedback_production_driver_auto_viz): subprocess to the sweep analysis
    # so each completed seed refreshes the ensemble figures + REPORT_ku35.md
    # from whatever seeds are present on disk. Race between parallel seeds is
    # benign — the analysis is idempotent and matplotlib savefig is atomic;
    # the last writer (typically the last seed to finish) wins. Best-effort:
    # a viz failure must NOT mask a valid seed result.
    try:
        import subprocess
        indir = (Path(args.out).resolve().parent if args.out
                 else PKG / "outputs" / "h3" / "production" / "ku35")
        subprocess.run(
            [sys.executable, str(PKG / "scripts" / "h3_ku35_sweep_analysis.py"),
             "--indir", str(indir)],
            check=True, cwd=str(PKG.parent),
        )
        print("FIGS auto-generated (visualize-at-closeout)", flush=True)
    except Exception as e:  # noqa: BLE001  (viz is non-critical to the seed result)
        print(f"WARN visualize-at-closeout failed (seed result still valid): {e}",
              flush=True)


if __name__ == "__main__":
    main()
