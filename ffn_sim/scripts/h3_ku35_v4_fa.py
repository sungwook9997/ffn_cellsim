#!/usr/bin/env python
"""KU-3.5 v4 — FA-anchored myosin-active cortex (full-fidelity tension driver).

The v3 driver (``h3_ku35_tension.py``) builds cortex + myosin + xlinks but
leaves the cortex a *free-floating shell*: with nothing anchoring it to a
substrate, myosin contraction merely shrinks the sphere and the measured
cortical tension collapses to a ~1.4e-4 mN/m floor (~3000× under the KU-3.5
[0.35, 0.65] mN/m band).

v4 closes that gap by building the **full mechanistic cell**:

    cortex (H.3) + myosin + crosslinkers
      + focal adhesions (H.4 ``resolve_h4`` → integrin clutches on a z=0
        substrate plane, with the internal substrate anchor that pulls the
        cortex shell onto the plane so clutch bonds can form)
      + enclosed-volume pressure (H.3 additive)
      + actin turnover (H.3 additive, cofilin sever + re-anneal)

built with ``reconcile_dt=True`` (unify every sub-action's dt to the
constrained production dt) and ``equilibrate=True`` (a pre-production settle
that brings the cell down onto the substrate so the FA clutch can engage).

Pipeline (mirrors v3 structure):

  * Phase 1 — warm-up (unconstrained BAOAB at CFL dt) to relax construction
    overlaps. Skipped on a fingerprint-matching checkpoint resume.
  * Phase 2 — constrained production (rigid actin backbone, M-SHAKE + Fixman,
    fast dt) with the FA layer + equilibration baked in by the builder, then
    the sampling loop measures γ_soft / γ_rigid / γ_total via the SAME
    method-of-planes used by v3 (imported, not re-implemented).

Honest-finding note (geometry): at the *physical* capture radius
(``capture_radius_R_FA`` ≈ 1.5 µm) the substrate (z=0) and the cortex shell
(r ~ R_cell ~ 10 µm) are disjoint at construction, so clutch bonds only form
after equilibration brings the cell onto the substrate. If they still don't
form at the physical radius, pass ``--fa-capture-radius`` to override (the
builder accepts the additive ``fa_clutch_capture_radius``). The smoke reports
the actual clutch-bond count either way — it does NOT fake clutch formation.

This is a SMOKE-capable driver. ``--smoke`` runs a tiny config that finishes
in a few minutes on CPU; it is too short to reach steady-state γ, so it reports
the *trend / magnitude* (is γ off the floor?), never a PASS.

Usage:
    python -m ffn_sim.scripts.h3_ku35_v4_fa --smoke --out /tmp/ku35_v4_smoke.json
    python -m ffn_sim.scripts.h3_ku35_v4_fa --n-fil 150 --n-sample 80 --device gpu
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
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume
from ffn_sim.cortex.turnover import resolve_turnover
from ffn_sim.bridge.fa import resolve_h4
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.common import checkpoint as _ckpt
from ffn_sim.common import integrity as _integrity
from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_production_device,
    validate_production_device_args,
)

# Reuse v3's measurement protocol verbatim — the γ method-of-planes (soft +
# rigid Lagrange) is the ratified KU-3.5 readout and MUST be identical so v3
# and v4 numbers are directly comparable.
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
    _tagpos,
    _restore_positions,
)

PKG = Path(__file__).resolve().parents[1]
H3_CFG = PKG / "configs" / "phase1_h3.yaml"
H4_CFG = PKG / "configs" / "phase1_h4.yaml"


def _run_fingerprint(
    *, n_fil: int, seed: int, dt_factor: float, equilibrate_steps: int,
    n_warmup: int, n_sample: int, interval: int, fa_capture_radius: float | None,
) -> str:
    """Resume-gating fingerprint for v4's physics + schedule parameters.

    Covers every parameter whose change must invalidate an existing
    checkpoint (filament count, seed, timestep, equilibration length, FA
    capture-radius override, sampling schedule). Output-only / operational
    flags are excluded so the same physical run resumes anywhere.
    """
    return _ckpt.compute_fingerprint(
        {
            "ku": "KU-3.5-v4-fa",
            "n_fil": int(n_fil),
            "seed": int(seed),
            "dt_factor": float(dt_factor),
            "equilibrate_steps": int(equilibrate_steps),
            "fa_capture_radius": (
                float(fa_capture_radius) if fa_capture_radius is not None else None
            ),
            "n_warmup": int(n_warmup),
            "n_sample": int(n_sample),
            "interval": int(interval),
        }
    )


def run(
    n_fil: int,
    seed: int,
    *,
    dt_factor: float = 0.001,
    n_warmup: int = 40_000,
    n_sample: int = 80,
    interval: int = 5_000,
    equilibrate_steps: int = 20_000,
    fa_capture_radius: float | None = None,
    device: str = "gpu",
    allow_cpu_dev: bool = False,
    out: str | None = None,
    checkpoint_every: int = 10,
    skip_integrity: bool = False,
    resume: bool = True,
) -> dict:
    """Build + run the KU-3.5 v4 FA-anchored cortex tension driver.

    Args:
        n_fil: Number of cortical filaments (mesoscopic ×40 scale).
        seed: RNG seed (construction + integrator).
        dt_factor: Production dt as a fraction of τ_bend.
        n_warmup: Unconstrained warm-up steps (relax overlaps).
        n_sample: Number of production samples.
        interval: Steps between samples.
        equilibrate_steps: Pre-production settle steps (cell onto substrate so
            the FA clutch can engage). Sized small for ``--smoke``.
        fa_capture_radius: Optional additive override for the FA clutch capture
            radius [m]. None = use the physical ``capture_radius_R_FA``.
        device: ``"cpu"`` or ``"gpu"``.
        allow_cpu_dev: Explicit CPU escape hatch for local smoke/dev runs.
        out: Output JSON path (``.npz`` frames derived alongside).
        checkpoint_every: Position-checkpoint cadence (requires ``out``).
        skip_integrity: Skip the pre-run module-integrity guard.
        resume: Resume from a fingerprint-matching checkpoint if present.

    Returns:
        Result dict (also written to ``out`` JSON).
    """
    require_production_device(
        device,
        allow_cpu_dev=allow_cpu_dev,
        hoomd_module=hoomd,
    )

    # Pre-run module-integrity guard (best-effort; never aborts a valid run).
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

    # ---- Resolve every mechanism (full-fidelity v4 cell) ----
    cfg = yaml.safe_load(open(H3_CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament
    F = p.n_filaments
    nca = F * N
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend

    p_myo = resolve_cortex_myosin(cfg, dt=dtc)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    p_ev = resolve_enclosed_volume(cfg, R_cell=p.R_cell)
    p_to = resolve_turnover(cfg, dt=dtc, rest_length=p.rest_length)

    # ---- FA layer scaled to the mesoscale cortex. ----
    # The governed phase1_h4.yaml resolves a FULL-cell FA layer
    # (35 nascent + 8 mature FAs × 50 integrins = 2150 integrins). That is
    # the native-scale physiological count, but the ×40 mesoscopic cortex
    # (CLAUDE.md, the ONLY sanctioned coarse-graining: ~1000 effective
    # filaments) has only n_fil·beads_per_filament cortex beads (~420 at
    # n_fil=60). Wiring 2150 integrin→cortex-bead clutch bonds onto ~420
    # beads piles dozens of HOOMD bonds onto single beads and trips HOOMD's
    # per-particle exclusion cap at nlist attach (the build raises). The FA
    # layer must therefore be coarse-grained to the SAME mesoscale as the
    # cortex — exactly as the verified FA-integration tests do (they use a
    # 3-FA × 4-integrin demo + L_box matched to the cortex box). We scale the
    # FA counts to the cortex bead budget (a derived, grid-invariant ratio,
    # NOT a tuned magic number) and match L_box so integrins land in-cell.
    cfg4 = yaml.safe_load(open(H4_CFG))
    b4 = cfg4["bridge"] if "bridge" in cfg4 else cfg4
    # Coarse-grain target: keep total integrins comfortably under the per-bead
    # exclusion budget. Empirically the build is stable up to ~100 integrins at
    # 420 cortex beads (n_int=200 already trips the cap); we target ≤ nca/8 so
    # the mean clutch fan-in per cortex bead stays O(1). Preserve the
    # nascent:mature ratio and a small per-FA cluster.
    n_int_budget = max(8, nca // 8)
    n_total_per_fa_meso = 4
    n_fa_meso = max(2, n_int_budget // n_total_per_fa_meso)
    ratio_mature = b4["fa"]["n_mature_per_cell"] / max(
        1, b4["fa"]["n_nascent_per_cell"] + b4["fa"]["n_mature_per_cell"]
    )
    n_mature_meso = max(1, int(round(n_fa_meso * ratio_mature)))
    n_nascent_meso = max(1, n_fa_meso - n_mature_meso)
    b4["fa"]["n_nascent_per_cell"] = int(n_nascent_meso)
    b4["fa"]["n_mature_per_cell"] = int(n_mature_meso)
    b4["fa"]["n_total_per_fa"] = int(n_total_per_fa_meso)
    # Match the FA substrate box to the cortex box so integrins/ligands land
    # inside the same periodic cell as the cortex shell (verified-test pattern).
    b4["fa"]["L_box"] = float(p.L_box)
    p_fa = resolve_h4(cfg4)
    n_fa_total_integrins = (
        (p_fa.n_nascent_per_cell + p_fa.n_mature_per_cell) * p_fa.n_total_per_fa
    )
    print(
        f"PROGRESS FA meso-scaled to cortex: {p_fa.n_nascent_per_cell} nascent + "
        f"{p_fa.n_mature_per_cell} mature FAs x {p_fa.n_total_per_fa} = "
        f"{n_fa_total_integrins} integrins (native phase1_h4 default is 2150; "
        f"coarse-grained to the {nca}-bead mesoscale cortex per CLAUDE.md x40 rule). "
        f"L_box matched to cortex ({p.L_box * 1e6:.1f} um).",
        flush=True,
    )

    # ---- FA co-tune: clamp the constrained dt to the global CFL bound. ----
    # The FA clutch's k_int_bare (1 pN/nm) is a SOFT (non-SHAKE-removed) bond,
    # so its binding CFL τ ≈ 6.5e-8 s sets a global soft-CFL bound dt_min ≈
    # 6.5e-9 s — ~10^5× smaller than the bare-cortex constrained dt
    # (0.001·τ_bend ≈ 7e-7 s). With reconcile_dt=True the builder hard-raises
    # in constrained mode if the requested constrained_dt exceeds that soft
    # bound (dt_reconcile.py "SHAKE timescale separation" contract), because
    # SHAKE removes the actin stretch but NOT the FA spring. The physically
    # correct co-tune (PI 2026-05-30 path B, optimization-and-viz memory) is to
    # run the constrained production at the FA-limited reconciled dt. We probe
    # the bound here via an UNCONSTRAINED compute_global_cfl_dt call (which only
    # REPORTS dt_min and never raises), then lower constrained_dt to it. No
    # frozen file is touched — only the dt scalar the driver requests.
    from ffn_sim.cell.dt_reconcile import compute_global_cfl_dt
    _bound = compute_global_cfl_dt(
        p, p_xlinks=p_xl, p_myosin=p_myo, p_fa=p_fa,
        p_enclosed_volume=p_ev, p_turnover=p_to,
        constrained=False, constrained_dt=None,
    )
    dt_constrained = min(dtc, float(_bound.dt_min))
    if dt_constrained < dtc:
        print(
            f"PROGRESS FA co-tune: lowering constrained dt from {dtc:.3e}s to "
            f"the FA-limited global CFL bound {dt_constrained:.3e}s "
            f"(binding element {getattr(_bound, 'limiting_element', 'fa.k_int_bare')}). "
            "Production runs at this finer dt so the soft FA clutch is stable.",
            flush=True,
        )

    dev = (
        hoomd.device.GPU()
        if device == "gpu"
        else hoomd.device.CPU(notice_level=0)
    )

    # ---- Resume probe (default ON) ----
    fingerprint = _run_fingerprint(
        n_fil=n_fil, seed=seed, dt_factor=dt_factor,
        equilibrate_steps=equilibrate_steps, n_warmup=n_warmup,
        n_sample=n_sample, interval=interval, fa_capture_radius=fa_capture_radius,
    )
    ckpt_state = _ckpt.load_checkpoint(out, fingerprint) if (resume and out) else None

    # ---- Phase 1 — unconstrained warm-up (relax overlaps). ----
    # The FA layer is NOT attached here: warm-up only relaxes the bare cortex
    # shell. The substrate anchor + FA clutch + equilibration are wired by the
    # builder in the constrained production build below (equilibrate=True),
    # which is where the cell is brought onto the substrate.
    pos_warm = None
    if ckpt_state is None:
        hw = build_cortex_full_simulation(
            p,
            p_xlinks=p_xl,
            p_myosin=p_myo,
            p_enclosed_volume=p_ev,
            p_turnover=p_to,
            device=dev,
            with_baoab=True,
            constrained=False,
            rng=np.random.default_rng(seed),
        )
        hw["sim"].run(0)
        hw["sim"].run(n_warmup)
        pos_warm = _tagpos(hw["sim"])
        # Warm-up positions are cortex-only ordering (no FA particles). The
        # production build appends FA integrins/ligands AFTER the cortex tags,
        # so the cortex block [0:nca] aligns; we restore only that block.
        del hw

    # ---- Phase 2 — FA-anchored constrained production build. ----
    # build_cortex_full_simulation with p_fa wires: integrins + ligands +
    # substrate plane, the dynamic Pereverzev clutch updater, the substrate
    # anchor, dt reconciliation, and (equilibrate=True) the cortex-onto-
    # substrate settle so clutch bonds can form before sampling.
    hc = build_cortex_full_simulation(
        p,
        p_xlinks=p_xl,
        p_myosin=p_myo,
        p_fa=p_fa,
        p_enclosed_volume=p_ev,
        p_turnover=p_to,
        device=dev,
        with_baoab=True,
        constrained=True,
        constrained_dt=dt_constrained,
        fa_clutch_capture_radius=fa_capture_radius,
        equilibrate=True,
        equilibrate_steps=equilibrate_steps,
        reconcile_dt=True,
        rng=np.random.default_rng(seed),
    )
    sim = hc["sim"]
    ma = hc["myosin_action"]
    act = hc["baoab_action"]
    n_fa_integrins = int(hc.get("n_fa_integrins", 0))
    n_fa_clutch_bonds = int(hc.get("n_fa_clutch_bonds", 0))
    cfl_result = hc.get("cfl_result")
    dt_used = hc.get("dt_used")
    eq_diag = hc.get("equilibrate_diagnostics")
    # The γ_rigid SHAKE-Lagrange conversion (T = λ·r0/Δt) needs the ACTUAL
    # integrator dt. Use the builder's dt_used (= dt_constrained, possibly
    # lowered further by B1 reconciliation) rather than the requested dtc.
    dt_meas = float(dt_used) if dt_used is not None else float(dt_constrained)

    # R1 — enable Lagrange-multiplier capture for the rigid-bond γ contribution
    # (identical to v3; no behavioural change, only the discarded λ is kept).
    if hasattr(act, "record_lambda"):
        act.record_lambda = True

    print(
        f"PROGRESS FA wired: n_fa_integrins={n_fa_integrins} "
        f"n_fa_clutch_bonds={n_fa_clutch_bonds} "
        f"fa_capture_radius={'physical' if fa_capture_radius is None else fa_capture_radius} "
        f"dt_used={dt_used} cfl_ok={getattr(cfl_result, 'ok', cfl_result)}",
        flush=True,
    )
    if n_fa_clutch_bonds == 0:
        print(
            "PROGRESS WARNING no FA clutch bonds formed after equilibration at "
            f"this capture radius ({'physical' if fa_capture_radius is None else fa_capture_radius}). "
            "Cortex is still effectively floating — γ will stay near the floor. "
            "Consider --fa-capture-radius (e.g. 5e-6) or longer "
            "--equilibrate-steps to settle the cell onto the substrate. "
            "NOTE: very large radii (>~8 µm at this mesoscale) over-fan the "
            "clutch onto single cortex beads and trip HOOMD's per-particle "
            "exclusion cap — keep the override modest.",
            flush=True,
        )

    # ---- Seed Phase-2 cortex positions (resume → checkpoint, else warm-up). ----
    # NOTE: with the FA layer the snapshot is LARGER than the cortex-only
    # warm-up (extra integrin/ligand particles after tag nca). We therefore
    # restore only the cortex block [0:nca] and leave FA particles where the
    # builder + equilibration placed them.
    if ckpt_state is not None:
        _restore_positions(sim, ckpt_state["positions_by_tag"])
    elif pos_warm is not None:
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            snap.particles.position[:nca] = pos_warm[:nca]
        sim.state.set_snapshot(snap)
    sim.run(0)

    frames = np.empty((n_sample, F, N, 3))
    diag = []
    start_sample = 0
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
        f"[start v4] device={device} n_fil={F} n_part={r_prev.shape[0]} "
        f"dt_s={dtc:.3e} dt_factor={dt_factor} n_warmup={n_warmup} "
        f"equilibrate_steps={equilibrate_steps} "
        f"n_sample={n_sample} interval={interval} "
        f"n_fa_integrins={n_fa_integrins} n_fa_clutch_bonds={n_fa_clutch_bonds} "
        f"hoomd={hoomd.version.version} gpu_build={hoomd.version.gpu_enabled}",
        flush=True,
    )
    for k in range(start_sample, n_sample):
        sim.run(interval)
        r = _tagpos(sim)
        frames[k] = r[:nca].reshape(F, N, 3)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        # LJ-CFL early-warning (raw vs min-image displacement; runaway ≫ box_L).
        box_L = float(sim.state.box.L[0])
        dvec_raw = r - r_prev
        dvec = dvec_raw - box_L * np.round(dvec_raw / box_L)
        max_disp = float(np.linalg.norm(dvec, axis=1).max())
        max_disp_raw = float(np.linalg.norm(dvec_raw, axis=1).max())
        r_prev = r
        # γ_soft (harmonic bonds incl. FA clutch) + γ_rigid (M-SHAKE Lagrange).
        gamma_soft = _tension_method_of_planes(sim, p.R_cell)
        gamma_rigid = _tension_method_of_planes_rigid(sim, act, p.R_cell, dt_meas)
        gamma_total = gamma_soft + gamma_rigid
        # Two distinct FA bond families (do NOT conflate them):
        #   * n_clutch_static = S2 fa_actin_clutch load-path bonds
        #     (integrin<->cortex-actin) formed at build. These carry the
        #     cortical tension gamma_soft picks up (S5 dynamic kinetics TODO).
        #   * n_integrin_bound = S1 dynamic Pereverzev catch bonds
        #     (integrin<->z=0 substrate-ligand) the IntegrinBondUpdater
        #     binds/unbinds (n_engaged). 0 until an integrin reaches a ligand.
        n_clutch_static = int(n_fa_clutch_bonds)
        n_integrin_bound = 0
        ia = hc.get("integrin_action")
        if ia is not None and hasattr(ia, "n_engaged"):
            try:
                n_integrin_bound = int(ia.n_engaged)
            except Exception:
                n_integrin_bound = 0
        n_clutch_live = n_clutch_static  # back-compat field (load-path bonds)
        diag.append(dict(
            step=int(sim.timestep), r_cortex_um=rmean * 1e6,
            r_over_r0=rmean / r0, myosin_engaged=int(ma.n_engaged),
            bind_total=int(ma.n_bind_total),
            step_advances=int(ma.n_step_advances_total),
            n_clutch_static=int(n_clutch_static),
            n_integrin_bound=int(n_integrin_bound),
            n_clutch_live=int(n_clutch_live),
            tension_soft_mN_per_m=gamma_soft * 1e3,
            tension_rigid_mN_per_m=gamma_rigid * 1e3,
            tension_total_mN_per_m=gamma_total * 1e3,
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
            f"clutch_static={int(n_clutch_static)} int_bound={int(n_integrin_bound)} "
            f"steps_adv={int(ma.n_step_advances_total)} "
            f"gamma_total_mN/m={gamma_total * 1e3:.3e} "
            f"gamma_soft={gamma_soft * 1e3:.3e} gamma_rigid={gamma_rigid * 1e3:.3e} "
            f"drift={float(act.max_constraint_drift):.2e} "
            f"max_disp_um={max_disp * 1e6:.3f}{warn}",
            flush=True,
        )

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

    outdir = PKG / "outputs" / "h3" / "production" / "ku35_v4_fa"
    outdir.mkdir(parents=True, exist_ok=True)
    if out is None:
        out_json = outdir / f"ku35_v4_n{n_fil}_s{seed}.json"
        frames_npz = outdir / f"ku35_v4_frames_n{n_fil}_s{seed}.npz"
    else:
        out_json = Path(out)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        frames_npz = out_json.with_suffix(".npz")
    np.savez_compressed(frames_npz, frames=frames)

    last_third = diag[max(0, len(diag) * 2 // 3):]
    g_plateau_total = float(np.mean([d["tension_total_mN_per_m"] for d in last_third]))
    g_plateau_soft = float(np.mean([d["tension_soft_mN_per_m"] for d in last_third]))
    g_plateau_rigid = float(np.mean([d["tension_rigid_mN_per_m"] for d in last_third]))

    # CFL-result is a dataclass-ish object; serialise defensively.
    def _ser(o):
        if o is None:
            return None
        if isinstance(o, (int, float, str, bool, list, dict)):
            return o
        d = {}
        for a in ("ok", "dt", "dt_cfl", "limiting", "ratio"):
            if hasattr(o, a):
                try:
                    d[a] = getattr(o, a)
                except Exception:
                    pass
        return d or str(o)

    result = dict(
        driver="h3_ku35_v4_fa",
        n_fil=int(F), seed=int(seed), dt_s=float(dt_meas), dt_factor=float(dt_factor),
        dt_requested_s=float(dtc),
        dt_speedup_vs_cfl=float(dt_meas / p.dt_cfl), r0_um=r0 * 1e6,
        n_fa_integrins=int(n_fa_integrins),
        n_fa_clutch_bonds_at_build=int(n_fa_clutch_bonds),
        n_clutch_static_final=int(diag[-1]["n_clutch_static"]) if diag else 0,
        n_integrin_bound_final=int(diag[-1]["n_integrin_bound"]) if diag else 0,
        n_clutch_live_final=int(diag[-1]["n_clutch_live"]) if diag else 0,
        fa_capture_radius_m=(float(fa_capture_radius) if fa_capture_radius is not None else None),
        equilibrate_steps=int(equilibrate_steps),
        dt_used=_ser(dt_used), cfl_result=_ser(cfl_result),
        equilibrate_diagnostics=_ser(eq_diag),
        r_final_over_r0=diag[-1]["r_over_r0"] if diag else None,
        myosin_engaged_final=diag[-1]["myosin_engaged"] if diag else 0,
        step_advances_final=diag[-1]["step_advances"] if diag else 0,
        tension_plateau_mN_per_m=g_plateau_total,
        tension_plateau_soft_mN_per_m=g_plateau_soft,
        tension_plateau_rigid_mN_per_m=g_plateau_rigid,
        tension_final_mN_per_m=diag[-1]["tension_total_mN_per_m"] if diag else None,
        tension_final_soft_mN_per_m=diag[-1]["tension_soft_mN_per_m"] if diag else None,
        tension_final_rigid_mN_per_m=diag[-1]["tension_rigid_mN_per_m"] if diag else None,
        ku35_target_mN_per_m=0.5, ku35_band_mN_per_m=[0.35, 0.65],
        n_motors=int(p_myo.n_motors_per_cell),
        max_drift=diag[-1]["max_drift"] if diag else None,
        note=("v4 FA-anchored: cortex + myosin + xlinks + FA (substrate clutch) "
              "+ enclosed-volume + turnover; equilibrate=True, reconcile_dt=True. "
              "γ_total = γ_soft + γ_rigid (method-of-planes, same as v3). "
              "A smoke run is too short for steady state — report trend only."),
        diag=diag,
    )
    out_json.write_text(json.dumps(result, indent=2))
    print(
        "RESULT " + json.dumps({k: v for k, v in result.items() if k != "diag"}),
        flush=True,
    )
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--n-warmup", type=int, default=40_000)
    ap.add_argument("--n-sample", type=int, default=80)
    ap.add_argument("--interval", type=int, default=5_000)
    ap.add_argument("--equilibrate-steps", type=int, default=20_000)
    ap.add_argument(
        "--fa-capture-radius", type=float, default=None,
        help="Additive override [m] for the FA clutch capture radius. "
             "None = physical capture_radius_R_FA (~1.5 µm). Pass e.g. 3e-6 if "
             "no clutch bonds form at the physical radius post-equilibration.",
    )
    add_production_device_args(ap, default="gpu")
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
    ap.add_argument(
        "--smoke", action="store_true",
        help="Tiny preset for a fast build+run sanity (n-fil 60, n-sample 3, "
             "interval 500, equilibrate-steps 500, n-warmup 2000). Reports γ "
             "trend off the floor; never a PASS (too short for steady state).",
    )
    args = ap.parse_args()
    validate_production_device_args(ap, args, hoomd_module=hoomd)

    if args.smoke:
        # Smoke preset overrides (only if the user left the default).
        if args.n_fil == 150:
            args.n_fil = 60
        if args.n_sample == 80:
            args.n_sample = 3
        if args.interval == 5_000:
            args.interval = 500
        if args.equilibrate_steps == 20_000:
            args.equilibrate_steps = 500
        if args.n_warmup == 40_000:
            args.n_warmup = 2_000

    run(
        args.n_fil, args.seed, dt_factor=args.dt_factor,
        n_warmup=args.n_warmup, n_sample=args.n_sample, interval=args.interval,
        equilibrate_steps=args.equilibrate_steps,
        fa_capture_radius=args.fa_capture_radius,
        device=args.device, out=args.out,
        allow_cpu_dev=args.allow_cpu_dev,
        checkpoint_every=args.checkpoint_every,
        skip_integrity=args.skip_integrity, resume=args.resume,
    )

    # Visualize-at-closeout (CLAUDE.md hard rule + production-driver auto-viz).
    # Best-effort subprocess to the existing sweep analysis so completed v4
    # seeds refresh the ensemble figures. A viz failure must NOT mask a valid
    # seed result. (--smoke skips viz: too short to be a meaningful figure.)
    if not args.smoke:
        try:
            import subprocess
            indir = (
                Path(args.out).resolve().parent if args.out
                else PKG / "outputs" / "h3" / "production" / "ku35_v4_fa"
            )
            subprocess.run(
                [sys.executable,
                 str(PKG / "scripts" / "h3_ku35_sweep_analysis.py"),
                 "--indir", str(indir)],
                check=True, cwd=str(PKG.parent),
            )
            print("FIGS auto-generated (visualize-at-closeout)", flush=True)
        except Exception as e:  # noqa: BLE001  (viz non-critical to the result)
            print(
                f"WARN visualize-at-closeout failed (seed result still valid): {e}",
                flush=True,
            )


if __name__ == "__main__":
    main()
