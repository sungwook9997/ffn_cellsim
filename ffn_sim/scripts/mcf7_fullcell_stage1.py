"""MCF7 FULL-CELL STAGE-1 test — cortex (STAGE-1 grip_walk) + ALL compartments ON.

The γ-floor was diagnosed in a bare-cortex sandbox. This driver builds the FULL
MCF7 cell (cortex + H.8 membrane + H.9 nucleus + H.10 cytoplasm + KU-3.1
enclosed-volume pressure), all additive modules turned ON with MCF7 cell-type
parameters, and runs the constrained grip_walk cortex (STAGE-1: CH1 r0≈0 + CH2
series-stall) so we can see γ_soft/γ_rigid + r/r0 + grip_s in the real
multi-compartment context (per the platform-test conventions).

MCF7 cell-type (PI-ratified 2026-06-02):
  R_cell = 7.5 µm (Wagner 2011) · cytoplasm η = 65.9 Pa·s (Hu 2024) ·
  nucleus E_nuc = 0.15 kPa chromatin small-strain (MCF7 microrheology PMC12173521),
  ratio_lamin = 5 (conservative vs the ~31× indentation/microrheology span) ·
  membrane γ_mem = 0.10 mN/m (KU-3.B1 mid-band) · cortex γ band [0.35,0.65] mN/m.

Compartments are md.force.Custom → INVISIBLE to the harmonic-bond method-of-planes
estimator; the cortex γ_soft/γ_rigid still measure the cortex channel. enclosed_volume
is the CHANGE-3c load-retention term (Young-Laplace ΔP=2γ/R) the cortex builds against.

Run (smoke first):  conda run -n ffn_sim python -m ffn_sim.scripts.mcf7_fullcell_stage1 --smoke
Real A/B (small):   conda run -n ffn_sim python -m ffn_sim.scripts.mcf7_fullcell_stage1 --n-fil 2000 --n-sample 6
"""
from __future__ import annotations

import argparse
import time
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cell.nucleus import resolve_nucleus
from ffn_sim.archive.hoomd_legacy.cell.membrane_surface import resolve_membrane_surface
from ffn_sim.archive.hoomd_legacy.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume import resolve_enclosed_volume
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)
from ffn_sim.scripts.h3_ku35_aggregation import (
    aggregation_ledger,
    format_ledger,
)
from ffn_sim.common import checkpoint as _ckpt
from ffn_sim.common.production_policy import (
    DEFAULT_MCF7_TURGOR_PA,
    add_production_device_args,
    require_full_cell_physiological_baseline,
    require_physiological_turgor,
    require_production_device,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# ---- MCF7 cell-type definition (PI-ratified) -------------------------------
MCF7 = dict(
    R_cell=7.5e-6,            # m   Wagner 2011 (15 µm diameter)
    eta_celltype="MCF7",      # cytoplasm.py table -> 65.9 Pa·s (Hu 2024)
    E_nuc=4700.0,             # Pa  4.7 kPa MCF7 indentation (in KU-3.B2.1 band [1,10] kPa);
    #                              the 0.15 kPa microrheology is a distinct small-strain quantity
    ratio_lamin=1.4,          # —   KU-3.B2.1 in-situ band minimum [1.4,5.0] (softest legal;
    #                              nucleus orthogonal to γ, softer => fewer CFL-beads => speed)
    gamma_mem=1.0e-4,         # N/m 0.10 mN/m, KU-3.B1 mid-band
    R_nuc_frac=0.25,          # nucleus:cell radius ratio
    n_nuc_beads=3000,         # CFL-safe at ratio_lamin=1.4 (k_hi ∝ 1/n_beads)
)

DEFAULT_TURGOR_PA = DEFAULT_MCF7_TURGOR_PA


def _require_production_device(device: str, *, allow_cpu_dev: bool) -> None:
    """Backward-compatible wrapper for tests/importers."""
    require_production_device(
        device,
        allow_cpu_dev=allow_cpu_dev,
        hoomd_module=hoomd,
    )


def _require_physiological_pressure(
    *,
    compartments_on: bool,
    only: str | None,
    turgor_pa: float | None,
    allow_unpressurized_dev: bool,
) -> None:
    """Backward-compatible wrapper for tests/importers."""
    if not compartments_on or only is not None:
        return
    probe = type("_TurgorProbe", (), {"turgor_dP0": turgor_pa})()
    require_physiological_turgor(
        probe,
        allow_unpressurized_dev=allow_unpressurized_dev,
    )


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def _resolve_compartments(p_cortex, *, turgor_pa=DEFAULT_TURGOR_PA):
    """Resolve the 4 MCF7 compartments (all ON).

    ``turgor_pa`` sets the BASELINE osmotic turgor Π₀ (``turgor_dP0``) — the
    resting intracellular pressure that PRE-TENSIONS the cortex at construction
    (ΔP = Π₀ at V = V0, cortex carries hoop tension γ = Π₀·R/2; physiological-
    baseline, CLAUDE.md).  Default Π₀ = 133 Pa (B3 PI-ratified 2026-06-06:
    band-implied, γ_band-centre 0.50 mN/m → Π₀ = 2γ/R at R = 7.5 µm; band
    [0.35,0.65] mN/m ⇒ Π₀ ≈ 93-173 Pa).  Stewart 2011 interphase 40 Pa is the
    lower alternative, superseded for the MCF7 baseline.  Π₀ = 0 / None leaves an
    UNPRESSURISED floppy shell and is a dev/attribution condition, not production.
    """
    R = p_cortex.R_cell
    p_nuc = resolve_nucleus(
        {"E_nuc": MCF7["E_nuc"], "ratio_lamin": MCF7["ratio_lamin"]},
        R_nuc=MCF7["R_nuc_frac"] * R, n_beads=MCF7["n_nuc_beads"],
    )
    ev_cfg = {} if turgor_pa is None else {"turgor_dP0": float(turgor_pa)}
    p_ev = resolve_enclosed_volume(ev_cfg, R_cell=R)
    p_mem = resolve_membrane_surface({"gamma_mem": MCF7["gamma_mem"]}, R_cell=R)
    p_cyto = resolve_cytoplasm(cell_type=MCF7["eta_celltype"])
    return dict(p_nucleus=p_nuc, p_enclosed_volume=p_ev,
                p_membrane_surface=p_mem, p_cytoplasm=p_cyto)


def _build(cfg, *, stepping_mode, force_scaling, constrained, compartments,
           dtc=None, seed=1, equilibrate=False, n_warmup=0, device="gpu", kon_scale=1.0,
           bind_scale=1.0, connected_mesh=False, cm_z_struct=3.7, cm_bundle_mult=2,
           v0_accel=1.0, motors_off=False, allow_cpu_dev: bool = False):
    require_production_device(
        device, allow_cpu_dev=allow_cpu_dev, hoomd_module=hoomd
    )

    from dataclasses import replace as _replace
    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dtc if dtc is not None else 0.001 * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    # FIX 2026-06-04: --v0-accel was a DEAD flag (threaded to run_arm but never
    # applied to p_myo). Apply it now (scales the NMIIA head sliding velocity).
    # NOTE per KU-3.5 gate structure: this is an accelerated-dynamics knob; the
    # ACTIVE gate must be read at F/F_stall ≤ 1 (Hill-valid), so v0_accel=1 is
    # the honest operating point and large v0_accel needs the F/F_stall check.
    if v0_accel != 1.0:
        p_myo = _replace(p_myo, v0_per_head=p_myo.v0_per_head * v0_accel)
    # motors_off: passive-baseline build (KU-3.5-active ON−OFF delta).
    if motors_off:
        p_myo = _replace(p_myo, n_motors_per_cell=0)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    if kon_scale != 1.0 and not connected_mesh:  # couple_accel: accelerate xlink
        p_xl = _replace(p_xl, k_on=p_xl.k_on * kon_scale)   # binding to PERCOLATE
    if bind_scale != 1.0 and not connected_mesh:  # mesoscale reach hack — REPLACED by
        # the connected-mesh path (which derives reach=√(A/n_fil) and SEEDS the
        # percolated bridge mesh directly; no kon/bind hacks needed).
        p_xl = _replace(p_xl, max_bind_dist=p_xl.max_bind_dist * bind_scale)
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    kw = dict(p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
              rng=np.random.default_rng(seed), **compartments)
    if connected_mesh:  # FAST HYBRID: seed bridge-different-filament + bundling +
        # adhered-baseline connected mesh on the fixed-N cortex (myosin unchanged).
        kw.update(connected_mesh=True, cm_z_struct=cm_z_struct,
                  cm_bundle_mult=cm_bundle_mult)
    if constrained:
        kw.update(constrained=True, constrained_dt=dtc)
    if equilibrate:
        kw.update(equilibrate=True, equilibrate_steps=n_warmup,
                  equilibrate_softstart_steps=max(300, n_warmup // 8))
    hw = build_cortex_full_simulation(p, **kw)
    return p, p_myo, p_xl, dtc, hw


def _cfg_for(n_fil, n_motors, n_xl, stepping_mode, force_scaling, backbone_nm=700):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["R_cell"] = MCF7["R_cell"]          # MCF7 7.5 µm
    # Minifilament backbone length override. Brief literal is 700 nm (14-bead grid
    # artifact); literature NMII bipolar minifilament ≈ 300 nm (Billington; bare
    # zone ~160 nm). 300 nm shrinks the placement exclusion (backbone+100 nm) from
    # 800→400 nm, so the NATIVE motor density fits without coarse-graining → denser,
    # more native-like motor field → better force TRANSMISSION (the γ wall).
    cfg["cortex"]["myosin"]["backbone_length"] = backbone_nm * 1e-9
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    cfg["cortex"]["myosin"]["stepping_mode"] = stepping_mode
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = bool(force_scaling)
    if n_xl is not None:
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = int(n_xl)
    return cfg


def _run_fingerprint(*, stepping_mode, n_fil, n_motors, n_xl, force_scaling,
                     v0_accel, n_warmup, n_sample_eff, interval_eff, backbone_nm,
                     kon_scale, bind_scale, connected_mesh, turgor_pa, motors_off,
                     compartments_on, only, aggregation, seed) -> str:
    """Resume-gating fingerprint over every run-defining parameter.

    Any change to physics or the sampling schedule invalidates an old
    checkpoint (so a resume can never splice an incompatible run). Device is
    deliberately excluded — CPU/GPU yield the same physics and the checkpoint
    restores positions exactly, so a GPU run resumes on CPU and vice-versa.
    """
    return _ckpt.compute_fingerprint({
        "driver": "mcf7_fullcell_stage1", "stepping_mode": stepping_mode,
        "n_fil": int(n_fil), "n_motors": int(n_motors), "n_xl": int(n_xl),
        "force_scaling": bool(force_scaling), "v0_accel": float(v0_accel),
        "n_warmup": int(n_warmup), "n_sample": int(n_sample_eff),
        "interval": int(interval_eff), "backbone_nm": float(backbone_nm),
        "kon_scale": float(kon_scale), "bind_scale": float(bind_scale),
        "connected_mesh": bool(connected_mesh),
        "turgor_pa": None if turgor_pa is None else float(turgor_pa),
        "motors_off": bool(motors_off), "compartments_on": bool(compartments_on),
        "only": only, "aggregation": bool(aggregation), "seed": int(seed),
    })


def run_arm(stepping_mode, *, n_fil, n_motors, n_xl, force_scaling, v0_accel,
            couple_accel, n_warmup, n_sample, interval, smoke, device="gpu",
            compartments_on=True, only=None, backbone_nm=700, kon_scale=1.0,
            bind_scale=1.0, connected_mesh=False, turgor_pa=DEFAULT_TURGOR_PA, motors_off=False,
            aggregation=False, resume=True, allow_unpressurized_dev=False,
            allow_cpu_dev: bool = False):
    cfg = _cfg_for(n_fil, n_motors, n_xl, stepping_mode, force_scaling, backbone_nm)
    # Resolve cortex once to get R_cell for the compartments.
    p0 = resolve_h3_derived(cfg)
    if compartments_on and only:
        full = _resolve_compartments(p0, turgor_pa=turgor_pa)
        comp = {only: full[only]}
        print(f"  [ONLY {only}] attribution run (other compartments OFF)", flush=True)
    elif compartments_on:
        comp = _resolve_compartments(p0, turgor_pa=turgor_pa)
        require_full_cell_physiological_baseline(
            comp,
            allow_unpressurized_dev=allow_unpressurized_dev,
        )
        print(f"  [compartments ON] nucleus k_chrom={comp['p_nucleus'].k_chrom:.2e} "
              f"E_nuc={MCF7['E_nuc']}Pa | enclosed dP_ref={comp['p_enclosed_volume'].dP_ref:.0f}Pa "
              f"turgor={comp['p_enclosed_volume'].turgor_dP0:.0f}Pa "
              f"| membrane gamma_mem={MCF7['gamma_mem']*1e3:.2f}mN/m "
              f"| cytoplasm eta={comp['p_cytoplasm'].eta_eff:.1f}Pa·s (MCF7)", flush=True)
    else:
        comp = {}
        print("  [compartments OFF] cortex-only at MCF7 R_cell=7.5µm (attribution baseline)",
              flush=True)

    # --- dt + resume fingerprint --------------------------------------------
    # checkpoint/resume: never redo the warm-up OR an already-completed sample
    # for the SAME config. A re-run (crash recovery, sample extension, or just
    # re-invoking an identical config) loads the last checkpoint, skips the
    # warm-up, and continues from the last completed sample. dt is computed up
    # front (not via the warm-up build) so the constrained build works on resume.
    seed = 1
    tau_bend = p0.gamma_b * p0.rest_length ** 3 / p0.bending_modulus
    dtc = 0.001 * tau_bend
    n = 2 if smoke else n_sample
    iv = 2000 if smoke else interval
    run_fp = _run_fingerprint(
        stepping_mode=stepping_mode, n_fil=n_fil, n_motors=n_motors, n_xl=n_xl,
        force_scaling=force_scaling, v0_accel=v0_accel, n_warmup=n_warmup,
        n_sample_eff=n, interval_eff=iv, backbone_nm=backbone_nm,
        kon_scale=kon_scale, bind_scale=bind_scale, connected_mesh=connected_mesh,
        turgor_pa=turgor_pa, motors_off=motors_off, compartments_on=compartments_on,
        only=only, aggregation=aggregation, seed=seed)
    ckpt_base = (PKG / "outputs" / "h3" / "production" / "ku35_active" / "ckpt"
                 / f"run_{run_fp[:16]}")
    ckpt_state = _ckpt.load_checkpoint(ckpt_base, run_fp) if resume else None
    tag = "motorsOFF" if motors_off else "motorsON"

    def _autoviz(ledgers):
        """Auto-viz the aggregation ledger (production-driver auto-viz rule)."""
        if not (aggregation and ledgers):
            return
        try:
            import json
            from ffn_sim.scripts.h3_ku35_aggregation import make_aggregation_figure
            outdir = PKG / "outputs" / "h3" / "production" / "ku35_active"
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / f"agg_ledger_{stepping_mode}_{tag}.json").write_text(
                json.dumps(ledgers, indent=2))
            make_aggregation_figure(
                ledgers,
                PKG / "outputs" / "h3" / "figs" / f"ku35_aggregation_{tag}.png",
                title=(f"KU-3.5-ACTIVE aggregation ledger ({tag}, n_fil={n_fil}) "
                       f"— {ledgers[-1]['verdict'].split(':')[0]}"))
        except Exception as exc:  # noqa: BLE001
            print(f"    [AGG] WARN auto-viz failed: {exc}", flush=True)

    # --- Restore from checkpoint, or warm up fresh --------------------------
    if ckpt_state is not None:
        seed_positions = ckpt_state["positions_by_tag"]
        _r = ckpt_state["diag_list"][0] if ckpt_state["diag_list"] else {}
        samples = [tuple(s) for s in _r.get("samples", [])]
        agg_ledgers = list(_r.get("agg_ledgers", []))
        r0 = _r.get("r0")
        start_sample = int(ckpt_state["sample_index"])
        print(f"  [RESUME] checkpoint hit (fp={run_fp[:12]}): skip warm-up, "
              f"resume at sample {start_sample}/{n}", flush=True)
    else:
        # Warm-up (unconstrained, literal v0) with all compartments.
        _, _, _, dtc, hw = _build(cfg, stepping_mode=stepping_mode,
                                  force_scaling=force_scaling, constrained=False,
                                  compartments=comp, equilibrate=True, n_warmup=n_warmup,
                                  device=device, kon_scale=kon_scale, bind_scale=bind_scale,
                                  connected_mesh=connected_mesh, dtc=dtc,
                                  v0_accel=v0_accel, motors_off=motors_off,
                                  allow_cpu_dev=allow_cpu_dev)
        hw["sim"].run(0)
        seed_positions = _tagpos(hw["sim"])
        del hw
        samples, agg_ledgers, r0, start_sample = [], [], None, 0

    # --- Fast path: every sample already done (re-run of a finished config) --
    if start_sample >= n:
        print(f"  [RESUME] all {n} samples already complete — returning cached "
              f"result (no warm-up, no production rebuild)", flush=True)
        _autoviz(agg_ledgers)
        return samples

    # --- Constrained production (rigid backbone) with all compartments ------
    p, p_myo_lit, p_xl, dtc, hc = _build(cfg, stepping_mode=stepping_mode,
                                         force_scaling=force_scaling,
                                         constrained=True, compartments=comp, dtc=dtc,
                                         device=device, kon_scale=kon_scale,
                                         bind_scale=bind_scale,
                                         connected_mesh=connected_mesh,
                                         v0_accel=v0_accel, motors_off=motors_off,
                                         allow_cpu_dev=allow_cpu_dev)
    sim = hc["sim"]
    act = hc["baoab_action"]
    myo_act = hc.get("myosin_action")
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    # Seed positions (the checkpoint's last state on resume, else the warm-up).
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = seed_positions   # get_snapshot is tag-ordered
    sim.state.set_snapshot(snap)
    sim.run(0)

    nca = p.n_filaments * p.beads_per_filament
    if r0 is None:
        r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    for k in range(start_sample, n):
        sim.run(iv)
        r = _tagpos(sim)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        g_soft = _tension_method_of_planes(sim, p.R_cell)
        g_rigid = _tension_method_of_planes_rigid(sim, act, p.R_cell, dtc)
        samples.append((rmean / r0, g_soft, g_rigid))
        print(f"  [{stepping_mode}] s={k+1}/{n} r/r0={rmean/r0:.5f} "
              f"g_soft={g_soft*1e3:.3e} g_rigid={g_rigid*1e3:.3e} "
              f"g_tot={(g_soft+g_rigid)*1e3:.3e} mN/m", flush=True)
        # KU-3.5-active force-AGGREGATION ledger: localise why per-head myosin
        # force does (not) become shell tension. Best-effort — a diagnostic
        # failure must never abort a valid tension sample.
        if aggregation:
            try:
                led = aggregation_ledger(
                    sim, R_cell=p.R_cell, p_myo=p_myo_lit,
                    myosin_action=myo_act, g_soft_gate=g_soft)
                agg_ledgers.append(led)
                print(format_ledger(led), flush=True)
                # Per-minifilament bipolar-stresslet readout (does local
                # head-walking become a COHERENT cortex-scale contractile
                # stresslet? — frac_complete_pairs + coherence).
                from ffn_sim.scripts.h3_ku35_stresslet import (
                    stresslet_ledger, format_stresslet)
                sled = stresslet_ledger(
                    sim, p_myo=p_myo_lit, myosin_action=myo_act,
                    beads_per_filament=p.beads_per_filament)
                print(format_stresslet(sled), flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"    [AGG] WARN ledger failed (sample still valid): {exc}",
                      flush=True)
        # Checkpoint after each completed sample (atomic). A re-run of this exact
        # config skips the warm-up + every completed sample (no redo-every-time).
        if resume:
            try:
                _ckpt.save_checkpoint(
                    ckpt_base, sample_index=k + 1, positions_by_tag=r,
                    diag_list=[{"samples": [list(s) for s in samples],
                                "agg_ledgers": agg_ledgers, "r0": r0}],
                    frames_so_far=np.empty((0,), dtype=np.float32),
                    params_fingerprint=run_fp, rng_seed=seed)
            except Exception as exc:  # noqa: BLE001
                print(f"    [CKPT] WARN save failed (sample still valid): {exc}",
                      flush=True)

    _autoviz(agg_ledgers)
    return samples


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true", help="tiny: n_fil=200, 2 samples")
    ap.add_argument("--n-fil", type=int, default=2000)
    ap.add_argument("--n-motors", type=int, default=200)
    ap.add_argument("--n-xl", type=int, default=2000)
    ap.add_argument("--force-scaling", action="store_true", default=True)
    ap.add_argument("--v0-accel", type=float, default=1.0,
                    help="scale NMIIA head v0 (accelerated-dynamics knob; was DEAD pre-6/4). "
                         "KU-3.5-active must be read at F/F_stall<=1, so 1.0=honest.")
    ap.add_argument("--couple-accel", action="store_true", default=True)
    ap.add_argument("--n-warmup", type=int, default=8000)
    ap.add_argument("--n-sample", type=int, default=6)
    ap.add_argument("--interval", type=int, default=20000)
    ap.add_argument("--arm", choices=["both", "binned_r0", "grip_walk"],
                    default="grip_walk")
    add_production_device_args(ap, default="gpu")
    ap.add_argument("--no-compartments", action="store_true",
                    help="cortex-only baseline (attribution: isolate the compartment effect)")
    ap.add_argument("--only", choices=["p_enclosed_volume", "p_cytoplasm",
                                       "p_nucleus", "p_membrane_surface"], default=None,
                    help="attribution: enable ONLY this one compartment")
    ap.add_argument("--backbone-nm", type=float, default=700,
                    help="myosin minifilament backbone length [nm] (700=brief, 300=literature NMII)")
    ap.add_argument("--kon-scale", type=float, default=1.0,
                    help="couple_accel: scale xlink k_on for percolation (z->[2,4]); native needs ~300")
    ap.add_argument("--bind-scale", type=float, default=1.0,
                    help="mesoscale-consistent xlink reach (scales max_bind_dist, NOT k); "
                         "bind×6 + n_xl/n_fil=1.5 → z=2.96 at mesoscale (fast STAGE-2)")
    ap.add_argument("--connected-mesh", action="store_true",
                    help="CORTEX REBUILD fast hybrid: seed the bridge-different-filament + "
                         "bundling + adhered-baseline CONNECTED mesh on the fixed-N cortex "
                         "(z→3.3, giant→0.99, 0 staples; replaces kon/bind hacks). The "
                         "γ-floor payoff test (connectivity → force transmission, Kadzik-Munro).")
    ap.add_argument("--motors-off", action="store_true",
                    help="passive-baseline build (n_motors=0) for the KU-3.5-active ON-OFF delta")
    ap.add_argument("--aggregation", action="store_true",
                    help="print the KU-3.5-active force-AGGREGATION ledger per sample "
                         "(Σ|F_head|, F/F_stall, η_agg/η_medium, branch verdict) — "
                         "diagnoses why per-head myosin force does not become shell tension")
    ap.add_argument("--no-resume", action="store_true",
                    help="disable checkpoint/resume (default: ON — a re-run of the same "
                         "config skips the warm-up + all completed samples; checkpoints "
                         "live under outputs/h3/production/ku35_active/ckpt/)")
    ap.add_argument("--turgor-pa", type=float, default=DEFAULT_TURGOR_PA,
                    help="override enclosed-volume intracellular pressure dP [Pa] "
                         "(KU-3.1 default 40 interphase; band [0.35,0.65] implies ~93-173 "
                         "metaphase). Dominant g_rigid lever (Young-Laplace γ=dP·R/2).")
    ap.add_argument("--allow-unpressurized-dev", action="store_true",
                    help="explicitly allow --turgor-pa 0/None for attribution/debug runs")
    args = ap.parse_args()

    if args.smoke:
        args.n_fil, args.n_motors, args.n_xl, args.n_warmup = 200, 20, 200, 1500

    try:
        _require_production_device(args.device, allow_cpu_dev=args.allow_cpu_dev)
        _require_physiological_pressure(
            compartments_on=not args.no_compartments,
            only=args.only,
            turgor_pa=args.turgor_pa,
            allow_unpressurized_dev=args.allow_unpressurized_dev,
        )
    except ValueError as exc:
        ap.error(str(exc))

    print(f"=== MCF7 FULL-CELL STAGE-1 (R_cell=7.5µm, all compartments ON, "
          f"n_fil={args.n_fil}, arm={args.arm}, smoke={args.smoke}) ===", flush=True)
    print(f"HOOMD {hoomd.version.version}", flush=True)

    arms = ["binned_r0", "grip_walk"] if args.arm == "both" else [args.arm]
    t0 = time.time()
    results = {}
    for arm in arms:
        print(f"\n--- arm: {arm} ---", flush=True)
        results[arm] = run_arm(
            arm, n_fil=args.n_fil, n_motors=args.n_motors, n_xl=args.n_xl,
            force_scaling=args.force_scaling, v0_accel=args.v0_accel,
            couple_accel=args.couple_accel, n_warmup=args.n_warmup,
            n_sample=args.n_sample, interval=args.interval, smoke=args.smoke,
            device=args.device, compartments_on=not args.no_compartments,
            only=args.only, backbone_nm=args.backbone_nm, kon_scale=args.kon_scale,
            bind_scale=args.bind_scale, connected_mesh=args.connected_mesh,
            turgor_pa=args.turgor_pa, motors_off=args.motors_off,
            aggregation=args.aggregation, resume=not args.no_resume,
            allow_unpressurized_dev=args.allow_unpressurized_dev,
            allow_cpu_dev=args.allow_cpu_dev,
        )
    dt = time.time() - t0
    print(f"\n=== DONE in {dt:.0f}s. Full cell (cortex+membrane+nucleus+cytoplasm"
          f"+enclosed-vol) ran stably with STAGE-1 grip_walk. ===", flush=True)
    for arm, s in results.items():
        if s:
            gs = np.mean([x[1] for x in s]) * 1e3
            gr = np.mean([x[2] for x in s]) * 1e3
            rr = np.mean([x[0] for x in s])
            print(f"  {arm}: <r/r0>={rr:.5f} <g_soft>={gs:.3e} <g_rigid>={gr:.3e} "
                  f"<g_tot>={gs+gr:.3e} mN/m  (band [0.35,0.65])", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
