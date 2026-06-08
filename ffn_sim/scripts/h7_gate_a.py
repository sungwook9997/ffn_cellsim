"""H.7 Gate-A experiment — does the active cortical tension RISE as contraction develops?

The σ_a/network investigation (H7_SIGMA_A_NETWORK_INVESTIGATION_2026-06-07) showed the active-γ
floor is GENERATION-limited, in two serial gates. Gate A: myosin heads bind+load but the grip
stretch s_grip stays ≈0 (heads don't WALK) → bond extensions never grow → γ_soft stays ~1e-4 mN/m,
~1000× below even the dipole ceiling. This driver runs the DECISIVE Gate-A test on the percolated
(connected_mesh) cortex: drive binding to equilibrium FAST via the numerical batch_steps decoupling
(the proven track1 trick — raise ONLY batch_steps, the D2 CFL gate auto-clamps; NO committed config
or kinetics changed), then run LONG enough for s_grip to grow toward ℓ₀, and watch γ vs s_grip.

  CONFIRM (Gate A was the wall): γ_soft / γ_soft_ik climb toward the band AS s_grip → ℓ₀, with
          mean F/F_stall ≤ 1 (Hill-valid). The floor was contraction-development.
  REFUTE  (a deeper wall): γ stays flat while s_grip → ℓ₀ and heads are organized → Gate B
          (buckling, M-SHAKE-forbidden) or a true generation deficit dominates → surface to PI.

This is the NO-PRESEED variant (default): binding develops from 0 via the kinetics, so it cannot be
accused of bypassing the kinetics. (A construction-time pre-equilibration variant is a later
refinement.) Only γ_soft / γ_soft_ik (the ACTIVE channel, F/F_stall ≤ 1) is the test — never fold in
γ_rigid / γ_passive (GATE-B channel-separation rule).

Usage:
    python -m ffn_sim.scripts.h7_gate_a --n-fil 300 --batch-steps 4000 --ticks 4000 \
        --measure-every 100 --device gpu
"""

from __future__ import annotations

import argparse
import json
import math
import time
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402
import hoomd  # noqa: E402

from ffn_sim.cortex.cortex import resolve_h3_derived  # noqa: E402
from ffn_sim.cortex.myosin import resolve_cortex_myosin  # noqa: E402
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers  # noqa: E402
from ffn_sim.cell.cell import build_cortex_full_simulation  # noqa: E402
from ffn_sim.cortex.cortical_tension import measure_cortical_tension  # noqa: E402
from ffn_sim.scripts.track1_gsoft_verify import _bound_metrics  # noqa: E402 (reuse, no copy)
from ffn_sim.common.production_policy import (  # noqa: E402
    add_production_device_args, validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
BAND = (0.35, 0.65)        # active cortical-tension band [mN/m]
DIPOLE_CEILING = 0.10      # independent-dipole generation ceiling [mN/m] (M1)
MN = 1.0e3                 # N/m -> mN/m


def _mean_s_grip_over_l0(ma, rest_length: float) -> tuple[float, float]:
    """Mean grip-stretch of BOUND heads, absolute [m] and normalised by ℓ₀."""
    bound = ma._head_bound_to_actin >= 0
    if not np.any(bound):
        return 0.0, 0.0
    s = float(np.mean(ma._head_grip_s[bound]))
    return s, s / rest_length


def _mean_backbone_r_over_r0(snap, rest_length: float) -> float:
    """Mean cortex-backbone bond length / r0 (buckling/condensation metric; <1 = condensed)."""
    types = list(snap.bonds.types)
    if "cortex-bond" not in types:
        return float("nan")
    cid = types.index("cortex-bond")
    tid = np.asarray(snap.bonds.typeid)
    grp = np.asarray(snap.bonds.group)[tid == cid]
    if grp.shape[0] == 0:
        return float("nan")
    pos = np.asarray(snap.particles.position)
    d = pos[grp[:, 1]] - pos[grp[:, 0]]
    L = np.linalg.norm(d, axis=1)
    return float(np.mean(L) / rest_length)


class SnapshotFreeMyosinGripClock(hoomd.custom.Action):
    """Advance the host-side myosin grip clock without touching HOOMD snapshots.

    This is a speed-diagnostic bridge, not a correctness replacement for
    :class:`MyosinStepUpdater`: it lets long GPU runs proceed after an initial
    warm-bind while keeping ``s_grip`` telemetry alive, but it does not retarget
    head-actin bonds or perform Bell-Evans rebinding/unbinding.
    """

    def __init__(self, myosin_action, p_myo, ell0_cortex: float) -> None:
        self.myosin_action = myosin_action
        self.p_myo = p_myo
        self.ell0_cortex = float(ell0_cortex)
        self.steps_run = 0

    def act(self, timestep: int) -> None:  # noqa: D401, ARG002
        ma = self.myosin_action
        bound = ma._head_bound_to_actin >= 0
        if not np.any(bound):
            self.steps_run += 1
            return
        s_cap = 2.0 * self.ell0_cortex * (1.0 - 1.0e-9)
        ds = self.p_myo.v0_per_head * self.p_myo.batch_dt
        ma._head_grip_s[bound] = np.minimum(s_cap, ma._head_grip_s[bound] + ds)
        self.steps_run += 1


def run(
    n_fil,
    batch_steps,
    ticks,
    measure_every,
    seed,
    device,
    connected_mesh,
    measure_mode,
    freeze_xlinks,
    snapshot_free_myosin,
    myosin_warm_bind_ticks,
    snapshot_free_grip_clock,
):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})
    # Route B production myosin (grip_walk + mesoscale force-scaling) — the physiological
    # generation config; the active-channel test must run it.
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    # Numerical batch_steps decoupling (track1): reach binding equilibrium fast; the §2 D2
    # CFL gate auto-clamps so this does NOT change kinetics. NOT a committed-config edit.
    cfg["cortex"]["myosin"]["batch_steps"] = int(batch_steps)

    p_cortex = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)

    cfl = p_myo.batch_dt * p_myo.k_off_max_at_zero_load
    if cfl > 1.0e-3 + 1e-12:
        raise RuntimeError(f"D2 CFL violated after clamp: batch_dt*k_off_max={cfl:.3e}")

    hw = build_cortex_full_simulation(
        p_cortex, p_xlinks=p_xl, p_myosin=p_myo, device=device, with_baoab=True,
        connected_mesh=connected_mesh, rng=np.random.default_rng(seed),
    )
    sim = hw["sim"]
    ma = hw["myosin_action"]
    if ma is None:
        raise RuntimeError("no MyosinStepUpdater on the sim")
    if freeze_xlinks:
        xu = hw.get("xlink_updater")
        if xu is not None:
            sim.operations.updaters.remove(xu)
    sim.run(0)

    R = p_cortex.R_cell
    ell0 = p_cortex.rest_length
    t0 = time.time()
    warm_bind_wall_s = 0.0
    grip_clock_action = None
    if snapshot_free_myosin:
        warm_ticks = int(myosin_warm_bind_ticks)
        if warm_ticks < 0:
            raise ValueError("myosin_warm_bind_ticks must be >= 0")
        if warm_ticks > 0:
            warm_t0 = time.time()
            sim.run(warm_ticks * p_myo.batch_steps)
            warm_bind_wall_s = time.time() - warm_t0

        myosin_updater = hw.get("myosin_updater")
        if myosin_updater is not None:
            sim.operations.updaters.remove(myosin_updater)
        if snapshot_free_grip_clock:
            grip_clock_action = SnapshotFreeMyosinGripClock(
                ma, p_myo, ell0_cortex=ell0
            )
            sim.operations.updaters.append(
                hoomd.update.CustomUpdater(
                    action=grip_clock_action,
                    trigger=hoomd.trigger.Periodic(p_myo.batch_steps),
                )
            )

    def _gamma():
        g = measure_cortical_tension(sim, R_cell=R, p_enclosed_volume=None)
        return g["gamma_soft"] * MN, g["gamma_soft_ik"] * MN

    rows = []
    n_meas = max(1, ticks // measure_every)
    tick_offset = int(myosin_warm_bind_ticks) if snapshot_free_myosin else 0
    for m in range(n_meas + 1):
        if m > 0:
            sim.run(measure_every * p_myo.batch_steps)
        tick = tick_offset + m * measure_every
        full_measure = measure_mode == "full" or (
            measure_mode == "final" and m == n_meas
        )
        if full_measure:
            snap = sim.state.get_snapshot()
            bf, fos, nb = _bound_metrics(ma, snap)
            gs, gik = _gamma()
            r_over_r0 = _mean_backbone_r_over_r0(snap, ell0)
        else:
            bound = ma._head_bound_to_actin >= 0
            nb = int(np.count_nonzero(bound))
            bf = float(nb / max(1, bound.size))
            fos = float("nan")
            gs = float("nan")
            gik = float("nan")
            r_over_r0 = float("nan")
        s_abs, s_rel = _mean_s_grip_over_l0(ma, ell0)
        phys_t = tick * p_myo.batch_dt
        rows.append(dict(tick=tick, phys_t_s=phys_t, bound_frac=bf, F_over_stall=fos,
                         n_bound=nb, gamma_soft=gs, gamma_soft_ik=gik,
                         s_grip_m=s_abs, s_grip_over_l0=s_rel, backbone_r_over_r0=r_over_r0))
        print(f"[tick {tick:6d} t={phys_t*1e3:7.2f}ms] bound={bf:.3f} F/Fs={fos:.2f} "
              f"γ_soft={gs:.4f} γ_ik={gik:.4f} mN/m  s_grip/ℓ₀={s_rel:.3f}  r/r0={r_over_r0:.4f}",
              flush=True)

    wall = time.time() - t0
    simulated_ticks = tick_offset + n_meas * measure_every
    simulated_baoab_steps = simulated_ticks * p_myo.batch_steps
    steps_per_s = simulated_baoab_steps / wall if wall > 0.0 else float("nan")
    # ---- verdict ----
    last = rows[-1]
    reached_band = last["gamma_soft_ik"] >= BAND[0]
    s_grown = last["s_grip_over_l0"] >= 0.5
    hill_valid = last["F_over_stall"] <= 1.0 + 1e-9
    if snapshot_free_myosin:
        verdict = (
            "DIAGNOSTIC ONLY: snapshot-free myosin fast path completed after "
            "warm-bind; long-run GPU speed is meaningful, but Gate-A physics "
            "correctness is not because dynamic bind/unbind/walk topology is frozen."
        )
    elif reached_band and hill_valid:
        verdict = "CONFIRM: active γ reached band as contraction developed (Gate A was the wall)."
    elif s_grown and not reached_band:
        verdict = ("REFUTE: s_grip developed but γ floored below band → deeper wall (Gate B "
                   "buckling / generation deficit). Surface to PI; do NOT chase by stiffening.")
    else:
        verdict = ("INCONCLUSIVE: s_grip did not develop enough in this run length — extend ticks "
                   f"(s_grip/ℓ₀={last['s_grip_over_l0']:.2f}).")
    return dict(rows=rows, wall_s=wall, verdict=verdict, n_fil=n_fil,
                batch_dt=p_myo.batch_dt, ell0=ell0, R_cell=R,
                F_stall_per_head=p_myo.F_stall_per_head,
                connected_mesh=connected_mesh,
                measure_mode=measure_mode,
                freeze_xlinks=freeze_xlinks,
                snapshot_free_myosin=snapshot_free_myosin,
                snapshot_free_grip_clock=snapshot_free_grip_clock,
                myosin_warm_bind_ticks=int(myosin_warm_bind_ticks),
                warm_bind_wall_s=warm_bind_wall_s,
                grip_clock_steps_run=(
                    int(grip_clock_action.steps_run)
                    if grip_clock_action is not None else 0
                ),
                simulated_ticks=simulated_ticks,
                simulated_baoab_steps=simulated_baoab_steps,
                steps_per_s=steps_per_s,
                physics_mode=(
                    "snapshot_free_grip_clock_static_topology"
                    if snapshot_free_myosin and snapshot_free_grip_clock
                    else (
                        "snapshot_free_static_topology"
                        if snapshot_free_myosin else "dynamic_myosin_snapshot_topology"
                    )
                ))


def _figure(res, out: Path):
    rows = res["rows"]
    t = np.array([r["phys_t_s"] for r in rows]) * 1e3  # ms
    gs = np.array([r["gamma_soft"] for r in rows])
    gik = np.array([r["gamma_soft_ik"] for r in rows])
    sg = np.array([r["s_grip_over_l0"] for r in rows])
    bf = np.array([r["bound_frac"] for r in rows])
    rr = np.array([r["backbone_r_over_r0"] for r in rows])
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    ax1.plot(t, gs, "-o", ms=3, color="#8a1f1f", label="γ_soft (MOP)")
    ax1.plot(t, gik, "-s", ms=3, color="#e8902a", label="γ_soft_ik (IK)")
    ax1.axhspan(BAND[0], BAND[1], color="#2f6fb0", alpha=0.15, label="band [0.35,0.65]")
    ax1.axhline(DIPOLE_CEILING, color="gray", ls="--", lw=1, label="dipole ceiling 0.10")
    ax1.set_xlabel("physical time (ms)"); ax1.set_ylabel("active γ (mN/m)")
    ax1.set_title("Gate A: active γ vs time"); ax1.legend(fontsize=8)
    ax2 = ax1.twinx()
    ax2.plot(t, sg, "-^", ms=3, color="green", alpha=0.7, label="s_grip/ℓ₀")
    ax2.set_ylabel("s_grip / ℓ₀", color="green"); ax2.set_ylim(0, 1.05)
    ax3.plot(t, bf, "-o", ms=3, label="bound fraction")
    ax3.plot(t, rr, "-s", ms=3, label="backbone r/r0 (buckling)")
    ax3.axhline(1.0, color="gray", ls=":", lw=1)
    ax3.set_xlabel("physical time (ms)"); ax3.set_title("engagement + buckling"); ax3.legend(fontsize=8)
    fig.suptitle(f"H.7 Gate-A (no-preseed, connected_mesh={res['connected_mesh']}) — {res['verdict'][:60]}",
                 fontweight="bold", fontsize=10)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150); plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--batch-steps", type=int, default=4000)
    ap.add_argument("--ticks", type=int, default=4000)
    ap.add_argument("--measure-every", type=int, default=100)
    ap.add_argument(
        "--measure-mode",
        choices=("full", "final"),
        default="full",
        help=(
            "full: snapshot+tension every progress point; final: progress uses "
            "myosin internal counters only and performs snapshot+tension once at the end."
        ),
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-connected-mesh", dest="connected_mesh", action="store_false")
    ap.set_defaults(connected_mesh=True)
    ap.add_argument(
        "--freeze-xlinks",
        action="store_true",
        help=(
            "Keep seeded connected-mesh crosslinks fixed by removing the dynamic "
            "xlink updater. This eliminates the xlink get_snapshot/set_snapshot "
            "path for Gate-A speed diagnostics; myosin kinetics remain active."
        ),
    )
    ap.add_argument(
        "--snapshot-free-myosin",
        action="store_true",
        help=(
            "Run a speed-diagnostic path: warm-bind with the dynamic myosin "
            "updater, then remove it so the long GPU segment has no myosin "
            "get_snapshot/set_snapshot. This freezes dynamic topology and is "
            "not a Gate-A correctness verdict."
        ),
    )
    ap.add_argument(
        "--myosin-warm-bind-ticks",
        type=int,
        default=1,
        help=(
            "Number of myosin batch ticks to run with the original dynamic "
            "updater before --snapshot-free-myosin removes it."
        ),
    )
    ap.add_argument(
        "--no-snapshot-free-grip-clock",
        dest="snapshot_free_grip_clock",
        action="store_false",
        help=(
            "Disable the cheap host-side s_grip clock used by "
            "--snapshot-free-myosin after topology is frozen."
        ),
    )
    ap.set_defaults(snapshot_free_grip_clock=True)
    ap.add_argument("--out", default=str(PKG / "outputs" / "h7" / "figs" / "h7_gate_a.png"))
    ap.add_argument("--json", default=str(PKG / "outputs" / "h7" / "production" / "h7_gate_a.json"))
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    res = run(args.n_fil, args.batch_steps, args.ticks, args.measure_every, args.seed,
              dev, args.connected_mesh, args.measure_mode, args.freeze_xlinks,
              args.snapshot_free_myosin, args.myosin_warm_bind_ticks,
              args.snapshot_free_grip_clock)
    _figure(res, Path(args.out))
    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.json, "w"), indent=1)
    print("=" * 70, flush=True)
    print(f"VERDICT: {res['verdict']}", flush=True)
    print(
        f"  wall={res['wall_s']:.0f}s  steps/s={res['steps_per_s']:.0f}  "
        f"mode={res['physics_mode']}  fig={args.out}",
        flush=True,
    )
    print("=" * 70, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
