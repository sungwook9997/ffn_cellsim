"""H.7 GATE-B smoke run — does permitting filament buckling/condensation (relaxed
M-SHAKE, compression side only) transmit Gate-A's real contraction into spanning
cortical tension?

Contract: ``docs/v2_audit/H7_GATE_B_CONTRACT_2026-06-08.md`` (LOCKED). Lever (from
design-investigation #2, ``H7_GATE_B_BUCKLING_INVESTIGATION_2_2026-06-08.md``):
~73% of load-bearing free spans are single M-SHAKE-rigid rods (no internal hinge),
so buckling is geometrically blocked by ANCHOR DENSITY irrespective of force. The
relaxed mode (``ConstrainedLeimkuhlerMatthewsBAOAB(compression_release=True)``)
makes the bond constraint UNILATERAL ``|s|≤ℓ₀``: tension side inextensible (the
physiological actin axis), compression side released → those single-segment rods
can yield axially and the network can condense.

This is the SMOKE science answer (cupy constrained integrator, CPU/GPU, ~160
filaments) — NOT the authoritative full-scale number (that is the gbook-GPU ×40
run [4], with the native relaxed mode). It exists to read the SIGN + ORDER of the
lever and to validate the controls before the expensive full run.

Observable (contract §4): the GATE criterion is on the MYOSIN-ATTRIBUTABLE cortical
tension, isolated by differencing over the myosin-OFF baseline at the SAME build /
turgor / mode:

    γ_active = [γ_soft + γ_rigid] |(myosin ON) − [γ_soft + γ_rigid] |(myosin OFF)

All channels (soft / rigid-Lagrange / passive-turgor / membrane) are reported
separately (B3/B4 discipline — turgor never folded into the active number).

Conditions run (per ``--mesoscale`` setting, both run by default):
  rigid  × {myo ON, myo OFF}  → γ_active(rigid)  = control §5.1 (≈ Gate-A plateau)
  relaxed × {myo ON, myo OFF}  → γ_active(relaxed) = the GATE number

Mesoscale-force is a PI-OPEN item (investigation #2 §margin: should a LOCAL buckling
span see the areal/aggregate-scaled per-head force ~8.48 pN, or the per-native
2.0 pN?). Not resolved here — run BOTH (``--mesoscale both``) and present the margin
two-sidedly to PI.

Pass/fail (contract §7, on γ_active in mN/m):
  CONFIRM  rises into [0.18, 0.40]  (Hosseini MCF7 interphase IQR)
  PARTIAL  ≥10× Gate-A plateau (≥~0.03) but < 0.18
  REFUTE   < ~0.03 (< 10× the Gate-A plateau)

Controls (contract §5), all reported:
  1. rigid-limit parity: relaxed-OFF reproduces the rigid γ within seed noise (the
     in-run superset check; the absolute 3.06e-3 reproduction is the full-scale run).
  2. tension-side inextensibility: max tension-side bond drift ≈ 0 (k_axial intact).
  3. condensation engaged: released-bond fraction + filament end-to-end shortening.
  4. myosin-OFF baseline: the differencing reference.

Usage (smoke, CPU dev):
    python -m ffn_sim.scripts.h7_gate_b_relaxed_smoke --n-filaments 160 \
        --warmup 600 --ticks 8 --interval 3000 --device cpu --allow-cpu-dev \
        --mesoscale both
Usage (smoke, gbook GPU — faster):
    python -u -m ffn_sim.scripts.h7_gate_b_relaxed_smoke --device gpu --mesoscale both
"""
from __future__ import annotations

import argparse
import json
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension
from ffn_sim.scripts.h7_native_fullcell_go import _tagpos

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"
_FIG_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
_MN_PER_M = 1.0e3  # N/m -> mN/m

# Gate-A plateau (full-scale native, H7_GATE_A_RESULT_2026-06-08): the rigid
# backbone active-soft contraction the lever must beat. Reported for context;
# the smoke gate uses its own in-run rigid baseline for the parity control.
GATE_A_PLATEAU_MN_M = 3.06e-3
CONFIRM_BAND = (0.18, 0.40)   # Hosseini MCF7 interphase IQR (contract §7)
PARTIAL_FLOOR = 0.03          # ≥10× Gate-A plateau


def _condensation_observable(sim, action, ell0: float) -> dict:
    """§5.3: did condensation engage? Filament end-to-end shortening + mean bond
    length, from the backbone chains. e2e/contour < 1 and mean-bond/ℓ₀ < 1 ⇒ the
    network condensed (rigid M-SHAKE holds mean-bond/ℓ₀ = 1 exactly)."""
    chains = getattr(action, "chains_tag_stacked", None)
    if chains is None:
        return {"available": False}
    chains = np.asarray(chains)
    pos = _tagpos(sim)
    F, npc = chains.shape
    contour = (npc - 1) * ell0
    # bond lengths along every chain
    seg = pos[chains[:, :-1]] - pos[chains[:, 1:]]
    blen = np.linalg.norm(seg, axis=2)                       # (F, m)
    e2e = np.linalg.norm(pos[chains[:, -1]] - pos[chains[:, 0]], axis=1)  # (F,)
    return {
        "available": True,
        "mean_bond_over_l0": float(np.mean(blen) / ell0),
        "frac_bonds_compressed": float(np.mean(blen < ell0)),
        "mean_e2e_over_contour": float(np.mean(e2e) / contour),
        "released_fraction_last_step": float(getattr(action, "released_fraction", 0.0)),
        "max_constraint_drift": float(getattr(action, "max_constraint_drift", 0.0)),
    }


def _mean_s_grip(ma, ell0: float) -> float:
    if ma is None:
        return 0.0
    bound = np.asarray(ma._head_bound_to_actin) >= 0
    if not np.any(bound):
        return 0.0
    return float(np.mean(np.asarray(ma._head_grip_s)[bound]) / ell0)


def _make_manifest(*, n_filaments, n_nuc_beads, myosin_on, mesoscale):
    """Suspended/rounded MCF7 (FA OFF, turgor ON; contract §6), myosin + mesoscale
    set per condition."""
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False  # suspended/rounded (§6)
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    if n_filaments is not None:
        co["n_filaments"] = int(n_filaments)
        co["demo_mode"] = True
    myo = co.setdefault("myosin", {})
    myo["mesoscale_force_scaling"] = bool(mesoscale)
    if not myosin_on:
        myo["n_motors_per_cell"] = 0  # §5.4 baseline: no motors (passive structure)
    if n_nuc_beads:  # 0/None → keep manifest default (3000, passes constrained CFL)
        m["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)
    return m


def _run_condition(*, manifest, device, seed, warmup, softstart, compression_release,
                   ticks, interval, measure_last, connected_mesh, label):
    """Two-phase build (unconstrained warm-up → constrained production seeded from
    warm positions), then sample the γ channels over the last ``measure_last`` ticks.
    Uses the cupy constrained Action (where ``compression_release`` lives)."""
    # ---- warm-up (unconstrained; binder engages, mesh relaxes) ----
    cell_w = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=False, with_baoab=True,
        equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=softstart, connected_mesh=connected_mesh,
    )
    cell_w.simulation.run(0)
    warm_pos = _tagpos(cell_w.simulation)
    del cell_w

    # ---- constrained production (cupy Action; relaxed flag set live) ----
    cell = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=True, equilibrate=False, connected_mesh=connected_mesh,
    )
    act = cell.baoab_action
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    act.compression_release = bool(compression_release)
    ell0 = float(act._chain_rest_length)
    # Euler buckling threshold (contract §8): F_crit = π²κ_B/ℓ₀² (single-segment,
    # DERIVED from the actin bending modulus κ_B = ℓ_p·k_BT — not a free knob). A
    # compressed backbone bond buckles (releases) only above this load; sub-
    # threshold compression stays rigid → the passive baseline is preserved.
    F_crit = None
    tau_bend = None
    if compression_release:
        kappa_B = float(cell.p_cortex.bending_modulus)
        F_crit = float(np.pi ** 2 * kappa_B / ell0 ** 2)
        act.release_load_crit = F_crit
        # τ_bend EMA window (PI 2026-06-08) — the segment bending-relaxation time
        # γ_b·ℓ₀³/κ_B; eligibility uses the SUSTAINED (τ_bend-averaged) load, not
        # the thermal per-step λ. Derived (cortex.py:450), not a free knob.
        tau_bend = float(cell.p_cortex.tau_bend)
        act.load_tau = tau_bend
    sim = cell.simulation
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = warm_pos
    sim.state.set_snapshot(snap)
    sim.run(0)

    R_cell = float(cell.p_cortex.R_cell)
    dtc = float(act.dt)
    pev = getattr(cell, "p_enclosed_volume", None)

    samples = []
    t0 = time.perf_counter()
    for k in range(ticks):
        sim.run(interval)
        # lightweight per-tick monitor (s_grip development + condensation engage)
        if k % max(1, ticks // 8) == 0 or k == ticks - 1:
            sg = _mean_s_grip(cell.myosin_action, ell0)
            rf = float(getattr(act, "released_fraction", 0.0))
            print(f"      · {label} tick {k+1}/{ticks} s_grip/ℓ0={sg:.3f} "
                  f"released={rf:.3f}", flush=True)
        if k >= ticks - measure_last:
            ct = measure_cortical_tension(
                sim, R_cell=R_cell, p_enclosed_volume=pev,
                lambda_accumulator=act, dt=dtc, n_planes=12,
            )
            cond = _condensation_observable(sim, act, ell0)
            samples.append({
                "tick": k + 1,
                "gamma_soft_mN_m": ct["gamma_soft"] * _MN_PER_M,
                "gamma_rigid_mN_m": ct["gamma_rigid"] * _MN_PER_M,
                "gamma_passive_mN_m": ct["gamma_passive"] * _MN_PER_M,
                "gamma_structural_mN_m": ct["gamma_structural"] * _MN_PER_M,
                "s_grip_over_l0": _mean_s_grip(cell.myosin_action, ell0),
                "condensation": cond,
            })
    sps = (ticks * interval) / max(time.perf_counter() - t0, 1e-9)

    def _mean(key):
        return float(np.mean([s[key] for s in samples])) if samples else 0.0

    last_cond = samples[-1]["condensation"] if samples else {}
    result = {
        "label": label,
        "compression_release": bool(compression_release),
        "gamma_soft_mN_m": _mean("gamma_soft_mN_m"),
        "gamma_rigid_mN_m": _mean("gamma_rigid_mN_m"),
        "gamma_passive_mN_m": _mean("gamma_passive_mN_m"),
        "gamma_structural_mN_m": _mean("gamma_structural_mN_m"),
        "s_grip_over_l0": _mean("s_grip_over_l0"),
        "condensation": last_cond,
        "n_samples": len(samples),
        "steps_per_s": sps,
        "R_cell_um": R_cell * 1e6,
        "ell0_nm": ell0 * 1e9,
        "F_crit_pN": (F_crit * 1e12) if F_crit is not None else None,
        "tau_bend_us": (tau_bend * 1e6) if tau_bend is not None else None,
        "n_cortex_actin": int(cell.n_cortex_actin),
    }
    print(f"    [{label}] γ_struct={result['gamma_structural_mN_m']:.4e} mN/m "
          f"(soft={result['gamma_soft_mN_m']:.3e} rigid={result['gamma_rigid_mN_m']:.3e}) "
          f"s_grip/ℓ0={result['s_grip_over_l0']:.3f} "
          f"cond(bond/ℓ0={last_cond.get('mean_bond_over_l0', float('nan')):.4f}, "
          f"released={last_cond.get('released_fraction_last_step', 0):.3f}) "
          f"{sps:.0f} steps/s", flush=True)
    del cell
    return result


def _verdict(gamma_active_relaxed_mN_m: float) -> str:
    g = gamma_active_relaxed_mN_m
    if CONFIRM_BAND[0] <= g <= CONFIRM_BAND[1]:
        return (f"CONFIRM: γ_active={g:.4f} mN/m ∈ [{CONFIRM_BAND[0]}, {CONFIRM_BAND[1]}] "
                "(Hosseini IQR) — buckling/condensation IS the transmission lever.")
    if g >= PARTIAL_FLOOR:
        return (f"PARTIAL: γ_active={g:.4f} mN/m ≥ {PARTIAL_FLOOR} (≥10× Gate-A plateau) "
                f"but < {CONFIRM_BAND[0]} — buckling contributes, insufficient; quantify "
                "the residual gap + next lever.")
    return (f"REFUTE: γ_active={g:.4f} mN/m < {PARTIAL_FLOOR} (< 10× Gate-A plateau) — "
            "buckling is NOT the lever; wall is soft long-range transmission or "
            "generation density.")


def _run_mesoscale(*, mesoscale, args, device):
    """Run the 4 conditions (rigid/relaxed × myo ON/OFF) at one mesoscale setting."""
    tag = "meso_on" if mesoscale else "meso_off"
    print(f"\n=== MESOSCALE-FORCE {'ON (areal-scaled ~8.48 pN/head)' if mesoscale else 'OFF (per-native 2.0 pN/head)'} ===",
          flush=True)
    conds = {}
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments
    for cr in (False, True):
        for myo_on in (True, False):
            name = f"{'relaxed' if cr else 'rigid'}_myo{'ON' if myo_on else 'OFF'}"
            manifest = _make_manifest(
                n_filaments=n_fil, n_nuc_beads=args.n_nuc_beads,
                myosin_on=myo_on, mesoscale=mesoscale)
            conds[name] = _run_condition(
                manifest=manifest, device=device, seed=args.seed,
                warmup=args.warmup, softstart=args.softstart,
                compression_release=cr, ticks=args.ticks, interval=args.interval,
                measure_last=args.measure_last, connected_mesh=not args.no_connected_mesh,
                label=f"{tag}/{name}")

    def _struct(name):
        return conds[name]["gamma_structural_mN_m"]

    gamma_active_rigid = _struct("rigid_myoON") - _struct("rigid_myoOFF")
    gamma_active_relaxed = _struct("relaxed_myoON") - _struct("relaxed_myoOFF")
    # control §5.1: relaxed-OFF (compression released, no myosin) vs rigid-OFF —
    # under turgor the cortex is in hoop tension, so few bonds compress ⇒ these
    # should agree (relaxed is a clean superset). Δ quantifies any drift.
    parity_off_delta = abs(_struct("relaxed_myoOFF") - _struct("rigid_myoOFF"))

    summary = {
        "mesoscale_force_scaling": mesoscale,
        "conditions": conds,
        "gamma_active_rigid_mN_m": gamma_active_rigid,
        "gamma_active_relaxed_mN_m": gamma_active_relaxed,
        "control_parity_off_delta_mN_m": parity_off_delta,
        "gate_a_plateau_mN_m": GATE_A_PLATEAU_MN_M,
        "verdict": _verdict(gamma_active_relaxed),
    }
    print(f"\n  --- mesoscale {'ON' if mesoscale else 'OFF'} SUMMARY ---", flush=True)
    print(f"    γ_active(rigid)   = {gamma_active_rigid:.4e} mN/m  "
          f"(control: should ≈ Gate-A plateau {GATE_A_PLATEAU_MN_M:.2e} at full scale)", flush=True)
    print(f"    γ_active(relaxed) = {gamma_active_relaxed:.4e} mN/m  ← GATE number", flush=True)
    print(f"    control §5.1 parity Δ(relaxed-OFF − rigid-OFF) = {parity_off_delta:.3e} mN/m", flush=True)
    rel_cond = conds["relaxed_myoON"]["condensation"]
    print(f"    §5.3 condensation (relaxed,myoON): mean-bond/ℓ0="
          f"{rel_cond.get('mean_bond_over_l0', float('nan')):.4f}, "
          f"released-frac={rel_cond.get('released_fraction_last_step', 0):.3f}, "
          f"e2e/contour={rel_cond.get('mean_e2e_over_contour', float('nan')):.4f}", flush=True)
    print(f"    VERDICT: {summary['verdict']}", flush=True)
    return summary


def _figure(results, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
    labels, ga_rigid, ga_relaxed = [], [], []
    for r in results:
        labels.append("meso ON" if r["mesoscale_force_scaling"] else "meso OFF")
        ga_rigid.append(r["gamma_active_rigid_mN_m"])
        ga_relaxed.append(r["gamma_active_relaxed_mN_m"])
    x = np.arange(len(labels))
    w = 0.35
    ax1.axhspan(CONFIRM_BAND[0], CONFIRM_BAND[1], color="#aed6f1", alpha=0.6,
                label=f"CONFIRM band {CONFIRM_BAND} (Hosseini IQR)")
    ax1.axhline(PARTIAL_FLOOR, color="orange", ls="--", lw=1.3, label=f"PARTIAL floor {PARTIAL_FLOOR}")
    ax1.axhline(GATE_A_PLATEAU_MN_M, color="grey", ls=":", lw=1.3,
                label=f"Gate-A plateau {GATE_A_PLATEAU_MN_M:.1e}")
    ax1.bar(x - w / 2, np.maximum(ga_rigid, 1e-6), w, color="#c0392b", label="γ_active rigid (control)")
    ax1.bar(x + w / 2, np.maximum(ga_relaxed, 1e-6), w, color="#2e86c1", label="γ_active relaxed (GATE)")
    ax1.set_yscale("log")
    ax1.set_xticks(x); ax1.set_xticklabels(labels)
    ax1.set_ylabel("γ_active = γ_struct(myoON − myoOFF)  [mN/m]")
    ax1.set_title("H.7 Gate-B smoke — γ_active vs CONFIRM band")
    ax1.legend(fontsize=7, loc="upper left")

    # condensation engagement (relaxed, myoON): mean-bond/ℓ0
    mb = [r["conditions"]["relaxed_myoON"]["condensation"].get("mean_bond_over_l0", np.nan)
          for r in results]
    rf = [r["conditions"]["relaxed_myoON"]["condensation"].get("released_fraction_last_step", 0)
          for r in results]
    ax2.bar(x - w / 2, mb, w, color="#16a085", label="mean bond / ℓ₀ (relaxed,myoON)")
    ax2.axhline(1.0, color="grey", ls=":", lw=1.3, label="rigid limit (=1)")
    ax2b = ax2.twinx()
    ax2b.bar(x + w / 2, rf, w, color="#8e44ad", alpha=0.7, label="released bond fraction")
    ax2b.set_ylabel("released bond fraction", color="#8e44ad")
    ax2.set_xticks(x); ax2.set_xticklabels(labels)
    ax2.set_ylabel("mean bond length / ℓ₀")
    ax2.set_title("§5.3 condensation engaged? (relaxed, myo ON)")
    ax2.legend(fontsize=7, loc="lower left")
    fig.suptitle("H.7 Gate-B smoke: does compression-release transmit contraction into γ?",
                 fontweight="bold")
    fig.savefig(out_png, dpi=140)
    print(f"  figure → {out_png}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160,
                    help="smoke scale; omit/0 for full ×40 scale")
    ap.add_argument("--n-nuc-beads", type=int, default=0,
                    help="0/omit → manifest default (3000; passes the constrained "
                         "nucleus CFL, as in Gate-A). A reduced count stiffens beads "
                         "and can violate the CFL under the constrained dt.")
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--softstart", type=int, default=150)
    ap.add_argument("--ticks", type=int, default=8, help="constrained production ticks")
    ap.add_argument("--interval", type=int, default=3000, help="steps per tick")
    ap.add_argument("--measure-last", type=int, default=3,
                    help="average γ over the last N ticks (plateau)")
    ap.add_argument("--mesoscale", choices=["on", "off", "both"], default="both",
                    help="mesoscale_force_scaling setting (PI-open; run both)")
    ap.add_argument("--no-connected-mesh", action="store_true",
                    help="disable connected_mesh (default ON — investigation #2 mesh)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tag", type=str, default="gateB_smoke")
    ap.add_argument("--out-json", type=str, default=None)
    ap.add_argument("--out-png", type=str, default=None)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    settings = {"on": [True], "off": [False], "both": [True, False]}[args.mesoscale]
    print("=" * 76, flush=True)
    print("H.7 GATE-B smoke — relaxed-constraint (unilateral M-SHAKE) buckling lever", flush=True)
    print(f"  scale: {'SMOKE n_filaments=%d' % args.n_filaments if args.n_filaments else 'FULL ×40'}"
          f"  connected_mesh={not args.no_connected_mesh}  device={args.device}", flush=True)
    print("=" * 76, flush=True)

    results = [_run_mesoscale(mesoscale=ms, args=args, device=dev) for ms in settings]

    out_json = Path(args.out_json) if args.out_json else _OUT_DIR / f"h7_{args.tag}.json"
    out_png = Path(args.out_png) if args.out_png else _FIG_DIR / f"h7_{args.tag}.png"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "scale": {"n_filaments": args.n_filaments,
                  "is_smoke": bool(args.n_filaments),
                  "connected_mesh": not args.no_connected_mesh},
        "contract": "H7_GATE_B_CONTRACT_2026-06-08.md",
        "confirm_band_mN_m": list(CONFIRM_BAND),
        "partial_floor_mN_m": PARTIAL_FLOOR,
        "results": results,
    }
    out_json.write_text(json.dumps(payload, indent=2))
    print(f"\n  json → {out_json}", flush=True)
    _figure(results, out_png)

    print("\n" + "=" * 76, flush=True)
    print("  GATE-B SMOKE VERDICTS:", flush=True)
    for r in results:
        ms = "meso ON " if r["mesoscale_force_scaling"] else "meso OFF"
        print(f"    [{ms}] {r['verdict']}", flush=True)
    print("  ⚠ SMOKE scale — authoritative γ is the full ×40 native run [4]. Mesoscale-force", flush=True)
    print("    choice is PI-open: both sides presented above.", flush=True)
    print("=" * 76, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
