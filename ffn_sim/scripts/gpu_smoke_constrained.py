"""GPU end-to-end smoke for the cupy-ported constrained BAOAB Action.

Validation-ladder step 2 (GPU_MAIN_PORT_2026-05-31.md §3b): run the SAME
constrained system on a CPU device and a GPU device and check that

  1. the GPU path runs end-to-end with no error (gpu_local_snapshot + cupy),
  2. SHAKE constraint drift stays ≤ tol on BOTH devices every block,
  3. the physics gates agree within RNG-stream noise (the cupy RNG is a
     different stream from numpy, so this is statistical, NOT bit-exact):
       (a) free rigid chain COM diffusion  D_com = kT/(n·γ),
       (b) rigid+Fixman filament bending equipartition ⟨E_bend⟩/angle ≈ the
           analytic flexible target ½k⟨(π-θ)²⟩.

Run on gbook (A5000):
    ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/gpu_smoke_constrained.py

It runs CPU then GPU automatically (skips GPU with a clear message if no CUDA
device / cupy is present, e.g. on the M1 dev box).
"""
from __future__ import annotations

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater  # noqa: F401  (parity import)
from ffn_sim.integrator.constrained_baoab import make_constrained_baoab_updater


def _make_device(kind: str):
    if kind == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


# ---------------------------------------------------------------------------
# Gate (a): free rigid chain COM diffusion
# ---------------------------------------------------------------------------
# A free (force-free) rigid chain of n_beads exercises the GPU-supported
# uniform-chain M-SHAKE + Fixman path (the cortex production path). The Fixman
# pseudo-force is internal (sums to zero per chain — Newton 3rd law), so the COM
# diffuses as D_com = kT/(n_beads·γ) regardless. (A 2-bead dimer would fall to
# the generic Gauss-Seidel shake_project, which is intentionally CPU-only —
# inherently serial per constraint, so not ported to the GPU.)
def freechain_diffusion(device_kind: str, *, n_beads=3, L0=1.0, kT=1.0,
                        gamma=1.0, dt=1e-3, n_blocks=400, block=200, seed=5):
    snap = hoomd.Snapshot()
    snap.particles.N = n_beads
    snap.particles.types = ["A"]
    pos = np.zeros((n_beads, 3))
    pos[:, 0] = (np.arange(n_beads) - (n_beads - 1) / 2.0) * L0
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = 0
    snap.configuration.box = [200, 200, 200, 0, 0, 0]
    snap.bonds.N = n_beads - 1
    snap.bonds.types = ["b"]
    snap.bonds.group[:] = [[i, i + 1] for i in range(n_beads - 1)]
    snap.bonds.typeid[:] = 0
    sim = hoomd.Simulation(device=_make_device(device_kind), seed=3)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=dt, forces=[], methods=[])
    pairs = np.array([[i, i + 1] for i in range(n_beads - 1)])
    act, upd = make_constrained_baoab_updater(
        kT=kT, gamma={"A": gamma}, dt=dt,
        constraint_pairs=pairs, constraint_lengths=np.full(n_beads - 1, L0),
        chains=[np.arange(n_beads)], seed=seed,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)

    coms, max_drift = [], 0.0
    for _ in range(n_blocks):
        sim.run(block)
        s = sim.state.get_snapshot()
        r = np.asarray(s.particles.position)
        coms.append(r.mean(0))
        max_drift = max(max_drift, act.max_constraint_drift)
    coms = np.array(coms)
    dcom = np.diff(coms, axis=0)
    msd = (dcom ** 2).sum(axis=1).mean()
    D_meas = msd / (6.0 * block * dt)
    D_expected = kT / (n_beads * gamma)
    return D_meas, D_expected, max_drift


# ---------------------------------------------------------------------------
# Gate (b): rigid+Fixman filament bending equipartition
# ---------------------------------------------------------------------------
def _analytic_flexible_E_bend(kth: float, kT: float = 1.0) -> float:
    trapz = getattr(np, "trapezoid", None) or np.trapz
    th = np.linspace(1e-6, np.pi - 1e-6, 400_000)
    w = np.sin(th) * np.exp(-0.5 * kth * (np.pi - th) ** 2 / kT)
    return float(trapz(0.5 * kth * (np.pi - th) ** 2 * w, th) / trapz(w, th))


def filament_bending(device_kind: str, *, n_beads=8, L0=1.0, kT=1.0, gamma=1.0,
                     kth=2.0, dt=1e-5, burn=100_000, samples=400, per=500,
                     seed=7):
    snap = hoomd.Snapshot()
    snap.particles.N = n_beads
    snap.particles.types = ["A"]
    pos = np.zeros((n_beads, 3))
    pos[:, 0] = (np.arange(n_beads) - (n_beads - 1) / 2.0) * L0
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = 0
    snap.configuration.box = [200, 200, 200, 0, 0, 0]
    snap.bonds.N = n_beads - 1
    snap.bonds.types = ["b"]
    snap.bonds.group[:] = [[i, i + 1] for i in range(n_beads - 1)]
    snap.bonds.typeid[:] = 0
    snap.angles.N = n_beads - 2
    snap.angles.types = ["a"]
    snap.angles.group[:] = [[i, i + 1, i + 2] for i in range(n_beads - 2)]
    snap.angles.typeid[:] = 0
    sim = hoomd.Simulation(device=_make_device(device_kind), seed=1)
    sim.create_state_from_snapshot(snap)
    angle = md.angle.Harmonic()
    angle.params["a"] = dict(k=kth, t0=np.pi)
    sim.operations.integrator = md.Integrator(dt=dt, forces=[angle], methods=[])
    pairs = np.array([[i, i + 1] for i in range(n_beads - 1)])
    act, upd = make_constrained_baoab_updater(
        kT=kT, gamma={"A": gamma}, dt=dt,
        constraint_pairs=pairs,
        constraint_lengths=np.full(n_beads - 1, L0),
        chains=[np.arange(n_beads)], seed=seed,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)
    sim.run(burn)

    e_bend, max_drift = [], 0.0
    for _ in range(samples):
        sim.run(per)
        s = sim.state.get_snapshot()
        r = np.asarray(s.particles.position)
        for i in range(n_beads - 2):
            b1 = r[i] - r[i + 1]
            b2 = r[i + 2] - r[i + 1]
            c = np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2))
            th = np.arccos(np.clip(c, -1.0, 1.0))
            e_bend.append(0.5 * kth * (np.pi - th) ** 2)
        max_drift = max(max_drift, act.max_constraint_drift)
    return float(np.mean(e_bend)), _analytic_flexible_E_bend(kth, kT), max_drift


def _run(device_kind: str, *, diff_kw=None, bend_kw=None):
    diff_kw = diff_kw or {}
    bend_kw = bend_kw or {}
    print(f"\n=== device: {device_kind.upper()} ===")
    Dm, De, drift_d = freechain_diffusion(device_kind, **diff_kw)
    rel_d = abs(Dm - De) / De
    print(f"[freechain] D_com meas={Dm:.4f} expected={De:.4f} rel={rel_d:.3f} "
          f"| max_drift={drift_d:.2e}")
    Em, Ea, drift_f = filament_bending(device_kind, **bend_kw)
    rel_f = abs(Em - Ea) / Ea
    print(f"[bending] <E_bend>/angle meas={Em:.4f} analytic={Ea:.4f} "
          f"rel={rel_f:.3f} | max_drift={drift_f:.2e}")
    return dict(D_meas=Dm, D_exp=De, rel_d=rel_d, drift_d=drift_d,
                E_meas=Em, E_an=Ea, rel_f=rel_f, drift_f=drift_f)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    # Reduced-step knobs: the per-step cupy dispatch makes the tiny-N GPU run
    # latency-bound, so allow shorter (still in-band) runs for the GPU gate.
    ap.add_argument("--burn", type=int, default=100_000)
    ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--per", type=int, default=500)
    ap.add_argument("--blocks", type=int, default=400)
    ap.add_argument("--block", type=int, default=200)
    args = ap.parse_args()
    diff_kw = dict(n_blocks=args.blocks, block=args.block)
    bend_kw = dict(burn=args.burn, samples=args.samples, per=args.per)

    cpu = _run("cpu", diff_kw=diff_kw, bend_kw=bend_kw)

    gpu_ok = False
    try:
        import cupy  # noqa: F401
        gpu = hoomd.device.GPU(notice_level=0)  # noqa: F841
        gpu_ok = True
    except Exception as exc:
        print(f"\n[GPU SKIP] no usable GPU device / cupy: {exc}")

    if not gpu_ok:
        return
    g = _run("gpu", diff_kw=diff_kw, bend_kw=bend_kw)

    print("\n=== CPU-vs-GPU agreement (within RNG-stream noise) ===")
    print(f"  dimer  D_com:   CPU {cpu['D_meas']:.4f}  GPU {g['D_meas']:.4f}")
    print(f"  bending E_bend: CPU {cpu['E_meas']:.4f}  GPU {g['E_meas']:.4f}")

    # Hard gates on BOTH devices: constraint drift + free-chain COM diffusion.
    assert cpu["drift_d"] < 1e-8 and g["drift_d"] < 1e-8, "freechain drift gate"
    assert cpu["drift_f"] < 1e-8 and g["drift_f"] < 1e-8, "bending drift gate"
    assert cpu["rel_d"] < 0.20 and g["rel_d"] < 0.20, "freechain D_com band"

    # Bending equipartition is only a VALID gate at full length: a single
    # filament's bending modes need ~1e5 burn steps to equilibrate, so a
    # reduced run under-samples and reads LOW on BOTH devices (not a bias —
    # under-equilibration). Hard-assert the 0.12 band only when burn is full;
    # otherwise it is informational with a gross-error bound (catches a real
    # temperature/sign bug while not penalising the deliberate short run).
    full_length = args.burn >= 80_000
    if full_length:
        assert cpu["rel_f"] < 0.12 and g["rel_f"] < 0.12, "bending equipartition band"
        print("\nALL GATES PASS (CPU + GPU)")
    else:
        for dev, r in (("CPU", cpu), ("GPU", g)):
            assert 0.30 < r["E_meas"] < 1.40, (
                f"{dev} bending E_bend={r['E_meas']:.3f} outside gross-error "
                "[0.30,1.40] — indicates a real bug, not under-sampling")
        print("\nHARD GATES PASS (drift + diffusion, CPU + GPU); "
              "equipartition INFORMATIONAL only (reduced burn < 80k → "
              "under-equilibrated; run full-length for the equipartition gate)")


if __name__ == "__main__":
    main()
