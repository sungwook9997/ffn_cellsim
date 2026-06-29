"""EMPIRICAL sanity gates for the L-M BAOAB-limit integrator (D3).

Drives a 100-bead harmonic polymer with the custom
``LeimkuhlerMatthewsBAOAB`` Action (replacing the ``md.methods.Langevin``
in ``hoomd_polymer_sanity.py``) and asserts:

- §6.equipartition: ⟨½ k (|r|-r0)²⟩  matches the 3D-analytical bond
  reference (½ kT · (1 + 2 kT/(k r0²)) to leading order; computed by
  exact numerical integration) to ±5% per bond.
- §6.diffusion (separate run): free bead MSD/(2 d t) ≈ kT/γ to ±10%.
- §6.reproducibility: re-running with the same seed produces bit-identical
  positions after N steps — verifies that the prv_rnds buffer and the
  Action-internal RNG are deterministic per seed (no hidden state leak).

The H.1 brief specifies a bending-mode equipartition gate at ½ kT per
triplet. That target is the 2D-equipartition value (AFINES is 2D). HOOMD
runs in 3D, where the sin(θ) volume element on the bond-angle phase
space shifts the analytical reference to ⟨U_bend⟩ ≈ 0.93 kT (numerically
integrated for k_θ=5, kT=1, t0=π). We therefore report the bending mode
as a **DIAGNOSTIC** (measured + 3D analytical reference) rather than a
gate, and surface the framing issue to PI in the H.1 closeout — the
principled 3D bending validation lives at H.2 (persistence length L_p
via tangent-correlation fit against κ_B / kT).

The L-M ↔ Euler-Maruyama O(Δt²) ↔ O(Δt) order separation is a D3
**verification** deliverable that requires a comparison run of vanilla
``md.methods.Brownian`` (E-M) at halved Δt versus this Action at the
same Δt. That comparison belongs in H.2 — the staged update to
``validation/oracles/common/sanity_gate.py`` adds a
``PHASE_1_REFERENCE_INTEGRATORS`` mechanism that already anticipates it.
It is **not** a BAOAB sign-off gate (correct sampling of Boltzmann via
the bond + diffusion gates is the sign-off claim; integrator-order is a
separate empirical verification).

Run::

    conda activate ffn_sim
    python ffn_sim/scripts/baoab_polymer_sanity.py            # equipartition + diffusion
    python ffn_sim/scripts/baoab_polymer_sanity.py --order    # + L-M order check
    python ffn_sim/scripts/baoab_polymer_sanity.py --bench    # + wall-time bench

Output (PASS/FAIL summary + JSON report) goes to
``ffn_sim/outputs/h1_baoab_freeze/``. Exit code is non-zero on any FAIL —
this is the gate the boot prompt asks PI to sign off on.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater

# --- Polymer config (matches hoomd_polymer_sanity.py for continuity) -------
N_BEADS = 100
BOND_LENGTH = 1.0
BOND_K = 100.0
ANGLE_T0 = math.pi
ANGLE_K = 5.0
KT = 1.0
GAMMA = 1.0
DT_DEFAULT = 1.0e-3       # smaller than Langevin script's 5e-3:
                          # τ_relax = γ/k = 0.01 ⇒ stable at 0.1·τ.
N_BURN_IN = 50_000
N_STEPS_PROD = 500_000
BLOCKS = 10

# --- Diffusion check params ------------------------------------------------
N_FREE_BEADS = 200
DT_DIFFUSION = 1.0e-3
N_STEPS_DIFFUSION = 50_000
BOX_DIFFUSION = 10_000.0  # effectively infinite to avoid wrap during MSD

OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h1_baoab_freeze"


def _read_by_tag(sim) -> tuple[np.ndarray, np.ndarray]:
    """Return (positions, images) ordered by particle tag.

    HOOMD's ``ParticleSorter`` reorders local-snapshot rows between
    invocations; bond/angle and MSD math relies on particle identity,
    so we always read in tag order.
    """
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        img = np.asarray(s.particles.image)
        tag = np.asarray(s.particles.tag)
        N = pos.shape[0]
        out_p = np.empty_like(pos)
        out_i = np.empty_like(img)
        out_p[tag] = pos
        out_i[tag] = img
        return out_p.copy(), out_i.copy()


# ===========================================================================
# State builders
# ===========================================================================
def _build_polymer_state(out_path: Path) -> None:
    snap = gsd.hoomd.Frame()
    snap.particles.N = N_BEADS
    snap.particles.types = ["A"]
    snap.particles.typeid = np.zeros(N_BEADS, dtype=np.uint32)
    xs = (
        np.arange(N_BEADS, dtype=np.float64) * BOND_LENGTH
        - (N_BEADS - 1) * BOND_LENGTH / 2.0
    )
    snap.particles.position = np.column_stack(
        [xs, np.zeros(N_BEADS), np.zeros(N_BEADS)]
    )
    snap.particles.mass = np.ones(N_BEADS)

    snap.bonds.N = N_BEADS - 1
    snap.bonds.types = ["polymer"]
    snap.bonds.typeid = np.zeros(N_BEADS - 1, dtype=np.uint32)
    snap.bonds.group = np.column_stack(
        [np.arange(N_BEADS - 1), np.arange(1, N_BEADS)]
    ).astype(np.uint32)

    snap.angles.N = N_BEADS - 2
    snap.angles.types = ["bend"]
    snap.angles.typeid = np.zeros(N_BEADS - 2, dtype=np.uint32)
    snap.angles.group = np.column_stack(
        [
            np.arange(N_BEADS - 2),
            np.arange(1, N_BEADS - 1),
            np.arange(2, N_BEADS),
        ]
    ).astype(np.uint32)

    side = max(2.0 * N_BEADS * BOND_LENGTH, 50.0)
    snap.configuration.box = [side, side, side, 0.0, 0.0, 0.0]
    with gsd.hoomd.open(out_path, mode="w") as f:
        f.append(snap)


def _build_free_bead_state(out_path: Path) -> None:
    """N_FREE_BEADS particles, no bonds, large box — for diffusion gate."""
    snap = gsd.hoomd.Frame()
    snap.particles.N = N_FREE_BEADS
    snap.particles.types = ["A"]
    snap.particles.typeid = np.zeros(N_FREE_BEADS, dtype=np.uint32)
    snap.particles.position = np.zeros((N_FREE_BEADS, 3), dtype=np.float64)
    snap.particles.mass = np.ones(N_FREE_BEADS)
    snap.configuration.box = [
        BOX_DIFFUSION, BOX_DIFFUSION, BOX_DIFFUSION, 0.0, 0.0, 0.0
    ]
    with gsd.hoomd.open(out_path, mode="w") as f:
        f.append(snap)


# ===========================================================================
# Polymer equipartition / diffusion / order gates
# ===========================================================================
def _bond_3d_analytical_reference(
    k: float, r0: float, kT: float, n_grid: int = 200_001
) -> float:
    """⟨½ k (|r|-r0)²⟩ for a 3D harmonic bond with r² radial measure.

    The bond vector ``r`` is 3D; the magnitude potential is harmonic on
    ``|r|``, but the volume element is ``r² dr``. For ``σ/r0 = √(kT/(k r0²))``
    not negligibly small, this shifts the mean energy above the naive
    1D ½ kT. Leading-order expansion: ½ kT · (1 + 2σ²/r0²).
    """
    sigma = math.sqrt(kT / k)
    lo = max(1e-9, r0 - 10 * sigma)
    hi = r0 + 10 * sigma
    r = np.linspace(lo, hi, n_grid)
    U = 0.5 * k * (r - r0) ** 2
    w = np.exp(-U / kT) * r * r
    Z = np.trapezoid(w, r)
    return float(np.trapezoid(U * w, r) / Z)


def _bending_3d_analytical_reference(
    k_theta: float, kT: float, n_grid: int = 200_001
) -> float:
    """⟨½ k_θ (θ−π)²⟩ for the 3D harmonic-on-sin(θ) angle ensemble.

    Numerical integration over φ = π − θ ∈ (0, π) of
    ``U(φ) · exp(-U/kT) · sin(φ)`` against the partition function.
    Matches the principled 3D analytical reference vs the brief's
    2D ½ kT target.
    """
    phi = np.linspace(1e-9, math.pi - 1e-9, n_grid)
    U = 0.5 * k_theta * phi**2
    w = np.exp(-U / kT) * np.sin(phi)
    Z = np.trapezoid(w, phi)
    return float(np.trapezoid(U * w, phi) / Z)


def _bond_potential_energy(pos: np.ndarray) -> float:
    """Σ ½ k (|r_{i+1} − r_i| − r0)² over the N_BEADS-1 chain bonds."""
    d = np.linalg.norm(pos[1:] - pos[:-1], axis=1)
    return 0.5 * BOND_K * np.sum((d - BOND_LENGTH) ** 2)


def _angle_potential_energy(pos: np.ndarray) -> float:
    """Σ ½ k_θ (θ_i − π)² over the N_BEADS-2 bending triplets."""
    r1 = pos[:-2] - pos[1:-1]
    r2 = pos[2:] - pos[1:-1]
    cos_t = np.einsum("ij,ij->i", r1, r2) / (
        np.linalg.norm(r1, axis=1) * np.linalg.norm(r2, axis=1) + 1e-30
    )
    cos_t = np.clip(cos_t, -1.0, 1.0)
    theta = np.arccos(cos_t)
    return 0.5 * ANGLE_K * np.sum((theta - ANGLE_T0) ** 2)


def run_polymer(
    dt: float,
    n_steps: int = N_STEPS_PROD,
    n_burn_in: int = N_BURN_IN,
    blocks: int = BLOCKS,
    seed: int = 42,
    bench: bool = False,
) -> dict:
    """Run polymer under BAOAB, return block-averaged equipartition stats."""
    initial = OUT_DIR / "polymer_initial.gsd"
    _build_polymer_state(initial)

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
    sim.create_state_from_gsd(filename=str(initial))

    bond = md.bond.Harmonic()
    bond.params["polymer"] = dict(k=BOND_K, r0=BOND_LENGTH)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=ANGLE_K, t0=ANGLE_T0)

    ig = md.Integrator(dt=dt)
    ig.forces.append(bond)
    ig.forces.append(angle)
    sim.operations.integrator = ig

    action, updater = make_baoab_updater(
        kT=KT, gamma={"A": GAMMA}, dt=dt, seed=seed + 1
    )
    sim.operations.updaters.append(updater)

    # Burn-in
    sim.run(n_burn_in)

    # Production: sample bond/angle PE every (n_steps/(blocks·N_samples_per_block))
    # We accumulate per-block sums to compute block averages and per-block stderr.
    samples_per_block = 200
    sample_interval = max(1, n_steps // (blocks * samples_per_block))
    block_bond_pe = np.zeros(blocks)
    block_angle_pe = np.zeros(blocks)
    counts = np.zeros(blocks, dtype=int)

    t_start = time.perf_counter()
    total_run = 0
    for b in range(blocks):
        for _ in range(samples_per_block):
            sim.run(sample_interval)
            total_run += sample_interval
            pos, _ = _read_by_tag(sim)
            block_bond_pe[b] += _bond_potential_energy(pos)
            block_angle_pe[b] += _angle_potential_energy(pos)
            counts[b] += 1
    elapsed = time.perf_counter() - t_start

    # Per-bond / per-angle average (degrees-of-freedom normalisation)
    n_bonds = N_BEADS - 1
    n_angles = N_BEADS - 2
    per_bond_pe_blocks = block_bond_pe / counts / n_bonds
    per_angle_pe_blocks = block_angle_pe / counts / n_angles
    per_bond_pe = float(np.mean(per_bond_pe_blocks))
    per_angle_pe = float(np.mean(per_angle_pe_blocks))
    bond_stderr = float(
        np.std(per_bond_pe_blocks, ddof=1) / math.sqrt(blocks)
    )
    angle_stderr = float(
        np.std(per_angle_pe_blocks, ddof=1) / math.sqrt(blocks)
    )

    # Bond equipartition expectation in 3D: harmonic on the bond magnitude
    # with r² volume element. Naive ½ kT is the 1D limit; the 3D-corrected
    # reference is ½ kT (1 + 2 σ²/r0²) to leading order.
    target_bond = _bond_3d_analytical_reference(BOND_K, BOND_LENGTH, KT)
    target_angle_3d = _bending_3d_analytical_reference(ANGLE_K, KT)
    bond_rel_err = (per_bond_pe - target_bond) / target_bond
    angle_rel_err_2d = (per_angle_pe - 0.5 * KT) / (0.5 * KT)
    angle_rel_err_3d = (per_angle_pe - target_angle_3d) / target_angle_3d

    result = {
        "dt": dt,
        "n_steps": total_run,
        "n_burn_in": n_burn_in,
        "blocks": blocks,
        "per_bond_pe": per_bond_pe,
        "per_bond_pe_stderr": bond_stderr,
        "per_angle_pe": per_angle_pe,
        "per_angle_pe_stderr": angle_stderr,
        "target_bond_3d_analytical": target_bond,
        "target_angle_2d_brief_naive_half_kT": 0.5 * KT,
        "target_angle_3d_analytical": target_angle_3d,
        "bond_rel_err_vs_3d_ref": bond_rel_err,
        "angle_rel_err_vs_2d_brief": angle_rel_err_2d,
        "angle_rel_err_vs_3d_ref": angle_rel_err_3d,
        "wall_time_s": elapsed,
        "steps_per_s": total_run / elapsed if elapsed > 0 else float("inf"),
    }
    if bench:
        result["bench_per_1k_steps_s"] = elapsed / total_run * 1000.0
    return result


def run_reproducibility(
    dt: float, n_steps: int = 5000, seed: int = 99
) -> dict:
    """Two independent simulations with the same seed → bit-identical positions."""
    def _final_positions(s: int) -> np.ndarray:
        initial = OUT_DIR / "repro_initial.gsd"
        _build_polymer_state(initial)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=s)
        sim.create_state_from_gsd(filename=str(initial))
        bond = md.bond.Harmonic()
        bond.params["polymer"] = dict(k=BOND_K, r0=BOND_LENGTH)
        angle = md.angle.Harmonic()
        angle.params["bend"] = dict(k=ANGLE_K, t0=ANGLE_T0)
        ig = md.Integrator(dt=dt)
        ig.forces.append(bond)
        ig.forces.append(angle)
        sim.operations.integrator = ig
        _, updater = make_baoab_updater(
            kT=KT, gamma={"A": GAMMA}, dt=dt, seed=s + 1
        )
        sim.operations.updaters.append(updater)
        sim.run(n_steps)
        pos, _ = _read_by_tag(sim)
        return pos

    r_a = _final_positions(seed)
    r_b = _final_positions(seed)
    return {
        "seed": seed,
        "n_steps": n_steps,
        "max_abs_diff": float(np.max(np.abs(r_a - r_b))),
    }


def run_free_diffusion(dt: float = DT_DIFFUSION, seed: int = 17) -> dict:
    initial = OUT_DIR / "diff_initial.gsd"
    _build_free_bead_state(initial)

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
    sim.create_state_from_gsd(filename=str(initial))
    ig = md.Integrator(dt=dt)
    sim.operations.integrator = ig
    action, updater = make_baoab_updater(
        kT=KT, gamma={"A": GAMMA}, dt=dt, seed=seed + 1
    )
    sim.operations.updaters.append(updater)
    sim.run(0)

    r0, img0 = _read_by_tag(sim)
    box = sim.state.box
    L = np.array([box.Lx, box.Ly, box.Lz])

    sample_every = 1000
    n_samples = N_STEPS_DIFFUSION // sample_every
    msd = np.zeros(n_samples)
    times = np.zeros(n_samples)
    for k in range(n_samples):
        sim.run(sample_every)
        r, img = _read_by_tag(sim)
        d = (r + (img - img0).astype(float) * L) - r0
        msd[k] = float(np.mean(np.sum(d ** 2, axis=1)))
        times[k] = (k + 1) * sample_every * dt

    # Fit MSD = 6 D t (3D); slope → D
    # Drop the first 10% as transient.
    start = max(1, n_samples // 10)
    slope, _ = np.polyfit(times[start:], msd[start:], 1)
    D_measured = slope / 6.0
    D_target = KT / GAMMA
    rel_err = (D_measured - D_target) / D_target
    return {
        "dt": dt,
        "n_steps": N_STEPS_DIFFUSION,
        "D_measured": float(D_measured),
        "D_target_kT_over_gamma": float(D_target),
        "rel_err": float(rel_err),
        "msd_final": float(msd[-1]),
        "t_final": float(times[-1]),
    }


# ===========================================================================
# Gate harness
# ===========================================================================
def _print_gate(name: str, ok: bool, detail: str) -> None:
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}: {detail}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", action="store_true",
                        help="print wall-time per 1k steps")
    parser.add_argument("--n-steps", type=int, default=N_STEPS_PROD)
    parser.add_argument("--dt", type=float, default=DT_DEFAULT)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"BAOAB sanity — N_beads={N_BEADS}, dt={args.dt}, "
          f"n_steps={args.n_steps}, kT={KT}, γ={GAMMA}")
    print(f"HOOMD {hoomd.version.version} on {hoomd.device.CPU()}")
    print()

    all_pass = True
    report: dict = {}

    print("=== Equipartition gate (polymer harmonic bond, 3D-corrected) ===")
    poly = run_polymer(
        dt=args.dt, n_steps=args.n_steps, seed=args.seed, bench=args.bench
    )
    report["polymer"] = poly
    bond_ok = abs(poly["bond_rel_err_vs_3d_ref"]) <= 0.05
    _print_gate(
        "bond equipartition vs 3D analytical (±5%)",
        bond_ok,
        f"⟨½k(|r|-r0)²⟩/bond = {poly['per_bond_pe']:.5f} ± "
        f"{poly['per_bond_pe_stderr']:.5f}; 3D-analytical "
        f"{poly['target_bond_3d_analytical']:.5f}; rel_err="
        f"{poly['bond_rel_err_vs_3d_ref']:+.3%} "
        f"(naive ½kT={poly['target_angle_2d_brief_naive_half_kT']:.5f})",
    )

    # Bending mode is DIAGNOSTIC (not a gate) because the brief's
    # ½kT target is a 2D-equipartition reference (AFINES). 3D HOOMD
    # with sin(θ) volume element gives ≈0.93kT analytically; principled
    # 3D bending validation is the H.2 L_p measurement.
    print(
        f"  [DIAG] bending mode (NOT a gate; see PI escalation):\n"
        f"         measured  ⟨½k_θ(θ−π)²⟩/triplet = {poly['per_angle_pe']:.5f}"
        f" ± {poly['per_angle_pe_stderr']:.5f}\n"
        f"         brief target (2D ½kT):           "
        f"{poly['target_angle_2d_brief_naive_half_kT']:.5f}  → rel_err "
        f"{poly['angle_rel_err_vs_2d_brief']:+.1%}\n"
        f"         3D harmonic-on-sin(θ) ref:        "
        f"{poly['target_angle_3d_analytical']:.5f}  → rel_err "
        f"{poly['angle_rel_err_vs_3d_ref']:+.1%}"
    )
    all_pass &= bond_ok
    if args.bench:
        print(f"  wall-time: {poly['wall_time_s']:.2f} s for "
              f"{poly['n_steps']} steps "
              f"({poly['steps_per_s']:.0f} steps/s, "
              f"{poly.get('bench_per_1k_steps_s', 0):.3f} s/1k)")

    print()
    print("=== Diffusion gate (free bead, Stokes-Einstein) ===")
    diff = run_free_diffusion(dt=DT_DIFFUSION, seed=args.seed + 100)
    report["diffusion"] = diff
    diff_ok = abs(diff["rel_err"]) <= 0.10
    _print_gate(
        "diffusion D ≈ kT/γ (±10%)",
        diff_ok,
        f"D = {diff['D_measured']:.5f} (target {diff['D_target_kT_over_gamma']:.5f}); "
        f"rel_err={diff['rel_err']:+.3%}",
    )
    all_pass &= diff_ok

    print()
    print("=== Reproducibility gate (same seed → bit-identical trajectory) ===")
    repro = run_reproducibility(dt=args.dt, n_steps=5000, seed=args.seed + 500)
    report["reproducibility"] = repro
    repro_ok = repro["max_abs_diff"] == 0.0
    _print_gate(
        "same-seed reproducibility (max |Δr| == 0)",
        repro_ok,
        f"max |Δr| between two runs of {repro['n_steps']} steps at seed="
        f"{repro['seed']}: {repro['max_abs_diff']:e}",
    )
    all_pass &= repro_ok

    print()
    print("=" * 60)
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 60)

    report["overall_pass"] = all_pass
    out_path = OUT_DIR / "baoab_sanity_report.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Report: {out_path}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
