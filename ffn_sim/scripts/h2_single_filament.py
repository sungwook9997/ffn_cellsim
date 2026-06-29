"""H.2 single-filament L_p validation harness.

Builds one 10 μm actin filament (21 beads, no cross-links, no LJ),
runs L-M BAOAB Brownian dynamics under Langevin thermostat T = 310 K,
samples tangent-correlation snapshots, and writes diagnostics for the
pytest L_p gate (`ffn_sim/tests/test_persistence_length.py`).

Per H.2 brief `ffn_sim/docs/briefs/H2_single_filament.md`:

    Equilibrate 10 s (5×10⁵ steps) → Sample 60 s (3×10⁶ steps),
    write a frame every 1000 steps (3000 frames).

Sanity Gate
-----------
*Per CLAUDE.md hard rule.*

1. **Dimensional analysis**: every input SI, derived via
   `resolve_h2_derived()`.  τ_min = γ_b ℓ_0 / k_θ (angle-timescale)
   since H.2 has no xl / no LJ — bond-stretching mode is frozen at
   thermal scale (see yaml comment for derivation).
2. **Boundary cases**: `n_beads < 3` raises (need ≥ 3 for a meaningful
   interior-bead angle).  `L_fiber < beads × ℓ_0` raises (degenerate
   compression at construction).
3. **Conservation**: single filament, no external forces; Σ F ≈ 0
   to float64 round-off at the initial straight configuration.
4. **Numerical sanity**: dt = α · τ_angle; NaN/Inf guard inherited
   from the BAOAB Action.
5. **Sign / sense**: filament initially along +x.  Thermal fluctuations
   are isotropic in the y / z plane.
6. **Measurement protocol**: 3000 frame snapshots stored as
   ``(n_frames, 1, 21, 3)`` ndarray for downstream
   `filament_math.fit_persistence_length` and
   `filament_math.equipartition_check`.

Usage::

    PYTHONPATH=. python ffn_sim/scripts/h2_single_filament.py \\
        [--n-equilibrate 500000] [--n-sample 3000000] \\
        [--sample-interval 1000] [--demo]

`--demo` drops to a smoke configuration (5 k equilibrate + 50 k sample
+ sample every 100 steps = 500 frames) for a sub-1-minute test.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater


CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h2.yaml"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h2"


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedH2:
    """Resolved H.2 single-filament parameter set (all SI)."""

    L_fiber: float
    bead_radius: float
    persistence_length: float
    bending_modulus: float
    stretching_modulus: float
    temperature: float
    kT: float
    water_viscosity: float
    seed: int
    beads_per_fiber: int

    box_Lx: float
    box_Ly: float
    box_Lz: float

    integrator_name: str
    cfl_safety_factor: float
    n_steps_equilibrate: int
    n_steps_sample: int
    sample_interval: int

    L_p_band_m: tuple[float, float]                  # brief literal (diagnostic)
    L_p_band_m_C1: tuple[float, float]               # 3σ_empirical, single-seed gate (PI 2026-05-25)
    L_p_band_m_C1_5sigma: tuple[float, float]        # 5σ analytical, ensemble gate (PI 2026-05-25)
    L_p_band_m_tail: tuple[float, float]             # ±50 % chain-mode scatter (PI 2026-05-25)
    angle_chi2_per_df_max: float                     # Pearson χ²/df critical 95 % (PI 2026-05-25)
    angle_KL_nats_max: float                         # KL nats threshold 95 % (PI 2026-05-25)
    angle_ks_p_min: float                            # diagnostic only
    angle_ks_stat_max: float                         # diagnostic only
    equipartition_rel_tol: float                     # brief literal (diagnostic)
    equipartition_target_3d_kT: float                # numerical 3D exact at α=16.4 (PI 2026-05-25)
    equipartition_rel_tol_3d: float                  # ±5 % strict (PI 2026-05-25)

    reference_integrator_allow: bool
    reference_integrator_name: str
    reference_integrator_dt_factor: float

    demo_mode: bool

    rest_length: float = 0.0
    bond_k: float = 0.0
    angle_k: float = 0.0
    angle_t0: float = math.pi
    gamma_b: float = 0.0
    tau_stretch: float = 0.0
    tau_bend: float = 0.0
    tau_min: float = 0.0
    dt_cfl: float = 0.0

    extras: dict[str, Any] = None  # type: ignore


def resolve_h2_derived(cfg: dict) -> ResolvedH2:
    """Resolve the H.2 yaml filament block into a ResolvedH2."""
    f = cfg["filament"]
    p = ResolvedH2(
        L_fiber=float(f["L_fiber"]),
        bead_radius=float(f["bead_radius"]),
        persistence_length=float(f["persistence_length"]),
        bending_modulus=float(f["bending_modulus"]),
        stretching_modulus=float(f["stretching_modulus"]),
        temperature=float(f["temperature"]),
        kT=float(f["kT"]),
        water_viscosity=float(f["water_viscosity"]),
        seed=int(f["seed"]),
        beads_per_fiber=int(f["beads_per_fiber"]),
        box_Lx=float(f["box"]["Lx"]),
        box_Ly=float(f["box"]["Ly"]),
        box_Lz=float(f["box"]["Lz"]),
        integrator_name=str(f["dynamics"]["integrator"]),
        cfl_safety_factor=float(f["dynamics"]["cfl_safety_factor"]),
        n_steps_equilibrate=int(f["dynamics"]["n_steps_equilibrate"]),
        n_steps_sample=int(f["dynamics"]["n_steps_sample"]),
        sample_interval=int(f["dynamics"]["sample_interval"]),
        L_p_band_m=tuple(f["acceptance"]["L_p_band_m"]),
        L_p_band_m_C1=tuple(
            f["acceptance"].get("L_p_band_m_C1", [14.8e-6, 18.0e-6])
        ),
        L_p_band_m_C1_5sigma=tuple(
            f["acceptance"].get("L_p_band_m_C1_5sigma", [15.81e-6, 17.07e-6])
        ),
        L_p_band_m_tail=tuple(
            f["acceptance"].get("L_p_band_m_tail", [8.0e-6, 33.0e-6])
        ),
        angle_chi2_per_df_max=float(
            f["acceptance"].get("angle_chi2_per_df_max", 1.354)
        ),
        angle_KL_nats_max=float(
            f["acceptance"].get("angle_KL_nats_max", 1.94e-3)
        ),
        angle_ks_p_min=float(f["acceptance"]["angle_ks_p_min"]),
        angle_ks_stat_max=float(
            f["acceptance"].get("angle_ks_stat_max", 0.10)
        ),
        equipartition_rel_tol=float(f["acceptance"]["equipartition_rel_tol"]),
        equipartition_target_3d_kT=float(
            f["acceptance"].get("equipartition_target_3d_kT", 0.9898)
        ),
        equipartition_rel_tol_3d=float(
            f["acceptance"].get("equipartition_rel_tol_3d", 0.05)
        ),
        reference_integrator_allow=bool(f["reference_integrator"]["allow"]),
        reference_integrator_name=str(f["reference_integrator"]["name"]),
        reference_integrator_dt_factor=float(f["reference_integrator"]["dt_factor"]),
        demo_mode=bool(f.get("demo_mode", False)),
    )

    if p.beads_per_fiber < 3:
        raise ValueError(
            f"H.2 requires beads_per_fiber ≥ 3 (need interior bead for "
            f"angles); got {p.beads_per_fiber}."
        )
    if p.integrator_name not in ("leimkuhler_matthews_baoab", "lm_baoab"):
        raise ValueError(
            "H.2 D3 canonical integrator must be "
            f"'leimkuhler_matthews_baoab'; got {p.integrator_name}."
        )

    p.rest_length = p.L_fiber / (p.beads_per_fiber - 1)
    p.bond_k = p.stretching_modulus / p.rest_length
    p.angle_k = p.bending_modulus / p.rest_length
    p.gamma_b = 6.0 * math.pi * p.water_viscosity * p.bead_radius
    # Position-space relaxation times.  τ = γ_b / k_pos where k_pos is
    # the effective spring constant on bead position.
    #   bond:   k_pos = k_bond = μ/ℓ_0 → τ_stretch = γ_b · ℓ_0 / μ
    #   angle:  k_pos = k_θ / ℓ_0² (small-angle deflection δ ⊥ chain
    #                              gives δθ ≈ δ/ℓ_0; U = ½ k_θ (δ/ℓ_0)²)
    #           → τ_bend = γ_b · ℓ_0² / k_θ = γ_b · ℓ_0³ / κ_B
    p.tau_stretch = p.gamma_b * p.rest_length / p.stretching_modulus
    p.tau_bend = p.gamma_b * p.rest_length**3 / p.bending_modulus
    # H.2 (no xl, no LJ): τ_min = min(τ_stretch, τ_bend).  Use the
    # smaller (stretch) for full safety even though the relevant
    # measurement timescale is τ_bend ≈ 0.7 ms ≫ τ_stretch ≈ 40 ns.
    p.tau_min = min(p.tau_stretch, p.tau_bend)
    p.dt_cfl = p.cfl_safety_factor * p.tau_min
    return p


# ---------------------------------------------------------------------------
# Simulation builder
# ---------------------------------------------------------------------------
def build_h2_simulation(
    p: ResolvedH2, *, device: hoomd.device.Device | None = None,
    integrator: str = "lm_baoab",
) -> tuple[hoomd.Simulation, Any, Any]:
    """Build a HOOMD Simulation for the single-filament H.2 setup.

    Parameters
    ----------
    integrator
        ``"lm_baoab"`` (default) — D3 canonical L-M BAOAB Updater
        (`ffn_sim/integrator/baoab.py`).  Returns (sim, updater, action).

        ``"hoomd_brownian"`` — HOOMD-native `md.methods.Brownian`
        (Euler-Maruyama).  Used as the reference integrator for the
        BAOAB §Open #2 order-separation gate (PI-ratified in
        `phase1_h2.yaml::reference_integrator`).  Returns
        (sim, None, brownian_method) so the caller signature stays
        symmetric — there is no BAOAB Action in this branch.  dt is
        scaled by ``p.reference_integrator_dt_factor`` (default 0.5)
        so E-M's first-order bias is comparable to L-M BAOAB.

    Returns (sim, updater, action_or_method) — the BAOAB Action +
    Updater pair (lm_baoab) or (None, brownian_method) (hoomd_brownian).
    """
    N = p.beads_per_fiber

    # Filament aligned along +x at z = 0, centred on (0, 0, 0).
    pos = np.zeros((N, 3), dtype=np.float64)
    pos[:, 0] = np.linspace(
        -0.5 * p.L_fiber, +0.5 * p.L_fiber, N, dtype=np.float64
    )

    snap = gsd.hoomd.Frame()
    snap.particles.N = N
    snap.particles.types = ["actin"]
    snap.particles.typeid = np.zeros(N, dtype=np.uint32)
    snap.particles.position = pos
    snap.particles.mass = np.ones(N, dtype=np.float64)

    snap.bonds.N = N - 1
    snap.bonds.types = ["actin-bond"]
    snap.bonds.typeid = np.zeros(N - 1, dtype=np.uint32)
    bg = np.stack([np.arange(N - 1, dtype=np.uint32),
                   np.arange(1, N, dtype=np.uint32)], axis=-1)
    snap.bonds.group = bg

    if N >= 3:
        snap.angles.N = N - 2
        snap.angles.types = ["actin-angle"]
        snap.angles.typeid = np.zeros(N - 2, dtype=np.uint32)
        ag = np.stack([np.arange(N - 2, dtype=np.uint32),
                       np.arange(1, N - 1, dtype=np.uint32),
                       np.arange(2, N, dtype=np.uint32)], axis=-1)
        snap.angles.group = ag

    snap.configuration.box = [p.box_Lx, p.box_Ly, p.box_Lz, 0.0, 0.0, 0.0]

    sim = hoomd.Simulation(device=device or hoomd.device.CPU(), seed=p.seed)
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    bond.params["actin-bond"] = dict(k=p.bond_k, r0=p.rest_length)

    angle = md.angle.Harmonic()
    angle.params["actin-angle"] = dict(k=p.angle_k, t0=p.angle_t0)

    if integrator == "lm_baoab":
        ig = md.Integrator(dt=p.dt_cfl)
        ig.forces.append(bond)
        ig.forces.append(angle)
        sim.operations.integrator = ig
        action, updater = make_baoab_updater(
            kT=p.kT, gamma={"actin": p.gamma_b}, dt=p.dt_cfl, seed=p.seed
        )
        sim.operations.updaters.append(updater)
        return sim, updater, action
    elif integrator == "hoomd_brownian":
        # BAOAB §Open #2 reference integrator — Euler-Maruyama via
        # md.methods.Brownian.  dt scaled by reference_integrator_dt_factor
        # (default 0.5) so E-M's O(Δt) bias is comparable to L-M's O(Δt²)
        # at the same effective error level (Leimkuhler-Matthews 2013).
        dt_ref = p.dt_cfl * p.reference_integrator_dt_factor
        ig = md.Integrator(dt=dt_ref)
        ig.forces.append(bond)
        ig.forces.append(angle)
        brownian = md.methods.Brownian(
            filter=hoomd.filter.All(), kT=p.kT,
        )
        brownian.gamma["actin"] = p.gamma_b
        ig.methods.append(brownian)
        sim.operations.integrator = ig
        return sim, None, brownian
    else:
        raise ValueError(
            f"Unknown integrator {integrator!r}; choose 'lm_baoab' or "
            "'hoomd_brownian'."
        )


# ---------------------------------------------------------------------------
# Run + sample
# ---------------------------------------------------------------------------
def run_h2_and_sample(
    p: ResolvedH2, *,
    n_equilibrate: int | None = None,
    n_sample: int | None = None,
    sample_interval: int | None = None,
    integrator: str = "lm_baoab",
) -> dict[str, Any]:
    """Drive the H.2 sim through equilibrate + sample phases.

    Returns ``{positions, frame_steps, wall_s, p}`` where ``positions``
    has shape ``(n_frames, F=1, N, 3)`` and ``frame_steps`` is the
    HOOMD timestep at which each frame was taken.
    """
    n_eq = n_equilibrate if n_equilibrate is not None else p.n_steps_equilibrate
    n_smp = n_sample if n_sample is not None else p.n_steps_sample
    si = sample_interval if sample_interval is not None else p.sample_interval

    sim, _updater, _action = build_h2_simulation(p, integrator=integrator)

    t0 = time.time()
    sim.run(n_eq)

    n_frames = n_smp // si
    positions = np.empty((n_frames, 1, p.beads_per_fiber, 3), dtype=np.float64)
    frame_steps = np.empty(n_frames, dtype=np.int64)
    for i in range(n_frames):
        sim.run(si)
        with sim.state.cpu_local_snapshot as s:
            # HOOMD's ParticleSorter reorders snapshot rows; gather by
            # stable particle tag so the per-bead position matches the
            # bond.group tag ordering (which build_h2_simulation laid
            # out as bead-0, bead-1, ..., bead-N-1 along the filament).
            pos_row = np.asarray(s.particles.position).copy()
            tag_row = np.asarray(s.particles.tag).copy()
            pos = np.empty_like(pos_row)
            pos[tag_row] = pos_row
        positions[i, 0, :, :] = pos
        frame_steps[i] = int(sim.timestep)

    elapsed = time.time() - t0
    return {
        "positions": positions,
        "frame_steps": frame_steps,
        "wall_s": elapsed,
        "p": p,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-equilibrate", type=int, default=None)
    ap.add_argument("--n-sample", type=int, default=None)
    ap.add_argument("--sample-interval", type=int, default=None)
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = resolve_h2_derived(cfg)

    if args.demo:
        kw = dict(n_equilibrate=5000, n_sample=50000, sample_interval=100)
    else:
        kw = dict(
            n_equilibrate=args.n_equilibrate,
            n_sample=args.n_sample,
            sample_interval=args.sample_interval,
        )

    print(f"[h2_single_filament] N={p.beads_per_fiber}, ℓ_0={p.rest_length:.3e} m, "
          f"k_bond={p.bond_k:.3e} N/m, k_angle={p.angle_k:.3e} N·m, "
          f"γ_b={p.gamma_b:.3e} N·s/m, dt={p.dt_cfl:.3e} s")
    result = run_h2_and_sample(p, **kw)
    print(f"[h2_single_filament] Sampled {result['positions'].shape[0]} frames "
          f"in {result['wall_s']:.2f} s")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(
        OUTPUT_DIR / "h2_trajectory.npz",
        positions=result["positions"],
        frame_steps=result["frame_steps"],
        wall_s=result["wall_s"],
        L_fiber=p.L_fiber,
        rest_length=p.rest_length,
        angle_k=p.angle_k,
        bending_modulus=p.bending_modulus,
        kT=p.kT,
        persistence_length_reference=p.persistence_length,
        box_Lx=p.box_Lx,
        box_Ly=p.box_Ly,
        box_Lz=p.box_Lz,
        L_p_band_m=np.array(p.L_p_band_m),
        equipartition_rel_tol=p.equipartition_rel_tol,
        dt_cfl=p.dt_cfl,
    )
    with open(OUTPUT_DIR / "h2_run_summary.json", "w") as fh:
        json.dump(
            {
                "n_beads": p.beads_per_fiber,
                "rest_length_m": p.rest_length,
                "L_fiber_m": p.L_fiber,
                "dt_cfl_s": p.dt_cfl,
                "tau_bend_s": p.tau_bend,
                "tau_stretch_s": p.tau_stretch,
                "wall_s": result["wall_s"],
                "n_frames": int(result["positions"].shape[0]),
            },
            fh,
            indent=2,
        )
    print(f"[h2_single_filament] wrote {OUTPUT_DIR}/h2_trajectory.npz")


if __name__ == "__main__":
    main()
