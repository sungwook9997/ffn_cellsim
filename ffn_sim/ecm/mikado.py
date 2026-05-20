"""H.1 ECM Mikado topology + HOOMD wiring.

Phase 1 Unit H.1 (PHASE_0_3_DECISIONS D3 + D4 + D7).

Builds a 3D-periodic HOOMD-blue simulation of a 2D Mikado fiber network
extruded to z=0 in a cubic 3D box (PI decision 2026-05-20: 3D-periodic +
2D projection, no z-pinning external field). Topology generation REUSES
``ffn_sim/validation/oracles/ecm/fiber_network.generate_2d_fiber_network``
for the geometry call only; the oracle's ``compute_energy`` /
``compute_forces`` are **not** imported here (CLAUDE.md hard rule: no v1
force kernel in v2 runtime).

Per-bead Hamiltonian wired into HOOMD:

    H = ½ k_bond Σ (|b_i| − ℓ₀)²              (md.bond.Harmonic, k=μ/ℓ₀)
      + ½ k_θ   Σ (θ_j − π)²                  (md.angle.Harmonic, t0=π)
      + Σ_{i<j} LJ_WCA(r_ij; ε, σ, r_cut)    (D7, repulsive only)

Mapping to oracle KU-1.24 Hamiltonian H_oracle:

    H_oracle = (μ/2ℓ₀) Σ (|b|−ℓ₀)² + (κ/ℓ₀) Σ (1 − cos θ_oracle)

where θ_oracle is the angle between consecutive bond vectors b_i, b_{i+1}
(0 at straight chain), while HOOMD uses θ_HOOMD between r_ji and r_jk at
central particle j (π at straight chain). The two forms are related by

    θ_oracle = π − θ_HOOMD,
    (1 − cos θ_oracle) = (1 − cos(π − θ_HOOMD)) = (1 + cos θ_HOOMD)
                       ≈ ½ (π − θ_HOOMD)²  as θ_HOOMD → π

so HOOMD ``k_θ = κ / ℓ₀`` matches the oracle in the small-bend limit
exactly, and the two energies agree to numerical precision at the
all-straight Mikado initial configuration (where both vanish). This
is the basis of the H.1 brief's energy-oracle gate at ≤ 1e-6 relative.

D3 integration is via the BAOAB-limit Action committed at
``ffn_sim/integrator/baoab.py`` (frozen 2026-05-20 by PI). The
``md.Integrator`` carries forces only (methods=[]); the Action advances
positions.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module." STATIC checks live in
``ffn_sim/tests/test_h1_mikado.py``; RUNTIME checks are enforced in
``resolve_derived``/``build_mikado_state``.*

1. **Dimensional analysis**

   - ``k_bond = μ / ℓ₀``  →  [N/m]. HOOMD ``md.bond.Harmonic.params['k']``
     is [E/length²] in HOOMD units; in SI, energy J = N·m so [k]
     = J/m² = N/m. ✓
   - ``k_θ = κ / ℓ₀``  →  [N·m² / m] = [N·m] = [J/rad²] since (θ−π)
     is in radians (dimensionless). HOOMD k is in [E/rad²] = [J]. ✓
   - ``γ_b = 6π η R``  →  [Pa·s · m] = [N·s/m]. ✓
   - ``ε_LJ = 0.5 · kT``  →  [J], σ_LJ = 2·R [m], r_cut = 2^(1/6)·σ.
     LJ U = 4ε[(σ/r)¹² − (σ/r)⁶], all dimensionless inside brackets. ✓
   - ``dt_CFL = α · min(τ_xl, τ_stretch, τ_bend)`` with τ in [s]. ✓
   - STATIC: ``test_h1_mikado.py::TestDimensional`` recomputes each
     derived value from the primary scales and asserts SI consistency.

2. **Boundary cases**

   - ``beads_per_fiber < 2``: forbidden (single-bead "fiber" has no
     bonds; oracle generator already raises). RUNTIME re-raised in
     ``resolve_derived``.
   - ``n_fibers == 0``: degenerate; the HOOMD state would have N=0
     particles. RUNTIME raises ``ValueError`` early.
   - ``L_z_over_L_box ≤ 0``: invalid box. RUNTIME raises.
   - ``target_segment_length`` outside ``acceptance.segment_length_range``
     OR ``expected_total_z`` outside ``acceptance.z_range``: surface
     to PI per CLAUDE.md "no gate-loosening" (raise unless
     ``demo_mode``).
   - ``n_fibers > acceptance.n_fibers_max``: cost ceiling violation;
     raise unless ``demo_mode``.

3. **Conservation invariants**

   - **At construction (all bonds at rest length, all chains straight)**:
     ``HOOMD potential_energy ≈ 0`` to numerical precision (each bond
     gives ½·k·0² = 0, each angle gives ½·k_θ·0² = 0), and the oracle's
     ``compute_energy`` also returns 0 on the same configuration.
     STATIC test ``TestEnergyOracle::test_initial_energy_is_zero``.
   - **Newton's 3rd law**: HOOMD's bond/angle force computes are
     conservative gradients by construction; the LJ pair list is
     symmetric. STATIC: ``net_force.sum()`` over all particles ≈ 0
     to float64 round-off in the initial state.
   - **Particle count = n_fibers · beads_per_fiber**, no spurious
     duplicates or omissions. STATIC.
   - **Bond count = n_fibers · (beads_per_fiber − 1)** exactly.
   - **Angle count = n_fibers · (beads_per_fiber − 2)** exactly.

4. **Numerical sanity**

   - ``dt_CFL`` must be > 0 and ≤ α·τ_min. RUNTIME assertion plus the
     KU-1.26 gate ``validation/oracles/common/sanity_gate.gate_unit1_2_dynamics``
     called when applicable.
   - All derived values finite. RUNTIME asserts ``np.isfinite``.
   - Position float64 (HOOMD CPU default). STATIC check.
   - Periodic wrap into [−L/2, L/2) (HOOMD convention). Oracle returns
     positions in [0, L); we shift by −L/2 before handing to GSD.
     RUNTIME assertion ``|x| ≤ L/2 + 1e-9`` per axis.
   - The BAOAB Action's own RUNTIME gate already catches NaN/Inf force
     and any `methods`-nonempty wiring error at attach.

5. **Sign / sense**

   - Stretched bond (|b| > ℓ₀): HOOMD harmonic force pulls bead i toward
     bead i+1 (attractive). Identical sign to oracle. STATIC test:
     stretch one bond by +δ, confirm HOOMD ``net_force`` on the two
     beads points along ±b̂ with magnitude k_bond·δ.
   - Bent triplet (θ_HOOMD < π): HOOMD harmonic-angle torque pushes the
     central bead to straighten the chain. STATIC: perturb one
     interior bead by +δ ŷ and confirm restoring force is along −ŷ.
   - LJ WCA (D7) is purely repulsive (cut at the LJ minimum); two beads
     at r < σ feel a force pushing them apart. STATIC: place two
     beads at r = σ/2, confirm both feel +force away from the pair midpoint.

6. **Measurement protocol**

   - **Topology smoke** (STATIC): generated arrays match the analytical
     counts in §3. Also checks the Mikado line-density relation
     ρ_L = N L_fiber / L_box² (KU-1.27).
   - **Energy oracle agreement** (STATIC): at the initial configuration,
     and at a small-perturbation configuration (Gaussian noise scale
     chosen so the small-bend approximation holds to ≤ 1e-6 rel),
     HOOMD ``ThermodynamicQuantities.potential_energy`` matches the
     oracle ``compute_energy`` (called on the same 2D-projected positions
     with the same parameters) to relative error ≤ 1e-6.
   - **Sanity-gate report**: ``resolve_derived`` returns a structured
     dict containing every derived value; the test suite asserts the
     full set of derivations and gate checks pass before any
     integration step runs.

References
----------
- Brief: ``ffn_sim/docs/briefs/H1_ecm_mikado.md`` §Topology, §Integration,
  §Excluded volume, §Validation acceptance.
- PHASE_0_3_DECISIONS.md §§D3, D4, D7.
- KU-1.1 (κ, ℓ_p), KU-1.2 (μ, R), KU-1.7 (ξ), KU-1.22 (L_fiber),
  KU-1.24 (oracle Hamiltonian), KU-1.26 (CFL τ_min), KU-1.27 (Mikado
  ρ_L → ℓ_c), KU-1.28 (xl_stiffness).
- ``ffn_sim/integrator/baoab.py`` (D3 frozen).
- ``ffn_sim/validation/oracles/ecm/fiber_network.py`` (geometry oracle —
  imported), ``ffn_sim/validation/oracles/ecm/fiber_mechanics.py`` (force
  oracle — **not** imported by this runtime module; test-only).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.validation.oracles.ecm.fiber_network import (
    FiberNetwork,
    generate_2d_fiber_network,
)


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedH1:
    """Fully-resolved H.1 parameter set (all SI).

    Populated by :func:`resolve_derived` from the raw YAML ``ecm:`` block.
    Every primary scale is preserved, and every ``derived`` value is
    computed from those primaries (no hand-tuning).
    """

    # Primary scales (verbatim from config)
    L_box: float
    L_fiber: float
    bead_radius: float
    persistence_length: float
    bending_modulus: float       # κ  [N·m²]
    stretching_modulus: float    # μ  [N]
    biological_mesh: float       # ξ  [m]
    target_segment_length: float # ℓ_c target  [m]
    xl_stiffness: float          # [N/m]
    temperature: float
    kT: float
    water_viscosity: float       # η  [Pa·s]
    seed: int
    beads_per_fiber: int
    dimensions: int
    L_z_over_L_box: float

    # D7 LJ
    lj_enabled: bool
    lj_epsilon_kT: float
    lj_sigma_factor: float

    # D3 dynamics
    integrator_name: str
    cfl_safety_factor: float

    # Acceptance
    segment_length_range: tuple[float, float]
    z_range: tuple[float, float]
    n_fibers_max: int
    demo_mode: bool

    # Derived
    L_z: float = 0.0
    rest_length: float = 0.0     # ℓ₀  [m]
    n_fibers: int = 0
    line_density: float = 0.0    # ρ_L  [1/m]
    mikado_ell_c_predicted: float = 0.0
    expected_xl_per_fiber: float = 0.0
    bond_k: float = 0.0          # μ / ℓ₀  [N/m]
    angle_k: float = 0.0         # κ / ℓ₀  [N·m / rad²]
    angle_t0: float = math.pi
    lj_epsilon: float = 0.0      # [J]
    lj_sigma: float = 0.0        # [m]
    lj_r_cut: float = 0.0        # [m]
    gamma_b: float = 0.0         # 6π η R  [N·s/m]
    tau_xl: float = 0.0
    tau_stretch: float = 0.0
    tau_bend: float = 0.0
    tau_min: float = 0.0
    dt_cfl: float = 0.0

    # Provenance
    extras: dict[str, Any] = field(default_factory=dict)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0, got {x!r}")


def resolve_derived(cfg: dict) -> ResolvedH1:
    """Resolve all derived parameters from a raw ``ecm:`` config block.

    The input ``cfg`` is the YAML ``ecm:`` sub-dict (as returned by
    ``yaml.safe_load``). The output is a :class:`ResolvedH1` with every
    ``derived:`` field populated. RUNTIME §1, §2, §4 gate checks are
    enforced here.
    """
    if "ecm" in cfg:
        cfg = cfg["ecm"]

    # ---- Primary ----
    p = ResolvedH1(
        L_box=float(cfg["L_box"]),
        L_fiber=float(cfg["L_fiber"]),
        bead_radius=float(cfg["bead_radius"]),
        persistence_length=float(cfg["persistence_length"]),
        bending_modulus=float(cfg["bending_modulus"]),
        stretching_modulus=float(cfg["stretching_modulus"]),
        biological_mesh=float(cfg["biological_mesh"]),
        target_segment_length=float(cfg["target_segment_length"]),
        xl_stiffness=float(cfg["xl_stiffness"]),
        temperature=float(cfg["temperature"]),
        kT=float(cfg["kT"]),
        water_viscosity=float(cfg["water_viscosity"]),
        seed=int(cfg["seed"]),
        beads_per_fiber=int(cfg["beads_per_fiber"]),
        dimensions=int(cfg["topology"]["dimensions"]),
        L_z_over_L_box=float(cfg["topology"]["L_z_over_L_box"]),
        lj_enabled=bool(cfg["excluded_volume"]["enabled"]),
        lj_epsilon_kT=float(cfg["excluded_volume"]["epsilon_kT"]),
        lj_sigma_factor=float(cfg["excluded_volume"]["sigma_factor"]),
        integrator_name=str(cfg["dynamics"]["integrator"]),
        cfl_safety_factor=float(cfg["dynamics"]["cfl_safety_factor"]),
        segment_length_range=tuple(cfg["acceptance"]["segment_length_range"]),
        z_range=tuple(cfg["acceptance"]["z_range"]),
        n_fibers_max=int(cfg["acceptance"]["n_fibers_max"]),
        demo_mode=bool(cfg.get("demo_mode", False)),
    )

    # ---- §2 boundary checks on primaries ----
    if p.beads_per_fiber < 2:
        raise ValueError(
            f"beads_per_fiber must be ≥ 2 (got {p.beads_per_fiber})"
        )
    if p.L_z_over_L_box <= 0:
        raise ValueError(
            f"topology.L_z_over_L_box must be > 0 (got {p.L_z_over_L_box})"
        )
    if p.dimensions != 3:
        raise ValueError(
            "topology.dimensions must be 3 (PI 2026-05-20: 3D-periodic + "
            f"2D projection); got {p.dimensions}."
        )
    if p.integrator_name not in ("leimkuhler_matthews_baoab", "lm_baoab"):
        raise ValueError(
            "dynamics.integrator must be the D3 canonical "
            f"'leimkuhler_matthews_baoab' or 'lm_baoab' (got {p.integrator_name})."
        )
    for name, x in [
        ("L_box", p.L_box), ("L_fiber", p.L_fiber),
        ("bead_radius", p.bead_radius),
        ("persistence_length", p.persistence_length),
        ("bending_modulus", p.bending_modulus),
        ("stretching_modulus", p.stretching_modulus),
        ("biological_mesh", p.biological_mesh),
        ("target_segment_length", p.target_segment_length),
        ("xl_stiffness", p.xl_stiffness),
        ("temperature", p.temperature), ("kT", p.kT),
        ("water_viscosity", p.water_viscosity),
        ("cfl_safety_factor", p.cfl_safety_factor),
    ]:
        _require_finite_positive(name, x)

    # ---- §1 derived: geometry + force constants ----
    p.L_z = p.L_box * p.L_z_over_L_box
    p.rest_length = p.L_fiber / (p.beads_per_fiber - 1)

    # Mikado segment-length inversion (KU-1.27):
    #   ℓ_c = π / (2 ρ_L),   ρ_L = N L_fiber / L_box²
    #   ⇒ N = π L_box² / (2 ℓ_c L_fiber)
    ell_c_target = p.target_segment_length
    n_fibers_real = math.pi * p.L_box**2 / (2.0 * ell_c_target * p.L_fiber)
    p.n_fibers = int(round(n_fibers_real))
    if p.n_fibers <= 0:
        raise ValueError(
            "Resolved n_fibers must be > 0; check L_box / L_fiber / "
            f"target_segment_length ({p.n_fibers_real=})."
        )
    p.line_density = p.n_fibers * p.L_fiber / (p.L_box**2)
    p.mikado_ell_c_predicted = math.pi / (2.0 * p.line_density)
    p.expected_xl_per_fiber = 2.0 * p.line_density * p.L_fiber / math.pi

    # ---- §1 derived: force constants ----
    p.bond_k = p.stretching_modulus / p.rest_length          # μ/ℓ₀
    p.angle_k = p.bending_modulus / p.rest_length            # κ/ℓ₀
    p.angle_t0 = math.pi                                     # HOOMD straight
    p.lj_sigma = p.lj_sigma_factor * p.bead_radius           # = 2 R
    p.lj_epsilon = p.lj_epsilon_kT * p.kT                    # 0.5 kT
    # WCA: r_cut at the LJ minimum 2^(1/6)·σ. Mathematical convention,
    # not a tunable — computed here at full double precision.
    p.lj_r_cut = (2.0 ** (1.0 / 6.0)) * p.lj_sigma
    p.gamma_b = 6.0 * math.pi * p.water_viscosity * p.bead_radius

    # ---- §1 derived: relaxation times + CFL ----
    p.tau_xl = p.gamma_b / p.xl_stiffness
    p.tau_stretch = p.gamma_b * p.rest_length / p.stretching_modulus
    p.tau_bend = p.gamma_b * p.rest_length**3 / p.bending_modulus
    p.tau_min = min(p.tau_xl, p.tau_stretch, p.tau_bend)
    p.dt_cfl = p.cfl_safety_factor * p.tau_min

    # ---- §2 acceptance checks (no gate-loosening) ----
    if not p.demo_mode:
        lo, hi = p.segment_length_range
        if not (lo <= p.mikado_ell_c_predicted <= hi):
            raise ValueError(
                f"Mikado ℓ_c predicted = {p.mikado_ell_c_predicted:.3e} m "
                f"outside KU-1.7 acceptance band [{lo:.3e}, {hi:.3e}]. "
                "Surface to PI (CLAUDE.md no-gate-loosening) or set demo_mode."
            )
        if p.n_fibers > p.n_fibers_max:
            raise ValueError(
                f"n_fibers = {p.n_fibers} exceeds cost ceiling "
                f"{p.n_fibers_max}; set demo_mode to override."
            )

    # ---- §4 finite checks on derived ----
    for name in [
        "L_z", "rest_length", "line_density", "mikado_ell_c_predicted",
        "expected_xl_per_fiber", "bond_k", "angle_k", "lj_epsilon",
        "lj_sigma", "lj_r_cut", "gamma_b", "tau_xl", "tau_stretch",
        "tau_bend", "tau_min", "dt_cfl",
    ]:
        v = getattr(p, name)
        if not (math.isfinite(v) and v > 0.0):
            raise RuntimeError(
                f"Derived {name}={v!r} is not finite-positive; "
                "config resolve produced a degenerate value."
            )

    return p


# ---------------------------------------------------------------------------
# Snapshot / Simulation builders
# ---------------------------------------------------------------------------
def build_mikado_state(
    p: ResolvedH1, *, with_cross_links: bool = True
) -> gsd.hoomd.Frame:
    """Build a HOOMD GSD frame from a resolved H.1 config.

    Calls the oracle's :func:`generate_2d_fiber_network` for the 2D
    Mikado geometry (REUSE only), pads to 3D with z=0, shifts positions
    to HOOMD's [-L/2, L/2) box convention, and assembles bonds + angles.

    Parameters
    ----------
    with_cross_links : bool, default True
        If True (M2 default), also call
        ``ffn_sim.ecm.cross_links.generate_xl_bonds`` and append the
        ``xl``-type bonds onto the frame. Set to False for M1
        topology-only tests.
    """
    if with_cross_links:
        # Generate xl bonds (which internally re-runs the geometry with
        # the same seed) and obtain the network it ran on, so the
        # subsequent ecm-bond frame uses the identical FiberNetwork.
        from ffn_sim.ecm.cross_links import add_xl_to_frame, generate_xl_bonds

        xl, net = generate_xl_bonds(p)
        snap = _network_to_frame(net, p)
        snap = add_xl_to_frame(snap, xl)
        return snap

    net = generate_2d_fiber_network(
        L_box=p.L_box,
        n_fibers=p.n_fibers,
        L_fiber=p.L_fiber,
        beads_per_fiber=p.beads_per_fiber,
        S_order=0.0,
        theta0=0.0,
        seed=p.seed,
        params={"resolved_from": "ffn_sim/configs/phase1_h1.yaml"},
    )
    return _network_to_frame(net, p)


def _network_to_frame(net: FiberNetwork, p: ResolvedH1) -> gsd.hoomd.Frame:
    F = p.n_fibers
    N = p.beads_per_fiber

    # (F, N, 2) → (F, N, 3) with z=0
    pos2d = net.bead_positions
    if pos2d.shape != (F, N, 2):
        raise RuntimeError(
            f"Oracle returned positions of shape {pos2d.shape}; "
            f"expected ({F}, {N}, 2). Did config drift?"
        )
    pos3d = np.zeros((F, N, 3), dtype=np.float64)
    # Oracle returns [0, L); HOOMD GSD expects [-L/2, L/2).
    pos3d[..., :2] = pos2d - 0.5 * p.L_box
    pos_flat = pos3d.reshape(F * N, 3)

    # §4 wrap sanity: all coords within [-L/2-eps, L/2+eps).
    L_half = 0.5 * p.L_box
    if (np.abs(pos_flat[:, 0]) > L_half + 1e-9).any() or \
       (np.abs(pos_flat[:, 1]) > L_half + 1e-9).any():
        raise RuntimeError(
            "Bead positions outside HOOMD box [-L/2, L/2) after shift; "
            "oracle wrap returned out-of-range coordinate."
        )

    n_part = F * N
    n_bonds = F * (N - 1)
    n_angles = F * (N - 2) if N >= 3 else 0

    # Bond groups: for each fiber f, bonds (f*N + i, f*N + i+1) for i in [0..N-2]
    f_idx = np.arange(F, dtype=np.int64).reshape(F, 1)
    i_idx = np.arange(N - 1, dtype=np.int64).reshape(1, N - 1)
    base = f_idx * N + i_idx                       # (F, N-1)
    bond_group = np.stack([base, base + 1], axis=-1).reshape(n_bonds, 2)

    if n_angles > 0:
        i_idx_a = np.arange(N - 2, dtype=np.int64).reshape(1, N - 2)
        base_a = f_idx * N + i_idx_a               # (F, N-2)
        angle_group = np.stack(
            [base_a, base_a + 1, base_a + 2], axis=-1
        ).reshape(n_angles, 3)
    else:
        angle_group = np.empty((0, 3), dtype=np.int64)

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_part
    snap.particles.types = ["actin_ecm"]
    snap.particles.typeid = np.zeros(n_part, dtype=np.uint32)
    snap.particles.position = pos_flat
    snap.particles.mass = np.ones(n_part, dtype=np.float64)

    snap.bonds.N = n_bonds
    snap.bonds.types = ["ecm-bond"]
    snap.bonds.typeid = np.zeros(n_bonds, dtype=np.uint32)
    snap.bonds.group = bond_group.astype(np.uint32)

    if n_angles > 0:
        snap.angles.N = n_angles
        snap.angles.types = ["ecm-angle"]
        snap.angles.typeid = np.zeros(n_angles, dtype=np.uint32)
        snap.angles.group = angle_group.astype(np.uint32)

    snap.configuration.box = [p.L_box, p.L_box, p.L_z, 0.0, 0.0, 0.0]
    return snap


def build_mikado_simulation(
    p: ResolvedH1,
    *,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
    with_cross_links: bool = True,
) -> tuple[hoomd.Simulation, Any, Any]:
    """Construct and wire a HOOMD Simulation for the H.1 Mikado.

    Parameters
    ----------
    with_baoab : bool, default True
        If True, attach the D3 BAOAB-limit Updater. If False, the
        Simulation contains forces only (bond + angle + LJ); no
        integration advance occurs on ``sim.run(n)`` for any n. The
        "False" mode is used by the energy-oracle test, which compares
        ``bond.energy + angle.energy`` against the oracle on a
        static (caller-perturbed) configuration. The first ``sim.run(0)``
        triggers an initial force evaluation; subsequent ``sim.run(1)``
        calls re-evaluate forces after the caller writes new positions
        via ``cpu_local_snapshot`` — ``run(0)`` after the first run
        does *not* trigger a fresh force eval (HOOMD 7 caches), so
        callers that mutate positions between reads must use
        ``sim.run(1)``.

    Returns
    -------
    sim : hoomd.Simulation
    updater : hoomd.update.CustomUpdater or None
        The BAOAB Updater wrapper if ``with_baoab=True``; else None.
    action : LeimkuhlerMatthewsBAOAB or None
        The Action instance if ``with_baoab=True``; else None.
    """
    snap = build_mikado_state(p, with_cross_links=with_cross_links)

    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p.seed
    )
    sim.create_state_from_snapshot(snap)

    # Bond force: U = ½ k_bond (|r| − r0)² for ecm-bond.
    # Cross-links (KU-1.28) use the same md.bond.Harmonic compute under
    # per-r0-bin types ("xl_b0" ... "xl_b{N-1}"), k=k_xl on every bin and
    # r0 = bin_center (PI 2026-05-20, option B — r0 type-binning so the
    # construction state is approximately force-free for cross-links and
    # production runs don't need a long xl-equilibration prelude).
    bond = md.bond.Harmonic()
    bond.params["ecm-bond"] = dict(k=p.bond_k, r0=p.rest_length)
    if with_cross_links:
        from ffn_sim.ecm.cross_links import (
            XL_N_BINS, xl_bin_rest_lengths, xl_bin_type_names,
        )
        bin_names = xl_bin_type_names()
        bin_r0 = xl_bin_rest_lengths()
        for i in range(XL_N_BINS):
            bond.params[bin_names[i]] = dict(
                k=p.xl_stiffness, r0=float(bin_r0[i])
            )

    # Angle force: U = ½ k_θ (θ − π)²  (small-bend match to oracle's
    # (κ/ℓ₀)(1 − cos θ_oracle); see module docstring.)
    angle = md.angle.Harmonic()
    angle.params["ecm-angle"] = dict(k=p.angle_k, t0=p.angle_t0)

    # D7 LJ WCA repulsive only.
    # We use md.nlist.Tree (BVH) instead of md.nlist.Cell here because the
    # ECM is extremely sparse: L_box ≈ 200 μm, LJ r_cut ≈ 112 nm so
    # L/r_cut ~ 1800. A uniform cell grid would pre-allocate ~5·10⁹
    # cells (the v2 trial run with md.nlist.Cell segfaulted on cell-list
    # allocation). BVH is O(N log N) and does not pre-allocate a cell
    # grid, so it handles the low-density Mikado without memory blowup.
    nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    lj.params[("actin_ecm", "actin_ecm")] = dict(
        epsilon=p.lj_epsilon, sigma=p.lj_sigma
    )
    lj.r_cut[("actin_ecm", "actin_ecm")] = p.lj_r_cut if p.lj_enabled else 0.0
    lj.mode = "shift"

    ig = md.Integrator(dt=p.dt_cfl)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    # methods=[] is REQUIRED for the L-M Action (it raises at attach
    # otherwise; BAOAB sign-off test confirms this contract).
    sim.operations.integrator = ig

    if with_baoab:
        action, updater = make_baoab_updater(
            kT=p.kT, gamma={"actin_ecm": p.gamma_b}, dt=p.dt_cfl, seed=p.seed
        )
        sim.operations.updaters.append(updater)
        return sim, updater, action

    return sim, None, None
