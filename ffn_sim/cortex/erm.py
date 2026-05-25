"""H.3 ERM tether — radial harmonic spring per cortex bead (KU-3.18).

Phase 1 H.3 brief §ERM tether. For each cortex bead i with position
``r_i``, applies a radial harmonic force toward the cell center
(origin) anchored at radius ``R_cell``:

    U_i = ½ k_ERM (|r_i| − R_cell)²
    F_i = −k_ERM (|r_i| − R_cell) · r̂_i

Implemented as a ``hoomd.md.force.Custom`` subclass so the force
contribution participates in HOOMD's ``net_force`` accumulator (which
the L-M BAOAB Updater reads each step).

Biology vs implementation note
------------------------------
The brief specifies ``k_ERM = 0.1 N/m`` per cortex bead. With
``kT = 4.28 · 10⁻²¹ J``, the per-bead radial thermal extent is
``σ_radial = √(kT/k_ERM) = 0.65 nm`` — three orders of magnitude
smaller than the 200 nm cortex thickness from KU-3.17. The brief value
therefore models the ERM tether as a TIGHT pinning at ``R_cell``, with
the biological 200 nm cortex thickness reflecting the membrane's
structural depth (cortex bundle radius + lipid bilayer), not the
thermal extent of single coarse-grained beads.

This is the brief's literal KU-3.18 anchor; we use it verbatim and
document the implication for downstream validation gates: cortex
beads do NOT thermally explore the 200 nm shell band — instead they
sit on the sphere ``r = R_cell`` to sub-nm tolerance, and rounding /
tension gates probe the EMERGENT large-scale curvature of the bead
ensemble, not the per-bead radial distribution.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_erm.py``.*

1. **Dimensional analysis**
   - ``k_ERM`` [N/m]; ``(|r| − R_cell)`` [m]; force ``= k_ERM · Δr``
     [N]. ✓
   - ``U = ½ k_ERM Δr²`` [J]. ✓
   - At construction (beads on sphere ``|r| = R_cell``), ``Δr = 0`` →
     ``F = 0`` and ``U = 0``. Force-free at construction.

2. **Boundary cases**
   - ``r = 0`` (bead at cell center): ``r̂`` is undefined; we set
     ``F = 0`` (per-bead) at that singular point. In practice cortex
     beads start at ``|r| = R_cell`` and the ``k_ERM = 0.1 N/m``
     pinning keeps them out of the singular region by many σ.
   - ``k_ERM ≤ 0``: caught by ``resolve_erm`` ValueError.
   - ``R_cell ≤ 0``: caught at resolve time.
   - Empty cortex (n_cortex_actin = 0): the custom force iterates zero
     particles; trivially correct.

3. **Conservation invariants**
   - Net force from ERM on the system: NOT zero (ERM is an EXTERNAL
     field anchored at the origin, breaking translational symmetry —
     this is intentional, modeling the cell-frame inertial reference).
   - Total ERM potential energy: ``Σ ½ k_ERM (|r_i| − R_cell)²`` —
     accumulated additively across cortex beads, monotonically
     non-negative.
   - The ERM force has NO contribution to xlink_head or any other
     non-actin_cortex particle (mask by tag range).

4. **Numerical sanity**
   - All positions read in float64.
   - Force compute matches the SI analytic expression to float64
     precision (STATIC test: stretch a single bead radially by δ,
     read net_force, confirm magnitude = ``k_ERM · δ`` to ≤ 1e-9
     relative).
   - Energy compute matches ``½ k_ERM · δ²``.
   - The custom force participates in HOOMD's net_force; the BAOAB
     Updater reads net_force per step (STATIC test: confirm
     net_force[cortex_bead] has ERM contribution after a sim.run(1)).

5. **Sign / sense**
   - Bead with ``|r| > R_cell`` (outside shell): force ``∝ −r̂`` —
     pulls INWARD toward cell center. STATIC test.
   - Bead with ``|r| < R_cell`` (inside shell): force ``∝ +r̂`` —
     pushes OUTWARD. STATIC test.
   - Bead exactly on shell (``|r| = R_cell``): F = 0.

6. **Measurement protocol**
   - Per-bead radial position at equilibrium is centered at R_cell with
     std-dev ``σ_radial = √(kT/k_ERM)``. STATIC test runs a few BAOAB
     steps and confirms ``|⟨r⟩ − R_cell| < 5 σ_radial`` (loose smoke
     gate; tight ±1 % gate is opt-in via H3_ERM_PRODUCTION=1).

References
----------
- Brief: ``ffn_sim/docs/briefs/H3_cortex.md`` §ERM tether.
- KU-3.18 (k_ERM = 0.1 N/m anchor — see brief; v1 yaml at
  ``validation/oracles/configs/phase1_unit3.yaml`` cites KU-3.18 for
  cortex active tension but not ERM directly; brief value adopted).
- HOOMD 7 ``md.force.Custom`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd
import hoomd.md as md


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedERM:
    """ERM tether parameters."""

    k_ERM: float                 # N/m  radial spring stiffness, KU-3.18 (0.1)
    R_cell: float                # m    anchor radius (matches cortex R_cell)
    cell_center: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Derived
    sigma_radial_thermal: float = 0.0   # √(kT/k_ERM)  [m]


def resolve_erm(cfg: dict, *, kT: float, R_cell: float) -> ResolvedERM:
    """Resolve ERM config block.

    ``cfg`` is the ``cortex.erm`` sub-dict OR the raw cortex config.
    ``kT`` and ``R_cell`` come from the host cortex resolved params.
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "erm" in cfg:
        cfg = cfg["erm"]

    p = ResolvedERM(
        k_ERM=float(cfg["k_ERM"]),
        R_cell=R_cell,
        cell_center=tuple(cfg.get("cell_center", [0.0, 0.0, 0.0])),
    )

    # §2 boundary checks
    if not (math.isfinite(p.k_ERM) and p.k_ERM > 0.0):
        raise ValueError(f"k_ERM must be finite and > 0; got {p.k_ERM}")
    if not (math.isfinite(p.R_cell) and p.R_cell > 0.0):
        raise ValueError(f"R_cell must be finite and > 0; got {p.R_cell}")

    # §1 derived
    p.sigma_radial_thermal = math.sqrt(kT / p.k_ERM)

    return p


# ---------------------------------------------------------------------------
# Custom force compute
# ---------------------------------------------------------------------------
class ERMHarmonic(md.force.Custom):
    """Per-bead radial harmonic spring to cell center (KU-3.18).

    Applies a radial spring force ``F = −k_ERM (|r| − R_cell) · r̂`` to
    every particle whose tag is in ``actin_cortex_tag_range`` (typically
    [0, n_cortex_actin)). Other particles (xlink_head, etc.) are left
    untouched.

    Parameters
    ----------
    p : ResolvedERM
        Resolved ERM config.
    actin_cortex_tag_range : tuple[int, int]
        Tag range [start, end) of cortex actin beads. Particles whose
        tag is in this range receive the ERM force.

    Notes
    -----
    HOOMD ``md.force.Custom`` participates in the integrator's
    ``net_force`` accumulator. The L-M BAOAB Updater reads net_force
    per step, so this force is felt by the dynamics.

    The custom force is computed in float64 from positions, applied to
    forces and energies via ``cpu_local_force_arrays``. The cell_center
    is fixed to the origin (per construction of the cortex topology in
    ``cortex.py``: filaments are placed on the sphere ``|r| = R_cell``
    at the box origin).
    """

    def __init__(
        self,
        p: ResolvedERM,
        actin_cortex_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(actin_cortex_tag_range[0])
        self.tag_end = int(actin_cortex_tag_range[1])
        if self.tag_end < self.tag_start:
            raise ValueError(
                f"actin_cortex_tag_range must satisfy end ≥ start; "
                f"got ({self.tag_start}, {self.tag_end})"
            )

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        # Read positions + tags from the read-only state snapshot,
        # compute forces analytically, then write into the per-row
        # force arrays. Both contexts use ROW indexing (not tag
        # indexing), and HOOMD guarantees row i in cpu_local_snapshot
        # is the same particle as row i in cpu_local_force_arrays on a
        # single rank.
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()

        # Mask: tags in [tag_start, tag_end).
        mask = (tag >= self.tag_start) & (tag < self.tag_end)
        # Center offset (typically zero).
        cx, cy, cz = self.p.cell_center
        dx = pos - np.array([cx, cy, cz], dtype=np.float64)
        r = np.linalg.norm(dx, axis=1)
        # Avoid division by zero at the singular center.
        r_safe = np.where(r > 0.0, r, 1.0)
        r_hat = dx / r_safe[:, None]
        F_radial_mag = -self.p.k_ERM * (r - self.p.R_cell)
        F_vec = F_radial_mag[:, None] * r_hat
        # Energy per particle: ½ k_ERM (|r| − R_cell)²
        U_per = 0.5 * self.p.k_ERM * (r - self.p.R_cell) ** 2
        # Zero out non-cortex beads.
        F_vec[~mask] = 0.0
        U_per[~mask] = 0.0

        with self.cpu_local_force_arrays as arrays:
            # Replace the force / energy arrays entirely (HOOMD
            # accumulates these into net_force / net_energy for
            # this force compute).
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per

    @staticmethod
    def predicted_sigma_radial(k_ERM: float, kT: float) -> float:
        """Analytic per-bead radial std-dev at thermal equilibrium."""
        return math.sqrt(kT / k_ERM)


# ---------------------------------------------------------------------------
# Public helper: attach ERM to an existing cortex simulation
# ---------------------------------------------------------------------------
def attach_erm_to_simulation(
    sim: hoomd.Simulation,
    p_erm: ResolvedERM,
    *,
    actin_cortex_tag_range: tuple[int, int],
    gamma_b: float | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> ERMHarmonic:
    """Append an ERMHarmonic custom force to an existing Integrator.

    Parameters
    ----------
    sim : hoomd.Simulation
        Already-built cortex simulation (must have an Integrator).
    p_erm : ResolvedERM
    actin_cortex_tag_range : tuple[int, int]
        Tag range of cortex actin beads.
    gamma_b : float or None
        Per-bead Stokes drag [N·s/m]. Required to gate the ERM CFL
        ``dt ≤ cfl_safety_factor · γ_b / k_ERM``. If None, the CFL
        gate is skipped (caller responsible).
    cfl_safety_factor : float, default 0.1
        Same convention as bond/angle CFL gates (D3 BAOAB).
    cfl_strict : bool, default True
        If True, raise on CFL violation. Set False for diagnostic-only
        gating (the caller knows what they're doing).

    Returns
    -------
    erm_force : ERMHarmonic
        The attached force compute (kept for introspection in tests).
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching ERM."
        )
    if gamma_b is not None:
        dt = float(ig.dt)
        tau_erm = gamma_b / p_erm.k_ERM
        dt_cfl_erm = cfl_safety_factor * tau_erm
        if dt > dt_cfl_erm and cfl_strict:
            raise RuntimeError(
                f"ERM CFL violated: dt = {dt:.3e} s > "
                f"{cfl_safety_factor:.2f} · τ_ERM = {dt_cfl_erm:.3e} s "
                f"(τ_ERM = γ_b / k_ERM = {tau_erm:.3e} s). "
                "Either reduce k_ERM (softer ERM, changes physics — "
                "requires PI sign-off vs brief literal KU-3.18 = 0.1 "
                "N/m) OR reduce dt (slower simulation; need dt ≤ "
                f"{dt_cfl_erm:.3e} s — {(dt / dt_cfl_erm):.1f}× smaller "
                "than current). Pass cfl_strict=False to skip this gate "
                "for diagnostic runs."
            )
    erm = ERMHarmonic(p_erm, actin_cortex_tag_range)
    ig.forces.append(erm)
    return erm
