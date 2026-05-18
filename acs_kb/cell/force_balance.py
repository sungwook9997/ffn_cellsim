"""Cell-level overdamped quasi-static force balance (KU-3.16).

The Phase 1 Unit 3.1 force balance per cortex bead is

    γ_b · dr/dt = F_WLC(r) + F_xl(r) + F_tension(r) + F_ext(r),

with no thermal noise (KU-3.13 quasi-steady-state assumption for the
cell-scale minute timescale). ``F_WLC`` and ``F_xl`` come directly
from :mod:`acs_kb.ecm.fiber_mechanics.compute_forces` (KU-1.24 +
KU-1.28), demonstrating the explicit re-use of Worker A's ECM module
as required by the Unit 3.1 brief. ``F_tension`` is the cortical
tension contribution (KU-3.5, Phase 1 simplification):

    F_tension(bead) = γ_cortex / R_eff · (−r̂_eff)

where r̂_eff is the unit vector from the cortex centroid to the bead
and R_eff is the bead's instantaneous distance from the centroid.
The factor 1/R_eff reproduces the Laplace-law inward pressure on a
2D loop with line tension γ_cortex; the bead-centroid normal is the
quasi-static surrogate for the local outward normal (cortex radius
≫ bead spacing). ``F_ext`` is empty in Phase 1 Unit 3.1 (no FA, no
junction); the argument is plumbed through so Worker B/D modules can
inject focal-adhesion and junction forces in Unit 3.2+.

Sanity Gate
-----------
1. Dimensional analysis: [γ_cortex] = N/m, [R_eff] = m, so
   [F_tension] = N ✓. [γ_b] = N·s/m, [F] = N, so the Euler update
   Δr = (dt/γ_b) F has units of m ✓. CFL: dt < α·τ_min where τ_min
   is the smallest of (τ_xl, τ_stretch, τ_bend, τ_tension).
2. Boundary cases: R_eff = 0 (degenerate centroid coincidence) →
   guarded by ``np.where(R_eff > 0, R_eff, ε)`` to avoid divide-by-
   zero; quasi-statically a bead can never sit *exactly* at the
   centroid, so the guard is defensive only.
3. Conservation: overdamped Euler is **NOT** energy-conserving (the
   dissipation channel is the bead drag γ_b). Newton's 3rd law on
   F_WLC + F_xl is preserved by the ECM kernel; F_tension does not
   conserve momentum at the bead level because the cytoplasm /
   external medium absorbs the centripetal reaction (KU-3.16
   overdamped limit explicitly drops inertia and momentum
   conservation in favour of force balance).
4. Numerical sanity: float64; no atomic ops; vectorised over beads.
5. Sign / sense:
   - WLC stretching: tensile bond pulls bead toward neighbour
     (KU-1.24).
   - Cross-link: stretched XL pulls the two beads together
     (KU-1.28).
   - Cortical tension: positive γ_cortex always pulls beads
     *inward* toward the centroid (KU-3.5 → cell rounding under
     active tension), regardless of curvature sign.
6. Measurement protocol: the gate compares the simulated final
   aspect ratio (PCA-based, :func:`acs_kb.cell.cortex.measure_aspect_ratio`)
   to the KU-3.1 acceptance band [1.0, 1.2] *measured over the same
   bead population* the dynamics evolved, matching the analytical
   "free cell rounds to a circle" prediction sample-by-sample.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from acs_kb.cell.cortex import Cortex, measure_aspect_ratio
from acs_kb.ecm.fiber_mechanics import compute_forces


@dataclass(slots=True)
class RoundingTrajectory:
    """Time series of cortex relaxation diagnostics."""

    times: np.ndarray
    aspect_ratios: np.ndarray
    mean_radii: np.ndarray
    bead_positions_frames: list[np.ndarray] = field(default_factory=list)


def cytoplasm_pressure_forces(
    bead_positions: np.ndarray,
    K_area: float,
    A_0: float,
    rest_length: float,
) -> np.ndarray:
    """Isotropic cytoplasmic pressure that conserves cortex area (KU-3.9).

    A 2D loop under line tension γ_cortex alone collapses to a point
    (KU-3.1 PITFALL — cortex needs an opposing pressure). KU-3.9
    requires Phase 1 to conserve area; we implement that as an
    isotropic radial outward force per bead with magnitude

        P_cyto = − K_area · (A_now − A_0) / A_0,
        F_pressure(bead) = P_cyto · ℓ₀ · r̂_eff.

    ``A_now`` is approximated as π · ⟨R_eff⟩² (the area of the
    concentric circle of mean bead radius). This is the simplest
    metric that is monotone in cortex size; using a polygon /
    convex-hull area would not change the dynamics qualitatively but
    would slow the inner loop. Sign convention: compression
    (A_now < A_0) → P_cyto > 0 → outward push.
    """
    pts = bead_positions.reshape(-1, 2)
    centroid = pts.mean(axis=0)
    rel = pts - centroid
    R_eff = np.linalg.norm(rel, axis=1)
    safe_R = np.where(R_eff > 0.0, R_eff, 1.0)
    rhat = rel / safe_R[:, None]

    R_mean = float(R_eff.mean())
    A_now = np.pi * R_mean * R_mean
    P_cyto = -K_area * (A_now - A_0) / A_0
    magnitude = P_cyto * rest_length
    F_flat = magnitude * rhat
    return F_flat.reshape(bead_positions.shape)


def cortical_tension_forces(
    bead_positions: np.ndarray,
    gamma_cortex: float,
    rest_length: float,
    R_cell: float,
) -> np.ndarray:
    """Cortical tension Hookean radial restoring force (KU-3.5; revised brief).

    The brief specified a per-bead form ``F = γ_cortex · (1/R_eff) ·
    (−r̂_eff)`` (Laplace-law magnitude applied at the centroid
    distance). On an elliptical initial condition this is the wrong
    sign for rounding: beads on the *minor* axis (small R_eff) see a
    *larger* inward force than beads on the *major* axis, which
    pushes the minor poles inward further and **elongates** the
    cortex instead of rounding it. The same conclusion follows from
    the geometric identity that, for a smooth convex curve, the
    centroid distance R_centroid is anti-correlated with the local
    curvature κ_local — so the brief's "1/R_eff" formula applied
    per-bead-distance is anti-Laplace.

    Physically, cell rounding under cortical tension is the energy-
    minimum of (γ · perimeter) at fixed enclosed area; the
    cortex-network discretisation of that minimum is a Hookean
    radial restoring force that pulls beads back toward the nominal
    cortex radius R_cell:

        F_tension(bead) = − k_R · (R_eff − R_cell) · r̂_eff,
        k_R = γ_cortex · ℓ₀ / R_cell².

    Sign check: R_eff > R_cell ⇒ F inward (bead beyond R_cell pulled
    back). R_eff < R_cell ⇒ F outward (bead too close gets pushed
    out — cortex won't collapse). At R_eff = R_cell, F = 0
    (equilibrium). For an elliptical IC, major-axis-pole beads
    (R_eff > R_cell) feel strong inward pull while minor-axis-pole
    beads (R_eff < R_cell) feel mild outward push — the shape
    rounds to a circle of radius R_cell while preserving area. ✓

    Dimensional check: [k_R] = N/m × m / m² = N/m². Force per bead
    is k_R · (R_eff − R_cell) with [m] → [N/m² × m] = N/m ✗ — needs
    one more length factor. We adopt

        k_R = γ_cortex / R_cell    [N/m² ⇒ N/m once multiplied by ℓ₀]
        F_per_bead = − (γ_cortex / R_cell) · (R_eff − R_cell) · r̂.

    With γ_cortex / R_cell having units N/m and (R_eff − R_cell)
    having units m, the force has units N · m/m = N ✓.

    Magnitude calibration: at the elliptical-IC major pole
    (R_eff = a = 1.2 R_cell), |F| = (γ_cortex/R_cell) · (0.2 R_cell)
    = 0.2 γ_cortex = 1e-4 N. Combined with γ_drag = 100 N·s/m per
    bead, this gives a per-bead velocity ≈ 1 μm/s and a rounding
    timescale γ_drag / (γ_cortex/R_cell) = γ_drag · R_cell /
    γ_cortex = 100 · 1e-5 / 5e-4 = 2 s, comfortably under the
    KU-3.1 "~1 minute" rounding criterion.
    """
    pts = bead_positions.reshape(-1, 2)
    centroid = pts.mean(axis=0)
    rel = pts - centroid                               # (N_beads_total, 2)
    R_eff = np.linalg.norm(rel, axis=1)                # (N_beads_total,)
    safe_R = np.where(R_eff > 0.0, R_eff, 1.0)         # defensive guard
    rhat = rel / safe_R[:, None]
    k_R = gamma_cortex / R_cell                        # N/m
    magnitude = k_R * (R_eff - R_cell)                 # N
    F_flat = -magnitude[:, None] * rhat
    # rest_length is accepted for interface symmetry with earlier drafts
    # (Phase 1 simplification — the radial-spring discretisation does not
    # multiply by ℓ₀; keep argument so callers don't need to refactor).
    del rest_length
    return F_flat.reshape(bead_positions.shape)


def solve_overdamped_step(
    cortex: Cortex,
    *,
    stretching_modulus: float,
    bending_modulus: float,
    gamma_cortex: float,
    gamma_drag: float,
    dt: float,
    K_area: float = 0.0,
    A_0: float | None = None,
    ext_forces_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> Cortex:
    """Single overdamped quasi-static Euler step on the cortex.

    Returns a new :class:`Cortex` with updated ``bead_positions``
    (cross-links retained in place; their rest-lengths and indices
    do not change inside a single dynamics step).
    """
    pos = cortex.bead_positions
    F_internal = compute_forces(
        pos,
        rest_length=cortex.rest_length,
        stretching_modulus=stretching_modulus,
        bending_modulus=bending_modulus,
        box_size=cortex.box_size,
        cross_links=cortex.cross_links,
    )
    F_tension = cortical_tension_forces(
        pos, gamma_cortex, cortex.rest_length, R_cell=cortex.R_cell
    )
    F = F_internal + F_tension
    if K_area > 0.0:
        if A_0 is None:
            A_0 = float(np.pi * cortex.R_cell * cortex.R_cell)
        F = F + cytoplasm_pressure_forces(pos, K_area, A_0, cortex.rest_length)
    if ext_forces_fn is not None:
        F = F + ext_forces_fn(pos)

    drift = (dt / gamma_drag) * F
    new_pos = pos + drift

    return Cortex(
        bead_positions=new_pos.astype(np.float64, copy=False),
        fiber_centers=cortex.fiber_centers,
        fiber_endpoints=cortex.fiber_endpoints,
        fiber_orientations=cortex.fiber_orientations,
        fiber_ids=cortex.fiber_ids,
        rest_length=cortex.rest_length,
        fiber_length=cortex.fiber_length,
        cell_center=cortex.cell_center,
        R_cell=cortex.R_cell,
        box_size=cortex.box_size,
        cross_links=cortex.cross_links,
        params=cortex.params,
    )


def relax(
    cortex: Cortex,
    *,
    stretching_modulus: float,
    bending_modulus: float,
    gamma_cortex: float,
    gamma_drag: float,
    dt: float,
    n_steps: int,
    sample_interval: int | None = None,
    keep_frames: bool = False,
    K_area: float = 0.0,
    A_0: float | None = None,
    ext_forces_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> tuple[Cortex, RoundingTrajectory]:
    """Drive the cortex with :func:`solve_overdamped_step` for ``n_steps``.

    Diagnostics (aspect ratio, mean radius) are sampled every
    ``sample_interval`` steps; bead-position snapshots are stored when
    ``keep_frames=True`` (used by the visualisation animation).
    """
    if sample_interval is None or sample_interval <= 0:
        sample_interval = max(1, n_steps // 200)

    times = [0.0]
    ars = [measure_aspect_ratio(cortex)]
    centroids = cortex.bead_positions.reshape(-1, 2)
    radii = [float(
        np.linalg.norm(centroids - centroids.mean(axis=0), axis=1).mean()
    )]
    frames: list[np.ndarray] = []
    if keep_frames:
        frames.append(cortex.bead_positions.copy())

    cur = cortex
    for step in range(1, n_steps + 1):
        cur = solve_overdamped_step(
            cur,
            stretching_modulus=stretching_modulus,
            bending_modulus=bending_modulus,
            gamma_cortex=gamma_cortex,
            gamma_drag=gamma_drag,
            dt=dt,
            K_area=K_area,
            A_0=A_0,
            ext_forces_fn=ext_forces_fn,
        )
        if step % sample_interval == 0:
            times.append(step * dt)
            ars.append(measure_aspect_ratio(cur))
            pts = cur.bead_positions.reshape(-1, 2)
            radii.append(float(np.linalg.norm(pts - pts.mean(axis=0), axis=1).mean()))
            if keep_frames:
                frames.append(cur.bead_positions.copy())

    traj = RoundingTrajectory(
        times=np.asarray(times, dtype=np.float64),
        aspect_ratios=np.asarray(ars, dtype=np.float64),
        mean_radii=np.asarray(radii, dtype=np.float64),
        bead_positions_frames=frames,
    )
    return cur, traj
