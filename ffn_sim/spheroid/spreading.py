"""L2.2 active spreading: add per-cell motility to the CBM and measure spreading-area growth.

Adds a per-cell active self-propulsion force (constant magnitude, persistent random
direction reassigned every tau steps — the Chen&Zou / Rosenbauer motility recipe) to the
cohesive CBM (``cbm.py``). The active "pressure" expands the cohesive aggregate to a larger
wetted area (active-wetting balance: motility spreads, cohesion + excluded volume resist);
``observables.projected_area`` tracks A(t).

⚠️ QUANTITATIVE STATUS (honest): the spreading AMOUNT is set by the traction-vs-cohesion
ratio. The cohesion D_e is currently a placeholder whose SCALE is wrong (cadherin-derived
D_e gives pN-scale holding; real cell-cell de-adhesion is nN-scale — Omidvar 2014, recovery
pending). So this module validates the MOTILITY MECHANISM (spreading happens, A grows then
plateaus, cluster stays connected) at a balanced REDUCED traction; the quantitative A/A0
sweep waits on the Omidvar-anchored nN cohesion + the scale-bridge traction. The active
force here is therefore a tunable mechanism parameter, NOT yet a literature anchor.

Implementation note (correctness): the per-cell force is indexed by particle order, so the
CBM particle SORTER is disabled (``operations.tuners.clear()``) to keep tag==local index on
a single CPU rank. Harmless for the order-independent G1 metrics.

Sanity Gate
-----------
- Dimensional: f_active [N], applied along a unit vector => force [N]. dt/gamma unchanged.
- Boundary: f_active = 0 reduces exactly to the G1 stable aggregate (no spreading).
- Sign-sense: larger f_active => larger spread area; larger cohesion D_e => smaller spread.
- Conservation: count conserved; with persistent-random (zero-mean) directions the COM
  drift is diffusive, not ballistic (no net propulsion of the whole cluster).
"""

from __future__ import annotations

from typing import Any

import hoomd
import numpy as np
import numpy.typing as npt
from hoomd import md

from ffn_sim.spheroid.cbm import build_cbm_simulation, get_positions
from ffn_sim.spheroid.observables import (
    detached_fraction,
    projected_area,
    radius_of_gyration,
)
from ffn_sim.spheroid.params import ResolvedL2


class ActiveMotility(md.force.Custom):
    """Per-cell active self-propulsion: constant |f| along a persistent random direction.

    Each cell carries a unit direction; all directions are resampled every ``tau_steps``
    (the finite persistence time). Assumes particle order == tag order (sorter disabled).
    """

    def __init__(
        self,
        *,
        f_mag: float,
        n_cells: int,
        tau_steps: int,
        seed: int,
        dim: int = 3,
    ) -> None:
        super().__init__()
        if f_mag < 0.0:
            raise ValueError("f_mag must be >= 0.")
        if tau_steps < 1:
            raise ValueError("tau_steps must be >= 1.")
        self._f = float(f_mag)
        self._tau = int(tau_steps)
        self._dim = int(dim)
        self._rng = np.random.default_rng(seed)
        self._dirs = self._sample(n_cells)
        self._next = self._tau

    def _sample(self, n: int) -> npt.NDArray[np.float64]:
        if self._dim == 2:
            th = self._rng.uniform(-np.pi, np.pi, n)
            return np.column_stack([np.cos(th), np.sin(th), np.zeros(n)])
        v = self._rng.normal(size=(n, 3))
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    def set_forces(self, timestep: int) -> None:
        if timestep >= self._next:
            self._dirs = self._sample(self._dirs.shape[0])
            self._next = timestep + self._tau
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = self._f * self._dirs


def build_cbm_spreading(
    resolved: ResolvedL2,
    n_cells: int,
    *,
    f_active: float,
    tau_steps: int,
    device: hoomd.device.Device | None = None,
    rng: np.random.Generator | None = None,
    dim: int = 3,
) -> tuple[hoomd.Simulation, Any, Any, ActiveMotility, float]:
    """Build a CBM sim with cohesion + excluded volume + per-cell active motility."""
    sim, action, updater, r_cut = build_cbm_simulation(
        resolved, n_cells, device=device, rng=rng
    )
    sim.operations.tuners.clear()  # tag==index for the per-cell active force
    active = ActiveMotility(
        f_mag=f_active, n_cells=n_cells, tau_steps=tau_steps,
        seed=resolved.seed + 1, dim=dim,
    )
    sim.operations.integrator.forces.append(active)
    return sim, action, updater, active, r_cut


def run_spreading(
    resolved: ResolvedL2,
    n_cells: int = 200,
    *,
    f_active: float,
    tau_steps: int = 200,
    settle_steps: int = 10_000,
    spread_steps: int = 40_000,
    n_samples: int = 9,
    device: hoomd.device.Device | None = None,
) -> dict[str, Any]:
    """Settle the aggregate (G1), then switch on motility and record A(t) while spreading.

    Returns the spread-area time series (normalised to the settled area = A/A0), the
    detached fraction, and the initial/settled/final positions.
    """
    # Settle with motility OFF (f_active is applied throughout; we settle first by running
    # a short pre-phase, but the force is constant — so we record A0 right after settle).
    sim, _a, _u, _active, _rc = build_cbm_spreading(
        resolved, n_cells, f_active=f_active, tau_steps=tau_steps, device=device
    )
    sim.run(0)
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos_settled = get_positions(sim)
    a0 = projected_area(pos_settled)

    ts, areas, rgs = [], [], []
    per = max(1, spread_steps // n_samples)
    done = 0
    while done < spread_steps:
        step = min(per, spread_steps - done)
        sim.run(step)
        done += step
        p = get_positions(sim)
        ts.append(done)
        areas.append(projected_area(p))
        rgs.append(radius_of_gyration(p))
    pos_final = get_positions(sim)

    det = detached_fraction(
        pos_final, d_crit=3.0 * resolved.morse_r0, neighbor_radius=1.5 * resolved.morse_r0
    )
    return {
        "n_cells": n_cells,
        "f_active": f_active,
        "tau_steps": tau_steps,
        "a0": a0,
        "t": np.array(ts),
        "area": np.array(areas),
        "area_over_a0": np.array(areas) / a0 if a0 > 0 else np.array(areas),
        "rg": np.array(rgs),
        "detached_fraction_final": det,
        "pos_init": pos_init,
        "pos_settled": pos_settled,
        "pos_final": pos_final,
    }
