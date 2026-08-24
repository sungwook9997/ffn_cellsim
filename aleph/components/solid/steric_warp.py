"""Warp-CUDA filament-filament excluded volume (I2b) — the WCA steric force over a device HashGrid.

Runtime primitive for the CLAUDE.md hard worked-example "LJ repulsive excluded volume ON from Phase 1"
and ENGINE_ARCHITECTURE_PLAN §3's living-cortex steric gap. Generalizes the existing MT-tip<->cortex
``ff/network_warp.py::soft_contact_kernel`` (a single hand-wired pair per MT tip) to ALL filament nodes
via a shared device ``wp.HashGrid``, so cortex / SF / MT filaments cannot interpenetrate.

Design (ported from ``dcm/dcm_neighbor_warp.py::cohesion_grid_kernel``, read-only):
  * Grid is built on a float32 copy of the node positions each solve (positions ~um => cell index is
    O(1), float32-safe); the FORCE math uses the f64 positions, so the grid only prunes far pairs that
    contribute exactly 0 -> the result equals the brute-force O(N^2) WCA sum to summation tolerance
    (the ``steric_reference.py`` oracle). One thread per node i, **own-row read-modify-write** accumulate
    ``force[i] += sum_j F_WCA`` — race-free (each node's row is written by exactly one thread), so no
    atomics and the sum is deterministic. Newton's third law: node j's thread adds the equal-opposite
    term, net internal force ~ 0.
  * Same-filament pairs are excluded via per-node ``fiber_id`` (a filament's own segment spacing is held
    by the bending/inextensibility laws, not steric self-repulsion) — the analogue of the ``edge_cell``
    same-cell skip in the DCM edge-edge kernel.

Force primitive contract (AC_PARALLEL_SESSIONS §1.4): ``StericForce.accumulate(state, out_force)`` adds
this primitive's per-node contribution to the shared inner-loop force array; the lead-owned integrator
sums ``MyosinForce`` (I3) + ``StericForce`` (I2b) + ``PressureCoupling`` (I1a). Field-INDEPENDENT: reads
only node positions, never another increment's live state.

Parameters: ``sigma`` (steric diameter, physical; params_i0b2b.yaml) and ``k_ev`` (contact stiffness, the
Magic-Number-Block numerical repulsion scale) -> ``epsilon = k_ev sigma^2 / (36 . 2^(2/3))``. ``force_cap``
bounds the near-core repulsion so the CFL stiffness stays finite (params_i0b2b.yaml).

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). No HOOMD, no CPU simulation path. This module is authored
SOURCE gated by the lead on the gbook A5000 against the ``steric_reference`` oracle (native-gate spec in
``INTEGRATION.md``); the dev Mac runs only the pure-NumPy oracles.
"""

from __future__ import annotations

import warp as wp

from aleph.components.solid.wca_analytic import epsilon_from_contact_stiffness, wca_cutoff

wp.init()

# WCA geometric constant, device-side (r_c / sigma = 2^(1/6)).
_TWO_POW_1_6 = wp.constant(wp.float64(2.0 ** (1.0 / 6.0)))


@wp.kernel
def pos_to_f32(pos: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3)):
    """Copy f64 node positions to a float32 array for the hash-grid build (grid indexing is O(1))."""
    i = wp.tid()
    p = pos[i]
    out[i] = wp.vec3(wp.float32(p[0]), wp.float32(p[1]), wp.float32(p[2]))


@wp.kernel
def wca_steric_kernel(
    grid: wp.uint64,                          # hash grid over node positions
    qpts: wp.array(dtype=wp.vec3),            # f32 node positions (grid query points)
    pos: wp.array(dtype=wp.vec3d),            # f64 node positions (force math)
    fiber_id: wp.array(dtype=wp.int32),       # per-node filament id (same-filament pairs skipped)
    active: wp.array(dtype=wp.int32),         # per-node active flag (>=0 live, <0 dormant/parked)
    radius: wp.float32,                       # query radius = r_c (f32)
    sigma: wp.float64,                        # steric diameter
    epsilon: wp.float64,                      # WCA energy scale
    f_cap: wp.float64,                        # max repulsive pair force (<=0 disables the cap)
    force: wp.array(dtype=wp.vec3d),          # (N,) own-row accumulate (force[i] += acc)
):
    """WCA excluded volume: node i is pushed away from every non-same-filament neighbour within ``r_c``.

    ``F(r) = (24 eps / r)(2 (sigma/r)^12 - (sigma/r)^6)`` for ``r < r_c``, purely repulsive; optionally
    capped at ``f_cap`` so the near-core stiffness (hence the CFL dt) stays bounded. Own-row write:
    each thread owns node i's row, sums its neighbours, and adds into ``force[i]`` (race-free, atomic-free).
    """
    i = wp.tid()
    if active[i] < wp.int32(0):
        return
    z = wp.float64(0.0)
    fi = fiber_id[i]
    pi = pos[i]
    r_c = _TWO_POW_1_6 * sigma
    acc = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, qpts[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and active[j] >= wp.int32(0) and fiber_id[j] != fi:
            r_vec = pi - pos[j]
            d = wp.length(r_vec)
            if d > wp.float64(1.0e-12) and d < r_c:
                sr6 = (sigma / d) * (sigma / d) * (sigma / d)
                sr6 = sr6 * sr6                               # (sigma/d)^6
                fmag = (wp.float64(24.0) * epsilon / d) * (wp.float64(2.0) * sr6 * sr6 - sr6)
                if f_cap > z and fmag > f_cap:
                    fmag = f_cap                              # bound near-core repulsion (CFL guard)
                acc = acc + r_vec * (fmag / d)                # push i AWAY from j
    force[i] = force[i] + acc


class StericForce:
    """Filament-filament excluded-volume force primitive (I2b) over a device ``wp.HashGrid``.

    Implements the ``accumulate(state, out_force)`` contract of AC_PARALLEL_SESSIONS §1.4: adds the WCA
    steric per-node force to a shared inner-loop force array. Owns its own hash grid + f32 position scratch;
    reads only node positions (field-independent).

    Args:
        fiber_id: Per-node filament id, host int array shape (N,). Same-filament pairs are excluded.
        sigma: Steric diameter [um] (physical excluded-volume shell; params_i0b2b.yaml).
        k_ev: Effective contact stiffness [pN/um] — the Magic-Number-Block numerical repulsion scale.
        device: Warp device (CUDA). No hard-coded id (I0-A hardware rule).
        active: Optional per-node active flag shape (N,) (>=0 live, <0 dormant); defaults to all live.
        force_cap: Optional max repulsive pair force [pN]; bounds near-core stiffness for the CFL.
        grid_dim: Hash-grid dimensions (nx, ny, nz); the cell size is the query radius r_c at build.
    """

    def __init__(
        self,
        fiber_id,
        sigma: float,
        k_ev: float,
        device,
        active=None,
        force_cap: float | None = None,
        grid_dim: tuple[int, int, int] = (128, 128, 128),
    ) -> None:
        import numpy as np

        self.sigma = float(sigma)
        self.k_ev = float(k_ev)
        self.epsilon = epsilon_from_contact_stiffness(k_ev, sigma)
        self.r_c = wca_cutoff(sigma)
        self.f_cap = float(force_cap) if force_cap is not None else -1.0
        self.device = device

        fid = np.ascontiguousarray(np.asarray(fiber_id, dtype=np.int32))
        self.n = int(fid.shape[0])
        self.fiber_id = wp.array(fid, dtype=wp.int32, device=device)
        if active is None:
            act = np.zeros(self.n, dtype=np.int32)            # all live (>=0)
        else:
            act = np.ascontiguousarray(np.asarray(active, dtype=np.int32))
        self.active = wp.array(act, dtype=wp.int32, device=device)
        self._qpts = wp.zeros(self.n, dtype=wp.vec3, device=device)
        self.grid = wp.HashGrid(grid_dim[0], grid_dim[1], grid_dim[2], device=device)

    def accumulate(self, state, out_force) -> None:
        """Add the WCA steric per-node force to ``out_force`` (the §1.4 inner-loop force primitive).

        Args:
            state: Live mechanical state exposing ``pos`` — a ``wp.array(dtype=wp.vec3d)`` of node
                positions [um], length N. (Duck-typed against the lead-owned integrator state.)
            out_force: Shared per-node force accumulator, ``wp.array(dtype=wp.vec3d)`` length N [pN].
        """
        pos = state.pos
        self._accumulate_pos(pos, out_force)

    def _accumulate_pos(self, pos, out_force) -> None:
        """Lower-level entry point taking the position array directly (used by the native parity gate)."""
        wp.launch(pos_to_f32, dim=self.n, inputs=[pos, self._qpts], device=self.device)
        self.grid.build(self._qpts, float(self.r_c))
        wp.launch(
            wca_steric_kernel,
            dim=self.n,
            inputs=[
                self.grid.id,
                self._qpts,
                pos,
                self.fiber_id,
                self.active,
                wp.float32(self.r_c),
                wp.float64(self.sigma),
                wp.float64(self.epsilon),
                wp.float64(self.f_cap),
                out_force,
            ],
            device=self.device,
        )
