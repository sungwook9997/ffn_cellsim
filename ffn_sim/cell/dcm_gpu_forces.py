"""GPU device-dispatch wiring for the DCM cell-cell tent contact (+ substrate).

This module is the WIRING that lets the large native-mesh + tent DCM spheroid run
GPU-resident on the lab RTX A5000, using the validated standalone kernels in
``ffn_sim/gpu_opt/{kernels_cpu,kernels_gpu}.py`` (all-parity, group_pair 827×; see
``gpu_opt/GPU_RESULT_A5000.md``). The kernels are standalone array math; here they
are dispatched through HOOMD's per-device local-array contexts inside
``md.force.Custom`` so positions/forces stay on whichever device the simulation
runs on:

  * on a ``hoomd.device.GPU``: read positions via ``gpu_local_snapshot`` as cupy,
    compute with ``kernels_gpu`` (uniform-grid neighbour list), write via
    ``gpu_local_force_arrays`` — NO host transfer in the hot loop.
  * on a ``hoomd.device.CPU``: read via ``cpu_local_snapshot`` as numpy, compute
    with ``kernels_cpu``, write via ``cpu_local_force_arrays``. This path is
    BIT-IDENTICAL to the existing ``cell/dcm_contact.py::DcmTentContact`` (same
    tent law, extracted as ``kernels_cpu.tent_contact_forces``).

The kernels and this dispatch are STANDALONE WIRING — they do not modify the live
DCM modules (``dcm_contact.py``, ``dcm_native_shell.py``, ``dcm_active.py``,
``integrator/baoab.py``). The GPU path is structurally complete and CPU-parity
tested on a CUDA-less Mac; the GPU speedup is validated separately on the gbook
A5000 (this Mac has no CUDA device).

Units SI.
"""

from __future__ import annotations

from contextlib import contextmanager

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.gpu_opt import kernels_cpu


# ---------------------------------------------------------------------------
# Device dispatch helper
# ---------------------------------------------------------------------------
def on_gpu(sim_or_state) -> bool:
    """True if the simulation/state runs on a ``hoomd.device.GPU``.

    Accepts either a ``hoomd.Simulation`` (has ``.device``) or a state object
    (has ``._simulation.device``) — the two contexts a Custom force sees.
    """
    dev = getattr(sim_or_state, "device", None)
    if dev is None:
        dev = getattr(getattr(sim_or_state, "_simulation", None), "device", None)
    return isinstance(dev, hoomd.device.GPU)


class DeviceDispatch:
    """Bundles the right snapshot/force-array contexts, array module, and kernel
    module for a Custom force, chosen by the simulation's device.

    The cupy import is GUARDED — it is only attempted on a GPU device, so this
    module imports and the CPU path runs on a machine with no CUDA / no cupy.

    Attributes (after construction):
        gpu (bool): True on a GPU device.
        xp: the array module (``cupy`` on GPU, ``numpy`` on CPU).
        kernels: the kernel module (``kernels_gpu`` on GPU, ``kernels_cpu`` else).
    """

    def __init__(self, force: md.force.Custom) -> None:
        self._force = force
        self._state = force._state
        self.gpu = on_gpu(self._state)
        if self.gpu:
            import cupy as cp                                # guarded GPU-only import
            from ffn_sim.gpu_opt import kernels_gpu
            self.xp = cp
            self.kernels = kernels_gpu
        else:
            self.xp = np
            self.kernels = kernels_cpu

    @contextmanager
    def snapshot(self):
        """Yield the device-appropriate local snapshot (gpu_* vs cpu_*)."""
        ctx = (self._state.gpu_local_snapshot if self.gpu
               else self._state.cpu_local_snapshot)
        with ctx as snap:
            yield snap

    @contextmanager
    def force_arrays(self):
        """Yield the device-appropriate local force-array context."""
        ctx = (self._force.gpu_local_force_arrays if self.gpu
               else self._force.cpu_local_force_arrays)
        with ctx as arr:
            yield arr


# ---------------------------------------------------------------------------
# GPU-capable cell-cell tent contact
# ---------------------------------------------------------------------------
class DcmTentContactGPU(md.force.Custom):
    """Device-dispatched SimuCell3D bilinear-tent cell-cell contact.

    Same construction signature and physics as
    ``cell/dcm_contact.py::DcmTentContact`` — the tent law lives in
    ``gpu_opt.kernels_{cpu,gpu}.tent_contact_forces`` (parity-matched). On a GPU
    device the force is computed in cupy via the K4-style uniform-grid neighbour
    list with NO host transfer; on a CPU device it is bit-identical to the
    original ``DcmTentContact`` (both call the numpy ``tent_contact_forces`` over
    a brute-force pair search → same FP ops).

    Args:
        cell_of_node: (N,) int, mutable; cell id of each node, −1 = dormant.
            (The proliferation updater can keep editing this in place; it is read
            fresh each ``set_forces``.)
        r_contact: m, effective contact radius (repulsion below this).
        c_adh: m, adhesion cutoff (tent peaks at c_adh/2).
        rep_strength: Pa/m, repulsion stiffness ξ.
        adh_strength: Pa/m, adhesion stiffness ω (0 ⇒ repulsion-only).
        patch_area: m², per-node contact patch area (= 4πR²/n_nodes).
        force_cap: N, per-pair force-magnitude cap (BAOAB int32 guard).
        cad_mult: optional (n_cells,) float; per-cell adhesion multiplier
            (mult = √(cad_mult_i·cad_mult_j)). None ⇒ all 1.
    """

    def __init__(self, *, cell_of_node: np.ndarray, r_contact: float, c_adh: float,
                 rep_strength: float, adh_strength: float, patch_area: float,
                 force_cap: float = 5.0e-8,
                 cad_mult: np.ndarray | None = None) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node            # (N,) int, −1 = dormant
        self.r_contact = float(r_contact)
        self.c_adh = float(c_adh)
        self.rep = float(rep_strength)
        self.omega = float(adh_strength)
        self.A = float(patch_area)
        self.force_cap = float(force_cap)
        self.cad_mult = cad_mult                    # (n_cells,) or None
        self.r_search = max(self.r_contact, self.c_adh)
        self._d: DeviceDispatch | None = None       # lazy (device known at run(0))

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            # local row → global (tag) order, then back. argsort is stable enough
            # for the unique-tag permutation (mirrors DcmTentContact exactly).
            perm = xp.argsort(tag)
            pos_g = pos[perm]
            cell_g = xp.asarray(self.cell_of_node)
            cad = None if self.cad_mult is None else xp.asarray(self.cad_mult)

            F_g = d.kernels.tent_contact_forces(
                pos_g, cell_g, self.r_contact, self.c_adh, self.rep,
                self.omega, self.A, self.force_cap, cad_mult=cad)

            # scatter global-order forces back to local rows: F[perm] = F_g
            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# (Optional) GPU-capable substrate force via K2 plane_well_forces
# ---------------------------------------------------------------------------
class DcmSubstrateForceGPU(md.force.Custom):
    """Device-dispatched adhesive substrate (capped-harmonic well at z = z0).

    Same physics as ``cell/dcm.py::DcmSubstrateForce`` but the per-node well is
    evaluated through ``gpu_opt.kernels_{cpu,gpu}.plane_well_forces`` (K2). The
    K2 well is ``|dz|≤w: F_z=−k·dz; dz<−w: F_z=+k·w; else 0`` which matches
    DcmSubstrateForce with ``k = k_well`` and ``w = adh_range``. Potential energy
    is left zero on this path (force-only; the original's U is diagnostic, not
    used by BAOAB). On GPU it is fully device-resident.

    Args mirror DcmSubstrateForce: z0, W_cs (J/node), adh_range (m), k_sub (N/m
    or None for a rigid dish).
    """

    def __init__(self, *, z0: float, W_cs: float, adh_range: float,
                 k_sub: float | None = None) -> None:
        super().__init__(aniso=False)
        self.z0 = float(z0)
        self.W_cs = float(W_cs)
        self.rng = float(adh_range)
        k_well = 2.0 * self.W_cs / (self.rng ** 2)
        self.k_well = k_well if k_sub is None else (k_well * k_sub) / (k_well + k_sub)
        self._d: DeviceDispatch | None = None

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            F = d.kernels.plane_well_forces(pos, self.z0, self.k_well, self.rng)
            U = xp.zeros(n, dtype=xp.float64)
        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# GPU-capable ACTIVE RIM TRACTION (coarse-grained lamellipodium + belt + clutch)
# ---------------------------------------------------------------------------
class DcmActiveRimTractionGPU(md.force.Custom):
    """Device-dispatched active OUTWARD traction on basal rim-cell nodes.

    Same force law, same construction signature, and same diagnostics as
    ``cell/dcm_active.py::ActiveRimTraction`` (the coarse-grained lamellipodium +
    actomyosin contraction-belt + FA-clutch engine at the DCM cell scale — see
    that class's docstring for the physics). The only difference is the device
    dispatch: this class routes its array math through :class:`DeviceDispatch`, so

      * on a ``hoomd.device.GPU`` the per-cell centroid / neighbour-crowding /
        rim-detection / per-node traction is computed in cupy on the
        ``gpu_local_snapshot`` (forces written to ``gpu_local_force_arrays``) —
        NO host transfer in the hot loop; and
      * on a ``hoomd.device.CPU`` the array module is numpy and the operations are
        the SAME ones (in the same order) as ``ActiveRimTraction.set_forces`` — so
        the CPU path is BIT-IDENTICAL (max abs force diff < 1e-12 N; see
        ``tests/test_dcm_active_gpu_parity.py``).

    The per-cell loop (rim detection, basal/apical node split, outward direction)
    stays a Python loop over the (small) active-cell count on both devices — the
    HOT cost is the all-node force write, which IS device-resident. The low-cadence
    STATE updaters (necrosis / pressure / junction / division) are NOT ported here;
    they remain the CPU ``hoomd.custom.Action`` updaters in ``dcm_active.py`` (they
    run batched at low cadence, off the per-step path — porting them buys nothing).

    Args: identical to ``ActiveRimTraction.__init__`` (cell_of_node, ranges,
    active, int_mult, R_cell, z0, f_act, f_cap, ramp_steps, contact_band,
    neighbour_factor, max_neighbours, integrin_switch_gain, belt_factor).
    """

    def __init__(self, *, cell_of_node: np.ndarray, ranges, active: np.ndarray,
                 int_mult: np.ndarray, R_cell: float, z0: float, f_act: float,
                 f_cap: float, ramp_steps: int, contact_band: float,
                 neighbour_factor: float, max_neighbours: int,
                 integrin_switch_gain: float, belt_factor: float,
                 arrest_radius_factor: float | None = None,
                 arrest_width: float = 0.18,
                 arrest_settle_steps: int = 4000) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node
        self.ranges = ranges
        self.active = active
        self.int_mult = int_mult
        self.R = float(R_cell)
        self.z0 = float(z0)
        self.f_act = float(f_act)
        self.f_cap = float(f_cap)
        self.ramp_steps = max(1, int(ramp_steps))
        self.contact_band = float(contact_band)
        self.r_neigh = float(neighbour_factor * R_cell)
        self.max_neigh = int(max_neighbours)
        self.switch_gain = float(integrin_switch_gain)
        self.belt = float(belt_factor)
        # SPREADING-ARREST (membrane-tension / contact-inhibition stall). The
        # outward rim traction is smoothly switched OFF as a rim cell's radial
        # distance r from the cluster centroid approaches a physiological maximum
        # r_max = arrest_radius_factor · R0_cluster (R0_cluster = the INITIAL
        # cluster radius, captured once at the first non-trivial call so the cap is
        # N-invariant — a multiplicative factor on the natural cluster scale, NOT an
        # empirical absolute length). The arrest gain is a SMOOTH ramp over a width
        # band so the stall is BAOAB-safe (no hard clip, no discontinuous force):
        #     g(r) = clip((r_max − r) / (arrest_width · r_max), 0, 1)
        # so traction is full deep inside (r ≤ (1−width)·r_max), ramps to 0 across
        # the band, and is exactly 0 beyond r_max. The mechanism is the
        # motor-clutch / membrane-tension feedback stall (Chan-Odde /
        # Elosegui-Artola; the same protrusion-stalls-as-tension-rises physics the
        # single-cell reference cell/spreading_drive.py models via the max_radius
        # footprint cap). arrest_radius_factor=None ⇒ NO arrest.
        #
        # DEFAULT = None (arrest OFF), and it is OPT-IN, on purpose (2026-06-11):
        #   1. The over-spread it was meant to cap is NOT active-traction-driven —
        #      the diagnostic sweep (outputs/h_dcm_gpu_lod/probe3_on.json) shows
        #      A/A₀ still BLOWS UP (gain→0.0, full arrest) because the runaway is a
        #      convex-hull footprint artifact of a few scattered cells + the
        #      now-fixed FLOAT artifact, not the rim push. So arrest does not solve
        #      the problem it was added for.
        #   2. Choosing the factor to land A/A₀ in a band (the old "1.7² ≈ 2.9 ∈
        #      [2,4]" note) is exactly the no-magic-number / no-fitting violation
        #      (CLAUDE.md): a tuned proxy, not a derived physiological setpoint.
        #   3. It is a LUMPED proxy for membrane-tension / contact-inhibition; the
        #      architectural rule prefers the mechanistic alternative. The correct
        #      cure for over-spread is the COHESIVE-start + real wetting baseline
        #      (cells adhered from t=0), not a traction cap.
        # The code + its dedicated test (test_dcm_spreading_arrest.py) are kept so a
        # future mechanistic membrane-tension force can be A/B'd, but production runs
        # OFF until a derived (non-tuned) cap is justified.
        self.arrest_factor = (None if arrest_radius_factor is None
                              else float(arrest_radius_factor))
        self.arrest_width = float(max(arrest_width, 1.0e-3))
        self._R0_cluster: float | None = None   # settled cluster radius (lazy)
        # R0 is captured AFTER the cluster has compacted from the gapped start (the
        # gapped t≈0 radius over-estimates the natural scale → an unreachable cap).
        # Captured at the first call with timestep ≥ arrest_settle_steps.
        self.arrest_settle_steps = max(0, int(arrest_settle_steps))
        # diagnostics (read by the driver) — kept on host (numpy) like the CPU class
        self.rim_cells: np.ndarray = np.empty(0, dtype=np.int64)
        self.f_per_cell: dict[int, float] = {}
        self._ramp = 0.0
        self.arrest_gain_mean: float = 1.0      # mean per-rim arrest gain (diag)
        self.rn_max: float = 0.0                # max rim radial distance [m] (diag)
        self._d: DeviceDispatch | None = None      # lazy (device known at run(0))

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def _capture_R0_cluster(self, rn, timestep: int = 0) -> None:
        """Record the SETTLED cluster radius once (max in-plane rim radius).

        ``rn`` is the (n_active,) array of in-plane radial distances of the active
        cell centroids from the cluster centroid. The max rim radius is also stored
        each call as ``rn_max`` (diagnostic). The arrest cap ``r_max`` is a
        multiplicative factor on this SETTLED cluster radius, so it is N-invariant
        (it scales with the cluster, not an absolute empirical length).

        Capture is DEFERRED until ``timestep ≥ arrest_settle_steps`` because the
        build starts cells GAPPED (spacing ~2.3·R); the gapped t≈0 radius
        over-estimates the natural scale (the cluster compacts before it spreads),
        which would put the cap out of reach and disable the arrest. Capturing after
        the settle anchors the cap to the real compacted cluster size.
        """
        try:
            self.rn_max = float(rn.max())
        except Exception:  # pragma: no cover — degenerate
            self.rn_max = 0.0
        if self.arrest_factor is None or self._R0_cluster is not None:
            return
        if timestep < self.arrest_settle_steps:
            return
        # guard a degenerate state (all centroids coincident) → defer until non-zero.
        if self.rn_max > 1.0e-9:
            self._R0_cluster = self.rn_max

    def _arrest_gain(self, xp, rn):
        """Per-active-cell smooth arrest gain g(r) ∈ [0, 1] (membrane-tension stall).

        ``g(r) = clip((r_max − r) / (arrest_width · r_max), 0, 1)`` with
        ``r_max = arrest_radius_factor · R0_cluster``. Returns all-ones (no arrest)
        when arrest is disabled or the initial cluster radius is not yet captured.

        Args:
            xp: the array module (numpy / cupy).
            rn: (n_active,) in-plane radial distances of active centroids.

        Returns:
            (n_active,) array of arrest gains in [0, 1].
        """
        if self.arrest_factor is None or self._R0_cluster is None:
            return xp.ones_like(rn)
        r_max = self.arrest_factor * self._R0_cluster
        if r_max <= 0.0:
            return xp.ones_like(rn)
        g = (r_max - rn) / (self.arrest_width * r_max)
        return xp.clip(g, 0.0, 1.0)

    def _footprint_radius(self, xp, act_blk, basal):
        """Robust colony FOOTPRINT radius — p90 in-plane radius of active basal nodes.

        ``act_blk`` is (na, nv, 3) active-cell node positions, ``basal`` (na, nv) the
        per-node basal mask. The footprint radius is the 90th-percentile in-plane
        distance of the active BASAL nodes from their collective centroid — a
        measure ROBUST to a few scattered outliers (which inflate a max / convex
        hull), so it tracks the bulk colony spread, the quantity that genuinely
        blows up. Returns a 0-d array (the array module's scalar).
        """
        b = basal.reshape(-1)                                  # (na*nv,)
        nodes_xy = act_blk[:, :, :2].reshape(-1, 2)            # (na*nv, 2)
        # basal-node centroid (mean of basal nodes only).
        w = b.astype(nodes_xy.dtype)
        denom = w.sum()
        # guard: if no basal node, fall back to all-node centroid (denom→nv*na).
        safe_denom = xp.where(denom > 0, denom, xp.asarray(nodes_xy.shape[0],
                                                           dtype=nodes_xy.dtype))
        cen_xy = (nodes_xy * w[:, None]).sum(axis=0) / safe_denom   # (2,)
        d = nodes_xy - cen_xy[None, :]
        r = xp.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2)              # (na*nv,)
        # p90 over the BASAL nodes ONLY, robust to a few scattered outliers (which
        # inflate a max/hull) → tracks the BULK footprint, not measurement noise.
        # cupy has no nanpercentile, so compute it via a sort (cupy-safe, no host
        # sync): set non-basal radii to 0 (basal radii are > 0), sort ascending → the
        # n_basal real radii occupy the TOP positions; the p90 is the basal radius at
        # rank ⌈0.9·n_basal⌉ from the bottom of that top block, i.e. ascending index
        # ``N − n_basal + ⌊0.9·(n_basal−1)⌋``.
        r_basal0 = xp.where(b, r, 0.0)                        # (N,), non-basal → 0
        rs = xp.sort(r_basal0)                                # ascending
        N = rs.shape[0]
        n_basal = b.sum()
        # guard: if no basal node, return 0 (caller treats as no spread).
        nb = xp.maximum(n_basal, 1)
        # 90th-percentile rank within the basal block (lower-interpolation index).
        idx_in_block = xp.floor(0.9 * (nb.astype(rs.dtype) - 1.0)).astype(xp.int64)
        idx = (N - nb) + idx_in_block
        idx = xp.clip(idx, 0, N - 1)
        return rs[idx]

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            pos_g = pos[perm]
            F_g = xp.zeros_like(pos_g)

            ramp = float(min(1.0, timestep / self.ramp_steps))
            self._ramp = ramp
            active_ids = xp.where(xp.asarray(self.active))[0]

            if int(active_ids.shape[0]) < 2 or ramp <= 0.0:
                F = xp.empty_like(pos)
                F[perm] = F_g
                U = xp.zeros(n, dtype=xp.float64)
                self.rim_cells = np.empty(0, dtype=np.int64)
                self.f_per_cell = {}
                with d.force_arrays() as arr:
                    arr.force[:] = F
                    arr.potential_energy[:] = U
                return

            # per-cell centroid (same stacking order as ActiveRimTraction).
            cents = xp.asarray(
                [pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                 for c in active_ids])
            cluster_cen = cents.mean(0)
            # neighbour count by centroid distance (rim = few neighbours)
            d2 = xp.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
            within = d2 < self.r_neigh ** 2
            xp.fill_diagonal(within, False)
            crowd = within.sum(axis=1)

            # SPREADING-ARREST keys on the COLONY FOOTPRINT RADIUS (same global gain
            # for every rim cell): the robust p90 in-plane radius of ALL active basal
            # nodes about their collective centroid. The over-spread is the basal
            # footprint exploding as cells DISPERSE (per-cell membrane radius and cell
            # centroids both stay bounded — only the footprint tracks the blow-up).
            # As the footprint → R_max=factor·R0_footprint the traction → 0 for every
            # rim cell → the colony stops. Same quantity + helper as the vec path
            # (uniform contiguous blocks here ⇒ act_blk reshape is valid → identical).
            zc = self.z0 + self.contact_band * self.R
            n_cells_all = len(self.ranges)
            nv_blk = self.ranges[0][1] - self.ranges[0][0]
            pos_blk_l = pos_g.reshape(n_cells_all, nv_blk, 3)
            act_blk_l = pos_blk_l[active_ids]
            basal_blk_l = act_blk_l[:, :, 2] < zc
            R_foot = self._footprint_radius(xp, act_blk_l, basal_blk_l)
            self._capture_R0_cluster(xp.reshape(R_foot, (1,)), timestep)
            arrest_scalar = float(self._arrest_gain(xp, R_foot))

            int_mult = xp.asarray(self.int_mult)
            rim_list: list[int] = []
            f_per_cell: dict[int, float] = {}
            for k in range(int(active_ids.shape[0])):
                c = int(active_ids[k])
                if int(crowd[k]) > self.max_neigh:
                    continue  # interior cell (well-coordinated) — not a rim cell
                lo, hi = self.ranges[c]
                cell_pos = pos_g[lo:hi]
                basal = cell_pos[:, 2] < zc
                if not bool(basal.any()):
                    continue  # not in substrate contact — cannot lamellipodiate
                rim_list.append(c)
                # in-plane outward direction (cluster centroid -> cell centroid)
                rxy = cents[k][:2] - cluster_cen[:2]
                rn = float(xp.hypot(rxy[0], rxy[1]))
                if rn < 1e-12:
                    continue  # cell at the very centre — no defined outward dir
                # cupy forbids building an array from a list of device scalars
                # (rxy[i]/rn are 0-d device arrays) -> stay in array ops: the
                # (2,) in-plane unit vector concatenated with a (1,) z=0.
                rhat = xp.concatenate(
                    [rxy / rn, xp.zeros(1, dtype=rxy.dtype)])
                gain = float(int_mult[c])
                # spreading-arrest: the single global colony-footprint arrest gain
                # reduces the outward traction toward 0 as the footprint → R_max.
                arrest = arrest_scalar
                fmag = ramp * self.f_act * gain * arrest
                fmag = min(fmag, self.f_cap)
                f_per_cell[c] = fmag
                # basal nodes pull OUTWARD (lamellipodium + clutch grip)
                idx = xp.where(basal)[0]
                F_g[lo:hi][idx] += fmag * rhat
                # apical nodes: weak INWARD contraction-belt tension
                if self.belt > 0.0:
                    apic = xp.where(~basal)[0]
                    if int(apic.shape[0]):
                        fb = min(self.belt * fmag, self.f_cap)
                        F_g[lo:hi][apic] += -fb * rhat

            # cap per node (BAOAB int32 guard) — never let a node exceed f_cap
            fn = xp.linalg.norm(F_g, axis=1)
            over = fn > self.f_cap
            if bool(over.any()):
                F_g[over] *= (self.f_cap / fn[over])[:, None]

            # diagnostics back to host (numpy) for the driver
            self.rim_cells = np.asarray(
                xp.asnumpy(xp.asarray(rim_list, dtype=xp.int64)) if d.gpu
                else np.array(rim_list, dtype=np.int64))
            self.f_per_cell = f_per_cell
            self.arrest_gain_mean = arrest_scalar
            self.rn_max = float(R_foot)

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# VECTORIZED GPU active rim traction — the per-cell Python loop removed
# ---------------------------------------------------------------------------
class DcmActiveRimTractionGPUVec(DcmActiveRimTractionGPU):
    """Drop-in, math-identical, FULLY-VECTORIZED active rim traction.

    WHY (gbook A5000 profile, 2026-06-11). The full-stack per-step profiler showed
    ``DcmActiveRimTractionGPU.set_forces`` is the single dominant cost at N=200 —
    ~22 ms/step, ~57% of the ~39 ms/step production total (passive total is only
    ~17 ms). The cost is NOT the ``gpu_local`` context (that floor is ~0.2 ms/force)
    and NOT cross-force dispatch (fusing 4 forces → 1 saves only ~0.7 ms): it is the
    PYTHON per-cell orchestration INSIDE this one force —

      * the centroid list-comprehension
        ``cents = xp.asarray([pos_g[lo:hi].mean(0) for c in active_ids])``
        (one device op + a host list build per active cell, ~10 ms at N=200), and
      * the ``for k in range(n_active)`` loop, each iteration doing ``int(c)`` /
        ``int(crowd[k])`` / ``bool(basal.any())`` / ``float(rn)`` / ``float(gain)``
        on 0-d cupy arrays — each a BLOCKING GPU→CPU scalar sync (~1000 syncs/step
        at N=200, ~11 ms) so the A5000 sits idle while Python dispatches.

    This subclass reproduces the EXACT same legacy basal-splay + apical contraction-
    belt law as the parent ``DcmActiveRimTractionGPU.set_forces`` (the law the GPU
    twin's frozen parity gate ``tests/test_dcm_active_gpu_parity.py`` pins to the
    legacy-mode ``ActiveRimTraction``), but with ZERO per-cell Python loop and ZERO
    per-cell scalar syncs: positions are reshaped to ``(n_cells, nv, 3)`` (the build
    uses uniform ``nv`` nodes per cell, contiguous tag blocks), centroids are one
    ``mean(axis=1)``, crowding is one O(n_active²) matrix, the per-cell mask /
    direction / magnitude are computed as whole arrays, and the per-node force is a
    single broadcast write. Output is bit-identical to the parent (max abs diff
    < 1e-10 N on the CPU path; see ``tests/test_dcm_active_vec_parity.py``).

    Construction signature is identical to ``DcmActiveRimTractionGPU`` — it is a true
    drop-in; ``dcm_gpu_build.build_gpu_dcm_simulation(..., fast_active=True)`` wires
    this class instead of the parent. The same mutable shared state (``active``,
    ``int_mult``, ``cell_of_node``) is read live each step.

    ASSUMPTION (asserted at first call): uniform per-cell node count and contiguous
    blocks, i.e. ``ranges[c] == (c·nv, (c+1)·nv)``. The GPU-friendly builds
    (``build_gpu_dcm_simulation`` / ``build_gpu_spheroid_prolif``) always satisfy
    this. If a build ever violated it, this subclass refuses (falls back is the
    caller's job) rather than silently miscomputing.
    """

    def _uniform_block(self) -> int:
        """Per-cell node count nv, asserting uniform contiguous ranges."""
        ranges = self.ranges
        nv = ranges[0][1] - ranges[0][0]
        # cheap structural check (host-side ints, no device sync)
        for c, (lo, hi) in enumerate(ranges):
            if lo != c * nv or hi != (c + 1) * nv:
                raise ValueError(
                    "DcmActiveRimTractionGPUVec requires uniform contiguous "
                    f"per-cell node blocks (ranges[{c}]=({lo},{hi}) != "
                    f"({c * nv},{(c + 1) * nv})).")
        return int(nv)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        nv = self._uniform_block()
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            pos_g = pos[perm]
            n_cells = len(self.ranges)

            ramp = float(min(1.0, timestep / self.ramp_steps))
            self._ramp = ramp
            active_mask = xp.asarray(self.active)
            active_ids = xp.where(active_mask)[0]
            na = int(active_ids.shape[0])

            if na < 2 or ramp <= 0.0:
                F = xp.empty_like(pos)
                F[perm] = xp.zeros_like(pos_g)
                U = xp.zeros(n, dtype=xp.float64)
                self.rim_cells = np.empty(0, dtype=np.int64)
                self.f_per_cell = {}
                with d.force_arrays() as arr:
                    arr.force[:] = F
                    arr.potential_energy[:] = U
                return

            # (n_cells, nv, 3) block view — uniform contiguous ranges (asserted).
            pos_blk = pos_g.reshape(n_cells, nv, 3)
            act_blk = pos_blk[active_ids]                 # (na, nv, 3)

            # per-active-cell centroid — ONE reduction, no list comprehension.
            cents = act_blk.mean(axis=1)                  # (na, 3)
            cluster_cen = cents.mean(axis=0)              # (3,)

            # crowding among active centroids (rim = few neighbours).
            diff = cents[:, None, :] - cents[None, :, :]
            d2 = xp.sum(diff * diff, axis=2)
            within = d2 < self.r_neigh ** 2
            xp.fill_diagonal(within, False)
            crowd = within.sum(axis=1)                    # (na,)

            zc = self.z0 + self.contact_band * self.R
            basal = act_blk[:, :, 2] < zc                 # (na, nv) bool
            has_basal = basal.any(axis=1)                 # (na,)

            # rim = crowd ≤ max_neigh AND has a basal node (the parent appends to
            # rim_list under EXACTLY these two conditions, BEFORE the rn check).
            is_rim = (crowd <= self.max_neigh) & has_basal   # (na,)

            # in-plane outward direction (cluster centroid → cell centroid).
            rxy = cents[:, :2] - cluster_cen[:2]          # (na, 2)
            rn = xp.sqrt(rxy[:, 0] ** 2 + rxy[:, 1] ** 2)  # (na,) == hypot
            # force contributes only for rim cells with a defined outward dir.
            contributes = is_rim & (rn >= 1e-12)          # (na,)
            rn_safe = xp.where(rn >= 1e-12, rn, 1.0)
            rhat = xp.zeros((na, 3), dtype=pos_g.dtype)   # (na, 3), z=0
            rhat[:, 0] = rxy[:, 0] / rn_safe
            rhat[:, 1] = rxy[:, 1] / rn_safe

            # SPREADING-ARREST keys on the COLONY FOOTPRINT RADIUS — the robust
            # (90th-percentile) in-plane radius of ALL active basal nodes about their
            # collective centroid. The over-spread is the basal FOOTPRINT exploding as
            # cells DISPERSE outward (the cluster scattering), NOT per-cell membrane
            # flattening nor cell-centroid drift (both stay bounded; measured directly
            # — only the footprint radius tracks the blow-up). This is the colony-scale
            # contact-inhibition-of-locomotion / membrane-tension stall: as the
            # footprint approaches its physiological max R_max = factor·R0_footprint,
            # the outward traction is smoothly switched off for EVERY rim cell (a
            # single global gain), so the colony stops and the footprint plateaus.
            R_foot = self._footprint_radius(xp, act_blk, basal)    # scalar (0-d)
            self._capture_R0_cluster(xp.reshape(R_foot, (1,)), timestep)
            arrest = self._arrest_gain(xp, R_foot)                 # scalar broadcast

            int_mult = xp.asarray(self.int_mult)
            gain = int_mult[active_ids]                   # (na,)
            fmag = xp.minimum(ramp * self.f_act * gain * arrest, self.f_cap)  # (na,)
            fmag = xp.where(contributes, fmag, 0.0)       # (na,) gate non-contrib

            # per-node force on each active cell block (na, nv, 3):
            #   basal nodes  += fmag · rhat
            #   apical nodes += −min(belt·fmag, f_cap) · rhat   (if belt > 0)
            fb = xp.minimum(self.belt * fmag, self.f_cap) if self.belt > 0.0 \
                else xp.zeros_like(fmag)
            # node weight per (cell, node): +fmag on basal, −fb on apical.
            node_w = xp.where(basal, fmag[:, None], -fb[:, None])  # (na, nv)
            # zero apical weight entirely when belt == 0 (fb is all-zero then).
            F_act_blk = node_w[:, :, None] * rhat[:, None, :]      # (na, nv, 3)

            # scatter the active blocks back into the full (n_cells, nv, 3) force.
            F_blk = xp.zeros((n_cells, nv, 3), dtype=pos_g.dtype)
            F_blk[active_ids] = F_act_blk
            F_g = F_blk.reshape(n, 3)

            # cap per node (BAOAB int32 guard) — never let a node exceed f_cap.
            fn = xp.linalg.norm(F_g, axis=1)
            over = fn > self.f_cap
            # always-vectorized cap (xp.where avoids the bool(any) host sync).
            scale = xp.where(over, self.f_cap / xp.where(fn > 0, fn, 1.0), 1.0)
            F_g = F_g * scale[:, None]

            # diagnostics (host): rim_cells = active cells with is_rim True.
            is_rim_host = (xp.asnumpy(is_rim) if d.gpu else np.asarray(is_rim))
            active_ids_host = (xp.asnumpy(active_ids) if d.gpu
                               else np.asarray(active_ids))
            self.rim_cells = active_ids_host[is_rim_host].astype(np.int64)
            # f_per_cell: cells that actually contributed (parent sets it after the
            # rn check) → contributes mask.
            contrib_host = (xp.asnumpy(contributes) if d.gpu
                            else np.asarray(contributes))
            fmag_host = (xp.asnumpy(fmag) if d.gpu else np.asarray(fmag))
            self.f_per_cell = {int(active_ids_host[i]): float(fmag_host[i])
                               for i in np.where(contrib_host)[0]}
            # arrest diagnostic: the single global colony-footprint arrest gain.
            self.arrest_gain_mean = float(arrest)
            self.rn_max = float(self._footprint_radius(xp, act_blk, basal))

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# GPU-capable PER-CELL TURGOR via the validated K1 mesh-pressure kernel
# ---------------------------------------------------------------------------
class DcmTurgorForceGPU(md.force.Custom):
    """Device-dispatched per-cell osmotic turgor via the K1 mesh-pressure kernel.

    Same construction signature and physics as ``cell/dcm.py::DcmTurgorForce`` —
    each cell's enclosed volume is the exact divergence-theorem volume of its
    triangulated shell, the osmotic pressure ``ΔP_c = turgor_dP0 + K_vol·(V0 −
    V_c)/V0`` is applied as an outward face-normal force split to the 3 facet
    nodes. The math lives in ``gpu_opt.kernels_{cpu,gpu}.mesh_pressure_forces``
    (K1, validated on the A5000 — 9.5×, parity OK), grouped by ``face_cell`` (=
    the cell id of each face). This is the GPU-FRIENDLY turgor path: ONE custom
    force, per-cell volume by face-grouping, NO per-cell mesh triangle types — so
    it does NOT trigger the native ``md.mesh.conservation.Volume`` many-triangle-
    types GPU stall (the A5000 0% util blocker at N=60).

    On a ``hoomd.device.GPU`` the force is computed in cupy via ``kernels_gpu``
    (NO host transfer in the hot loop); on a ``hoomd.device.CPU`` it routes
    through ``kernels_cpu`` and is BIT-IDENTICAL to ``DcmTurgorForce`` (both call
    the same divergence-theorem volume + face-normal pressure over the same
    tag-ordered positions → same FP ops; max abs force diff < 1e-12 N, see
    ``tests/test_dcm_turgor_gpu_parity.py``).

    Args: identical to ``DcmTurgorForce.__init__`` (faces, face_cell, n_cells,
    V0, turgor_dP0, K_vol).
    """

    def __init__(self, *, faces: np.ndarray, face_cell: np.ndarray, n_cells: int,
                 V0: float, turgor_dP0: float, K_vol: float) -> None:
        super().__init__(aniso=False)
        self.faces = np.asarray(faces, dtype=np.int64)          # (F,3) global idx
        self.face_cell = np.asarray(face_cell, dtype=np.int64)  # (F,)
        self.n_cells = int(n_cells)
        self.V0 = float(V0)
        self.turgor_dP0 = float(turgor_dP0)
        self.K_vol = float(K_vol)
        self._d: DeviceDispatch | None = None       # lazy (device known at run(0))

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            # local-row order -> tag(global) order, then scatter back (mirrors the
            # CPU DcmTurgorForce permutation exactly so faces index the same nodes).
            perm = xp.argsort(tag)
            pos_g = pos[perm]
            faces = xp.asarray(self.faces)
            face_cell = xp.asarray(self.face_cell)

            F_g = d.kernels.mesh_pressure_forces(
                pos_g, faces, face_cell, self.n_cells,
                self.V0, self.turgor_dP0, self.K_vol)

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# FUSED single-callback force — turgor + tent + substrate + active in ONE pass
# ---------------------------------------------------------------------------
class DcmFusedForceGPU(md.force.Custom):
    """All per-step DCM custom forces summed in ONE ``set_forces`` callback.

    Computes TURGOR (K1 mesh-pressure) + CELL-CELL TENT contact + adhesive
    SUBSTRATE well + (optional) vectorized ACTIVE rim traction in a SINGLE
    ``md.force.Custom`` — one ``gpu_local_snapshot`` enter, one position read, one
    ``argsort``, one ``gpu_local_force_arrays`` write — instead of 3-4 separate
    Custom forces each paying their own context enter/exit + position read each
    step. The per-pair / per-cell MATH is IDENTICAL to the separate forces (it
    calls the same kernels and the same vectorized active law), so the fused total
    equals the sum of the separate forces to floating-point round-off (max abs diff
    < 1e-10 N on the CPU path; see ``tests/test_dcm_fused_parity.py``).

    HONEST SCOPE (gbook A5000 profile, 2026-06-11). Fusing removes the redundant
    per-force context/Python overhead, which the profile measured at only ~0.2 ms
    PER FORCE (so ~0.6-0.7 ms/step saved going 4 → 1) — a SMALL win (~2% of the
    ~39 ms/step production total). The DOMINANT cost is INSIDE the active traction's
    per-cell Python loop, addressed by :class:`DcmActiveRimTractionGPUVec` (which
    this fused force uses for its active term). Fusing is provided for completeness
    and composes with the vectorized active; the big lever is the vectorization.

    Construction takes the SAME parameters as the four separate forces. The active
    block is optional (``with_active=False`` ⇒ passive turgor + tent + substrate).
    Mutable shared state (``cell_of_node``, ``cad_mult``, ``active``, ``int_mult``)
    is read live each step exactly as the separate forces do.
    """

    def __init__(self, *, n_cells: int,
                 # turgor (K1)
                 faces: np.ndarray, face_cell: np.ndarray, V0: float,
                 turgor_dP0: float, K_vol: float,
                 # tent contact (K5)
                 cell_of_node: np.ndarray, r_contact: float, c_adh: float,
                 rep_strength: float, adh_strength: float, patch_area: float,
                 contact_force_cap: float = 5.0e-8,
                 cad_mult: np.ndarray | None = None,
                 # substrate (K2)
                 z0: float = 0.0, W_cs: float = 0.0, adh_range: float = 1.0,
                 k_sub: float | None = None,
                 # active rim traction (vectorized) — optional
                 with_active: bool = False, ranges=None,
                 active: np.ndarray | None = None,
                 int_mult: np.ndarray | None = None, R_cell: float = 7.5e-6,
                 f_act: float = 0.0, f_cap: float = 0.0, ramp_steps: int = 1,
                 contact_band: float = 0.5, neighbour_factor: float = 2.6,
                 max_neighbours: int = 9, belt_factor: float = 0.0,
                 arrest_radius_factor: float | None = None,
                 arrest_width: float = 0.18,
                 arrest_settle_steps: int = 4000) -> None:
        super().__init__(aniso=False)
        self.n_cells = int(n_cells)
        # turgor
        self.faces = np.asarray(faces, dtype=np.int64)
        self.face_cell = np.asarray(face_cell, dtype=np.int64)
        self.V0 = float(V0)
        self.turgor_dP0 = float(turgor_dP0)
        self.K_vol = float(K_vol)
        # tent
        self.cell_of_node = cell_of_node
        self.r_contact = float(r_contact)
        self.c_adh = float(c_adh)
        self.rep = float(rep_strength)
        self.omega = float(adh_strength)
        self.A = float(patch_area)
        self.contact_force_cap = float(contact_force_cap)
        self.cad_mult = cad_mult
        # substrate (well stiffness derived as in DcmSubstrateForceGPU)
        k_well = 2.0 * float(W_cs) / (float(adh_range) ** 2) if adh_range else 0.0
        self.sub_k = (k_well if k_sub is None
                      else (k_well * k_sub) / (k_well + k_sub)) if k_well else 0.0
        self.z0 = float(z0)
        self.sub_range = float(adh_range)
        # active (vectorized law params)
        self.with_active = bool(with_active)
        self.ranges = ranges
        self.active = active
        self.int_mult = int_mult
        self.R = float(R_cell)
        self.f_act = float(f_act)
        self.f_cap = float(f_cap)
        self.ramp_steps = max(1, int(ramp_steps))
        self.contact_band = float(contact_band)
        self.r_neigh = float(neighbour_factor * R_cell)
        self.max_neigh = int(max_neighbours)
        self.belt = float(belt_factor)
        # spreading-arrest (same law as DcmActiveRimTractionGPUVec)
        self.arrest_factor = (None if arrest_radius_factor is None
                              else float(arrest_radius_factor))
        self.arrest_width = float(max(arrest_width, 1.0e-3))
        self.arrest_settle_steps = max(0, int(arrest_settle_steps))
        self._R0_cluster: float | None = None
        # active diagnostics (mirror the standalone vec force)
        self.rim_cells: np.ndarray = np.empty(0, dtype=np.int64)
        self.f_per_cell: dict[int, float] = {}
        self._ramp = 0.0
        self.arrest_gain_mean: float = 1.0
        self.rn_max: float = 0.0
        self._d: DeviceDispatch | None = None

    # arrest helpers (identical law to DcmActiveRimTractionGPU)
    _capture_R0_cluster = DcmActiveRimTractionGPU._capture_R0_cluster
    _arrest_gain = DcmActiveRimTractionGPU._arrest_gain

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def _active_term(self, d, pos_g, timestep: int):
        """Vectorized active rim traction force in GLOBAL (tag) order → (N,3).

        Bit-identical to :meth:`DcmActiveRimTractionGPUVec.set_forces`'s F_g (same
        legacy basal-splay + apical belt law); written here so the fused force pays
        the snapshot/argsort once. Updates ``rim_cells`` / ``f_per_cell``.
        """
        xp = d.xp
        n = int(pos_g.shape[0])
        ranges = self.ranges
        nv = ranges[0][1] - ranges[0][0]
        n_cells = len(ranges)
        ramp = float(min(1.0, timestep / self.ramp_steps))
        self._ramp = ramp
        active_ids = xp.where(xp.asarray(self.active))[0]
        na = int(active_ids.shape[0])
        if na < 2 or ramp <= 0.0:
            self.rim_cells = np.empty(0, dtype=np.int64)
            self.f_per_cell = {}
            return xp.zeros_like(pos_g)

        pos_blk = pos_g.reshape(n_cells, nv, 3)
        act_blk = pos_blk[active_ids]
        cents = act_blk.mean(axis=1)
        cluster_cen = cents.mean(axis=0)
        diff = cents[:, None, :] - cents[None, :, :]
        d2 = xp.sum(diff * diff, axis=2)
        within = d2 < self.r_neigh ** 2
        xp.fill_diagonal(within, False)
        crowd = within.sum(axis=1)
        zc = self.z0 + self.contact_band * self.R
        basal = act_blk[:, :, 2] < zc
        has_basal = basal.any(axis=1)
        is_rim = (crowd <= self.max_neigh) & has_basal
        rxy = cents[:, :2] - cluster_cen[:2]
        rn = xp.sqrt(rxy[:, 0] ** 2 + rxy[:, 1] ** 2)
        contributes = is_rim & (rn >= 1e-12)
        rn_safe = xp.where(rn >= 1e-12, rn, 1.0)
        rhat = xp.zeros((na, 3), dtype=pos_g.dtype)
        rhat[:, 0] = rxy[:, 0] / rn_safe
        rhat[:, 1] = rxy[:, 1] / rn_safe
        # spreading-arrest gain on the per-cell OWN BASAL MEMBRANE RADIUS (same law
        # as the vec force): the global colony-footprint arrest gain (p90 radius of
        # all active basal nodes about their collective centroid → R_max stall).
        R_foot = self._footprint_radius(xp, act_blk, basal)
        self._capture_R0_cluster(xp.reshape(R_foot, (1,)), timestep)
        arrest = self._arrest_gain(xp, R_foot)
        gain = xp.asarray(self.int_mult)[active_ids]
        fmag = xp.minimum(ramp * self.f_act * gain * arrest, self.f_cap)
        fmag = xp.where(contributes, fmag, 0.0)
        fb = xp.minimum(self.belt * fmag, self.f_cap) if self.belt > 0.0 \
            else xp.zeros_like(fmag)
        node_w = xp.where(basal, fmag[:, None], -fb[:, None])
        F_act_blk = node_w[:, :, None] * rhat[:, None, :]
        F_blk = xp.zeros((n_cells, nv, 3), dtype=pos_g.dtype)
        F_blk[active_ids] = F_act_blk
        F_act = F_blk.reshape(n, 3)
        fn = xp.linalg.norm(F_act, axis=1)
        over = fn > self.f_cap
        scale = xp.where(over, self.f_cap / xp.where(fn > 0, fn, 1.0), 1.0)
        F_act = F_act * scale[:, None]
        # diagnostics
        is_rim_h = xp.asnumpy(is_rim) if d.gpu else np.asarray(is_rim)
        ids_h = xp.asnumpy(active_ids) if d.gpu else np.asarray(active_ids)
        contrib_h = xp.asnumpy(contributes) if d.gpu else np.asarray(contributes)
        fmag_h = xp.asnumpy(fmag) if d.gpu else np.asarray(fmag)
        self.rim_cells = ids_h[is_rim_h].astype(np.int64)
        self.f_per_cell = {int(ids_h[i]): float(fmag_h[i])
                           for i in np.where(contrib_h)[0]}
        self.arrest_gain_mean = float(arrest)
        self.rn_max = float(R_foot)
        return F_act

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            pos_g = pos[perm]

            # --- TURGOR (K1) ---
            F_g = d.kernels.mesh_pressure_forces(
                pos_g, xp.asarray(self.faces), xp.asarray(self.face_cell),
                self.n_cells, self.V0, self.turgor_dP0, self.K_vol)

            # --- CELL-CELL TENT (K5) ---
            cad = None if self.cad_mult is None else xp.asarray(self.cad_mult)
            F_g = F_g + d.kernels.tent_contact_forces(
                pos_g, xp.asarray(self.cell_of_node), self.r_contact, self.c_adh,
                self.rep, self.omega, self.A, self.contact_force_cap,
                cad_mult=cad)

            # --- SUBSTRATE (K2) — z-only well, position already global order ---
            if self.sub_k > 0.0:
                F_g = F_g + d.kernels.plane_well_forces(
                    pos_g, self.z0, self.sub_k, self.sub_range)

            # --- ACTIVE rim traction (vectorized) ---
            if self.with_active:
                F_g = F_g + self._active_term(d, pos_g, timestep)

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U
