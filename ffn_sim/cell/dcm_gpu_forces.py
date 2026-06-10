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
                 integrin_switch_gain: float, belt_factor: float) -> None:
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
        # diagnostics (read by the driver) — kept on host (numpy) like the CPU class
        self.rim_cells: np.ndarray = np.empty(0, dtype=np.int64)
        self.f_per_cell: dict[int, float] = {}
        self._ramp = 0.0
        self._d: DeviceDispatch | None = None      # lazy (device known at run(0))

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

            int_mult = xp.asarray(self.int_mult)
            rim_list: list[int] = []
            f_per_cell: dict[int, float] = {}
            zc = self.z0 + self.contact_band * self.R
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
                rhat = xp.asarray([rxy[0] / rn, rxy[1] / rn, 0.0])
                gain = float(int_mult[c])
                fmag = ramp * self.f_act * gain
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

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U
