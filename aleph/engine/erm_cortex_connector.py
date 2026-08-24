r"""Graph-owned ``membrane_erm_cortex`` ERM connector — the SECOND dynamic-runtime KMC participant.

The first dynamic slice (:mod:`aleph.engine.cortex_motor_slice`) drives the NMII actuator + its
``nmii_cortex_motor`` connector through :class:`~aleph.engine.transaction.CellTransaction`, so active
cortical tension EMERGES from head-binding EVENTS.  This module is the next participant by direct analogy: the
ERM (ezrin/radixin/moesin) membrane--cortex clutch.  It is the smallest, highest-value next slice because its
two endpoints — the ``membrane`` and ``cortex`` — are ALREADY native ``SurfaceComponentStateOwner``s
(:mod:`aleph.engine.surface_body`), so it adds NO new native dependency, and it is SIMPLER than NMII:
there is no point-to-segment attach query.  ERM pairs are FIXED by radial pairing at build; a shed tether
rebinds geometrically within a capture radius INSIDE the KMC, so ``propose_events`` is a near-no-op.

:class:`ErmCortexConnector` implements the :class:`~aleph.engine.surface_body.SurfacePairConnector`
protocol.  It OWNS the ERM SoA (the pair indices, the ``bound`` flag, the per-tether ``rest`` length, the RNG
epoch) plus the snapshot twins, and drives the two-array split kernels landed in
:mod:`aleph.components.incumbent.erm_tether`:

* :meth:`accumulate_pair` launches ``erm_bell_force_pair_kernel`` — the unilateral tensile Hookean load
  scattered ``-f`` into the membrane-owned force array and ``+f`` into the cortex-owned force array (Newton's
  3rd law across the two NEVER-merged owner arrays);
* :meth:`snapshot_candidate` / :meth:`rollback` snapshot and reject-gated bit-restore the ``bound`` (int32) +
  ``rest`` (float64) + ``rng_epoch`` (int32) SoA, so a rejected outer step advances no ERM state;
* :meth:`commit_irreversible` launches the accepted-gated ``erm_bell_kmc_pair_kernel`` (Bell-SLIP detach under
  the live tensile load + capture-gated geometric rebind) and the epoch kernel;
* :meth:`propose_events` is a near-no-op (pairs are fixed; rebind lives inside the KMC).

Every CUDA op flows through an injected launcher/copy, so the module is CPU-importable and the structural
gates drive it with recording doubles; the native CUDA validation on the full 70,686-filament cortex is the
Lead's gbook gate (:mod:`aleph.scripts.ac_gate_b_erm_cortex_native`).

Sanity Gate:
    * ownership: the ERM state (pairs / bound / rest / epoch) is OWNED here, never in a surface owner and never
      co-located in a shared array; the tether is the ONLY membrane<->cortex mechanical coupling.
    * boundary/sign: the two-array force is momentum-conserving (membrane ``-f`` + cortex ``+f`` = 0); a
      compressed / broken / resting tether is force-free (a molecular linker is a tether, not a strut).
    * conservation: ``bound``/``rest``/``epoch`` are snapshot and reject-gated bit-restored; the KMC and the
      epoch advance are device-gated on the accepted predicate — a rejected step is a bit-exact no-op with no
      host ``bool()``/``numpy()`` on the predicate.
    * kinetics: ERM is a pure Bell SLIP bond (``bell_kinetics_analytic`` §1) — the off-rate RISES with tensile
      load, so under load the Bell shedding IS the physical bleb-onset (it replaces a hard ``f_rupt`` threshold
      with the mechanistic force-dependent detachment).
"""

from __future__ import annotations

import math
from typing import Callable

import warp as wp

from aleph.components.incumbent.erm_tether import (
    erm_bell_force_pair_kernel,
    erm_bell_kmc_pair_kernel,
    increment_erm_epoch_if_accepted_kernel,
)
from aleph.engine.surface_body import (
    CORTEX_COMPONENT,
    ERM_CONNECTOR,
    MEMBRANE_COMPONENT,
)

__all__ = [
    "ErmCortexConnector",
    "_restore_float64_if_rejected_kernel",
    "_restore_int32_if_rejected_kernel",
]

# ── Injected CUDA seam: one launch/copy indirection keeps the wiring exercisable on a CUDA-free host. ──
LaunchFn = Callable[..., None]
CopyFn = Callable[..., None]


def _default_launch(kernel: object, *, dim: int, inputs: list, outputs: list | None = None,
                    device: object | None = None) -> None:
    """Production launcher: forward to ``wp.launch`` on the owning CUDA device."""
    if outputs is None:
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    else:
        wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


def _default_copy(dst: object, src: object) -> None:
    """Production D2D snapshot/restore: forward to ``wp.copy``."""
    wp.copy(dst, src)


# ── Reject-gated restore kernels for the connector-owned reversible ERM SoA (select the REJECTED branch
# in-device: ``accepted[0]==0`` restores the pre-candidate value; an accepted step keeps the candidate so the
# accepted-gated KMC commit is the only thing that mutates it).  The connector's SoA is not touched by the
# force pass (a pure accumulator), so this reject-gated restore is a bit-exact no-op on an accepted step. ────
@wp.kernel
def _restore_int32_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.int32),
    snap: wp.array(dtype=wp.int32),
) -> None:
    """Restore an int32 array ``live <- snap`` iff the outer step was rejected (device-only, no host branch)."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _restore_float64_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.float64),
    snap: wp.array(dtype=wp.float64),
) -> None:
    """Restore a float64 array ``live <- snap`` iff the outer step was rejected (device-only, no host branch)."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


class ErmCortexConnector:
    """``membrane_erm_cortex`` ERM connector: the two-array Bell tether + accepted-step KMC.

    Implements :class:`~aleph.engine.surface_body.SurfacePairConnector`.  The membrane is
    ``component_a`` and the cortex ``component_b`` (matching the graph edge and the argument order
    :class:`~aleph.engine.surface_body.SurfaceBody.accumulate_mechanics` uses).

    Args:
        membrane_idx_d / cortex_idx_d: ``(Ne,)`` int32 endpoint indices (fixed by radial pairing).  In the
            split-ownership world these index the membrane-owned and cortex-owned arrays respectively; in the
            native driver both surfaces alias one global array and these carry GLOBAL indices (logical split).
        bound_d: ``(Ne,)`` int32 per-tether bound flag (1=intact, 0=shed).
        rest_d: ``(Ne,)`` float64 per-tether rest length [µm] (formation length; force-free at rest).
        rng_epoch_d: ``(1,)`` int32 KMC RNG epoch; advances only on an accepted step.
        bound_snap_d / rest_snap_d / rng_epoch_snap_d: the accepted-step snapshot twins.
        detach_events_d / attach_events_d: ``(1,)`` int32 cumulative accepted Bell detach/attach counters.
        membrane_pos_d / cortex_pos_d: the LIVE membrane/cortex position arrays the accepted KMC reads to
            recompute the tensile load (mutated in place by the inner solve, so always the converged geometry).
        k_erm: ERM linker stiffness [pN/µm] (Braunger 2014, PI-ratified 4.6e3; KB-3.B1.6 pending).
        k_on / k_off0 / bell_force_pn / capture_radius_um: the source-gated Bell on/off contract (all PI-GAP —
            supplied by the caller, never invented here).
        launch / copy / device: injected CUDA seam.
    """

    #: fidelity markers (this is a pure force + Bell-slip KMC clutch — no aggregate/lumped substitute).
    individual_tether_state = True
    bell_slip_kinetics = True
    aggregate_or_two_anchor = False

    def __init__(
        self,
        *,
        membrane_idx_d: wp.array,
        cortex_idx_d: wp.array,
        bound_d: wp.array,
        rest_d: wp.array,
        rng_epoch_d: wp.array,
        bound_snap_d: wp.array,
        rest_snap_d: wp.array,
        rng_epoch_snap_d: wp.array,
        detach_events_d: wp.array,
        attach_events_d: wp.array,
        membrane_pos_d: wp.array,
        cortex_pos_d: wp.array,
        k_erm: float,
        k_on: float,
        k_off0: float,
        bell_force_pn: float,
        capture_radius_um: float,
        name: str = ERM_CONNECTOR,
        component_a: str = MEMBRANE_COMPONENT,
        component_b: str = CORTEX_COMPONENT,
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        if name != ERM_CONNECTOR:
            raise ValueError(f"ERM connector must register as {ERM_CONNECTOR!r}")
        if {component_a, component_b} != {MEMBRANE_COMPONENT, CORTEX_COMPONENT}:
            raise ValueError("ERM connector must join membrane and cortex")
        for label, value in (("k_erm", k_erm), ("k_on", k_on), ("k_off0", k_off0),
                             ("bell_force_pn", bell_force_pn), ("capture_radius_um", capture_radius_um)):
            if value is None or not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"ERM connector {label} must be a supplied positive-finite value; got {value!r}")
        self.name = name
        self.component_a = component_a
        self.component_b = component_b
        self.n_erm = int(bound_d.shape[0])
        self.membrane_idx_d = membrane_idx_d
        self.cortex_idx_d = cortex_idx_d
        self.bound_d = bound_d
        self.rest_d = rest_d
        self.rng_epoch_d = rng_epoch_d
        self.bound_snap_d = bound_snap_d
        self.rest_snap_d = rest_snap_d
        self.rng_epoch_snap_d = rng_epoch_snap_d
        self.detach_events_d = detach_events_d
        self.attach_events_d = attach_events_d
        self._membrane_pos_d = membrane_pos_d
        self._cortex_pos_d = cortex_pos_d
        self.k_erm = float(k_erm)
        self.k_on = float(k_on)
        self.k_off0 = float(k_off0)
        self.bell_force_pn = float(bell_force_pn)
        self.capture_radius_um = float(capture_radius_um)
        self._launch = launch
        self._copy = copy
        self._device = device
        self.launched_kernels: list[object] = []
        self.ledger_calls: list[object] = []

    def accumulate_pair(
        self,
        pos_a: wp.array,
        force_a: wp.array,
        pos_b: wp.array,
        force_b: wp.array,
    ) -> None:
        """Scatter the two-array Bell tether load: membrane node ``-f`` / cortex node ``+f`` (Newton-3rd).

        ``pos_a``/``force_a`` are the membrane-owned arrays and ``pos_b``/``force_b`` the cortex-owned arrays
        (the order :class:`SurfaceBody.accumulate_mechanics` passes).  The kernel mutates neither ``bound`` nor
        ``rest`` — only the accepted-step KMC commit may (solver iterations are not biological time).
        """
        self.launched_kernels = []
        self._launch(
            erm_bell_force_pair_kernel, dim=self.n_erm,
            inputs=[pos_a, force_a, pos_b, force_b, self.membrane_idx_d, self.cortex_idx_d,
                    self.bound_d, self.rest_d, wp.float64(self.k_erm)], device=self._device,
        )
        self.launched_kernels.append(erm_bell_force_pair_kernel)

    def snapshot_candidate(self) -> None:
        """D2D-snapshot the reversible ERM SoA (``bound`` + ``rest`` + ``rng_epoch``): ``snap <- live``."""
        self._copy(self.bound_snap_d, self.bound_d)
        self._copy(self.rest_snap_d, self.rest_d)
        self._copy(self.rng_epoch_snap_d, self.rng_epoch_d)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated bit-restore of the ERM SoA; the accepted branch keeps the candidate for the KMC commit.

        The force pass never mutates the SoA, so on an accepted step this is a bit-exact no-op and
        :meth:`commit_irreversible` is the only writer; on a rejected step the ``bound``/``rest``/``epoch``
        state is restored to the pre-step snapshot, so no shed/rebind or epoch tick survives a rejection.
        """
        self._launch(_restore_int32_if_rejected_kernel, dim=self.n_erm,
                     inputs=[accepted, self.bound_d, self.bound_snap_d], device=self._device)
        self._launch(_restore_float64_if_rejected_kernel, dim=self.n_erm,
                     inputs=[accepted, self.rest_d, self.rest_snap_d], device=self._device)
        self._launch(_restore_int32_if_rejected_kernel, dim=1,
                     inputs=[accepted, self.rng_epoch_d, self.rng_epoch_snap_d], device=self._device)

    def propose_events(self, rates: object, dt_phys: float, rng_seed: int, neighbors: object) -> None:
        """Near-no-op: ERM pairs are FIXED by radial pairing; a shed tether rebinds geometrically-within-capture
        INSIDE the accepted KMC.  Present so the ERM connector is a valid (kinetic) event proposer; it writes no
        force, no host state, and no candidate buffer (there is no point-to-segment attach query, unlike NMII).
        """
        return None

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Accepted-gated ERM Bell KMC on the converged geometry: SLIP detach under load + capture-gated rebind.

        Reads the LIVE (converged) membrane/cortex positions to recompute each tether's tensile load, so the
        Bell off-rate sees the true resting/loaded strain.  Under tensile load the slip off-rate RISES, so the
        most-loaded tethers shed — the mechanistic bleb-onset that replaces a hard ``f_rupt``.  The epoch
        advances only on acceptance.
        """
        self.launched_kernels = []
        seed = wp.int32(int(rng_seed) & 0x7FFFFFFF)
        self._launch(
            erm_bell_kmc_pair_kernel, dim=self.n_erm,
            inputs=[self._membrane_pos_d, self._cortex_pos_d, self.membrane_idx_d, self.cortex_idx_d,
                    self.bound_d, wp.float64(self.k_erm), self.rest_d, wp.float64(self.k_on),
                    wp.float64(self.k_off0), wp.float64(self.bell_force_pn), wp.float64(self.capture_radius_um),
                    wp.float64(dt_phys), seed, self.rng_epoch_d, accepted,
                    self.detach_events_d, self.attach_events_d], device=self._device,
        )
        self.launched_kernels.append(erm_bell_kmc_pair_kernel)
        self._launch(increment_erm_epoch_if_accepted_kernel, dim=1,
                     inputs=[accepted, self.rng_epoch_d], device=self._device)
        self.launched_kernels.append(increment_erm_epoch_if_accepted_kernel)

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle without deciding acceptance (force/topology reduction lands device-side)."""
        self.ledger_calls.append(ledger)
