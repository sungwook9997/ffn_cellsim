"""Commit-safe ERM membrane-cortex tethers for the Active Cell runtime.

Force evaluation is deliberately pure: an inner quasi-static iteration may discover a loaded tether, but it
may not advance stochastic biology.  Detachment/rebinding is a separate CUDA KMC kernel called once only after
the outer physical transaction is accepted.  This mirrors the nuclear-envelope rupture split and prevents a
rejected mechanical candidate from leaving a bleb behind.

Two explicitly separated modes exist:

* ``erm_tether_force_kernel`` + ``erm_tether_rupture_kernel`` preserve the historical membrane-energy-derived
  hard threshold as a DIAGNOSTIC compatibility path.  That threshold is not a single-molecule ERM kinetic
  datum and therefore cannot authorize production.
* ``erm_bell_force_kernel`` + ``erm_bell_kmc_kernel`` implement the production-form mechanism: a unilateral
  Hookean tether carries its full tensile load during the accepted interval, then detaches with the Bell slip
  probability ``1-exp(-dt*k_off0*exp(F/F0))`` or rebinds within a capture radius with
  ``1-exp(-dt*k_on)``.  A successful rebind sets the rest length to the formation length, so binding never
  injects an unrecorded prestress.

The Bell parameters intentionally have NO biological defaults.  :class:`ERMBellKinetics` validates a complete,
source-labelled contract; the cell production latch refuses partial or absent kinetics.  This lets the full
mechanism land without inventing the currently absent MCF7 rates/capture distance.

Sanity Gate:
    * Force: every intact, extended tether applies equal-and-opposite tensile Hookean forces; a compressed,
      resting, or broken tether applies zero.  A membrane--actin linker is a molecular tether, not a strut.
    * Time: force kernels never write ``bound`` or ``rest``.  Only the accepted outer-step KMC kernel may
      change them; solver iterations are not biological time.
    * Kinetics: Bell off-rate is monotone in tensile force; rejected transactions are exact no-ops; rebind is
      capture-gated and force-free at formation because ``rest := current length``.
    * Legacy threshold: ``f_rupt`` is retained only for diagnostic compatibility and is never relabelled as a
      molecular Bell force or a production rupture law.
    * Residency: both kernels are Warp-CUDA operations; no host state participates in the update.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.motor.hand import attach_prob, bell_off_rate, detach_prob

__all__ = [
    "ERMBellKinetics",
    "erm_tether_force_kernel",
    "erm_tether_rupture_kernel",
    "erm_bell_force_kernel",
    "erm_bell_kmc_kernel",
    "erm_bell_force_pair_kernel",
    "erm_bell_kmc_pair_kernel",
    "erm_bell_force_pair_reference",
    "increment_erm_epoch_if_accepted_kernel",
]


@dataclass(frozen=True, slots=True)
class ERMBellKinetics:
    """Complete source-gated Bell on/off contract for one ERM tether population.

    No field has a biological default.  ``formation_length`` is the only implemented rebind policy: it is
    explicit in the config and its provenance must be covered by ``source`` before production use.
    """

    k_on_s: float
    k_off0_s: float
    bell_force_pn: float
    capture_radius_um: float
    source: str
    rebind_rest_policy: str = "formation_length"

    def __post_init__(self) -> None:
        values = {
            "k_on_s": self.k_on_s,
            "k_off0_s": self.k_off0_s,
            "bell_force_pn": self.bell_force_pn,
            "capture_radius_um": self.capture_radius_um,
        }
        for name, value in values.items():
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"ERM Bell {name} must be finite and positive")
        if not self.source.strip():
            raise ValueError("ERM Bell kinetics require an explicit literature/provenance source")
        if self.rebind_rest_policy != "formation_length":
            raise ValueError("only the explicit ERM rebind_rest_policy='formation_length' is implemented")


@wp.kernel
def erm_tether_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_erm: wp.float64,
    rest: wp.array(dtype=wp.float64),
    f_rupt: wp.float64,
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Accumulate a sub-threshold ERM spring without mutating its irreversible bound state.

    A super-threshold tether contributes no candidate force: the accepted commit will latch its detached
    state, so retaining its Hookean force during convergence would certify a balance that commit invalidates.
    Compression is likewise force-free: ezrin resists membrane--actin separation but cannot push the two
    surfaces apart as a fictitious molecular strut.
    """
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    delta = pos[mi] - pos[ci]
    length = wp.length(delta)
    if length < wp.float64(1.0e-12):
        return
    tension = k_erm * (length - rest[k])
    if tension <= wp.float64(0.0):
        return
    if tension > f_rupt:
        return
    tether_force = tension * delta / length
    wp.atomic_add(force, mi, -tether_force)
    wp.atomic_add(force, ci, tether_force)


@wp.kernel
def erm_bell_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_erm: wp.float64,
    rest: wp.array(dtype=wp.float64),
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Accumulate the full unilateral tensile load for a Bell-kinetic tether without mutating state."""
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    delta = pos[mi] - pos[ci]
    length = wp.length(delta)
    if length < wp.float64(1.0e-12):
        return
    tension = k_erm * (length - rest[k])
    if tension <= wp.float64(0.0):
        return
    tether_force = tension * delta / length
    wp.atomic_add(force, mi, -tether_force)
    wp.atomic_add(force, ci, tether_force)


@wp.kernel
def erm_tether_rupture_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_erm: wp.float64,
    rest: wp.array(dtype=wp.float64),
    f_rupt: wp.float64,
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Commit irreversible ERM rupture once at an accepted outer physical-time boundary."""
    k = wp.tid()
    if accepted[0] == 0 or bound[k] == 0:
        return
    length = wp.length(pos[membrane_idx[k]] - pos[cortex_idx[k]])
    tension = k_erm * (length - rest[k])
    if tension > f_rupt:
        bound[k] = 0


@wp.kernel
def erm_bell_kmc_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_erm: wp.float64,
    rest: wp.array(dtype=wp.float64),
    k_on: wp.float64,
    k_off0: wp.float64,
    bell_force: wp.float64,
    capture_radius: wp.float64,
    tau: wp.float64,
    rng_seed: wp.int32,
    rng_epoch: wp.array(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
    detach_events: wp.array(dtype=wp.int32),
    attach_events: wp.array(dtype=wp.int32),
) -> None:
    """Advance one ERM Bell on/off KMC tick at an accepted outer physical-time boundary.

    Bound tethers detach from their current tensile load.  Compression carries zero force and therefore uses
    the zero-force Bell rate.  Unbound fixed endpoint pairs may rebind only within ``capture_radius``; formation
    is force-free because the current separation becomes the new rest length.  The device ``accepted`` latch
    makes a rejected transaction a bit-exact no-op without a host decision.
    """
    k = wp.tid()
    if accepted[0] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    length = wp.length(pos[mi] - pos[ci])
    state = wp.rand_init(rng_seed ^ rng_epoch[0], k)
    draw = wp.float64(wp.randf(state))
    if bound[k] == 1:
        tension = wp.max(k_erm * (length - rest[k]), wp.float64(0.0))
        p_off = bell_off_rate(tension, k_off0, bell_force)
        if draw < detach_prob(tau, p_off):
            bound[k] = 0
            wp.atomic_add(detach_events, 0, 1)
        return
    if length <= capture_radius and draw < attach_prob(tau, k_on):
        bound[k] = 1
        rest[k] = length
        wp.atomic_add(attach_events, 0, 1)


@wp.kernel
def increment_erm_epoch_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
) -> None:
    """Advance the ERM RNG epoch only when the outer physical transaction commits."""
    if accepted[0] != 0:
        rng_epoch[0] = rng_epoch[0] + 1


# ── Two-array SPLIT-OWNERSHIP ERM (membrane and cortex are NEVER-merged state owners) ─────────────────
# These mirror ``erm_bell_force_kernel`` / ``erm_bell_kmc_kernel`` bit-for-formula, but the membrane node is
# addressed in the membrane-owned ``pos_m``/``force_m`` arrays and the cortex node in the cortex-owned
# ``pos_c``/``force_c`` arrays — the split-ownership relayout the ``ac/engine`` component/connector graph
# requires (co-location in one array is NEVER a connection; the ERM tether is the ONLY coupling).  The
# single-array kernels above are UNCHANGED — the legacy ``MembraneCompartment`` still owns them.  Newton's 3rd
# law holds ACROSS the two arrays: the membrane node gets ``-f`` and the cortex node ``+f`` (their resultant is
# zero), exactly as the single-array kernel loads the two nodes of one combined array.


@wp.kernel
def erm_bell_force_pair_kernel(
    pos_m: wp.array(dtype=wp.vec3d),
    force_m: wp.array(dtype=wp.vec3d),
    pos_c: wp.array(dtype=wp.vec3d),
    force_c: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    k_erm: wp.float64,
) -> None:
    """Two-array unilateral Bell tether force; membrane node ``-f`` / cortex node ``+f`` (Newton-3rd, no merge).

    Identical unilateral Hookean tensile law as :func:`erm_bell_force_kernel` (compression and a broken tether
    are force-free; a molecular linker is a tether, not a strut), but the membrane endpoint indexes the
    membrane-owned ``pos_m``/``force_m`` while the cortex endpoint indexes the cortex-owned ``pos_c``/
    ``force_c``.  The two arrays are addressed independently and never concatenated.  Mutates no ``bound``/
    ``rest`` state — only the accepted-step KMC may (Sanity Gate: solver iterations are not biological time).
    """
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    delta = pos_m[mi] - pos_c[ci]
    length = wp.length(delta)
    if length < wp.float64(1.0e-12):
        return
    tension = k_erm * (length - rest[k])
    if tension <= wp.float64(0.0):
        return
    tether_force = tension * delta / length
    wp.atomic_add(force_m, mi, -tether_force)
    wp.atomic_add(force_c, ci, tether_force)


@wp.kernel
def erm_bell_kmc_pair_kernel(
    pos_m: wp.array(dtype=wp.vec3d),
    pos_c: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_erm: wp.float64,
    rest: wp.array(dtype=wp.float64),
    k_on: wp.float64,
    k_off0: wp.float64,
    bell_force: wp.float64,
    capture_radius: wp.float64,
    tau: wp.float64,
    rng_seed: wp.int32,
    rng_epoch: wp.array(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
    detach_events: wp.array(dtype=wp.int32),
    attach_events: wp.array(dtype=wp.int32),
) -> None:
    """Two-array ERM Bell on/off KMC tick, accepted-gated (same law as :func:`erm_bell_kmc_kernel`, split pos).

    A bound tether detaches from its current tensile load with the Bell SLIP probability
    ``1-exp(-tau*k_off0*exp(F/F0))`` — the correct passive-linker law (ERM is a pure slip bond,
    ``bell_kinetics_analytic`` §1).  An unbound fixed pair may rebind within ``capture_radius`` with
    ``1-exp(-tau*k_on)``; formation is force-free because ``rest := current length`` (binding injects no
    unrecorded prestress).  The membrane endpoint reads ``pos_m``, the cortex endpoint ``pos_c`` — the two
    arrays are never merged.  The device ``accepted`` latch makes a rejected transaction a bit-exact no-op
    without a host decision; the RNG folds ``rng_epoch`` so a re-proposed rejected step re-draws identically.
    """
    k = wp.tid()
    if accepted[0] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    length = wp.length(pos_m[mi] - pos_c[ci])
    state = wp.rand_init(rng_seed ^ rng_epoch[0], k)
    draw = wp.float64(wp.randf(state))
    if bound[k] == 1:
        tension = wp.max(k_erm * (length - rest[k]), wp.float64(0.0))
        p_off = bell_off_rate(tension, k_off0, bell_force)
        if draw < detach_prob(tau, p_off):
            bound[k] = 0
            wp.atomic_add(detach_events, 0, 1)
        return
    if length <= capture_radius and draw < attach_prob(tau, k_on):
        bound[k] = 1
        rest[k] = length
        wp.atomic_add(attach_events, 0, 1)


def erm_bell_force_pair_reference(
    pos_m: npt.NDArray,
    pos_c: npt.NDArray,
    membrane_idx: npt.NDArray,
    cortex_idx: npt.NDArray,
    bound: npt.NDArray,
    rest: npt.NDArray,
    k_erm: float,
) -> dict[str, npt.NDArray]:
    """NumPy host oracle (bit-for-formula twin of :func:`erm_bell_force_pair_kernel`).

    Returns ``{force_m, force_c, tension}``: the two never-merged owner-force arrays each ERM tether scatters
    into (membrane node ``-f``, cortex node ``+f``) and the per-tether tensile magnitude.  The acceptance the
    CPU gate checks: the two scatters are equal-and-opposite (``sum(force_m)+sum(force_c)=0`` — Newton's 3rd
    law across the split arrays) and a resting tether (``rest == length``) is force-free at formation.
    """
    pos_m = np.asarray(pos_m, dtype=np.float64)
    pos_c = np.asarray(pos_c, dtype=np.float64)
    mi = np.asarray(membrane_idx, dtype=np.int64)
    ci = np.asarray(cortex_idx, dtype=np.int64)
    bound = np.asarray(bound, dtype=np.int64)
    rest = np.asarray(rest, dtype=np.float64)
    force_m = np.zeros_like(pos_m, dtype=np.float64)
    force_c = np.zeros_like(pos_c, dtype=np.float64)
    tension = np.zeros(mi.shape[0], dtype=np.float64)
    for k in range(mi.shape[0]):
        if bound[k] == 0:
            continue
        delta = pos_m[mi[k]] - pos_c[ci[k]]
        length = float(np.linalg.norm(delta))
        if length < 1.0e-12:
            continue
        t = float(k_erm) * (length - float(rest[k]))
        if t <= 0.0:
            continue  # unilateral: compression / rest carries no force
        tension[k] = t
        f = t * delta / length
        force_m[mi[k]] += -f
        force_c[ci[k]] += f
    return {"force_m": force_m, "force_c": force_c, "tension": tension}
