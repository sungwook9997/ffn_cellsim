r"""Topology-reforming crosslink KMC — Warp CUDA kernel SOURCE (I4).

Device realization of the reattach KMC (``ac.weave.crosslink_kmc`` is the host oracle). Each crosslinker end is
a device state (bound flag + partner node + load); per tick a bound end Bell-detaches and a free end rebinds to
a NEW nearest cross-fiber node found through the device neighbour search — so the network topology REFORMS on
the GPU with no authoritative host state (I0-A). The bound/detach probabilities REUSE the I3 hand device
functions (``ac.motor.hand.attach_prob`` / ``detach_prob`` / ``bell_off_rate``) — consumed READ-ONLY per the
§1.4 frozen interface; this track never edits the KMC core.

The nearest-cross-fiber-partner query here is a PLACEHOLDER signature (``nearest_dist`` / ``nearest_idx`` are
filled by a device ``wp.HashGrid`` pass the lead owns, mirroring ``ac.solid`` / ``dcm_neighbor_warp``); the
same-fiber exclusion is enforced by the caller when building the candidate list. Authored on the dev Mac,
NEVER launched here (no CUDA); the lead launches it on the gbook A5000.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). HOOMD is never imported.

Sanity Gate (native, lead): device bound fraction -> k_on/(k_on+k_off0) unloaded (matches host oracle); the
bound-pair set reforms while the count is steady; crosslinker count conserved; OFF (k_on=k_off0=0) is
bit-identical to a static-crosslink regression.
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.motor.hand import attach_prob, bell_off_rate, detach_prob  # I3 API, read-only

__all__ = ["crosslink_detach_kernel", "crosslink_reattach_kernel"]


@wp.kernel
def crosslink_detach_kernel(
    bound: wp.array(dtype=wp.int32),
    partner: wp.array(dtype=wp.int32),
    loads: wp.array(dtype=wp.float64),
    k_off0: wp.float64,
    f0: wp.float64,
    tau: wp.float64,
    rng_seed: wp.int32,
) -> None:
    """Bell-slip detach pass: each bound crosslinker end detaches with prob ``1 - exp(-tau*p_off(load))``.

    Uses the I3 hand device functions bit-for-formula (``bell_off_rate`` + ``detach_prob``). Device RNG per
    end (no host draw); a detached end sets ``partner = -1`` and waits for the reattach pass.
    """
    e = wp.tid()
    if bound[e] == 0:
        return
    p_off = bell_off_rate(loads[e], k_off0, f0)
    pdet = detach_prob(tau, p_off)
    state = wp.rand_init(rng_seed, e)
    if wp.float64(wp.randf(state)) < pdet:
        bound[e] = 0
        partner[e] = -1


@wp.kernel
def crosslink_reattach_kernel(
    bound: wp.array(dtype=wp.int32),
    partner: wp.array(dtype=wp.int32),
    nearest_dist: wp.array(dtype=wp.float64),   # (E,) nearest cross-fiber node distance (device HashGrid pass)
    nearest_idx: wp.array(dtype=wp.int32),      # (E,) that node index (already same-fiber-excluded)
    reach: wp.float64,
    k_on: wp.float64,
    tau: wp.float64,
    rng_seed: wp.int32,
) -> None:
    """Reattach pass: each free end within ``reach`` of a NEW cross-fiber node rebinds with ``attach_prob``.

    The rebind partner (``nearest_idx``) may DIFFER from the pre-detach partner -> the topology reforms.
    Crosslinker ends are never created or destroyed (count conserved); the same-fiber exclusion is applied
    when the neighbour list is built (mirrors the host oracle's ``fiber_of_node`` filter).
    """
    e = wp.tid()
    if bound[e] == 1:
        return
    if nearest_dist[e] > reach:
        return
    state = wp.rand_init(rng_seed + wp.int32(1), e)
    if wp.float64(wp.randf(state)) < attach_prob(tau, k_on):
        bound[e] = 1
        partner[e] = nearest_idx[e]
