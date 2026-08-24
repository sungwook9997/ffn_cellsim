r"""I4 weave presets — construct NMII hands via the I3 API + the ``walk_dir`` attach hand-off — I4.

Downstream tracks CONSTRUCT hands via the frozen I3 interface (``ac.motor.hand``) in their OWN ``presets.py``
and never edit the KMC core (§1.4). This module is I4-weave's preset layer:

  1. :func:`nmii_hand_params` — build an ``ac.motor.hand.NMIIHandParams`` (device struct) from the I0-B4/I0-B3
     ledger. Every magnitude is a GAP/provisional value passed IN by the lead at the native gate (never chosen
     here to lift a force band); this factory only assembles them.
  2. :func:`fill_walk_dir_kernel` — the Warp SOURCE of the I4 attach hand-off: after the I3 ``attach_kernel``
     binds a head to an actin node, OVERWRITE ``walk_dir[h]`` from that filament's BARBED-end polarity so the
     power stroke is a directed contraction (INTEGRATION.md §1, HARD). Authored, NEVER launched on the Mac.
  3. :func:`host_fill_walk_dir` — the host-numpy mirror (``ac.weave.walk_dir.fill_walk_dir``) that the CPU gate
     checks bit-for-behaviour against the kernel.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). HOOMD is never imported.
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.motor.hand import NMIIHandParams, allocate_hand_state  # I3 API, read-only
from aleph.components.weave.walk_dir import fill_walk_dir as host_fill_walk_dir

__all__ = ["nmii_hand_params", "allocate_hand_state", "fill_walk_dir_kernel", "host_fill_walk_dir"]


def nmii_hand_params(
    *,
    k_on: float,
    k_off0: float,
    f0: float,
    v0: float,
    f_stall: float,
    kappa: float,
    k_xb: float,
    r0_head: float,
    r0_xb: float = 0.0,
    capture_radius: float,
) -> NMIIHandParams:
    """Assemble an :class:`ac.motor.hand.NMIIHandParams` from lead-supplied I0-B3/I0-B4 values (engine units).

    This is a pure assembler — it CHOOSES nothing. Each argument is an I0-B3 evidence decision (see
    ``ac/motor/params_i0b3.yaml``); passing a GAP value that was tuned to hit a gate violates the hard rule.
    """
    p = NMIIHandParams()
    p.k_on = wp.float64(k_on)
    p.k_off0 = wp.float64(k_off0)
    p.f0 = wp.float64(f0)
    p.v0 = wp.float64(v0)
    p.f_stall = wp.float64(f_stall)
    p.kappa = wp.float64(kappa)
    p.k_xb = wp.float64(k_xb)
    p.r0_head = wp.float64(r0_head)
    p.r0_xb = wp.float64(r0_xb)
    p.capture_radius = wp.float64(capture_radius)
    return p


@wp.kernel
def fill_walk_dir_kernel(
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),           # actin node each bound head grabbed (-1 if free)
    node_fiber: wp.array(dtype=wp.int32),       # (N,) fiber id of every actin node
    barbed_node: wp.array(dtype=wp.int32),      # (F,) global barbed-end node of each fiber
    pos: wp.array(dtype=wp.vec3d),              # (N,3) actin node positions
    walk_dir: wp.array(dtype=wp.vec3d),         # (H,) per-head walk direction (overwritten here)
) -> None:
    """The I4 attach hand-off (device): set a bound head's ``walk_dir`` from its actin filament's barbed end.

    Myosin walks toward the BARBED (+) end, so ``walk_dir = unit(pos[barbed] - pos[anchor])``. A free head
    (bound == 0) is zeroed (passive, the OFF/regression default). This is launched by the lead right AFTER the
    I3 ``attach_kernel`` each KMC tick so the power stroke walks along the ACTUAL actin the head grabbed (not
    the frozen minifilament-axis default). Bit-identical to ``ac.weave.walk_dir.fill_walk_dir`` on CPU.
    """
    h = wp.tid()
    if bound[h] == 0 or anchor[h] < 0:
        walk_dir[h] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        return
    a = anchor[h]
    b = barbed_node[node_fiber[a]]
    v = pos[b] - pos[a]
    n = wp.length(v)
    if n < wp.float64(1e-12):
        walk_dir[h] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        return
    walk_dir[h] = v / n
