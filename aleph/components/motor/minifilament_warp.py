r"""Head-resolved Stam-Hocky NMII minifilament — Warp CUDA force primitive (``MyosinForce.accumulate``).

The production actomyosin motor (P2): an explicit bipolar minifilament of backbone beads + individual heads,
each head a crossbridge spring to its bound actin site, assembled ENTIRELY on the GPU with no per-step host
state (I0-A). This REPLACES the lumped ``ff/network_warp.myosin_kernel`` (constant ``f_myo``) and the aggregate
``ff/myosin_linear.minifilament_kernel`` (a single two-anchor force) — both demoted to Warp diagnostics only
(see INTEGRATION.md; the deletion is a same-commit-as-fine-motor guard against the f_myo+fine-motor
double-count, hard-truth #6).

Topology (ported read-only from the retired HOOMD sources `cortex/myosin.py` + `bridge/motor.py`, which
lived under `archive/hoomd_legacy/` until they were DELETED at `1a9ded66` — read `git show 1a9ded66^:...`
if the derivation is ever challenged; they were never executed;
see minifilament_topology.MinifilamentTopology): ``[backbone_0..n_bb-1, head+_0..H-1, head-_0..H-1]`` with
three bond families:
  * backbone-backbone: stiff harmonic rigid-rod (k_backbone = k_backbone_factor * k_xb), rest = L_bb/(n_bb-1);
  * head-backbone ARM: harmonic (k_head_spring), rest = ``r0_head`` (= topology head_offset_um); rigid rod;
  * head-actin:        DYNAMIC power-stroke crossbridge (k_xb), rest = ``r0_xb`` (~0, SEPARATE from r0_head!);
                       created by the hand.py KMC attach; only for BOUND heads.

⭐ ACTIVE-FORCE ROOT (the P2 fix): the crossbridge is NOT a spring to a fixed anchor — its attachment point on
actin ADVANCES along the head's barbed-end direction ``walk_dir`` by the walked ``abscissa`` (the myosin power
stroke), ``x_att = pos[anchor] + abscissa * walk_dir``. A stepping head therefore stretches its crossbridge and
generates a DIRECTED contractile force that slides actin relative to the backbone. The load fed back to the KMC
is the tangential tension ``k_xb (dot(x_att - pos[head], walk_dir) - r0_xb)`` (:func:`compute_head_loads_kernel`),
which RISES with abscissa; ``step_detach_kernel`` advances abscissa by ``tau * hill_velocity(load)``, so the
walk self-limits when the tension reaches ``f_stall`` (v -> 0). Force-velocity and the ensemble stall EMERGE
from this coupled mechanics<->kinetics loop; nothing is imposed (no ``N_side * F_head``). With ``abscissa = 0``
(or ``walk_dir = 0``) the crossbridge reduces to the passive compliance spring toward the grabbed node.

Force assembly (``@wp.func``/kernels below): each bond adds ``F = k (|r| - r0) r_hat`` to its two endpoints via
``wp.atomic_add`` (Newton's 3rd law inherited). The bipolar minifilament PULLS its two actin anchors together
(contractile, never extensile). ``MyosinForce.accumulate(pos, out_force)`` sums all three families per §1.4
(the lead-owned integrator sums MyosinForce + StericForce + PressureCoupling).

⚠ Runtime: Warp CUDA GPU only (I0-A). Authored on the dev Mac (no CUDA) but NOT launched here; the native
topology/count, per-head Newton, GPU-residency, and ensemble-stall gates run on the lead's gbook A5000 against
the CPU-green analytic oracles in this package.
"""

from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.motor.backbone_warp import (
    angle_harmonic_kernel,
    arm_bending_cfl_stiffness,
    backbone_bending_cfl_stiffness,
)
from aleph.components.motor.hand import NMIIHandParams, allocate_hand_state
from aleph.components.motor.minifilament_topology import MinifilamentTopology

# head-arm orientation rest angle: the lever arm stands perpendicular to the backbone axis (see
# minifilament_topology.head_arm_angle_triples / backbone_warp F6 fix).
ARM_REST_ANGLE = float(np.pi / 2.0)
BACKBONE_REST_ANGLE = float(np.pi)


# ── Internal-bond force kernels (float64, engine units um-pN-s) ──────────────────────────────────
@wp.kernel
def harmonic_bond_kernel(
    pos: wp.array(dtype=wp.vec3d),
    bonds: wp.array(dtype=wp.int32, ndim=2),
    k: wp.float64,
    r0: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Harmonic bond ``F = k (|r| - r0) r_hat`` on each pair ``bonds[t] = (i, j)`` (backbone or head-backbone)."""
    t = wp.tid()
    i = bonds[t, 0]
    j = bonds[t, 1]
    d = pos[j] - pos[i]
    length = wp.length(d)
    if length < wp.float64(1.0e-12):
        return
    fmag = k * (length - r0)
    f = (fmag / length) * d
    wp.atomic_add(force, i, f)     # i pulled toward j when stretched
    wp.atomic_add(force, j, -f)


@wp.func
def attachment_point(
    anchor_pos: wp.vec3d,
    abscissa: wp.float64,
    walk_dir: wp.vec3d,
) -> wp.vec3d:
    """The crossbridge attachment point on actin, ADVANCED by the power stroke.

    A bound head has walked a displacement ``abscissa`` (>= 0) along its actin filament since it grabbed the
    node at ``anchor_pos``; its zero-strain attachment reference therefore sits ``abscissa`` further along the
    barbed-end direction ``walk_dir``::

        x_att = anchor_pos + abscissa * walk_dir

    This is the wire that carries the myosin power stroke into the mechanics: as the head steps (abscissa grows)
    the attachment point moves away from the head, the crossbridge spring stretches, and a DIRECTED contractile
    force appears (:func:`crossbridge_kernel`). At ``abscissa = 0`` (or ``walk_dir = 0``) it collapses to the
    grabbed node — no active component, just the passive compliance spring.
    """
    return anchor_pos + abscissa * walk_dir


@wp.kernel
def crossbridge_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Dynamic head-actin crossbridge WITH the power stroke: ``F = k_xb (|r| - r0_xb) r_hat`` where ``r`` runs
    from the head to the power-stroke-ADVANCED attachment point ``x_att = pos[anchor] + abscissa * walk_dir``.

    ``head_node[h]`` = the head particle's node index; ``anchor[h]`` = its bound actin node index (-1 if free);
    ``abscissa[h]`` = the walked displacement along actin [µm]; ``walk_dir[h]`` = its barbed-end unit direction.
    As the head steps, ``x_att`` advances along ``walk_dir`` and the spring stretches, so the minifilament PULLS
    the head toward the barbed end (+f) and drags the actin node the opposite way (-f) — the DIRECTED, active
    contraction. With ``abscissa = 0`` this reduces to the passive compliance spring toward the grabbed node
    (r0_xb ~ 0, so a freshly-bound head sitting on its actin bears ~0 force — no construction pre-stress).
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    a = anchor[h]
    if a < 0:
        return
    hn = head_node[h]
    x_att = attachment_point(pos[a], abscissa[h], walk_dir[h])
    d = x_att - pos[hn]
    length = wp.length(d)
    if length < wp.float64(1.0e-12):
        return
    fmag = k_xb * (length - r0_xb)
    f = (fmag / length) * d
    wp.atomic_add(force, hn, f)    # head pulled toward the advanced attachment (barbed / walk direction)
    wp.atomic_add(force, a, -f)    # actin dragged the opposite way (Newton 3rd law) — the contractile reaction


@wp.kernel
def compute_head_loads_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    loads: wp.array(dtype=wp.float64),
    loads_full: wp.array(dtype=wp.float64),
):
    r"""Fill the two per-head crossbridge loads the KMC needs — the F4/F5 split.

    ``loads[h]`` (Hill / step): the TANGENTIAL along-actin (walk-direction) component of the crossbridge
    strain — the tension that actually resists the WALK::

        loads[h] = k_xb * ( dot(x_att - pos[head], walk_dir) - r0_xb )
                 = k_xb * ( dot(pos[anchor] - pos[head], walk_dir) + abscissa - r0_xb )

    so the walked ``abscissa`` enters DIRECTLY (d(load)/d(abscissa) = k_xb): the more the head steps, the higher
    the resisting tension. ``step_detach_kernel`` advances abscissa by ``tau·hill_velocity(loads)``, lifting it
    toward ``f_stall`` where ``hill_velocity → 0`` — force-velocity self-limits (nothing imposed).

    ``loads_full[h]`` (Bell / detach): the FULL crossbridge tension magnitude ``k_xb·(|x_att − pos[head]| −
    r0_xb)`` — a crossbridge slips under its TOTAL load, not just the tangential component (Bell 1978 reads the
    bond force magnitude). Splitting the two (F4) stops the walk-tangent load from being (wrongly) reused as the
    Bell force; at collinear isometric stall the transverse component ~0 so the two coincide, but off-axis they
    differ and the Bell slip must see the whole force.
    """
    h = wp.tid()
    if bound[h] == 0:
        loads[h] = wp.float64(0.0)
        loads_full[h] = wp.float64(0.0)
        return
    a = anchor[h]
    if a < 0:
        loads[h] = wp.float64(0.0)
        loads_full[h] = wp.float64(0.0)
        return
    x_att = attachment_point(pos[a], abscissa[h], walk_dir[h])
    d = x_att - pos[head_node[h]]
    loads[h] = k_xb * (wp.dot(d, walk_dir[h]) - r0_xb)      # tangential → Hill
    loads_full[h] = k_xb * (wp.length(d) - r0_xb)           # full |F| → Bell slip


# ── Host-side topology construction (numpy; uploaded to device once) ─────────────────────────────
def _as_int32_2d(arr: wp.array | np.ndarray | None, device: str | None):
    """Return ``arr`` as a device int32 (N,2)/(N,3) Warp array (or ``None``); pass Warp arrays through."""
    if arr is None:
        return None
    if isinstance(arr, wp.array):
        return arr
    a = np.ascontiguousarray(np.asarray(arr, dtype=np.int32))
    if a.shape[0] == 0:
        return None
    return wp.array(a, dtype=wp.int32, device=device)


def build_minifilament_nodes(
    topo: MinifilamentTopology,
    centre: np.ndarray,
    axis: np.ndarray,
) -> dict[str, np.ndarray]:
    """Build ONE minifilament's node positions + internal bond index arrays (host numpy, for device upload).

    Args:
        topo: the bipolar topology (counts + geometry; counts are I0-B3 GAPs, invariants are not).
        centre: minifilament centre ``(3,)`` [µm] in the cell frame.
        axis: unit backbone axis ``(3,)``.

    Returns:
        dict of numpy arrays: ``positions`` (n_particles, 3), ``backbone_bonds`` (n_bb-1, 2),
        ``head_bonds`` (2H, 2) head<->backbone, ``head_node`` (2H,) the head particle indices,
        ``walk_dir`` (2H, 3) the DEFAULT per-head barbed-end/walk direction (unit),
        ``backbone_angles`` (n_bb-2, 3) consecutive backbone bending triples (F6),
        ``head_arm_angles`` (2H, 3) head-arm orientation triples (F6). All bond/triple indices are LOCAL to
        this minifilament (0-based); the caller offsets them to global node indices.

    ``walk_dir`` default is derived purely from the bipolar topology geometry: a myosin head walks toward the
    actin barbed end, which in a bipolar minifilament sits at the periphery (the two anti-parallel arms pull
    their actins toward the minifilament centre — contraction). So each head's walk direction points OUTWARD
    along the backbone axis away from the centre: the ``+`` arm (first H heads, on the ``-axis`` half) walks
    ``-ax``, the ``-`` arm walks ``+ax``; the resulting crossbridge reaction (``-walk_dir``) drags both actins
    inward -> contractile. This is a geometry default only — the production engine OVERWRITES ``walk_dir`` per
    head from the bound actin filament's true barbed-end polarity (I4-weave hand-off, INTEGRATION.md §1).
    """
    local = topo.positions()                      # axis along +x, centred at origin
    ax = axis / (np.linalg.norm(axis) + 1e-30)
    # rotate local +x onto `ax` (Rodrigues); for the source spec a simple basis map suffices
    x_hat = np.array([1.0, 0.0, 0.0])
    v = np.cross(x_hat, ax)
    s = np.linalg.norm(v)
    c = float(np.dot(x_hat, ax))
    if s < 1e-12:
        rot = np.eye(3) if c > 0 else np.diag([-1.0, -1.0, 1.0])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rot = np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))
    positions = (local @ rot.T) + centre
    backbone_bonds = np.column_stack([np.arange(topo.n_bb - 1), np.arange(1, topo.n_bb)]).astype(np.int32)
    head_bonds = topo.head_backbone_bonds().astype(np.int32)
    head_node = np.arange(topo.n_bb, topo.n_particles, dtype=np.int32)
    h_per_side = topo.n_heads_per_side
    walk_dir = np.zeros((topo.n_heads, 3), dtype=np.float64)
    walk_dir[:h_per_side] = -ax                   # + arm (first H, -axis half): barbed end is outward = -ax
    walk_dir[h_per_side:] = +ax                   # - arm (next H, +axis half):  barbed end is outward = +ax
    return {
        "positions": positions,
        "backbone_bonds": backbone_bonds,
        "head_bonds": head_bonds,
        "head_node": head_node,
        "walk_dir": walk_dir,
        "backbone_angles": topo.backbone_angle_triples().astype(np.int32),
        "head_arm_angles": topo.head_arm_angle_triples().astype(np.int32),
    }


class MyosinForce:
    """Head-resolved NMII force primitive — owns ``accumulate(pos, out_force)`` per §1.4.

    Holds the device topology (backbone + head-backbone bonds, head node indices), the per-head hand state
    (from ``hand.allocate_hand_state``), and the ``NMIIHandParams`` + rigid-rod stiffnesses. The lead-owned
    integrator calls :meth:`accumulate` to add the minifilament force to the global node force, then
    :meth:`step_kinetics` to advance the KMC once the mechanics are converged.

    This is a device-runtime SOURCE (authored, not launched on the dev Mac). Construction inputs (bond arrays,
    head-node map, params) are built host-side then uploaded once; there is no authoritative per-step host
    state thereafter.
    """

    def __init__(
        self,
        backbone_bonds: wp.array,
        head_bonds: wp.array,
        head_node: wp.array,
        params: NMIIHandParams,
        k_backbone: float,
        r0_backbone: float,
        k_head_spring: float,
        n_hands: int,
        walk_dir: wp.array | np.ndarray | None = None,
        device: str | None = None,
        backbone_angles: wp.array | np.ndarray | None = None,
        head_arm_angles: wp.array | np.ndarray | None = None,
        k_theta_backbone: float = 0.0,
        k_theta_arm: float = 0.0,
        segment_len_um: float = 0.0,
    ) -> None:
        self.backbone_bonds = backbone_bonds
        self.head_bonds = head_bonds
        self.head_node = head_node
        self.params = params
        self.k_backbone = wp.float64(k_backbone)
        self.r0_backbone = wp.float64(r0_backbone)
        self.k_head_spring = wp.float64(k_head_spring)   # head<->backbone ARM stiffness; rest = params.r0_head
        self.state = allocate_hand_state(n_hands, device=device)
        self.loads = wp.zeros(n_hands, dtype=wp.float64, device=device)         # tangential → Hill (step)
        self.loads_full = wp.zeros(n_hands, dtype=wp.float64, device=device)    # full |F| → Bell (detach), F4/F5
        self.segment_runtime = None
        if walk_dir is not None:
            # per-head barbed-end/walk direction (default geometry from build_minifilament_nodes; the production
            # engine overwrites it from actin polarity). Zero (allocate default) => passive head, no directed force.
            self.state["walk_dir"] = wp.array(np.asarray(walk_dir, dtype=np.float64),
                                              dtype=wp.vec3d, device=device)

        # ── F6 angle-harmonic bending (backbone rigid rod + head-arm orientation). Optional: when the triples
        # and k_theta are supplied the head can no longer swing free / the rod cannot bow, so the crossbridge
        # tension transmits collinearly (transmission → 1). Absent (defaults) => the pre-fix distance-only
        # topology, kept for backward compatibility. ──────────────────────────────────────────────────────
        self.k_theta_backbone = wp.float64(k_theta_backbone)
        self.k_theta_arm = wp.float64(k_theta_arm)
        self.backbone_angles = _as_int32_2d(backbone_angles, device)
        self.head_arm_angles = _as_int32_2d(head_arm_angles, device)
        # CFL: bending adds an effective LINEAR stiffness folded into the composed kmax (Round-2 boot note).
        cfl = max(float(k_backbone), float(k_head_spring), float(params.k_xb))
        if self.backbone_angles is not None and k_theta_backbone > 0.0 and segment_len_um > 0.0:
            cfl = max(cfl, backbone_bending_cfl_stiffness(float(k_theta_backbone), float(segment_len_um)))
        if self.head_arm_angles is not None and k_theta_arm > 0.0:
            r0h = float(params.r0_head) if float(params.r0_head) > 0.0 else 1.0
            cfl = max(cfl, arm_bending_cfl_stiffness(float(k_theta_arm), r0h))
        self.cfl_stiffness = cfl

    def enable_segment_runtime(
        self,
        seg_node_a: wp.array,
        seg_node_b: wp.array,
        seg_polarity: wp.array,
        max_segment_length_um: float,
        *,
        grid_dim: int = 128,
        device: str | None = None,
    ) -> None:
        """Switch this composed motor from diagnostic node anchors to production segment anchors."""
        from aleph.components.motor.segment_motor import SegmentMotorRuntime

        runtime = SegmentMotorRuntime(
            self.head_node,
            self.params,
            seg_node_a,
            seg_node_b,
            seg_polarity,
            max_segment_length_um,
            grid_dim=grid_dim,
            device=device,
        )
        self.segment_runtime = runtime
        # Preserve the public state/load handles used by diagnostics and ledgers, but point them at the
        # authoritative production segment state once the composed cell opts in.
        self.state = runtime.state
        self.loads = runtime.loads
        self.loads_full = runtime.loads_full

    def accumulate(self, pos: wp.array, out_force: wp.array) -> None:
        """Add the minifilament's backbone + head-backbone-arm + power-stroke crossbridge + F6 bending forces.

        The crossbridge is the ACTIVE term: its attachment point on actin is advanced along ``walk_dir`` by the
        walked ``abscissa`` (:func:`crossbridge_kernel`), so a stepping head generates a directed contractile
        force. Backbone-backbone (rigid rod) and head<->backbone (arm, rest ``params.r0_head`` == head offset)
        are the passive structure; the crossbridge rest is the SEPARATE ``params.r0_xb`` (~0). The F6 angle-
        harmonic bending (backbone rest π, head-arm rest π/2) gives the rod and arms the transverse stiffness
        the distance springs lack, so the tension transmits collinearly (see backbone_warp). (§1.4 contract.)
        """
        wp.launch(harmonic_bond_kernel, dim=self.backbone_bonds.shape[0],
                  inputs=[pos, self.backbone_bonds, self.k_backbone, self.r0_backbone], outputs=[out_force])
        wp.launch(harmonic_bond_kernel, dim=self.head_bonds.shape[0],
                  inputs=[pos, self.head_bonds, self.k_head_spring, self.params.r0_head], outputs=[out_force])
        if self.segment_runtime is None:
            wp.launch(crossbridge_kernel, dim=self.state["bound"].shape[0],
                      inputs=[pos, self.head_node, self.state["bound"], self.state["anchor"],
                              self.state["abscissa"], self.state["walk_dir"],
                              self.params.k_xb, self.params.r0_xb], outputs=[out_force])
        else:
            self.segment_runtime.accumulate(pos, out_force)
        if self.backbone_angles is not None and float(self.k_theta_backbone) > 0.0:
            wp.launch(angle_harmonic_kernel, dim=self.backbone_angles.shape[0],
                      inputs=[pos, self.backbone_angles, self.k_theta_backbone,
                              wp.float64(BACKBONE_REST_ANGLE)], outputs=[out_force])
        if self.head_arm_angles is not None and float(self.k_theta_arm) > 0.0:
            wp.launch(angle_harmonic_kernel, dim=self.head_arm_angles.shape[0],
                      inputs=[pos, self.head_arm_angles, self.k_theta_arm,
                              wp.float64(ARM_REST_ANGLE)], outputs=[out_force])

    def compute_loads(self, pos: wp.array) -> None:
        """Refresh the per-head crossbridge loads from the current (converged) geometry — the F4/F5 split.

        Fills ``loads`` (tangential, → Hill/step) AND ``loads_full`` (full |F|, → Bell/detach). Reads the walked
        ``abscissa`` + ``walk_dir`` so the tangential load rises as the head steps — the mechanics<->kinetics
        coupling that makes force-velocity / ensemble stall emergent (see compute_head_loads_kernel).
        """
        if self.segment_runtime is None:
            wp.launch(compute_head_loads_kernel, dim=self.state["bound"].shape[0],
                      inputs=[pos, self.head_node, self.state["bound"], self.state["anchor"],
                              self.state["abscissa"], self.state["walk_dir"],
                              self.params.k_xb, self.params.r0_xb], outputs=[self.loads, self.loads_full])
        else:
            self.segment_runtime.compute_loads(pos)

    def commit_segment_kinetics(
        self,
        pos: wp.array,
        accepted: wp.array,
        tau: float,
        rng_seed: int,
    ) -> None:
        """Commit one production segment KMC tick using the scheduler-owned acceptance scalar."""
        if self.segment_runtime is None:
            raise RuntimeError("production segment runtime is not configured")
        self.segment_runtime.commit_kinetics(pos, accepted, tau, rng_seed)

    def step_kinetics(
        self,
        nearest_dist: wp.array,
        nearest_idx: wp.array,
        tau: float,
        rng_seed: int,
    ) -> None:
        """Advance one KMC tick on device: attach (free heads in range) then step+detach (bound heads).

        Uses the freshly computed :attr:`loads` (call :meth:`compute_loads` after the mechanical solve first).
        Offsets the RNG seed between the two sub-steps so their per-head streams are independent.
        """
        if self.segment_runtime is not None:
            raise RuntimeError("segment motor kinetics must commit through commit_segment_kinetics")
        from aleph.components.motor.hand import attach_kernel, step_detach_kernel

        n = self.state["bound"].shape[0]
        seed_attach = wp.int32(rng_seed & 0x7FFFFFFF)
        seed_detach = wp.int32((rng_seed ^ 0x5BD1E995) & 0x7FFFFFFF)  # decorrelate the two sub-step streams
        wp.launch(attach_kernel, dim=n,
                  inputs=[self.state["bound"], self.state["anchor"], self.state["abscissa"],
                          nearest_dist, nearest_idx, self.params, wp.float64(tau), seed_attach])
        # F4/F5 split: Hill/step reads the tangential load, Bell/detach reads the full |F| load.
        wp.launch(step_detach_kernel, dim=n,
                  inputs=[self.state["bound"], self.state["anchor"], self.state["abscissa"],
                          self.loads, self.loads_full, self.params, wp.float64(tau), seed_detach])


__all__ = [
    "harmonic_bond_kernel",
    "attachment_point",
    "crossbridge_kernel",
    "compute_head_loads_kernel",
    "build_minifilament_nodes",
    "MyosinForce",
    "ARM_REST_ANGLE",
    "BACKBONE_REST_ANGLE",
]
