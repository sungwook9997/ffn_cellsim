r"""Implicit preconditioned-CG inner solve for the ``sf_arc`` + ``nmii`` motor slice.

WHAT THIS CLOSES.  The GATE-B SF-motor slice (:mod:`aleph.engine.sf_motor_slice`) established the
MECHANISM natively — myosin heads bind by ``k_on`` events and the SF axial tension plus the inward FA traction
emerge from zero — but its QUANTITATIVE claim was BLOCKED because the driver's inner solve was an explicit
overdamped relax that does not converge: the free-node residual sat at ~15% of the reported tension, and a 10x
longer relax moved ``T_max`` by 3.4x at an identical bound population (``outputs/ac/gate_b_sf_motor/REPORT.md``).
Every tension magnitude was therefore a relaxation TRANSIENT.

The cause is structural, not a tuning matter.  The SOURCED α-actinin dorsal↔arc crosslink (4.6e5 pN/µm, several
of them meeting at one arc-apex node) dominates the Gershgorin bound, so λ_max ≈ 3.7e6 pN/µm and the CFL-stable
explicit step is ~1.4e-7 µm/pN, while the SF axial mode relaxes orders of magnitude slower.  No number of
explicit iterations fixes that ratio — the fix is an implicit solve, exactly as the cortex needed at GATE A.

WHAT IT REUSES, AND WHY THAT IS THE POINT.  :class:`aleph.components.incumbent.implicit_mechanics.ProjectedAnalyticCG` is
the project's audited Warp-CUDA PCG: its device-resident recurrence (dot products, axpy updates, the
convergence/curvature/breakdown latches, the fixed-launch schedule with no per-iteration host physics readback)
touches ``cell`` NOWHERE — ``cell`` enters only through ``_project`` / ``_stiffness`` /
``_build_preconditioner`` and the flag-guarded fiber-block / multigrid / fiber-quotient branches.  This subclass
therefore overrides exactly those three plus ``__init__``, turns every cortex-only branch off, and INHERITS
``solve`` / ``_operator`` / ``_precondition`` verbatim.  ``ac/cell`` is not modified (it is feature-frozen);
reusing its kernels from the engine is the documented design (CLAUDE.md: the engine "OWNS/orchestrates and
progressively binds the Warp physics kernels that live in ``ac/cell/`` … ``ac/motor/`` … and ``ff/``").

**No new physics is introduced here.**  Every stiffness term is an existing kernel and each is the tangent of a
force the slice already launches:

===========================================  =========================================================
force the slice launches                     tangent this solver launches
===========================================  =========================================================
``ff.network_warp.link_spring_kernel``        ``add_pair_stiffness_kernel`` (SF axial links)
   (same kernel, on ``arc_joints``)           ``add_pair_stiffness_kernel`` (α-actinin dorsal↔arc)
``ff.forces_warp.cytosim_bending_kernel``     ``add_bending_stiffness_kernel`` (exact; NF2007 bending is linear)
``ac.motor.minifilament_warp.harmonic_bond``  ``add_uniform_pair_stiffness_kernel`` (NMII backbone rod, head arm)
``ac.motor.backbone_warp.angle_harmonic``     ``add_angle_gauss_newton_stiffness_kernel`` (backbone π, arm π/2)
``crossbridge_segment_split_r0bind_kernel``   ``add_segment_crossbridge_stiffness_kernel``
===========================================  =========================================================

The crossbridge tangent is EXACT for the ``r0bind`` force variant: that force is
``f = k_xb·(⟨x_att − x_head, ŵ⟩ − r0_bind[h])·ŵ``, so ``∂f/∂x = k_xb·ŵŵᵀ`` independently of ``r0_bind`` — the
per-head zero-strain reference shifts the force but not its derivative.  The two angle terms are PSD
Gauss-Newton approximations rather than exact Hessians; that is safe because the solve is ``A dx = F`` and at
the fixed point ``F = 0 ⇒ dx = 0`` for ANY non-singular SPD ``A``.  An inexact ``A`` changes the convergence
RATE, never the equilibrium the gate measures.

THE PROJECTOR IS A DIRICHLET MASK, NOT INEXTENSIBILITY — the one thing that must not be inherited.
``ProjectedAnalyticCG``'s ``P`` is the exact NF2007 inextensibility projector ``P = I − Jᵀ(JJᵀ)⁻¹J`` over
per-segment squared-length constraints (``ac/cell/inner_mechanics.py::project_constraint_forces_kernel``).
**The SF passive model has no such constraint**: it carries a Hookean axial spring whose stiffness
``k_axial_pn_per_um`` is an explicitly documented modelling GAP precisely BECAUSE NF2007 treats the backbone as
inextensible and rejects an axial penalty spring (:mod:`aleph.engine.sf_mechanics`).  Inheriting ``P``
would annihilate the axial component of the residual while leaving the axial stiffness inside ``K`` — the solve
would converge to a state whose SF axial tension is never equilibrated, i.e. WORSE than the explicit relax it
replaces, and the tension diagnostic (which reads ``k·(L−r0)``) would be measuring a frozen length.

What SF does need projected is the **FA anchor Dirichlet condition**.  A ventral stress fiber is a closed
contractile dipole only because its outer ends are held; the driver currently enforces that by zeroing the
FORCE on pinned nodes, which is sufficient for an explicit step (``x += step·F`` with ``F = 0`` leaves ``x``
fixed) but NOT for a coupled implicit solve: with a masked right-hand side alone, ``dx`` at a pinned node of
``(aI + K)dx = b`` is still non-zero through the stiffness coupling to its neighbours, so the anchors would
drift and the FA-traction sign/magnitude diagnostics would become meaningless.  Substituting the 0/1 mask for
``P`` makes the inherited operator ``M(aI+K)M + a(I−M)``, which is SPD, non-singular on the pinned block, and
solves to exactly ``dx_pinned = 0``.  The "projected" architecture composes with a mask unchanged — only the
projector implementation is swapped.

THE REGULARIZER IS DERIVED, NOT TUNED (Magic-Number Block).  ``a = regularization[0]`` is set by the caller to
``1/(dt/γ)`` — the reciprocal of the mobility step the driver already derives from the assembled stiffnesses
via its CFL bound.  Then ``aI + K = (γ/dt)I + K`` is literally the IMEX overdamped step, it reduces to the
current explicit update as ``K → 0``, and it is what makes the operator non-singular for the parts of this
slice that are genuinely unanchored: ``TRANSVERSE_ARC`` bundles carry no FA/LINC anchor at all, the
``PERINUCLEAR_CAP`` LINC sites are not pinned by the driver, and every unbound minifilament is a free rigid
body.  It is a mobility the gate already computed, not a constant chosen to make anything pass.  The parent's
``omitted_regularization_base`` is deliberately NOT used — it reads cortex/membrane/nucleus fields this slice
does not have.

NO MULTIGRID.  The fiber-arclength V-cycle that unblocked the cortex is deliberately not used here.  Its
line smoother approximates each fiber by ``D_external + D2ᵀ diag(α) D2`` — external diagonal plus bending only
— which presumes the intra-fiber AXIAL operator is absent because it is a hard constraint.  That presumption is
exactly what SF violates, so the property the V-cycle relies on (near-zero row sums annihilating the fiber
translation mode) does not hold, and at this population (a few hundred nodes, ≤9 nodes per fiber) the hierarchy
would collapse after one level anyway.  Node-Jacobi preconditioning is the right tool at this size, and it is
also the direct fix for the stiffness contrast that broke the explicit path: each node gets its own stiffness
scale, so the 4.6e5 α-actinin node no longer dictates the step for the 1e3 axial nodes.

SPLIT OWNERSHIP IS PRESERVED.  ``sf_arc`` and ``nmii`` own two never-merged device arrays, and the landed
gates assert that (``sf_motor_slice`` Sanity Gate: "the minifilament particle array is NOT the SF node array").
The crossbridge tangent kernel, like its force counterpart, needs head and segment nodes addressable in one
index space, so this solver allocates its OWN concatenated scratch ``[SF nodes | NMII particles]`` and
gathers/scatters into it.  The two authoritative state arrays are never aliased or merged — only solver scratch
is contiguous, the same discipline the split crossbridge force already follows.

engine units: length µm, force pN, stiffness pN/µm.  Runtime: NVIDIA Warp on CUDA only (I0-A).  The module is
CPU-importable (kernels JIT lazily); :meth:`SFImplicitCG.__init__` allocates device memory and is a gbook lane.

Sanity Gate (self-tested in tests/ac/engine/test_sf_implicit.py):
  * operator/force consistency: for each family, the launched tangent matches a finite-difference of the force
    the slice launches (pure-NumPy references, no CPU-Warp launch).
  * SPD: ``vᵀ(aI + MKM + a(I−M))v > 0`` for random ``v`` with ``a > 0``; every family's tangent is PSD.
  * Dirichlet: a pinned node's ``dx`` is exactly zero, and masking the RHS alone provably is NOT sufficient
    (the reference solve shows non-zero pinned displacement without the mask in the operator).
  * ownership: the two state arrays are never aliased; only solver scratch is concatenated (pointer check).
  * derived-not-tuned: the regularizer equals the caller's CFL mobility reciprocal; no literal appears here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

# Reused audited tangents + preconditioner diagonals + the CG machinery (ac/cell is NOT modified; CLAUDE.md
# designates these kernels as engine-bindable physics).  Import order mirrors the force launch order.
from aleph.components.incumbent.implicit_mechanics import (
    ProjectedAnalyticCG,
    add_angle_gauss_newton_preconditioner_kernel,
    add_angle_gauss_newton_stiffness_kernel,
    add_bending_preconditioner_kernel,
    add_bending_stiffness_kernel,
    add_pair_preconditioner_kernel,
    add_pair_stiffness_kernel,
    add_segment_crossbridge_preconditioner_kernel,
    add_segment_crossbridge_stiffness_kernel,
    add_uniform_pair_preconditioner_kernel,
    add_uniform_pair_stiffness_kernel,
    set_mass_action_kernel,
    set_preconditioner_kernel,
)

__all__ = [
    "SFImplicitCG",
    "SFImplicitTopology",
    "gather_two_kernel",
    "mask_vec3_kernel",
    "scatter_two_kernel",
]

#: Rest angles of the NMII internal angle-harmonics (backbone straight rod, head arm perpendicular).  These are
#: the SAME constants the force kernels use — imported rather than re-literal'd so they cannot drift.
from aleph.components.motor.minifilament_warp import ARM_REST_ANGLE, BACKBONE_REST_ANGLE


# ── scratch gather/scatter + the Dirichlet mask (the only new kernels; no physics in them) ─────────────────
@wp.kernel
def gather_two_kernel(
    first: wp.array(dtype=wp.vec3d),
    second: wp.array(dtype=wp.vec3d),
    n_first: wp.int32,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Gather two never-merged owner arrays into the solver's contiguous scratch ``[first | second]``."""
    i = wp.tid()
    if i < n_first:
        out[i] = first[i]
    else:
        out[i] = second[i - n_first]


@wp.kernel
def scatter_two_kernel(
    source: wp.array(dtype=wp.vec3d),
    n_first: wp.int32,
    first: wp.array(dtype=wp.vec3d),
    second: wp.array(dtype=wp.vec3d),
) -> None:
    """Scatter a contiguous scratch vector back onto the two owner arrays (``+=`` is the caller's business)."""
    i = wp.tid()
    if i < n_first:
        first[i] = source[i]
    else:
        second[i - n_first] = source[i]


@wp.kernel
def add_scatter_two_kernel(
    source: wp.array(dtype=wp.vec3d),
    n_first: wp.int32,
    first: wp.array(dtype=wp.vec3d),
    second: wp.array(dtype=wp.vec3d),
) -> None:
    """Add a contiguous scratch displacement onto the two owner position arrays (the implicit step apply)."""
    i = wp.tid()
    if i < n_first:
        first[i] = first[i] + source[i]
    else:
        second[i - n_first] = second[i - n_first] + source[i]


@wp.kernel
def mask_vec3_kernel(
    source: wp.array(dtype=wp.vec3d),
    pinned: wp.array(dtype=wp.int32),
    destination: wp.array(dtype=wp.vec3d),
) -> None:
    """Dirichlet projector ``M``: copy ``source`` but zero every pinned (held) degree of freedom.

    This replaces the parent's NF2007 inextensibility projector.  Applied to the operator input AND output (the
    inherited ``_operator`` does both), it turns the operator into ``M(aI+K)M + a(I−M)`` — SPD, non-singular on
    the pinned block, and with a masked residual it solves to exactly ``dx = 0`` there.
    """
    i = wp.tid()
    if pinned[i] != wp.int32(0):
        destination[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    else:
        destination[i] = source[i]


@dataclass(frozen=True, slots=True)
class SFImplicitIndexMap:
    """The concatenated index space as pure host tables — no device, no Warp allocation.

    Split out of :class:`SFImplicitTopology` so the index arithmetic that defines the solver's
    address space can be checked without a GPU and without a Warp-CPU allocation.  The repo's
    contract is that Warp executes on CUDA only (I0-A), and the static guard
    ``test_ac_tests_do_not_launch_production_kernels_on_cpu`` enforces it over test SOURCE — so a
    test that reached this arithmetic by passing ``device="cpu"`` was violating the contract to
    check something that never needed a device in the first place.

    Attributes:
        n_sf: Count of ``sf_arc`` nodes, occupying ``[0, n_sf)``.
        n_nmii: Count of ``nmii`` particles, occupying ``[n_sf, n_sf + n_nmii)``.
        n: Total width of the concatenated space.
        backbone_bonds: NMII backbone bond pairs, shifted into the concatenated space.
        head_bonds: NMII head bond pairs, shifted.
        backbone_angles: NMII backbone angle triples, shifted.
        head_arm_angles: NMII head-arm angle triples, shifted.
        head_node_shifted: Per-head NMII particle index, shifted.
        pinned: ``(n,)`` Dirichlet mask; SF anchors held, NMII particles never held.
        n_pinned: Number of held nodes.
    """

    n_sf: int
    n_nmii: int
    n: int
    backbone_bonds: np.ndarray
    head_bonds: np.ndarray
    backbone_angles: np.ndarray
    head_arm_angles: np.ndarray
    head_node_shifted: np.ndarray
    pinned: np.ndarray
    n_pinned: int


def build_sf_implicit_index_map(
    *,
    n_sf: int,
    n_nmii: int,
    backbone_bonds: np.ndarray,
    head_bonds: np.ndarray,
    backbone_angles: np.ndarray,
    head_arm_angles: np.ndarray,
    head_node: np.ndarray,
    pinned_sf_nodes: np.ndarray,
) -> SFImplicitIndexMap:
    """Shift the NMII tables into the concatenated space and build the Dirichlet mask.

    SF-side tables are already SF-local and are not shifted; NMII-side tables shift by ``+n_sf``
    exactly once.  The mask holds only SF anchors — an NMII particle is never a Dirichlet site,
    because a held minifilament would be a reaction the motor never generated.

    Args:
        n_sf: Count of ``sf_arc`` nodes.
        n_nmii: Count of ``nmii`` particles.
        backbone_bonds: NMII-local backbone bond pairs, shape ``(m, 2)``.
        head_bonds: NMII-local head bond pairs, shape ``(m, 2)``.
        backbone_angles: NMII-local backbone angle triples, shape ``(m, 3)``.
        head_arm_angles: NMII-local head-arm angle triples, shape ``(m, 3)``.
        head_node: NMII-local particle index per head.
        pinned_sf_nodes: SF-local indices of the held FA anchors.

    Returns:
        The host-side :class:`SFImplicitIndexMap`.

    Raises:
        ValueError: If any pin index falls outside the ``sf_arc`` node block.
    """
    n_sf = int(n_sf)
    n_nmii = int(n_nmii)
    n = n_sf + n_nmii
    shift = np.int32(n_sf)

    pinned = np.zeros(n, np.int32)
    sites = np.asarray(pinned_sf_nodes, np.int64).reshape(-1)
    if sites.size:
        if int(sites.min()) < 0 or int(sites.max()) >= n_sf:
            raise ValueError("pinned_sf_nodes must index the sf_arc node block")
        pinned[sites] = 1

    return SFImplicitIndexMap(
        n_sf=n_sf,
        n_nmii=n_nmii,
        n=n,
        backbone_bonds=backbone_bonds.reshape(-1, 2).astype(np.int32) + shift,
        head_bonds=head_bonds.reshape(-1, 2).astype(np.int32) + shift,
        backbone_angles=backbone_angles.reshape(-1, 3).astype(np.int32) + shift,
        head_arm_angles=head_arm_angles.reshape(-1, 3).astype(np.int32) + shift,
        head_node_shifted=np.asarray(head_node).astype(np.int32) + shift,
        pinned=pinned,
        n_pinned=int(pinned.sum()),
    )


class SFImplicitTopology:
    """Device topology + parameters of the SF-motor operator, in the solver's concatenated index space.

    Built once from the slice's already-existing host/device tables.  Index convention: entries ``[0, n_sf)``
    are ``sf_arc`` nodes and ``[n_sf, n_sf + n_nmii)`` are ``nmii`` particles.  SF-side tables (links, arc
    joints, bending triples, and the crossbridge's ``seg_a``/``seg_b``) are already in SF-local indices and need
    no shift; NMII-side tables (backbone bonds, head bonds, both angle triple sets, ``head_node``) are shifted
    by ``+n_sf`` ONCE at build time into solver-owned copies.  The authoritative state arrays are untouched.
    """

    def __init__(
        self,
        *,
        device: str,
        n_sf: int,
        n_nmii: int,
        sf_topology: object,
        actuator: object,
        connector_state: object,
        k_xb_pn_per_um: float,
        r0_xb_um: float,
        pinned_sf_nodes: np.ndarray,
    ) -> None:
        self.device = str(device)
        self.n_sf = int(n_sf)
        self.n_nmii = int(n_nmii)
        self.n = self.n_sf + self.n_nmii
        self.k_xb = float(k_xb_pn_per_um)
        self.r0_xb = float(r0_xb_um)

        def up(array: np.ndarray, dtype: object) -> wp.array:
            return wp.array(np.ascontiguousarray(array), dtype=dtype, device=self.device)

        # ── SF side (already SF-local indices) ────────────────────────────────────────────────────────────
        self.links_d = up(sf_topology.links.reshape(-1, 2).astype(np.int32), wp.int32)
        self.link_k_d = up(sf_topology.link_k.astype(np.float64), wp.float64)
        self.link_r0_d = up(sf_topology.link_r0.astype(np.float64), wp.float64)
        self.n_links = int(sf_topology.n_links)

        self.triples_d = up(sf_topology.bend_triples.reshape(-1, 3).astype(np.int32), wp.int32)
        self.alpha_d = up(sf_topology.bend_alpha.astype(np.float64), wp.float64)
        self.n_triples = int(sf_topology.n_triples)

        self.arc_d = up(sf_topology.arc_joints.reshape(-1, 2).astype(np.int32), wp.int32)
        self.arc_k_d = up(sf_topology.arc_k.astype(np.float64), wp.float64)
        self.arc_r0_d = up(sf_topology.arc_r0.astype(np.float64), wp.float64)
        self.n_arc = int(sf_topology.n_arc_joints)

        # ── NMII side: shift into the concatenated space ONCE (solver-owned copies) ───────────────────────
        # The arithmetic lives in build_sf_implicit_index_map so it has exactly one definition and can be
        # checked without a device; this method only uploads what that returns.
        mech = actuator.mechanics
        index_map = build_sf_implicit_index_map(
            n_sf=self.n_sf,
            n_nmii=self.n_nmii,
            backbone_bonds=mech.backbone_bonds_d.numpy(),
            head_bonds=mech.head_bonds_d.numpy(),
            backbone_angles=mech.backbone_angles_d.numpy(),
            head_arm_angles=mech.head_arm_angles_d.numpy(),
            head_node=actuator.head_node_d.numpy(),
            pinned_sf_nodes=pinned_sf_nodes,
        )
        self.index_map = index_map
        self.backbone_bonds_d = up(index_map.backbone_bonds, wp.int32)
        self.head_bonds_d = up(index_map.head_bonds, wp.int32)
        self.backbone_angles_d = up(index_map.backbone_angles, wp.int32)
        self.head_arm_angles_d = up(index_map.head_arm_angles, wp.int32)
        self.n_backbone_bonds = int(self.backbone_bonds_d.shape[0])
        self.n_head_bonds = int(self.head_bonds_d.shape[0])
        self.n_backbone_angles = int(self.backbone_angles_d.shape[0])
        self.n_head_arm_angles = int(self.head_arm_angles_d.shape[0])
        self.k_backbone = float(mech.k_backbone)
        self.r0_backbone = float(mech.r0_backbone)
        self.k_head_spring = float(mech.k_head_spring)
        self.r0_head = float(mech.r0_head)
        self.k_theta_backbone = float(mech.bending.k_theta_backbone)
        self.k_theta_arm = float(mech.bending.k_theta_arm)

        # crossbridge: head_node is NMII-local (shift); seg_a/seg_b are SF-local (no shift).  The binding SoA
        # itself is the connector's LIVE state and is referenced, never copied — the tangent must see the same
        # bound set the force saw this candidate.
        self.head_node_shifted_d = up(index_map.head_node_shifted, wp.int32)
        self.state = connector_state
        self.n_heads = int(connector_state.n_heads)

        # ── Dirichlet mask over the concatenated space (SF FA anchors held; NMII always free) ─────────────
        self.pinned_d = up(index_map.pinned, wp.int32)
        self.n_pinned = index_map.n_pinned


class SFImplicitCG(ProjectedAnalyticCG):
    """Preconditioned-CG solve of ``(aI + M K M + a(I−M)) dx = M F`` for the SF-motor slice.

    Overrides exactly the four members that read the incumbent ``cell`` (``__init__``, :meth:`_project`,
    :meth:`_stiffness`, :meth:`_build_preconditioner`) and inherits ``solve`` / ``_operator`` /
    ``_precondition`` — the audited device-resident recurrence and its convergence/curvature/breakdown latches
    — unchanged.  ``self.cell`` is deliberately set to ``None`` so any inherited path that still reached for it
    would raise loudly instead of degrading silently.

    Args:
        topology: the built :class:`SFImplicitTopology` (concatenated index space + Dirichlet mask).
        max_iterations: PCG launch budget.  A fixed compute budget, never interpreted as physical time; the
            caller must not authorize a candidate from an unconverged solve (see the driver's line search).
        cg_check_every: optional host early-exit cadence for the outer PCG loop.  ``0`` runs the full fixed
            budget (the device latch already no-ops post-convergence updates); a positive value reads only the
            device ``active`` flag — a control-flow readback, not an authoritative physics roundtrip — and the
            returned ``dx`` is bit-identical either way.
    """

    def __init__(self, topology: SFImplicitTopology, *, max_iterations: int = 200,
                 cg_check_every: int = 0) -> None:
        if max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        if cg_check_every < 0:
            raise ValueError("cg_check_every must be nonnegative")
        self.topology = topology
        self.cell = None                     # loud failure if an inherited path reaches for it
        self.device = topology.device
        self.n = int(topology.n)
        self.max_iterations = int(max_iterations)
        self.cg_check_every = int(cg_check_every)

        # Every cortex-only branch of the inherited solve/_precondition is switched OFF here.
        self.use_fiber_block = False
        self.use_erm_augmented_block = False
        self.erm_tension_side_precond = False
        self.multigrid = False
        self.mg = None
        self._mg_reg = None
        self.fq_coarse = None
        self.fq_mode = "pathA"
        self._fq_reg = None
        self.coarse_modes = 0
        self.coarse_basis_d = None
        self.coarse_iterations = 0
        self.max_fiber_nodes = 0
        self.mem_edges_d = None

        with wp.ScopedDevice(self.device):
            self.r = wp.zeros(self.n, dtype=wp.vec3d)
            self.z = wp.zeros(self.n, dtype=wp.vec3d)
            self.projected_z = wp.zeros(self.n, dtype=wp.vec3d)
            self.p = wp.zeros(self.n, dtype=wp.vec3d)
            self.ap = wp.zeros(self.n, dtype=wp.vec3d)
            self.dx = wp.zeros(self.n, dtype=wp.vec3d)
            self.k_input = wp.zeros(self.n, dtype=wp.vec3d)
            self.k_output = wp.zeros(self.n, dtype=wp.vec3d)
            self.projected_k = wp.zeros(self.n, dtype=wp.vec3d)
            self.preconditioner = wp.zeros(self.n, dtype=wp.vec3d)
            self.coarse_external_diagonal = wp.zeros(self.n, dtype=wp.vec3d)
            # solver-owned concatenated scratch: the two OWNER arrays are never merged or aliased.
            self.cat_pos = wp.zeros(self.n, dtype=wp.vec3d)
            self.cat_rhs = wp.zeros(self.n, dtype=wp.vec3d)
            self.rr = wp.zeros(1, dtype=wp.float64)
            self.rr_initial = wp.zeros(1, dtype=wp.float64)
            self.rz = wp.zeros(1, dtype=wp.float64)
            self.rz_old = wp.zeros(1, dtype=wp.float64)
            self.p_ap = wp.zeros(1, dtype=wp.float64)
            self.alpha = wp.zeros(1, dtype=wp.float64)
            self.beta = wp.zeros(1, dtype=wp.float64)
            self.active = wp.zeros(1, dtype=wp.int32)
            self.converged = wp.zeros(1, dtype=wp.int32)
            self.iterations = wp.zeros(1, dtype=wp.int32)
            self.regularization = wp.zeros(1, dtype=wp.float64)
            self.finite = wp.ones(1, dtype=wp.int32)

    # ── the three cell-dependent overrides ────────────────────────────────────────────────────────────────
    def _project(self, pos: wp.array, source: wp.array, destination: wp.array,
                 finite: wp.array) -> None:
        """Apply the FA-anchor Dirichlet mask ``M`` (NOT the parent's NF2007 inextensibility projector)."""
        wp.launch(mask_vec3_kernel, dim=self.n,
                  inputs=[source, self.topology.pinned_d, destination], device=self.device)

    def _stiffness(self, pos: wp.array, vector: wp.array, out: wp.array,
                   regularization: wp.array) -> None:
        """Launch ``out = a·vector + K·vector`` with one reused tangent per force the slice launches."""
        t = self.topology
        d = self.device
        wp.launch(set_mass_action_kernel, dim=self.n, inputs=[vector, regularization, out], device=d)
        if t.n_links:
            wp.launch(add_pair_stiffness_kernel, dim=t.n_links,
                      inputs=[pos, vector, t.links_d, t.link_k_d, t.link_r0_d, out], device=d)
        if t.n_arc:
            wp.launch(add_pair_stiffness_kernel, dim=t.n_arc,
                      inputs=[pos, vector, t.arc_d, t.arc_k_d, t.arc_r0_d, out], device=d)
        if t.n_triples:
            wp.launch(add_bending_stiffness_kernel, dim=t.n_triples,
                      inputs=[vector, t.triples_d, t.alpha_d, out], device=d)
        if t.n_backbone_bonds:
            wp.launch(add_uniform_pair_stiffness_kernel, dim=t.n_backbone_bonds,
                      inputs=[pos, vector, t.backbone_bonds_d, wp.float64(t.k_backbone),
                              wp.float64(t.r0_backbone), out], device=d)
        if t.n_head_bonds:
            wp.launch(add_uniform_pair_stiffness_kernel, dim=t.n_head_bonds,
                      inputs=[pos, vector, t.head_bonds_d, wp.float64(t.k_head_spring),
                              wp.float64(t.r0_head), out], device=d)
        if t.n_backbone_angles and t.k_theta_backbone > 0.0:
            wp.launch(add_angle_gauss_newton_stiffness_kernel, dim=t.n_backbone_angles,
                      inputs=[pos, vector, t.backbone_angles_d, wp.float64(t.k_theta_backbone),
                              wp.float64(BACKBONE_REST_ANGLE), out], device=d)
        if t.n_head_arm_angles and t.k_theta_arm > 0.0:
            wp.launch(add_angle_gauss_newton_stiffness_kernel, dim=t.n_head_arm_angles,
                      inputs=[pos, vector, t.head_arm_angles_d, wp.float64(t.k_theta_arm),
                              wp.float64(ARM_REST_ANGLE), out], device=d)
        if t.n_heads:
            wp.launch(add_segment_crossbridge_stiffness_kernel, dim=t.n_heads,
                      inputs=[pos, vector, t.head_node_shifted_d, t.state.bound_d, t.state.seg_a_d,
                              t.state.seg_b_d, t.state.bary_t_d, t.state.abscissa_d, t.state.walk_dir_d,
                              wp.float64(t.k_xb), wp.float64(t.r0_xb), out], device=d)

    def _build_preconditioner(self, pos: wp.array, regularization: wp.array,
                              finite: wp.array) -> None:
        """Node-Jacobi diagonal of the same families — the fix for the 4.6e5 vs 1e3 stiffness contrast."""
        t = self.topology
        d = self.device
        diag = self.preconditioner
        wp.launch(set_preconditioner_kernel, dim=self.n, inputs=[regularization, diag], device=d)
        if t.n_links:
            wp.launch(add_pair_preconditioner_kernel, dim=t.n_links,
                      inputs=[pos, t.links_d, t.link_k_d, t.link_r0_d, diag], device=d)
        if t.n_arc:
            wp.launch(add_pair_preconditioner_kernel, dim=t.n_arc,
                      inputs=[pos, t.arc_d, t.arc_k_d, t.arc_r0_d, diag], device=d)
        if t.n_triples:
            wp.launch(add_bending_preconditioner_kernel, dim=t.n_triples,
                      inputs=[t.triples_d, t.alpha_d, diag], device=d)
        if t.n_backbone_bonds:
            wp.launch(add_uniform_pair_preconditioner_kernel, dim=t.n_backbone_bonds,
                      inputs=[pos, t.backbone_bonds_d, wp.float64(t.k_backbone),
                              wp.float64(t.r0_backbone), diag], device=d)
        if t.n_head_bonds:
            wp.launch(add_uniform_pair_preconditioner_kernel, dim=t.n_head_bonds,
                      inputs=[pos, t.head_bonds_d, wp.float64(t.k_head_spring),
                              wp.float64(t.r0_head), diag], device=d)
        if t.n_backbone_angles and t.k_theta_backbone > 0.0:
            wp.launch(add_angle_gauss_newton_preconditioner_kernel, dim=t.n_backbone_angles,
                      inputs=[pos, t.backbone_angles_d, wp.float64(t.k_theta_backbone),
                              wp.float64(BACKBONE_REST_ANGLE), diag], device=d)
        if t.n_head_arm_angles and t.k_theta_arm > 0.0:
            wp.launch(add_angle_gauss_newton_preconditioner_kernel, dim=t.n_head_arm_angles,
                      inputs=[pos, t.head_arm_angles_d, wp.float64(t.k_theta_arm),
                              wp.float64(ARM_REST_ANGLE), diag], device=d)
        if t.n_heads:
            wp.launch(add_segment_crossbridge_preconditioner_kernel, dim=t.n_heads,
                      inputs=[pos, t.head_node_shifted_d, t.state.bound_d, t.state.seg_a_d,
                              t.state.seg_b_d, t.state.bary_t_d, t.state.abscissa_d, t.state.walk_dir_d,
                              wp.float64(t.k_xb), wp.float64(t.r0_xb), diag], device=d)

    # ── the caller-facing implicit step ───────────────────────────────────────────────────────────────────
    def step(
        self,
        *,
        sf_position_d: wp.array,
        sf_force_d: wp.array,
        nmii_position_d: wp.array,
        nmii_force_d: wp.array,
        mobility_step: float,
    ) -> wp.array:
        r"""Take ONE implicit overdamped step from the current force, returning the applied displacement.

        Solves ``(aI + M K M + a(I−M)) dx = M F`` with ``a = 1/mobility_step`` and adds ``dx`` to both owner
        position arrays.  ``mobility_step`` is the caller's ``dt/γ`` — the SAME CFL mobility the explicit path
        uses — so this reduces to the explicit update when ``K → 0`` and is the IMEX overdamped step otherwise.
        Nothing here is tuned: the only scalar is the caller's derived mobility.

        The caller owns force zeroing + accumulation before each call (the slice's ``zero_forces`` +
        ``accumulate``), exactly as with the explicit path.  Positions are updated in place on BOTH owner
        arrays; the arrays are never merged (only this solver's scratch is contiguous).

        Returns:
            The device displacement ``dx`` in the concatenated space (for diagnostics; already applied).
        """
        if not (mobility_step > 0.0 and np.isfinite(mobility_step)):
            raise ValueError(f"mobility_step must be positive-finite, got {mobility_step!r}")
        d = self.device
        n_sf = wp.int32(self.topology.n_sf)
        self.regularization.fill_(wp.float64(1.0 / float(mobility_step)))
        self.finite.fill_(wp.int32(1))
        wp.launch(gather_two_kernel, dim=self.n,
                  inputs=[sf_position_d, nmii_position_d, n_sf, self.cat_pos], device=d)
        wp.launch(gather_two_kernel, dim=self.n,
                  inputs=[sf_force_d, nmii_force_d, n_sf, self.cat_rhs], device=d)
        # Mask the residual so the pinned block's right-hand side is zero; the operator masking (via _project,
        # applied by the inherited _operator on input AND output) is what makes dx_pinned exactly zero.
        self._project(self.cat_pos, self.cat_rhs, self.cat_rhs, self.finite)
        dx = self.solve(self.cat_pos, self.cat_rhs, self.regularization, self.finite)
        wp.launch(add_scatter_two_kernel, dim=self.n,
                  inputs=[dx, n_sf, sf_position_d, nmii_position_d], device=d)
        return dx

    def converged_flag(self) -> bool:
        """Host read of the PCG convergence latch, for OUT-OF-LOOP reporting only (never a physics decision)."""
        return bool(int(self.converged.numpy()[0]))

    def iteration_count(self) -> int:
        """Host read of the PCG iteration count reached, for OUT-OF-LOOP reporting only."""
        return int(self.iterations.numpy()[0])
