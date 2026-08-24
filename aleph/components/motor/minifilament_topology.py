r"""Bipolar minifilament topology/count invariants + per-head Newton force-balance closure (host NumPy).

Host-side acceptance oracle (pure NumPy, NO Warp) for two I3 gates: "topology/count invariants" and
"per-head Newton closure" (§3 line 243). It also carries the ``k_xb`` MASTER-force-knob working-stroke-strain
arbiter (the "1 pN/µm breaks stall" check from §3 Session D).

TOPOLOGY (ported read-only from the retired HOOMD `cortex/myosin.py` + `bridge/motor.py`, deleted at
`1a9ded66`; recoverable with `git show 1a9ded66^:aleph/archive/hoomd_legacy/...` — NEVER executed):
one bipolar minifilament is ``[backbone_0 .. backbone_{n_bb-1}, head+_0 .. head+_{H-1}, head-_0 .. head-_{H-1}]``
with ``H = n_heads_per_side`` heads on each anti-parallel side. Bonds (static, at construction):
  * backbone-backbone: ``n_bb - 1`` consecutive bonds (a connected chain), stiff rigid-rod (no bending potential);
  * head-backbone:     ``2H`` bonds (each head anchors to exactly one backbone bead), rest length r0_head;
  * head-actin:        DYNAMIC (0 at construction; created by the KMC attach at runtime) — NOT static topology.
Invariant totals: ``n_particles = n_bb + 2H``; ``n_static_bonds = (n_bb - 1) + 2H``. Bipolar symmetry: exactly
``H`` heads per side. These hold for ANY (n_bb >= 2, H >= 1) — the count-invariant gate is magnitude-free
(the specific n_bb=14, H in {10, 28} numbers are I0-B3 GAPs; the invariants are not).

PER-HEAD NEWTON CLOSURE: a bound head is a crossbridge spring (stiffness ``k_xb``, rest ``r0_head``) between a
backbone bead and its actin site; the actin site is held by an anchor spring (``k_anchor``, rest ``x0_a``). The
inner mechanical solve must drive the net nodal force to zero. For this linear system the equilibrium is
closed-form and Newton converges in ONE step (quadratic); the gate checks the residual ``|dE/dx| -> 0`` to
machine precision AND that the crossbridge force is equal-and-opposite on its two endpoints (Newton's 3rd law).

k_xb MASTER-KNOB arbiter: at equilibrium the head bears ``F_head`` only if its strain ``F_head / k_xb`` is
physical. Working stroke is ~5-20 nm, so ``k_xb`` must be ~100-1000 pN/µm. The archived cortex value
``k_head_actin = 1 pN/µm`` gives strain ``F_head / 1`` = 0.5-2 µm — larger than the whole ~0.3 µm minifilament,
so the stall geometry collapses ("1 pN/µm breaks stall"). :func:`working_stroke_strain` exposes this; the gate
flags the unphysical branch. ``k_xb`` is an I0-B3 GAP — this oracle certifies physicality, it does not choose.

Sanity Gate (of the oracle itself — tests/ac/motor/test_topology_oracle.py):
  * count invariants: n_particles = n_bb + 2H; n_static_bonds = (n_bb-1) + 2H; H heads per side (bipolar);
  * connectivity: backbone bonds form a single connected chain; every head bonds exactly one backbone bead;
  * geometry: the two head sets sit on opposite sides of the backbone axis (net zero transverse offset);
  * Newton closure: residual -> 0 (< 1e-12) within 1-2 steps; converges to the closed-form equilibrium;
  * Newton 3rd law: crossbridge force on backbone bead = -(force on actin node) exactly;
  * k_xb arbiter: strain(F_head, k_xb=100-1000) in [5, 20] nm band; strain(F_head, k_xb=1) > L_bb (unphysical).

References:
  retired HOOMD cortex/myosin.py (topology/bonds), bridge/motor.py (axial bipolar split) — both deleted
  at `1a9ded66`, recoverable from git history; Stam-Hocky
  2015; Billington 2013 (n_bb, H, 301 nm backbone). Newton-Raphson closure of a linear spring lattice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

__all__ = [
    "MinifilamentTopology",
    "head_newton_residual",
    "straddle_frame",
    "working_stroke_strain",
]


@dataclass(frozen=True, slots=True)
class MinifilamentTopology:
    """One bipolar Stam-Hocky minifilament's discrete topology (counts + geometry), engine units µm.

    Attributes:
        n_bb: explicit backbone beads (>= 2). I0-B3 GAP reference value 14.
        n_heads_per_side: heads ``H`` on each anti-parallel side (>= 1). I0-B3 GAP: 10 (AFINES) / 28 (Billington).
        backbone_length_um: backbone contour length ``L_bb`` [µm]. I0-B3 GAP: 0.301 µm (Billington 2013).
        head_offset_um: perpendicular head<->backbone offset ``r0_head`` [µm]. Archived 0.200 µm.
    """

    n_bb: int
    n_heads_per_side: int
    backbone_length_um: float
    head_offset_um: float

    def __post_init__(self) -> None:
        if self.n_bb < 2:
            raise ValueError("n_bb must be >= 2 (a backbone is a chain)")
        if self.n_heads_per_side < 1:
            raise ValueError("n_heads_per_side must be >= 1")
        if self.backbone_length_um <= 0.0 or self.head_offset_um < 0.0:
            raise ValueError("backbone_length_um must be > 0 and head_offset_um >= 0")

    @property
    def n_heads(self) -> int:
        """Total heads (both sides) = ``2 * n_heads_per_side``."""
        return 2 * self.n_heads_per_side

    @property
    def n_particles(self) -> int:
        """Total particles = backbone beads + heads = ``n_bb + 2H``."""
        return self.n_bb + self.n_heads

    @property
    def n_backbone_bonds(self) -> int:
        """Consecutive backbone-backbone bonds = ``n_bb - 1`` (a connected chain)."""
        return self.n_bb - 1

    @property
    def n_head_bonds(self) -> int:
        """Head<->backbone bonds = ``2H`` (one per head)."""
        return self.n_heads

    @property
    def n_static_bonds(self) -> int:
        """Static bonds at construction = ``(n_bb - 1) + 2H`` (head<->actin is dynamic, not counted)."""
        return self.n_backbone_bonds + self.n_head_bonds

    @property
    def segment_length_um(self) -> float:
        """Backbone bead spacing = ``L_bb / (n_bb - 1)`` [µm] (the backbone-backbone rest length)."""
        return self.backbone_length_um / (self.n_bb - 1)

    def positions(self) -> npt.NDArray[np.float64]:
        """Reference bead positions ``(n_particles, 3)`` [µm] for the axial-bipolar layout (motor.py style).

        Backbone beads on the x-axis over ``[-L/2, +L/2]``; ``+`` heads flank the first half at ``+offset`` in
        +y, ``-`` heads flank the second half at ``-offset`` in -y (opposite transverse sides = bipolar).
        Ordering matches the port spec: ``[backbone..., head+..., head-...]``.
        """
        xs = np.linspace(-0.5 * self.backbone_length_um, 0.5 * self.backbone_length_um, self.n_bb)
        bb = np.column_stack([xs, np.zeros(self.n_bb), np.zeros(self.n_bb)])
        half = max(self.n_bb // 2, 1)
        h = self.n_heads_per_side
        # + side maps to backbone beads [0, half), - side to [n_bb-half, n_bb) (axial half-split, motor.py:321-349)
        bb_plus = (np.arange(h) * half) // h
        bb_minus = self.n_bb - 1 - (np.arange(h) * half) // h
        head_plus = bb[bb_plus] + np.array([0.0, self.head_offset_um, 0.0])
        head_minus = bb[bb_minus] + np.array([0.0, -self.head_offset_um, 0.0])
        return np.vstack([bb, head_plus, head_minus])

    def placed_positions(
        self,
        centre: npt.NDArray[np.float64],
        e_x: npt.NDArray[np.float64],
        e_y: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        r"""Reference bead positions ``(n_particles, 3)`` [µm] mapped into a caller-chosen ORTHONORMAL frame.

        Same local axial-bipolar layout as :meth:`positions` (backbone on local ``+x``, ``+`` heads at ``+y``
        offset, ``−`` heads at ``−y``), but rotated so local ``x̂ → e_x`` and local ``ŷ → e_y`` (``e_z = e_x × e_y``)
        and translated to ``centre``. This is the *placement/orientation* degree of freedom :meth:`positions`
        (and :func:`aleph.components.motor.minifilament_warp.build_minifilament_nodes`, which only aligns ``x̂ → axis``
        and leaves the transverse ``ŷ`` arbitrary) does not expose: to make a bipolar minifilament STRADDLE two
        anti-parallel actin filaments, the backbone axis must run along the shared filament line (``e_x``) AND
        the head-offset direction must point across at the filaments (``e_y``) so the ``+`` heads land on one
        filament and the ``−`` heads on the anti-parallel other. ``e_y`` is re-orthogonalised against ``e_x``;
        the rigid geometry (backbone straight, arms perpendicular at ``head_offset_um``) is preserved exactly, so
        the backbone-bend and head-arm angle springs stay force-free at rest just as in :meth:`positions`.

        Args:
            centre: minifilament centre ``(3,)`` [µm] in the cell frame.
            e_x: backbone-axis direction ``(3,)`` (need not be unit; normalised here).
            e_y: transverse head-offset direction ``(3,)`` (re-orthogonalised against ``e_x``); the ``+`` heads
                sit on the ``+e_y`` side.

        Returns:
            ``(n_particles, 3)`` world-frame bead positions, ordering ``[backbone…, head+…, head−…]``.
        """
        local = self.positions()
        ex = np.asarray(e_x, dtype=np.float64)
        ex = ex / (np.linalg.norm(ex) + 1e-30)
        ey = np.asarray(e_y, dtype=np.float64)
        ey = ey - np.dot(ey, ex) * ex                 # Gram-Schmidt: strip the e_x component
        ey = ey / (np.linalg.norm(ey) + 1e-30)
        ez = np.cross(ex, ey)
        rot = np.column_stack([ex, ey, ez])           # local axes as columns → local→world map
        return (local @ rot.T) + np.asarray(centre, dtype=np.float64)

    def head_backbone_bonds(self) -> npt.NDArray[np.int64]:
        """Head<->backbone bond pairs ``(2H, 2)`` as ``[head_particle_index, backbone_bead_index]``."""
        half = max(self.n_bb // 2, 1)
        h = self.n_heads_per_side
        bb_plus = (np.arange(h) * half) // h
        bb_minus = self.n_bb - 1 - (np.arange(h) * half) // h
        plus_idx = self.n_bb + np.arange(h)
        minus_idx = self.n_bb + h + np.arange(h)
        return np.vstack([
            np.column_stack([plus_idx, bb_plus]),
            np.column_stack([minus_idx, bb_minus]),
        ]).astype(np.int64)

    def backbone_angle_triples(self) -> npt.NDArray[np.int64]:
        r"""Consecutive backbone bending triples ``(n_bb−2, 3)`` as ``[i−1, i, i+1]`` (vertex = the middle bead).

        The F6 backbone-bending topology: an angle-harmonic on each interior bead (rest ``θ₀ = π``, straight)
        turns the distance-spring chain into a rigid rod that TRANSMITS the axial crossbridge tension instead
        of bowing (see :mod:`aleph.components.motor.backbone_warp`). Empty ``(0, 3)`` for ``n_bb = 2`` (no interior
        bead). The count-invariant ``n_backbone_angles = max(n_bb − 2, 0)`` holds for any ``n_bb ≥ 2``.
        """
        if self.n_bb < 3:
            return np.zeros((0, 3), dtype=np.int64)
        mid = np.arange(1, self.n_bb - 1)
        return np.column_stack([mid - 1, mid, mid + 1]).astype(np.int64)

    def head_arm_angle_triples(self) -> npt.NDArray[np.int64]:
        r"""Head-arm orientation triples ``(2H, 3)`` as ``[backbone_neighbor, backbone_bead, head]`` (vertex =
        the backbone bead the head anchors to).

        The F6 head-arm topology: an angle-harmonic at the anchoring backbone bead between its backbone
        neighbor and the head (rest ``θ₀ = π/2`` — the lever arm stands perpendicular to the backbone) gives
        the arm the transverse stiffness a single distance spring lacks, so a head pulled along ``walk_dir``
        LOADS the backbone instead of swinging (see :mod:`aleph.components.motor.backbone_warp`). The neighbor is the
        bead one step toward the backbone CENTRE (``bead+1`` on the ``−x`` half, ``bead−1`` on the ``+x`` half)
        — a mirror-symmetric choice so the ``+`` and ``−`` arms are stiffened identically and the two clamp
        reactions stay equal-and-opposite (an asymmetric neighbor rule biases ``|reaction_A| ≠ |reaction_B|``).
        Ordering matches :meth:`head_backbone_bonds` / :meth:`positions` (``[+ heads, − heads]``).
        """
        half = max(self.n_bb // 2, 1)
        h = self.n_heads_per_side
        bb_plus = (np.arange(h) * half) // h
        bb_minus = self.n_bb - 1 - (np.arange(h) * half) // h
        beads = np.concatenate([bb_plus, bb_minus])
        head_idx = np.concatenate([self.n_bb + np.arange(h), self.n_bb + h + np.arange(h)])
        centre = 0.5 * (self.n_bb - 1)
        # step toward the backbone centre (symmetric); clamp so an end bead still has a valid in-range neighbor
        neighbor = np.where(beads < centre, beads + 1, beads - 1)
        neighbor = np.clip(neighbor, 0, self.n_bb - 1)
        return np.column_stack([neighbor, beads, head_idx]).astype(np.int64)

    @property
    def n_backbone_angles(self) -> int:
        """Backbone bending triples = ``max(n_bb − 2, 0)`` (one per interior bead)."""
        return max(self.n_bb - 2, 0)

    @property
    def n_head_arm_angles(self) -> int:
        """Head-arm orientation triples = ``2H`` (one per head)."""
        return self.n_heads

    def transverse_imbalance_um(self) -> float:
        """Net transverse (y) offset of all heads [µm]; ~0 for a symmetric bipolar minifilament."""
        pos = self.positions()
        return float(pos[self.n_bb:, 1].sum())


def straddle_frame(
    a: npt.NDArray[np.float64],
    t_a: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64],
    t_b: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    r"""Placement frame for a bipolar minifilament straddling two anti-parallel actin filaments.

    Given actin filament A (a point ``a`` on it, unit tangent ``t_a`` toward its barbed end) and an
    anti-parallel neighbour filament B (point ``b``, tangent ``t_b ≈ −t_a``), returns the frame the
    minifilament is placed in so its heads TOUCH actin (defect#3 fix):

      * ``e_x`` (backbone axis) runs along the SHARED filament line ``t_a − t_b`` (``≈ 2 t_a`` when the two are
        anti-parallel) — the heads then walk along the actin, and each half-backbone hugs one filament;
      * ``e_y`` (head-offset direction) points from the minifilament centre toward A, so the ``+`` heads
        (local ``+y`` = ``+head_offset``) land on A and the ``−`` heads on B;
      * ``centre`` = the midpoint ``½(a + b)``.

    When the two filaments are separated by ≈ ``2·head_offset_um`` the ``±head_offset`` heads land ON the two
    filament lines (perpendicular residual → 0), so the head↔actin crossbridge starts at ~nm extension and the
    arm rests force-free — a balanced, self-anchored contractile dipole (the anti-parallel pair contracts
    inward). Degenerate inputs (parallel ``t_a‖t_b`` or ``b`` on A's line) fall back to a valid orthonormal
    frame rather than raising. All vectors are engine units [µm]; ``e_x``/``e_y`` are returned unit + orthogonal.

    Args:
        a: a point on filament A ``(3,)`` [µm] (the ``+``-head anchor).
        t_a: unit tangent of filament A ``(3,)`` (toward its barbed end).
        b: a point on the anti-parallel filament B ``(3,)`` [µm] (the ``−``-head anchor).
        t_b: unit tangent of filament B ``(3,)``.

    Returns:
        ``(centre, e_x, e_y)`` — the arguments :meth:`MinifilamentTopology.placed_positions` consumes.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ta = np.asarray(t_a, dtype=np.float64)
    tb = np.asarray(t_b, dtype=np.float64)
    axis = ta - tb                                    # anti-parallel → ≈ 2·t_a (the shared filament line)
    if np.linalg.norm(axis) < 1e-9:                   # t_a ‖ t_b (parallel, not anti-) → fall back to t_a
        axis = ta.copy()
    if np.linalg.norm(axis) < 1e-9:                   # both tangents degenerate → fall back to a→b
        axis = b - a
    e_x = axis / (np.linalg.norm(axis) + 1e-30)
    to = b - a
    n_perp = to - np.dot(to, e_x) * e_x               # component of a→b perpendicular to the backbone axis
    if np.linalg.norm(n_perp) < 1e-9:                 # b on A's line → any perpendicular works
        seed = np.array([1.0, 0.0, 0.0]) if abs(e_x[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        n_perp = seed - np.dot(seed, e_x) * e_x
    n_hat = n_perp / (np.linalg.norm(n_perp) + 1e-30)  # points A → B
    e_y = -n_hat                                       # local +y (the + heads) points back toward A
    centre = 0.5 * (a + b)
    return centre, e_x, e_y


def working_stroke_strain(f_stall_head: float, k_xb: float) -> float:
    """Crossbridge strain ``F_head / k_xb`` [µm] a head must take on to bear its stall force.

    The k_xb MASTER-knob arbiter. Physical working stroke ~5-20 nm needs ``k_xb`` ~100-1000 pN/µm; the archived
    cortex ``k_xb = 1 pN/µm`` gives strain ~0.5-2 µm (larger than the whole minifilament -> stall collapses).

    Args:
        f_stall_head: per-head stall force [pN].
        k_xb: crossbridge stiffness [pN/µm].

    Returns:
        Equilibrium crossbridge strain [µm].
    """
    if k_xb <= 0.0:
        raise ValueError("k_xb must be positive")
    return f_stall_head / k_xb


def head_newton_residual(
    x_init: float,
    k_xb: float,
    k_anchor: float,
    r0_head: float,
    x0_a: float,
    n_iter: int = 3,
) -> dict[str, object]:
    """Newton-Raphson closure of one bound head's linear spring system; residual must vanish.

    The bound head couples a backbone bead (fixed at 0) to an actin node at ``x`` via a crossbridge spring
    (``k_xb``, rest ``r0_head``); the actin node is held by an anchor spring (``k_anchor``, rest ``x0_a``).
    Total energy ``E(x) = 1/2 k_xb (x - r0_head)^2 + 1/2 k_anchor (x - x0_a)^2``; the inner solve minimizes E.
    For this linear system Newton reaches the closed-form equilibrium
    ``x* = (k_xb r0_head + k_anchor x0_a)/(k_xb + k_anchor)`` in one step.

    Args:
        x_init: initial actin-node position [µm].
        k_xb: crossbridge stiffness [pN/µm].
        k_anchor: actin anchor stiffness [pN/µm].
        r0_head: crossbridge rest length [µm] (the head↔actin crossbridge rest — post-2026-07-16 split this is
            ``NMIIHandParams.r0_xb`` ~ 0, NOT the head↔backbone arm rest; the tests pass 0.0).
        x0_a: actin anchor rest position [µm].
        n_iter: Newton iterations to run.

    Returns:
        dict with ``residuals`` (|dE/dx| after each step), ``x`` (final position), ``x_star`` (closed form),
        ``f_head`` (crossbridge force on the actin node at convergence [pN], negative = pulls toward backbone),
        ``newton3rd_residual`` (|f_on_backbone + f_on_actin|, must be 0), and ``converged`` (bool).
    """
    if k_xb <= 0.0 or k_anchor <= 0.0:
        raise ValueError("k_xb and k_anchor must be positive")
    if n_iter < 1:
        raise ValueError("n_iter must be >= 1")
    k_tot = k_xb + k_anchor
    x_star = (k_xb * r0_head + k_anchor * x0_a) / k_tot

    def dE(x: float) -> float:
        return k_xb * (x - r0_head) + k_anchor * (x - x0_a)

    x = float(x_init)
    residuals: list[float] = []
    for _ in range(n_iter):
        grad = dE(x)
        residuals.append(abs(grad))
        x = x - grad / k_tot  # Hessian = k_tot (constant); one step is exact for this linear system
    residuals.append(abs(dE(x)))

    # crossbridge force on the actin node and (Newton's 3rd law) its equal-and-opposite on the backbone bead
    f_on_actin = -k_xb * (x - r0_head)
    f_on_backbone = +k_xb * (x - r0_head)
    return {
        "residuals": residuals,
        "x": x,
        "x_star": x_star,
        "f_head": f_on_actin,
        "newton3rd_residual": abs(f_on_actin + f_on_backbone),
        "converged": residuals[-1] < 1e-12,
    }
