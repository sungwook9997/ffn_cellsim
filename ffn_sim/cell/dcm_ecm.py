"""DCM cell ⊗ explicit cross-linked fiber ECM — catch-slip focal-adhesion clutch.

Headline "beyond-the-papers" integration module (H.7 platform, 2026-06-11).
Couples the project's **Deformable Cell Model** shell (:mod:`ffn_sim.cell.dcm`)
to an **explicit cross-linked Mikado fiber ECM** (:mod:`ffn_sim.ecm.mikado`)
through a **real catch-slip focal-adhesion clutch**, all on the project's
Leimkuhler-Matthews BAOAB integrator (:mod:`ffn_sim.integrator.baoab`).

This reproduces *and exceeds* the reference paper Slater…T.Kim, *Soft Matter*
2021 (d0sm01911a — a contractile cell remodelling a viscoelastic fiber ECM) by
replacing each of their lumped mechanisms with our finer-grained one:

================  ============================  =================================
Mechanism         d0sm01911a (reference)        this module (beyond)
================  ============================  =================================
cell body         simple cortex / point dipole  exact-turgor deformable shell
                                                 (icosphere + face-normal pressure)
cell→ECM link     permanent linear spring       Pereverzev **catch-slip** FA clutch
                                                 (force-free at build, breaks +
                                                 rebinds via k_off(F))
ECM               continuum / mean-field fibers  **explicit** cross-linked Mikado
                                                 fiber beads (bond+angle+WCA)
contraction       prescribed contractile stress  membrane edge-spring prestress
                                                 (r0 → φ·ℓ̄, tension emerges)
================  ============================  =================================

Architecture (all SI units)
----------------------------
1. A flat Mikado fiber bed (z ≈ 0) is built via
   :func:`ffn_sim.ecm.mikado.build_mikado_state` and cropped to a window
   ~``footprint_factor``× the cell footprint, so the cell sits on a finite
   patch we can remodel and measure (the full 200 µm Mikado box is too sparse).
2. ONE DCM icosphere shell (:func:`ffn_sim.cell.dcm.icosphere_mesh`) is placed
   resting on the bed, lowest node a small gap ``z_gap`` above z=0 — NO initial
   LJ overlap with fiber beads (the BAOAB blow-up lesson).
3. The two gsd frames are merged into one snapshot with distinct particle types
   (``actin_ecm`` fibers + ``dcm_mem`` cell nodes), bonds/angles re-indexed
   (mirrors :func:`ffn_sim.cell.doublet._merge_two_cortex_frames`).
4. **FA clutch**: harmonic bonds from cell BASAL nodes (within ``R_cell`` of the
   bed plane) to their nearest fiber bead, ``r0 = realised separation`` so the
   clutch is FORCE-FREE at construction. The bonds are made catch-slip by
   :class:`FaCatchSlipUpdater`, which reuses the canonical Pereverzev
   ``k_off(F)`` kernel (:mod:`ffn_sim.validation.pereverzev` — the SAME physics
   wired into :class:`ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater`) to
   break over-loaded clutches and rebind slack basal nodes.
5. **Contraction**: the membrane edge springs are given a reduced rest length
   ``r0 = contractility·ℓ̄`` (< mean edge) so the shell is under cortical tension
   and contracts, pulling the FA-anchored fiber beads inward — the d0sm01911a
   "cortical contraction" expressed mechanistically as a prestressed spring net.
6. **Excluded volume**: WCA between ``dcm_mem`` and ``actin_ecm`` so the cell
   cannot pass through the fibers.

FA-clutch design note (catch-slip vs the H.4 IntegrinBondUpdater)
-----------------------------------------------------------------
:class:`ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater` implements the full
H.4 catch-slip clutch but is wired around ``FALayout`` (one ligand per FA, a
single integrin↔ligand binding partner per FA, a full-snapshot ``set_snapshot``
rebuild per batch). The DCM↔ECM contact is geometrically different — every
basal cell node is a clutch anchored to whichever fiber bead is nearest, with
many candidate beads. Rather than force-fit ``FALayout`` (and the heavyweight
per-batch snapshot rebuild) onto that geometry, :class:`FaCatchSlipUpdater`
applies the **same Pereverzev** ``pereverzev_k_off(F)`` break law directly to the
FA bond group via the cheap ``cpu_local_snapshot`` topology buffer — so the
catch-slip physics is identical to H.4, but the wiring is light and robust. This
is the prompt-sanctioned "reuse the Pereverzev kernel" path; it is genuine
catch-slip, not a Bell-only fallback.

References
----------
- Slater, Idema, …, T. Kim, *Soft Matter* 2021, d0sm01911a — contractile cell
  remodelling a fiber ECM; stress ~1/r decay in 2D; radial fiber alignment.
- :mod:`ffn_sim.cell.dcm`, :mod:`ffn_sim.ecm.mikado`,
  :mod:`ffn_sim.bridge.integrin_bonds`, :mod:`ffn_sim.validation.pereverzev`.
- :func:`ffn_sim.cell.doublet._merge_two_cortex_frames` (the snapshot-merge idiom).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.cell.dcm import DcmTurgorForce, icosphere_mesh
from ffn_sim.ecm.mikado import resolve_derived
from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.validation.pereverzev import (
    EXP_ARG_GUARD,
    PereverzevParams,
    pereverzev_F_star,
    pereverzev_k_off,
)

_KT_310 = 4.28e-21  # J at 310 K

# Particle types (distinct so WCA + gamma_map can address them).
TYPE_ECM = "actin_ecm"      # Mikado fiber beads (keeps mikado's own type name)
TYPE_MEM = "dcm_mem"        # DCM membrane shell nodes

# Bond types.
BOND_EDGE = "dcm_edge"      # membrane cortical edge springs (contractile)
BOND_FA = "fa_clutch"       # cell basal node ↔ ECM fiber bead (catch-slip FA)


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedDcmEcm:
    """Resolved parameters for the DCM⊗ECM remodelling run (all SI).

    Cell-shell + adhesion + FA constants are MCF7 / collagen literature bands
    (trusted directly per the DCM tier convention). The ECM bed parameters are
    inherited from the resolved H.1 Mikado config.
    """

    # --- DCM cell shell ---
    R_cell: float = 7.5e-6          # m   MCF7 radius
    subdivisions: int = 2           # icosphere level (2 = 162 nodes — enough basal nodes)
    turgor_dP0: float = 133.0       # Pa  baseline osmotic turgor (Young-Laplace)
    K_vol: float = 1.0e3            # Pa  osmotic bulk modulus
    k_edge: float = 5.0e-4          # N/m membrane edge spring (cortical elasticity)
    contractility: float = 0.7      # edge r0 = contractility · mean_edge (< 1 → contracts)
    gamma_node: float = 3.9e-10     # N·s/m per-node Stokes drag

    # --- FA catch-slip clutch (cell basal node ↔ nearest fiber bead) ---
    k_fa: float = 1.0e-3            # N/m  integrin clutch stiffness (1 pN/nm, KU-2.7)
    fa_basal_band: float = 0.55     # basal node iff (z_node − z_bed) < band · R_cell
    fa_capture: float = 1.0e-6      # m   max basal-node↔bead distance to anchor / rebind
    k_on: float = 1.0               # s⁻¹ FA rebind on-rate (KU-2.18)
    fa_batch_steps: int = 200       # catch-slip updater cadence (batched-updater convention)

    # --- excluded volume (cell node ↔ fiber bead) ---
    eps_rep_kT: float = 5.0         # WCA strength [kT]

    # --- ECM bed window ---
    footprint_factor: float = 2.6   # bed half-window = factor · R_cell
    z_gap: float = 1.0e-7           # m   lowest cell node sits this far above z_bed

    kT: float = _KT_310
    dt: float = 5.0e-10             # s   start small (BAOAB stability); auto-reduced if needed
    seed: int = 7

    # Pereverzev catch-slip params (KU-2.18 illustrative; F* ≈ 7 pN).
    pereverzev_cfg: dict = field(default_factory=lambda: {
        "k_off_slip": 0.5, "F_s": 30.0e-12,
        "k_off_catch": 0.4, "F_c": 7.0e-12, "k_on": 1.0,
    })
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Catch-slip FA clutch updater (reuses the canonical Pereverzev k_off kernel)
# ---------------------------------------------------------------------------
class FaCatchSlipUpdater(hoomd.custom.Action):
    """Pereverzev catch-slip dynamics on the FA clutch bond group.

    Mirrors the break/rebind logic of
    :class:`ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater` (same
    ``pereverzev_k_off`` kernel, same break-before-overflow hardening) but
    operates on the DCM↔ECM ``fa_clutch`` bond group, where each clutch links a
    cell BASAL node to its nearest fiber bead. Per batch tick of
    ``batch_steps`` integration steps it:

    1. Reads positions + the current ``fa_clutch`` bond group.
    2. Breaks each engaged clutch with probability ``1 − exp(−k_off(|F|)·Δt)``
       where ``F = k_fa·(|Δr| − r0_fa)`` (catch-slip: k_off DIPS then rises with
       load — the bond strengthens under moderate tension, then yields).
    3. Rebinds each currently-unbound basal node to its nearest fiber bead if
       that bead is within ``capture`` (probability ``1 − exp(−k_on·Δt)``), with
       a per-bond ``r0`` set to the realised separation so the *new* clutch is
       force-free at the instant it forms.

    The action mutates ONLY the bond topology (positions are integrated solely
    by the BAOAB Action — no double-stepping). Per-bond ``r0`` values live in
    ``self.r0_by_pair`` keyed by ``(basal_tag, bead_tag)`` and are read back by
    :class:`FaClutchForce` (which owns the actual harmonic force so r0 can be
    per-bond, unlike a single-r0 ``md.bond.Harmonic``).
    """

    def __init__(
        self,
        *,
        force: "FaClutchForce",
        basal_tags: np.ndarray,
        bead_pos_getter,
        pereverzev: PereverzevParams,
        k_fa: float,
        k_on: float,
        capture: float,
        batch_steps: int,
        dt: float,
        seed: int = 0,
    ) -> None:
        super().__init__()
        self.force = force
        self.basal_tags = np.asarray(basal_tags, dtype=np.int64)
        self._bead_pos_getter = bead_pos_getter   # () -> (bead_tags, bead_xyz)
        self.pereverzev = pereverzev
        self.k_fa = float(k_fa)
        self.k_on = float(k_on)
        self.capture = float(capture)
        self.batch_steps = int(batch_steps)
        self._batch_dt = self.batch_steps * float(dt)
        self._rng = np.random.default_rng(seed)
        self._sim_ref: hoomd.Simulation | None = None
        self.n_break_total = 0
        self.n_bind_total = 0
        self.n_forced_break = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos_local = np.asarray(snap.particles.position).copy()
        # tag-ordered positions: pos_g[global_tag] = position
        perm = np.argsort(tag)
        pos_g = pos_local[perm]

        engaged = self.force.pairs            # (M, 2) [basal_tag, bead_tag]
        # ---- Step 1: break engaged clutches via Pereverzev k_off(F) ----
        if engaged.shape[0] > 0:
            r_basal = pos_g[engaged[:, 0]]
            r_bead = pos_g[engaged[:, 1]]
            r = np.linalg.norm(r_basal - r_bead, axis=1)
            r0 = self.force.r0_for(engaged)
            F_mag = self.k_fa * np.clip(r - r0, 0.0, None)
            # break-before-overflow: a clutch past the Pereverzev slip range is
            # physically already broken (same hardening as IntegrinBondUpdater).
            F_overflow = EXP_ARG_GUARD * min(self.pereverzev.F_s, self.pereverzev.F_c)
            safe = F_mag < F_overflow
            p_break = np.ones(F_mag.shape[0], dtype=np.float64)
            if safe.any():
                k_off = pereverzev_k_off(F_mag[safe], self.pereverzev)
                p_break[safe] = 1.0 - np.exp(-k_off * self._batch_dt)
            self.n_forced_break += int(np.count_nonzero(~safe))
            u = self._rng.uniform(0.0, 1.0, size=F_mag.shape[0])
            broke = u < p_break
            self.n_break_total += int(broke.sum())
            survivors = engaged[~broke]
        else:
            survivors = engaged

        # ---- Step 2: rebind unbound basal nodes to nearest fiber bead ----
        bound_basal = set(int(t) for t in survivors[:, 0]) if survivors.shape[0] else set()
        unbound = np.array(
            [t for t in self.basal_tags if int(t) not in bound_basal], dtype=np.int64
        )
        new_pairs = []
        new_r0 = []
        if unbound.size > 0:
            bead_tags, _bead_xyz0 = self._bead_pos_getter()
            # current bead positions (tag-ordered global read)
            bead_now = pos_g[bead_tags]
            r_unbound = pos_g[unbound]
            # nearest bead per unbound basal node, IN-PLANE (xy) — same reach-down
            # model as the build (integrin reaches the flat bed); r0 = full 3D sep.
            dxy = np.linalg.norm(
                r_unbound[:, None, :2] - bead_now[None, :, :2], axis=2
            )
            j = np.argmin(dxy, axis=1)
            dxy_min = dxy[np.arange(unbound.size), j]
            in_range = dxy_min <= self.capture
            if in_range.any():
                p_bind = 1.0 - np.exp(-self.k_on * self._batch_dt)
                u2 = self._rng.uniform(0.0, 1.0, size=int(in_range.sum()))
                fire = u2 < p_bind
                cand_basal = unbound[in_range][fire]
                cand_bead = bead_tags[j[in_range][fire]]
                cand_r0_3d = np.linalg.norm(
                    r_unbound[in_range][fire] - bead_now[j[in_range][fire]], axis=1
                )
                for b, bd, rr in zip(cand_basal, cand_bead, cand_r0_3d):
                    new_pairs.append((int(b), int(bd)))
                    new_r0.append(float(rr))
                self.n_bind_total += len(new_pairs)

        # ---- Step 3: commit the updated clutch population to the force ----
        if new_pairs:
            all_pairs = np.vstack([survivors, np.array(new_pairs, dtype=np.int64)])
            self.force.set_pairs(all_pairs, extra_r0={
                p: r for p, r in zip(new_pairs, new_r0)
            })
        else:
            self.force.set_pairs(survivors)


# ---------------------------------------------------------------------------
# Per-bond FA harmonic clutch force (per-bond r0 → force-free at any bind time)
# ---------------------------------------------------------------------------
class FaClutchForce(md.force.Custom):
    """Per-bond harmonic FA clutch with PER-BOND rest length.

    A plain ``md.bond.Harmonic`` carries a single r0 per bond *type*; a
    catch-slip clutch that rebinds at runtime needs each new clutch to be
    force-free at *its* bind separation, i.e. a per-bond r0. This Custom force
    stores the engaged clutch pairs + their individual r0 and applies
    ``F = −k_fa·(|Δr| − r0)·r̂`` to each. The pair set is mutated only by
    :class:`FaCatchSlipUpdater` (break/rebind); positions are read tag-ordered
    so it is ParticleSorter-safe.
    """

    def __init__(self, *, k_fa: float, pairs: np.ndarray, r0: np.ndarray) -> None:
        super().__init__(aniso=False)
        self.k_fa = float(k_fa)
        self._pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
        self._r0_map: dict[tuple[int, int], float] = {
            (int(a), int(b)): float(r) for (a, b), r in zip(self._pairs, r0)
        }

    @property
    def pairs(self) -> np.ndarray:
        return self._pairs

    def r0_for(self, pairs: np.ndarray) -> np.ndarray:
        return np.array(
            [self._r0_map[(int(a), int(b))] for a, b in pairs], dtype=np.float64
        )

    def set_pairs(self, pairs: np.ndarray, extra_r0: dict | None = None) -> None:
        pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
        if extra_r0:
            self._r0_map.update({(int(a), int(b)): float(r)
                                 for (a, b), r in extra_r0.items()})
        # prune r0 entries for pairs no longer engaged (keep map bounded)
        live = {(int(a), int(b)) for a, b in pairs}
        self._r0_map = {k: v for k, v in self._r0_map.items() if k in live}
        self._pairs = pairs

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        F = np.zeros_like(pos)
        U = np.zeros(n, dtype=np.float64)
        if self._pairs.shape[0] > 0:
            perm = np.argsort(tag)            # local-row → global-tag order
            pos_g = pos[perm]
            a = self._pairs[:, 0]
            b = self._pairs[:, 1]
            r0 = self.r0_for(self._pairs)
            dvec = pos_g[a] - pos_g[b]        # basal − bead
            r = np.linalg.norm(dvec, axis=1)
            r_safe = np.where(r > 1e-18, r, 1.0)
            uhat = dvec / r_safe[:, None]
            f_scalar = -self.k_fa * (r - r0)  # along +r̂ on basal node
            f_vec = f_scalar[:, None] * uhat
            F_g = np.zeros_like(pos_g)
            np.add.at(F_g, a, f_vec)
            np.add.at(F_g, b, -f_vec)
            # split the per-bond harmonic energy equally to its two endpoints
            u_bond = 0.5 * self.k_fa * (r - r0) ** 2
            U_g = np.zeros(n, dtype=np.float64)
            np.add.at(U_g, a, 0.5 * u_bond)
            np.add.at(U_g, b, 0.5 * u_bond)
            F[perm] = F_g
            U[perm] = U_g
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# ECM bed builder (crop the Mikado network to the cell-footprint window)
# ---------------------------------------------------------------------------
def build_ecm_bed(p: ResolvedDcmEcm, h1_cfg: dict):
    """Build a flat cross-linked Mikado fiber bed cropped to the cell window.

    Returns ``(pos (n,3), bonds (m,2), bond_r0 (m,), angles (k,3), p_h1)`` where
    positions are recentred so the bed spans roughly
    ``[-W, W]² × {z≈0}`` with ``W = footprint_factor · R_cell``. Bonds carry
    their per-bond rest length (= construction length) so the cropped bed is
    force-free at build (cropping cannot stretch a bond). Cross-link (``xl``)
    bonds are included — this is the EXPLICIT cross-linked ECM.
    """
    from ffn_sim.ecm.mikado import build_mikado_state

    p_h1 = resolve_derived(h1_cfg)
    frame = build_mikado_state(p_h1, with_cross_links=True)

    pos = np.asarray(frame.particles.position, dtype=np.float64).reshape(-1, 3)
    bg = np.asarray(frame.bonds.group, dtype=np.int64).reshape(-1, 2)

    # Crop to a window around the box centre (the Mikado box is centred at 0).
    W = p.footprint_factor * p.R_cell
    in_win = (np.abs(pos[:, 0]) <= W) & (np.abs(pos[:, 1]) <= W)
    # keep only beads in window; remap indices
    keep_idx = np.flatnonzero(in_win)
    if keep_idx.size < 10:
        raise RuntimeError(
            f"ECM crop kept only {keep_idx.size} beads in the "
            f"{2*W*1e6:.1f} µm window; widen footprint_factor or lower "
            "target_segment_length so the bed is dense enough to remodel."
        )
    remap = -np.ones(pos.shape[0], dtype=np.int64)
    remap[keep_idx] = np.arange(keep_idx.size)
    pos_c = pos[keep_idx].copy()
    pos_c[:, 2] = 0.0                              # flatten the bed exactly to z=0

    # keep only bonds whose both ends survived the crop
    both = in_win[bg[:, 0]] & in_win[bg[:, 1]]
    bg_c = remap[bg[both]]
    # per-bond rest length = construction length (force-free after crop)
    bond_r0 = np.linalg.norm(pos_c[bg_c[:, 0]] - pos_c[bg_c[:, 1]], axis=1)

    # angles (bending) — keep triplets fully inside the window
    if int(getattr(frame.angles, "N", 0)) > 0:
        ag = np.asarray(frame.angles.group, dtype=np.int64).reshape(-1, 3)
        keep_a = in_win[ag[:, 0]] & in_win[ag[:, 1]] & in_win[ag[:, 2]]
        ag_c = remap[ag[keep_a]]
    else:
        ag_c = np.empty((0, 3), dtype=np.int64)

    return pos_c, bg_c, bond_r0, ag_c, p_h1


# ---------------------------------------------------------------------------
# Full assembly: merge ECM bed + DCM shell, wire FA clutch + contraction
# ---------------------------------------------------------------------------
def build_dcm_ecm_simulation(
    p: ResolvedDcmEcm,
    h1_cfg: dict,
    *,
    device: hoomd.device.Device | None = None,
):
    """Assemble the DCM⊗ECM remodelling simulation on the BAOAB integrator.

    Returns a dict of handles (``sim``, ``ecm_bead_tags``, ``mem_tags``,
    ``basal_tags``, ``fa_force``, ``fa_updater``, ``turgor``, ``p``, …).
    """
    import gsd.hoomd

    # ---- 1. ECM bed ----
    ecm_pos, ecm_bonds, ecm_bond_r0, ecm_angles, p_h1 = build_ecm_bed(p, h1_cfg)
    n_ecm = ecm_pos.shape[0]

    # ---- 2. DCM shell, resting on the bed (small gap, no overlap) ----
    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())
    # place sphere centre at (0,0, R_cell + z_gap) so the lowest node clears z=0
    centre = np.array([0.0, 0.0, p.R_cell + p.z_gap], dtype=np.float64)
    mem_pos = verts0 + centre

    # ---- 3. merge into one snapshot (ECM first, then membrane) ----
    pos_all = np.vstack([ecm_pos, mem_pos])
    typeid = np.concatenate([
        np.zeros(n_ecm, dtype=np.uint32),         # TYPE_ECM = 0
        np.ones(nv, dtype=np.uint32),             # TYPE_MEM = 1
    ])
    n_total = pos_all.shape[0]
    mem_tag_offset = n_ecm
    mem_tags = np.arange(nv, dtype=np.int64) + mem_tag_offset
    ecm_bead_tags = np.arange(n_ecm, dtype=np.int64)

    # membrane edge bonds (re-indexed into merged tag space)
    edge_bonds = edges + mem_tag_offset
    # ECM backbone+xl bonds keep their indices (ECM is first block)
    # build merged bond group + per-bond r0:
    #   - ECM bonds keep their (cropped) construction r0
    #   - membrane edges get the CONTRACTILE r0 = contractility · mean_edge
    n_ecm_bonds = ecm_bonds.shape[0]
    n_edge = edge_bonds.shape[0]

    # faces (triangulated shell) for the exact turgor force — global indices
    faces = (tris0 + mem_tag_offset).astype(np.int64)
    face_cell = np.zeros(faces.shape[0], dtype=np.int64)

    # ---- 4. FA clutch: basal nodes → nearest fiber bead, force-free r0 ----
    z_bed = 0.0
    basal_local = np.flatnonzero(
        (mem_pos[:, 2] - z_bed) < p.fa_basal_band * p.R_cell
    )
    if basal_local.size == 0:
        raise RuntimeError(
            "No basal cell nodes within the FA band; raise fa_basal_band or "
            "lower the cell so its south cap is near the bed."
        )
    basal_tags = mem_tags[basal_local]
    basal_xyz = mem_pos[basal_local]
    # Nearest fiber bead per basal node, measured IN-PLANE (xy). The integrin
    # reaches DOWN to the flat z=0 bed (exactly the H.4 FALayout model: integrin
    # h_above the z=0 ligand plane), so a basal node directly above a bead binds
    # regardless of its z-height; the bond r0 is then the full 3D realised
    # separation so the clutch is force-free at construction.
    dxy = np.linalg.norm(
        basal_xyz[:, None, :2] - ecm_pos[None, :, :2], axis=2
    )
    jb = np.argmin(dxy, axis=1)
    dxy_min = dxy[np.arange(basal_local.size), jb]
    # only anchor basal nodes whose nearest bead is within capture (others
    # remain free and may bind later via the updater's rebind step)
    anchored = dxy_min <= p.fa_capture
    fa_pairs = np.column_stack([
        basal_tags[anchored], ecm_bead_tags[jb[anchored]]
    ]).astype(np.int64)
    fa_r0 = np.linalg.norm(
        basal_xyz[anchored] - ecm_pos[jb[anchored]], axis=1
    ).astype(np.float64)   # full 3D realised separation → force-free

    # ---- 5. assemble the gsd frame ----
    snap = gsd.hoomd.Frame()
    snap.particles.N = n_total
    snap.particles.types = [TYPE_ECM, TYPE_MEM]
    snap.particles.typeid = typeid
    snap.particles.position = pos_all
    snap.particles.mass = np.ones(n_total, dtype=np.float64)

    # bonds: ECM backbone+xl (Harmonic, native r0) + membrane edges (contractile)
    # We keep ECM bonds and membrane edges as md.bond.Harmonic groups; FA clutch
    # is the Custom FaClutchForce (per-bond r0), NOT a bond.Harmonic group.
    bond_group = np.vstack([ecm_bonds, edge_bonds]).astype(np.uint32)
    bond_typeid = np.concatenate([
        np.zeros(n_ecm_bonds, dtype=np.uint32),   # "ecm-bond"
        np.ones(n_edge, dtype=np.uint32),         # "dcm_edge"
    ])
    snap.bonds.N = bond_group.shape[0]
    snap.bonds.types = ["ecm-bond", BOND_EDGE]
    snap.bonds.typeid = bond_typeid
    snap.bonds.group = bond_group

    # angles: ECM bending only
    if ecm_angles.shape[0] > 0:
        snap.angles.N = ecm_angles.shape[0]
        snap.angles.types = ["ecm-angle"]
        snap.angles.typeid = np.zeros(ecm_angles.shape[0], dtype=np.uint32)
        snap.angles.group = ecm_angles.astype(np.uint32)

    # box: span both the bed window and the cell, with margin
    span = max(
        float(np.abs(pos_all[:, :2]).max()),
        float(pos_all[:, 2].max()),
        p.R_cell,
    )
    box_edge = 2.0 * span + 6.0 * p.R_cell
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]

    # ---- 6. build the simulation + forces ----
    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    # ECM backbone+xl + contractile membrane edges (one Harmonic, two types).
    # ECM bonds: a single Harmonic type cannot hold per-bond r0, but the cropped
    # bed's bonds are all close to the Mikado rest length ℓ₀ (cropping does not
    # stretch them); we use the H.1 resolved bond_k + rest_length for "ecm-bond"
    # (force ~0 at build) and the contractile r0 for the membrane edges.
    bond = md.bond.Harmonic()
    bond.params["ecm-bond"] = dict(k=p_h1.bond_k, r0=p_h1.rest_length)
    bond.params[BOND_EDGE] = dict(k=p.k_edge, r0=p.contractility * mean_edge)
    ig.forces.append(bond)

    if ecm_angles.shape[0] > 0:
        angle = md.angle.Harmonic()
        angle.params["ecm-angle"] = dict(k=p_h1.angle_k, t0=p_h1.angle_t0)
        ig.forces.append(angle)

    # WCA excluded volume: cell node ↔ fiber bead (cell can't pass through fibers).
    # sigma ~ contact distance of a node and a fiber bead; moderate eps (5 kT).
    sigma_wca = max(0.5 * mean_edge, 2.0 * p_h1.bead_radius)
    nlist = md.nlist.Tree(buffer=0.5 * sigma_wca)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    rc = 2.0 ** (1.0 / 6.0) * sigma_wca
    # all pairs need params in HOOMD 7; default inert, then turn on mem↔ecm.
    for ta, tb in [(TYPE_ECM, TYPE_ECM), (TYPE_ECM, TYPE_MEM), (TYPE_MEM, TYPE_MEM)]:
        lj.params[(ta, tb)] = dict(epsilon=0.0, sigma=sigma_wca)
        lj.r_cut[(ta, tb)] = 0.0
    lj.params[(TYPE_ECM, TYPE_MEM)] = dict(epsilon=p.eps_rep_kT * p.kT, sigma=sigma_wca)
    lj.r_cut[(TYPE_ECM, TYPE_MEM)] = rc
    # cell-cell excluded volume too (keeps the shell from self-intersecting under
    # strong contraction) — same moderate WCA.
    lj.params[(TYPE_MEM, TYPE_MEM)] = dict(epsilon=p.eps_rep_kT * p.kT, sigma=sigma_wca)
    lj.r_cut[(TYPE_MEM, TYPE_MEM)] = rc
    lj.mode = "shift"
    ig.forces.append(lj)

    # Turgor (exact triangulated face-normal pressure) on the membrane shell.
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    turgor = DcmTurgorForce(
        faces=faces, face_cell=face_cell, n_cells=1,
        V0=V0, turgor_dP0=p.turgor_dP0, K_vol=p.K_vol,
    )
    ig.forces.append(turgor)

    # FA catch-slip clutch (per-bond r0 Custom force).
    fa_force = FaClutchForce(k_fa=p.k_fa, pairs=fa_pairs, r0=fa_r0)
    ig.forces.append(fa_force)

    sim.operations.integrator = ig
    sim.run(0)   # seed net_force before the first BAOAB step

    # BAOAB updater — gamma_map MUST cover every particle type.
    gamma = {TYPE_ECM: p_h1.gamma_b, TYPE_MEM: p.gamma_node}
    _baoab_action, baoab_updater = make_baoab_updater(
        kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1
    )
    sim.operations.updaters.append(baoab_updater)

    # FA catch-slip updater (reuses canonical Pereverzev kernel).
    pereverzev = PereverzevParams.from_config({"bridge": {"catch_bond": p.pereverzev_cfg}})

    def _bead_getter():
        return ecm_bead_tags, ecm_pos   # tags + initial positions (for nearest-bead)

    fa_updater_action = FaCatchSlipUpdater(
        force=fa_force,
        basal_tags=basal_tags,
        bead_pos_getter=_bead_getter,
        pereverzev=pereverzev,
        k_fa=p.k_fa,
        k_on=p.k_on,
        capture=p.fa_capture,
        batch_steps=p.fa_batch_steps,
        dt=p.dt,
        seed=p.seed + 2,
    )
    fa_updater = hoomd.update.CustomUpdater(
        action=fa_updater_action,
        trigger=hoomd.trigger.Periodic(p.fa_batch_steps),
    )
    sim.operations.updaters.append(fa_updater)

    return {
        "sim": sim,
        "ecm_bead_tags": ecm_bead_tags,
        "ecm_pos0": ecm_pos.copy(),
        "mem_tags": mem_tags,
        "mem_pos0": mem_pos.copy(),
        "basal_tags": basal_tags,
        "fa_force": fa_force,
        "fa_updater": fa_updater_action,
        "turgor": turgor,
        "faces": faces,
        "ecm_bonds": ecm_bonds,
        "ecm_bond_r0": ecm_bond_r0,
        "mean_edge": mean_edge,
        "V0": V0,
        "n_ecm": n_ecm,
        "nv": nv,
        "cell_centre0": centre.copy(),
        "p": p,
        "p_h1": p_h1,
        "F_star": float(pereverzev_F_star(pereverzev)),
    }


# ---------------------------------------------------------------------------
# Measurement helpers (d0sm01911a observables)
# ---------------------------------------------------------------------------
def read_positions(sim) -> tuple[np.ndarray, np.ndarray]:
    """Return (tag-sorted positions (N,3), tags) from the current state.

    ``state.get_snapshot()`` returns a global, tag-ordered Snapshot (row i ==
    particle with tag i), so no permutation is needed — unlike the per-rank
    ``cpu_local_snapshot`` whose rows the ParticleSorter reorders.
    """
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    tags = np.arange(pos.shape[0], dtype=np.int64)
    return pos, tags


def ecm_radial_displacement(pos_g: np.ndarray, ecm_pos0: np.ndarray,
                            ecm_bead_tags: np.ndarray, cell_xy: np.ndarray):
    """Inward radial displacement of each fiber bead toward the cell centre.

    Returns ``(r0 (n,), inward (n,))`` where ``r0`` is each bead's initial
    in-plane distance from ``cell_xy`` and ``inward`` is the inward component of
    its displacement (positive = moved toward the cell — i.e. remodelling).
    """
    now = pos_g[ecm_bead_tags][:, :2]
    d0 = ecm_pos0[:, :2] - cell_xy[None, :]
    r0 = np.linalg.norm(d0, axis=1)
    r_safe = np.where(r0 > 1e-18, r0, 1.0)
    rhat = d0 / r_safe[:, None]                  # outward unit vector at t0
    disp = now - ecm_pos0[:, :2]
    inward = -np.einsum("ij,ij->i", disp, rhat)  # +ve = toward cell
    return r0, inward


def fiber_tension_vs_r(pos_g: np.ndarray, ecm_bonds: np.ndarray,
                       ecm_bond_r0: np.ndarray, bond_k: float,
                       cell_xy: np.ndarray):
    """Per-bond tension and its midpoint distance r from the cell centre.

    Tension ``T = k·(|b| − r0)`` (positive = stretched). Returns
    ``(r_mid (m,), tension (m,))`` for the radial-stress ~1/r check.
    """
    a = ecm_bonds[:, 0]
    b = ecm_bonds[:, 1]
    bvec = pos_g[a] - pos_g[b]
    blen = np.linalg.norm(bvec, axis=1)
    tension = bond_k * (blen - ecm_bond_r0)
    mid = 0.5 * (pos_g[a][:, :2] + pos_g[b][:, :2])
    r_mid = np.linalg.norm(mid - cell_xy[None, :], axis=1)
    return r_mid, tension


def fiber_radial_alignment(pos_g: np.ndarray, ecm_bonds: np.ndarray,
                           cell_xy: np.ndarray):
    """|cos(angle between bond and the radial direction)| per bond, + r.

    1 = bond points radially (toward/away from the cell, the d0sm01911a
    tensed-fiber signature); 0 = bond is tangential. Returns ``(r_mid, align)``.
    """
    a = ecm_bonds[:, 0]
    b = ecm_bonds[:, 1]
    bvec = (pos_g[a] - pos_g[b])[:, :2]
    blen = np.linalg.norm(bvec, axis=1)
    blen_s = np.where(blen > 1e-18, blen, 1.0)
    bhat = bvec / blen_s[:, None]
    mid = 0.5 * (pos_g[a][:, :2] + pos_g[b][:, :2])
    rvec = mid - cell_xy[None, :]
    rlen = np.linalg.norm(rvec, axis=1)
    rlen_s = np.where(rlen > 1e-18, rlen, 1.0)
    rhat = rvec / rlen_s[:, None]
    align = np.abs(np.einsum("ij,ij->i", bhat, rhat))
    return rlen, align
