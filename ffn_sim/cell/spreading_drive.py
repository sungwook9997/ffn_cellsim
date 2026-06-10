"""H.7 — membrane-tracked leading-edge nucleation + FA molecular-clutch ratchet
for the DYNAMIC single-cell spreading run.

Two live ``hoomd.custom.Action`` updaters that close the gaps which made the
H.5 lamellipodium stall instead of spread:

1. :class:`LeadingEdgeNucleationUpdater` — membrane-recruited Arp2/3 nucleation
   at the ADVANCING front.

   The H.5 lamellipodium pins its WAVE / NPF reservoir at the construction
   contact ring, and :class:`ffn_sim.cell.lamellipodium.ArpBranchingUpdater`
   only nucleates daughters within ``r_branch_eff`` (~100 nm, the Arp2/3
   physical reach) of a WAVE bead. Once the protrusive front advances more than
   one reach beyond the fixed ring, no new barbed ends nucleate, the existing
   ends cap out, and spreading stalls (measured A/A₀ ≈ 1.56 plateau at the
   physiological cap/elong balance). In vivo the Arp2/3 / WAVE machinery is
   continuously recruited to the advancing membrane edge, so nucleation TRACKS
   the front. This updater is the mesoscale proxy for that membrane-recruited
   nucleation: each batch it promotes the outermost basal lamellipodial actin
   beads back to barbed-end status with an outward-radial tangent, so the
   existing :class:`~ffn_sim.cell.lamellipodium.BarbedEndElongationUpdater`
   keeps advancing the front. It mutates ONLY the Python
   :class:`~ffn_sim.cell.lamellipodium.LamellipodiumState` bookkeeping
   (``barbed_end_tags`` / ``tangent_of`` / ``capped_tags``) — never particle
   positions — so the BAOAB integrator remains the sole stepper (no
   double-stepping; CLAUDE.md hard rule).

2. :class:`FrontClutchRatchetUpdater` — the FA molecular clutch at the front.

   The Explore audit found the spreading loop OPEN: lamellipodial actin grows
   freely along its tangents and the FA integrin↔ligand clutches engage the
   substrate, but the two are never coupled — the clutch does not grip the
   advancing actin, so retrograde flow is never converted to anchored
   protrusion (the defining step of the Chan-Odde / Elosegui-Artola motor-clutch
   model). This updater closes the loop: each batch, for an advancing front
   actin bead within ``grip_radius`` of a substrate ligand, it adds a
   FORCE-FREE harmonic grip bond (``front_clutch_grip``) anchoring that bead to
   the substrate (the engaged molecular clutch). The gripped front is thereby
   ratcheted onto the substrate — the adhesion footprint. Toggling this updater
   off is the clutch-OFF control that isolates the clutch's mechanical role.

Both are batched D2-style Updaters (fire every ``batch_steps`` integration
steps), mirroring the lamellipodium / integrin updater convention.
"""

from __future__ import annotations

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md


# ---------------------------------------------------------------------------
# Basal adhesion tether (the FA molecular-clutch ensemble, as a substrate force)
# ---------------------------------------------------------------------------
class BasalAdhesionTether(md.force.Custom):
    """Hold the lamellipodial sheet on the substrate — the FA clutch ensemble.

    The lamellipodium is a thin BASAL sheet that a spreading cell keeps pressed
    onto the substrate through its integrin focal adhesions (the molecular
    clutch). In the explicit build the advancing lamellipodial actin is NOT yet
    bond-coupled to a substrate ligand field, so without a holding force the
    cortex / membrane line tension lifts and retracts the protruded front out of
    the basal plane and spreading stalls (the clutch-OFF failure mode, observed
    as the A/A₀≈2.1 staircase plateau). This Custom force is the mesoscale
    representation of the engaged FA clutch ensemble: a harmonic z-tether pulling
    every ``actin_lamel`` bead toward the basal contact plane ``z_basal``,

        ``U_i = ½ k_adh (z_i − z_basal)²``   (z-component only),

    with ``k_adh`` derived from the FA clutch stiffness ``k_int`` (the engaged
    integrin-clutch spring, KU-2.4) — a force-free tether for a bead already in
    the plane, so it injects no energy at construction. Toggling it off is the
    clutch-OFF control that isolates the clutch's mechanical role in spreading.
    Force is z-only, so it does NOT bias the in-plane footprint growth (the
    spreading observable) — it only keeps the sheet basal.
    """

    def __init__(self, *, z_basal: float, k_adh: float, actin_typeid: int):
        super().__init__(aniso=False)
        self.z_basal = float(z_basal)
        self.k_adh = float(k_adh)
        self.actin_typeid = int(actin_typeid)  # global type index for actin_lamel

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        # The cpu_local_snapshot exposes per-rank typeid (global type indices) +
        # positions, but NOT the global type-name list — so we compare against a
        # precomputed actin_lamel typeid resolved at attach time.
        with self._state.cpu_local_snapshot as snap:
            tid = np.asarray(snap.particles.typeid).copy()
            pos = np.asarray(snap.particles.position).copy()
        mask = tid == self.actin_typeid
        dz = pos[:, 2] - self.z_basal
        F = np.zeros_like(pos)
        F[:, 2] = -self.k_adh * dz
        U = 0.5 * self.k_adh * dz * dz
        F[~mask] = 0.0
        U[~mask] = 0.0
        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F
            arrays.potential_energy[:] = U


# ---------------------------------------------------------------------------
# Leading-edge nucleation (membrane-tracked Arp2/3 proxy)
# ---------------------------------------------------------------------------
class LeadingEdgeNucleationUpdater(hoomd.custom.Action):
    """Membrane-tracked protrusion engine: advance the lamellipodial front by
    appending actin beads outward at the leading edge.

    Why this is needed (and not a mock): the H.5 Bell-Evans elongation updater is
    batch-CFL-throttled to ≤ ~10⁻³ monomer-addition events per barbed end per
    tick (the D2 accuracy bound), so advancing a leading filament ~5 µm needs
    ~10⁴ ticks — and the Arp2/3 nucleation is pinned at the construction WAVE
    ring, so the front caps out and stalls at A/A₀ ≈ 1.5. Brute-forcing the
    physiological spread to A/A₀ ≈ 3 is ~10⁸ BAOAB steps (infeasible on CPU in
    the window). Per PI authorization (2026-06-10, band-match literature),
    this updater drives the protrusion at a band-matched effective velocity
    ``v_front`` (literature single-cell spreading 3–12 µm/min): each tick, for
    each leading-edge tip, it appends a new ``actin_lamel`` bead one ℓ₀ outward
    along the radial tangent with probability ``p = v_front · Δt_tick · S / ℓ₀``
    (S = the run's kinetic acceleration). The appended bead is FORCE-FREE (rest
    length along the tangent) and is thereafter a full BAOAB degree of freedom —
    it relaxes under the cortex / turgor / membrane / FA-clutch forces like any
    other particle. So the FRONT VELOCITY is the single band-matched constitutive
    input (exactly as in continuum spreading models, Betorz 2023); everything
    downstream (network mechanics, clutch grip, traction, footprint shape) is
    emergent BAOAB physics. Mutates topology via the lamellipodium's own
    ``_extend_snapshot_with_new_actins`` helper + the Python state bookkeeping;
    never moves an existing particle (no double-stepping).

    Parameters
    ----------
    lamel_state : LamellipodiumState
        The live lamellipodium bookkeeping (mutated in place).
    z_basal, band : float
        Basal contact plane [m] and contact-band half-width [m]; only tips
        within ``band`` of ``z_basal`` advance (the contact lamella).
    rest_length : float
        ℓ₀ [m] — the appended monomer step (the lamellipodium backbone bond r0).
    p_advance : float
        Per-tip append probability per tick = v_front·Δt_tick·S/ℓ₀, clamped to
        [0, 1]. The harness computes it from the band-matched v_front.
    advance_margin : float
        A tip qualifies as "leading" if its in-plane radius is within this [m]
        of the current max front radius (advance the genuine edge, not the bulk).
    max_advance : int
        Max tips advanced per tick (bounds density / cost).
    """

    def __init__(
        self, *, lamel_state, z_basal: float, band: float, rest_length: float,
        p_advance: float = 1.0, advance_margin: float = 1.0e-6,
        max_advance: int = 60, max_radius: float = float("inf"), seed: int = 7,
    ) -> None:
        super().__init__()
        self.state = lamel_state
        self.z_basal = float(z_basal)
        self.band = float(band)
        self.rest_length = float(rest_length)
        self.p_advance = float(min(max(p_advance, 0.0), 1.0))
        self.advance_margin = float(advance_margin)
        self.max_advance = int(max_advance)
        self.max_radius = float(max_radius)  # cap the footprint (physiological
        # maximal spread + keep the front inside the box; no wall artefacts)
        self._rng = np.random.default_rng(seed)
        self._sim_ref: hoomd.Simulation | None = None
        self._n_promoted = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        from ffn_sim.cell.lamellipodium import _extend_snapshot_with_new_actins

        sim = self._sim_ref
        assert sim is not None
        snap = sim.state.get_snapshot()
        if snap.communicator.rank != 0:
            return
        types = list(snap.particles.types)
        if "actin_lamel" not in types:
            return
        tid = np.asarray(snap.particles.typeid)
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        actin_rows = np.where(tid == types.index("actin_lamel"))[0]  # row == tag
        if actin_rows.size == 0:
            return
        ap = pos[actin_rows]
        near = np.abs(ap[:, 2] - self.z_basal) <= self.band
        if not near.any():
            return
        rows = actin_rows[near]
        rad = np.hypot(ap[near, 0], ap[near, 1])

        # Leading-edge tips = genuine chain ends (NOT a parent of any bead) in
        # the outer radial shell. Appending outward from a true tip lands in open
        # space (force-free); appending from an interior bead would overlap the
        # network → BAOAB blow-up, so interior beads are excluded.
        parents = set(self.state.parent_of.values())
        r_lead = float(rad.max())
        tip = (rad >= r_lead - self.advance_margin) & np.array(
            [int(t) not in parents for t in rows]
        )
        if not tip.any():
            return
        cand_rows = rows[tip]
        cand_rad = rad[tip]
        order = np.argsort(cand_rad)[::-1][: self.max_advance]

        new_positions, new_bonds = [], []
        next_tag = int(snap.particles.N)
        for j in order:
            if self._rng.uniform() >= self.p_advance:
                continue
            tag = int(cand_rows[j])
            rxy = pos[tag, :2]
            rnorm = float(np.hypot(rxy[0], rxy[1]))
            if rnorm < 1.0e-12:
                continue
            if rnorm + self.rest_length > self.max_radius:
                continue  # footprint cap reached (maximal spread / box guard)
            tangent = np.array([rxy[0] / rnorm, rxy[1] / rnorm, 0.0])
            r_new = pos[tag] + self.rest_length * tangent
            new_tag = next_tag + len(new_positions)
            new_positions.append(r_new)
            new_bonds.append((tag, new_tag))
            # Bookkeeping: old tip becomes interior parent; new bead is the tip.
            self.state.parent_of[new_tag] = tag
            self.state.tangent_of[new_tag] = tangent
            self.state.capped_tags.discard(tag)
            try:
                self.state.barbed_end_tags.remove(tag)
            except ValueError:
                pass
            self.state.barbed_end_tags.append(new_tag)
        if not new_positions:
            return
        write = _extend_snapshot_with_new_actins(
            snap, new_positions, new_bonds,
            actin_type_name="actin_lamel", bond_type_name="lamel_actin_bond",
        )
        sim.state.set_snapshot(write)
        self._n_promoted += len(new_positions)

    @property
    def n_promoted(self) -> int:
        return self._n_promoted


# ---------------------------------------------------------------------------
# Front molecular-clutch ratchet (closes the actin↔clutch↔substrate loop)
# ---------------------------------------------------------------------------
class FrontClutchRatchetUpdater(hoomd.custom.Action):
    """Grip advancing front actin to the substrate via force-free clutch bonds.

    Each batch tick, for an advancing lamellipodial actin bead in the contact
    band whose nearest substrate ligand is within ``grip_radius`` and that is
    not yet gripped, add a force-free harmonic ``bond_type_name`` bond
    (integrin/ligand-style clutch) anchoring the bead to that ligand. This is
    the engaged molecular clutch: the gripped front is ratcheted onto the
    substrate, converting protrusion into a stable adhesion footprint.

    The bond ``r0`` is the current bead↔ligand separation (EXACT-r0, the
    platform's force-free-at-construction convention) so adding it never injects
    energy or trips the BAOAB guard. ``bond_type_name`` must already be
    registered on the simulation's ``md.bond.Harmonic`` force with a finite k
    (wired by the harness before the first run).
    """

    def __init__(
        self, *, z_basal: float, band: float, grip_radius: float,
        bond_type_name: str = "front_clutch_grip", max_grip_per_tick: int = 30,
    ) -> None:
        super().__init__()
        self.z_basal = float(z_basal)
        self.band = float(band)
        self.grip_radius = float(grip_radius)
        self.bond_type_name = str(bond_type_name)
        self.max_grip_per_tick = int(max_grip_per_tick)
        self._sim_ref: hoomd.Simulation | None = None
        self._gripped_actin: set[int] = set()
        self._n_grips = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    @property
    def n_grips(self) -> int:
        return self._n_grips

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        snap = sim.state.get_snapshot()
        if snap.communicator.rank != 0:
            return
        types = list(snap.particles.types)
        btypes = list(snap.bonds.types)
        if "actin_lamel" not in types or "ligand" not in types:
            return
        if self.bond_type_name not in btypes:
            raise RuntimeError(
                f"bond type {self.bond_type_name!r} not registered; got {btypes}"
            )
        tid = np.asarray(snap.particles.typeid)
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        actin_rows = np.where(tid == types.index("actin_lamel"))[0]
        lig_rows = np.where(tid == types.index("ligand"))[0]
        if actin_rows.size == 0 or lig_rows.size == 0:
            return
        ap = pos[actin_rows]
        near = np.abs(ap[:, 2] - self.z_basal) <= self.band
        cand_actin = actin_rows[near]
        if cand_actin.size == 0:
            return
        cand_pos = pos[cand_actin]
        lig_pos = pos[lig_rows]

        new_bonds = []
        new_r0 = []
        n_added = 0
        for ai, arow in enumerate(cand_actin):
            atag = int(arow)
            if atag in self._gripped_actin:
                continue
            d = np.linalg.norm(lig_pos - cand_pos[ai], axis=1)
            jmin = int(np.argmin(d))
            if d[jmin] <= self.grip_radius:
                new_bonds.append((atag, int(lig_rows[jmin])))
                new_r0.append(float(d[jmin]))
                self._gripped_actin.add(atag)
                n_added += 1
                if n_added >= self.max_grip_per_tick:
                    break
        if not new_bonds:
            return

        # Append the new force-free grip bonds to the topology.
        bg = np.asarray(snap.bonds.group, dtype=np.int64)
        bt = np.asarray(snap.bonds.typeid, dtype=np.uint32)
        gid = btypes.index(self.bond_type_name)
        add_bg = np.array(new_bonds, dtype=np.int64)
        add_bt = np.full(add_bg.shape[0], gid, dtype=np.uint32)
        merged_bg = (np.concatenate([bg, add_bg], axis=0) if bg.size
                     else add_bg).astype(np.uint32)
        merged_bt = (np.concatenate([bt, add_bt]) if bt.size else add_bt).astype(np.uint32)

        write = hoomd.Snapshot()
        write.particles.N = int(snap.particles.N)
        write.particles.types = list(types)
        write.particles.typeid[:] = np.asarray(snap.particles.typeid)
        write.particles.position[:] = pos
        write.particles.velocity[:] = np.asarray(snap.particles.velocity)
        write.particles.mass[:] = np.asarray(snap.particles.mass)
        write.particles.image[:] = np.asarray(snap.particles.image)
        write.configuration.box = list(snap.configuration.box)
        write.bonds.N = int(merged_bg.shape[0])
        write.bonds.types = list(btypes)
        write.bonds.group[:] = merged_bg
        write.bonds.typeid[:] = merged_bt
        for grp in ("angles", "dihedrals", "impropers"):
            src = getattr(snap, grp)
            dst = getattr(write, grp)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)
        sim.state.set_snapshot(write)
        self._n_grips += n_added
