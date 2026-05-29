"""D2 Pereverzev catch-slip dynamic integrin↔ligand bond updater (KU-2.5).

Phase 1 H.4 brief §Integrin off-rate. The updater is a
``hoomd.custom.Action`` wrapped in a periodic ``CustomUpdater`` that runs
every ``integrin_batch_steps`` integration steps (D2 batched-updater
convention; gated by ``gate_bell_evans_batch_cfl``).

Per batch tick the action:

1. Reads the current particle positions via ``sim.state.get_snapshot()``.
2. For every **currently engaged** integrin bond:
   - Computes the bond force ``F = k_int^eff · (|r_int − r_lig| − r0)``
     directly from the snapshot (HOOMD's per-bond force accessor returns
     the same value but requires a force re-eval; reading positions is
     cheaper and matches the analytic harmonic).
   - Samples the Pereverzev break probability
     ``p_break = 1 − exp(−k_off(|F|) · Δt_batch)``.
   - On fire: marks the integrin unbound and removes the bond from the
     topology buffer.
3. For every **currently unbound** integrin within ``capture_radius_R_FA``
   of a free ligand: samples ``p_on = 1 − exp(−k_on · Δt_batch)``;
   on fire, registers the new bond.
4. Rebuilds ``bonds.group / bonds.typeid`` from the updated engaged mask
   and calls ``sim.state.set_snapshot(snap)``.

D2 batch CFL contract
---------------------
The brief D2 specification (also ``gate_bell_evans_batch_cfl``) requires
``batch_steps · Δt · k_off_max ≤ 10⁻³``. ``ResolvedH4.resolve_h4`` enforces
this at config-resolve time by shrinking ``integrin_batch_steps`` if the
bound is violated. The updater's ``__init__`` re-asserts the bound so a
caller-constructed updater (outside ``build_h4_simulation``) cannot
silently violate D2.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks live in
``ffn_sim/tests/test_h4_topology.py``; the KU-2.5 emergent-vs-oracle gate
lives in ``ffn_sim/tests/validation/test_ku25_catch_peak.py``.*

1. **Dimensional analysis**
   - ``Δt_batch = batch_steps · dt`` [s].
   - ``k_off(F) · Δt_batch`` [dimensionless] → break probability.
   - ``k_on · Δt_batch`` [dimensionless] → bind probability.
   - All position math in metres, force in Newtons.

2. **Boundary cases**
   - No engaged bonds at startup: handled — engaged-loop is empty, only
     the binding loop runs.
   - All integrins out of capture range (FA centre too far from ligand
     plane): no binding occurs; the updater idles. RUNTIME pass.
   - ``k_off · Δt_batch ≥ 1`` (D2 violation): the break probability
     saturates at 1 (every engaged bond breaks every batch). Detected
     by the resolve_h4 re-shrink of batch_steps OR (if a caller
     bypasses resolve_h4) by ``__init__``'s assertion.

3. **Conservation invariants**
   - Particle count is invariant — only ``bonds.group`` mutates.
   - Engaged-clutch population is bounded by ``n_total_per_fa`` per FA.
   - Force on the ligand-anchored side equals minus the integrin-side
     force (Newton 3rd law; inherited from md.bond.Harmonic).

4. **Numerical sanity**
   - Per-bond positions read in float64; bond.group rebuilt as uint32
     (HOOMD GSD requirement).
   - RNG isolated per Action (``np.random.default_rng`` with the
     simulation seed offset by 1, so the BAOAB Updater's seed and the
     IntegrinBondUpdater's seed are independent — RNG-contamination
     check is a STATIC gate).
   - When ``set_snapshot`` rebuilds the neighbor list, any in-flight
     force computes are invalidated; we structure the act() so the next
     ``sim.run(N)`` re-evaluates forces before the next bond-event tick.

5. **Sign / sense**
   - Bond force magnitude ``F = k_int · max(0, |Δr| − r0)``; with r0=0
     the magnitude is just ``k_int · |Δr|``. Always ≥ 0 — fed to
     Pereverzev k_off(F ≥ 0) which is well-defined.
   - Break probability monotonically increases with F in the slip
     regime (F > F*) and decreases in the catch regime (F < F*); this
     IS the catch-slip signature, validated by KU-2.5 oracle gate.

6. **Measurement-protocol consistency**
   - The KU-2.5 catch-peak validation test (``test_ku25_catch_peak.py``)
     instruments the updater with a fixed-force harness: a single
     integrin held at distance ``F_target / k_int`` from a single
     ligand, run for ``n_steps`` batches, count break events,
     ``k_off_emergent = n_break / (n_steps · Δt_batch · n_alive)``.
     This matches the oracle's ``pereverzev_k_off(F_target)`` within
     ± 5 % per ``acceptance.pereverzev_rel_tolerance``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.validation.pereverzev import pereverzev_k_off

if TYPE_CHECKING:
    from ffn_sim.bridge.fa import FALayout, ResolvedH4


class IntegrinBondUpdater(hoomd.custom.Action):
    """Pereverzev catch-slip dynamic integrin↔ligand bond updater (D2).

    Parameters
    ----------
    p : ResolvedH4
        Resolved H.4 parameters. Provides batch_steps, k_int_bare,
        capture_radius_R_FA, k_on, Pereverzev params.
    layouts : list[FALayout]
        Per-FA tag bookkeeping (integrin tag ranges, ligand tags).
    bond : md.bond.Harmonic
        The HOOMD bond force; needed only so the Action can assert the
        per-type k matches ``p.k_int_bare`` (sanity check; the bond force
        owns the canonical k value during the run).
    seed_offset : int, default 1
        Offset added to the simulation seed when initialising this
        Action's RNG, so it is independent of the BAOAB Action's RNG.

    Notes
    -----
    The action does NOT mutate particle positions; it only mutates the
    bond topology. The L-M BAOAB Action remains the sole position
    integrator (CLAUDE.md hard rule: no double-stepping).
    """

    def __init__(
        self,
        *,
        p: "ResolvedH4",
        layouts: "list[FALayout]",
        bond: md.bond.Harmonic,
        seed_offset: int = 1,
    ) -> None:
        super().__init__()
        self.p = p
        self.layouts = list(layouts)
        self.bond = bond
        self._rng = np.random.default_rng(p.seed + seed_offset)

        batch_dt = p.integrin_batch_steps * p.dt
        events_max_per_batch = batch_dt * p.bond_event_rate_max
        if events_max_per_batch > 1.0e-3:
            raise RuntimeError(
                f"IntegrinBondUpdater violates D2 batch CFL: "
                f"batch_dt·k_off_max = {batch_dt:.3e}·{p.bond_event_rate_max:.3e}"
                f" = {events_max_per_batch:.3e} > 1e-3 ceiling. "
                f"Shrink integrin_batch_steps or pre-resolve via resolve_h4."
            )
        self._batch_dt = batch_dt
        self._n_FAs = len(layouts)

        # --- S5 tag-space unification (2026-05-30) -------------------------
        # Per-integrin runtime state (_engaged / _ligand_for_integrin /
        # _fa_for_integrin) is indexed by a CONTIGUOUS LOCAL index in
        # [0, n_int_total), NOT by the integrin's global HOOMD tag. Previously
        # these arrays were indexed directly by global tag (lo = L.
        # integrin_tag_start; self._engaged[lo:hi] = ...), which only works if
        # integrins occupy global tags [0, n_int_total). The whole-cell builder
        # places the focal-adhesion block in natural append order (AFTER cortex
        # / myosin / xlink / lamellipodium), so integrin tags are an arbitrary
        # contiguous-per-FA set, not [0, n_int).
        #
        # We therefore build an explicit global-tag -> local-index map (and its
        # inverse, local-index -> global-tag, as an array). Snapshot POSITIONS
        # are still read by global tag (pos[tag]) because a HOOMD snapshot is
        # tag-ordered (row i == particle with tag i). When integrins do occupy
        # [0, n_int) the local index equals the global tag, so this is
        # bit-for-bit identical to the previous behavior (the H.4 standalone
        # build_h4_simulation path is unaffected).
        n_int_total = sum(L.n_total for L in layouts)
        self._engaged = np.zeros(n_int_total, dtype=bool)
        # Map each integrin (LOCAL index) to its FA's ligand GLOBAL tag.
        self._ligand_for_integrin = np.full(n_int_total, -1, dtype=np.int64)
        # And to its FA index.
        self._fa_for_integrin = np.full(n_int_total, -1, dtype=np.int64)
        # LOCAL index -> GLOBAL integrin tag (used for snapshot position reads).
        self._integrin_global_tag = np.full(n_int_total, -1, dtype=np.int64)
        # GLOBAL integrin tag -> LOCAL state index (used to resolve bond cols).
        self._integrin_tag_to_local: dict[int, int] = {}
        local = 0
        for fa_idx, L in enumerate(layouts):
            for k in range(L.n_total):
                gtag = int(L.integrin_tag_start) + k
                self._engaged[local] = bool(L.initial_engaged[k])
                self._ligand_for_integrin[local] = int(L.ligand_tag)
                self._fa_for_integrin[local] = fa_idx
                self._integrin_global_tag[local] = gtag
                self._integrin_tag_to_local[gtag] = local
                local += 1

        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run: int = 0
        self._n_break_total: int = 0
        self._n_bind_total: int = 0

        # k_int contract: must match what the bond.Harmonic was wired with.
        try:
            k_wired = float(bond.params["integrin_ligand"]["k"])
        except (KeyError, AttributeError):
            k_wired = float("nan")
        if np.isfinite(k_wired) and not np.isclose(k_wired, p.k_int_bare):
            raise RuntimeError(
                f"bond.Harmonic.params['integrin_ligand']['k']={k_wired:.3e} "
                f"!= ResolvedH4.k_int_bare={p.k_int_bare:.3e}; the wiring "
                "is inconsistent. The Updater uses k_int_bare for analytic "
                "force computation."
            )

    # ------------------------------------------------------------------
    # HOOMD Action lifecycle
    # ------------------------------------------------------------------
    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None

        # HOOMD 7 ``state.get_snapshot()`` returns a read-only Snapshot
        # view. To mutate topology we must construct a *fresh*
        # ``hoomd.Snapshot()`` carrying the same particle data but with
        # an updated bonds block, then call ``state.set_snapshot()``.
        # This costs an O(N) particle-data copy per batch; with
        # batch_steps=100 and N ≈ 2.2k particles the per-tick cost is
        # tens of microseconds — negligible vs the integration cost of
        # the 100 BAOAB steps between ticks.
        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32).copy()
        bond_type_names = list(read_snap.bonds.types)

        try:
            integrin_type = bond_type_names.index("integrin_ligand")
        except ValueError:
            raise RuntimeError(
                "snap.bonds.types must include 'integrin_ligand'; got "
                f"{bond_type_names}."
            )

        # Partition current bonds: integrin_ligand vs others (motor etc.).
        is_int_bond = bt == integrin_type
        int_bonds = bg[is_int_bond]
        other_bg = bg[~is_int_bond]
        other_bt = bt[~is_int_bond]

        # ---- Step 1: break engaged bonds via Pereverzev k_off ----
        # int_bonds columns: [integrin_tag, ligand_tag] by construction
        # (we always insert in this order when binding below). The
        # bond.group convention in HOOMD GSD is (tag_a, tag_b) per bond;
        # both columns are stable across ParticleSorter reorderings
        # because they store tags.
        if int_bonds.shape[0] > 0:
            int_tags = int_bonds[:, 0]
            lig_tags = int_bonds[:, 1]
            r_int = pos[int_tags]
            r_lig = pos[lig_tags]
            dr = r_int - r_lig
            r = np.linalg.norm(dr, axis=1)
            F_mag = self.p.k_int_bare * np.clip(r - self.p.integrin_r0, 0.0, None)
            k_off = pereverzev_k_off(F_mag, self.p.pereverzev)
            p_break = 1.0 - np.exp(-k_off * self._batch_dt)
            u = self._rng.uniform(0.0, 1.0, size=int_bonds.shape[0])
            broke = u < p_break
            n_broke = int(broke.sum())
            self._n_break_total += n_broke

            if n_broke > 0:
                # int_tags are GLOBAL integrin tags (bond col 0). Map to LOCAL
                # state indices before clearing _engaged (S5 tag-space
                # unification; identity map when integrins occupy [0, n_int)).
                broke_local = np.fromiter(
                    (self._integrin_tag_to_local[int(t)]
                     for t in int_tags[broke]),
                    dtype=np.int64, count=int(n_broke),
                )
                self._engaged[broke_local] = False
                # Keep only the survivors.
                int_bonds = int_bonds[~broke]

        # ---- Step 2: bind unbound integrins within capture radius ----
        # _engaged is LOCAL-indexed (S5), so flatnonzero gives LOCAL indices.
        # Positions are read by GLOBAL tag (snapshot is tag-ordered);
        # _ligand_for_integrin is LOCAL-indexed and stores GLOBAL ligand tags.
        unbound_local = np.flatnonzero(~self._engaged)
        if unbound_local.size > 0:
            unbound_global = self._integrin_global_tag[unbound_local]
            r_int_u = pos[unbound_global]
            lig_tags_u = self._ligand_for_integrin[unbound_local]
            r_lig_u = pos[lig_tags_u]
            d = np.linalg.norm(r_int_u - r_lig_u, axis=1)
            in_range = d <= self.p.capture_radius_R_FA
            n_candidates = int(in_range.sum())
            if n_candidates > 0:
                p_bind = 1.0 - np.exp(-self.p.k_on * self._batch_dt)
                u2 = self._rng.uniform(0.0, 1.0, size=n_candidates)
                bind_mask = u2 < p_bind
                n_bound = int(bind_mask.sum())
                self._n_bind_total += n_bound
                if n_bound > 0:
                    bind_local = unbound_local[in_range][bind_mask]
                    bind_int_tags = self._integrin_global_tag[bind_local]
                    bind_lig_tags = lig_tags_u[in_range][bind_mask]
                    self._engaged[bind_local] = True
                    # Bond col 0 = GLOBAL integrin tag, col 1 = GLOBAL ligand
                    # tag (so the break loop's tag->local lookup round-trips).
                    new_bonds = np.stack(
                        [bind_int_tags, bind_lig_tags], axis=1
                    ).astype(np.int64)
                    int_bonds = np.concatenate([int_bonds, new_bonds], axis=0)

        # ---- Step 3: build fresh Snapshot with updated bond topology ----
        new_int_typeid = np.full(int_bonds.shape[0], integrin_type, dtype=np.uint32)
        if other_bg.shape[0] > 0:
            new_bg = np.concatenate([int_bonds, other_bg], axis=0).astype(np.uint32)
            new_bt = np.concatenate([new_int_typeid, other_bt]).astype(np.uint32)
        else:
            new_bg = int_bonds.astype(np.uint32)
            new_bt = new_int_typeid

        write_snap = hoomd.Snapshot()
        # Copy particle data verbatim from read_snap.
        N_part = int(read_snap.particles.N)
        write_snap.particles.N = N_part
        write_snap.particles.types = list(read_snap.particles.types)
        write_snap.particles.typeid[:] = np.asarray(read_snap.particles.typeid)
        write_snap.particles.position[:] = pos
        write_snap.particles.velocity[:] = np.asarray(read_snap.particles.velocity)
        write_snap.particles.mass[:] = np.asarray(read_snap.particles.mass)
        write_snap.particles.image[:] = np.asarray(read_snap.particles.image)
        # Box copied as-is.
        box = read_snap.configuration.box
        write_snap.configuration.box = list(box)

        # Bonds block: updated.
        write_snap.bonds.N = int(new_bg.shape[0])
        write_snap.bonds.types = list(bond_type_names)
        if new_bg.shape[0] > 0:
            write_snap.bonds.group[:] = new_bg
            write_snap.bonds.typeid[:] = new_bt

        # Pass-through angles / dihedrals / impropers / pairs if any
        # (week-1 H.4 has none, but the H.5 cortex hookup will introduce
        # cortex angles; pre-empt by carrying them across).
        for grp_name in ("angles", "dihedrals", "impropers"):
            src = getattr(read_snap, grp_name)
            dst = getattr(write_snap, grp_name)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)

        sim.state.set_snapshot(write_snap)

        self._steps_run += 1

    # ------------------------------------------------------------------
    # Introspection (used by tests + REPORT)
    # ------------------------------------------------------------------
    @property
    def n_engaged(self) -> int:
        return int(self._engaged.sum())

    @property
    def n_engaged_per_fa(self) -> np.ndarray:
        """Shape (n_FAs,) — engaged count per FA. Used by KU-2.17 test."""
        out = np.zeros(self._n_FAs, dtype=np.int64)
        for fa_idx in range(self._n_FAs):
            mask = self._fa_for_integrin == fa_idx
            out[fa_idx] = int(self._engaged[mask].sum())
        return out

    @property
    def n_break_total(self) -> int:
        return self._n_break_total

    @property
    def n_bind_total(self) -> int:
        return self._n_bind_total

    @property
    def steps_run(self) -> int:
        return self._steps_run

    @property
    def batch_dt(self) -> float:
        return self._batch_dt


def make_integrin_updater(
    *,
    sim: hoomd.Simulation,
    p: "ResolvedH4",
    layouts: "list[FALayout]",
    bond: md.bond.Harmonic,
    seed_offset: int = 1,
) -> tuple[IntegrinBondUpdater, hoomd.update.CustomUpdater]:
    """Build the IntegrinBondUpdater + wrap in a Periodic CustomUpdater.

    Fires every ``p.integrin_batch_steps`` steps.
    """
    action = IntegrinBondUpdater(
        p=p, layouts=layouts, bond=bond, seed_offset=seed_offset,
    )
    updater = hoomd.update.CustomUpdater(
        action=action,
        trigger=hoomd.trigger.Periodic(p.integrin_batch_steps),
    )
    return action, updater
