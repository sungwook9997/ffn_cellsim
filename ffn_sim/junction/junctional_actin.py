"""Junctional actin coupling — the α-catenin/vinculin belt linking the
cadherin junction to each cell's actin cortex (STUB; KU-4.x junction track).

One-line purpose
----------------
Mechanically couple every E-cadherin trans-junction particle (``cadherin``,
built by ``ffn_sim.junction.cadherin``) to the nearby actin cortex of the
SAME cell through an explicit α-catenin/vinculin clutch — a *catch* coupling:
α-catenin unfolds under tension and recruits vinculin, reinforcing the
cadherin→cortex link (Yonemura 2010; le Duc 2010; Buckley 2014; Yao 2014).

This is the third leg of the cell–cell adhesion mechanism: cadherin trans-
dimers carry the cell-to-cell load (``ffn_sim.junction.cadherin``), and the
junctional-actin belt carries that load from the cadherin tail into each
cell's cortical actin network. Without this leg the cadherin junction floats
mechanically uncoupled from the cytoskeleton — the cortex never feels the
junction and the junction never feels cortical tension.

Explicit mechanism (full-fidelity per CLAUDE.md architectural rule)
-------------------------------------------------------------------
NO lumped line tension, NO mesh-as-physics, NO proxy spring standing in for
the catenin complex. Each coupling is an EXPLICIT dynamic harmonic bond:

* particle types added: ``junc_actin`` — an α-catenin/vinculin *coupling
  head*, one per cadherin particle that engages the cortex. (The cadherin
  particles themselves are owned by ``ffn_sim.junction.cadherin``; this
  module adds only the coupling-head cloud and its bonds.)
* bonds added: a dynamic catch clutch
  ``junc_actin_couple_b{i}`` (per-r0-binned harmonic, exactly the H.3
  ``xlink_attach_b{i}`` binning pattern) joining a ``junc_actin`` head to one
  nearby ``actin_cortex`` bead of the SAME cell, plus a permanent harmonic
  ``junc_actin_anchor`` tethering each ``junc_actin`` head to its parent
  ``cadherin`` particle (the α-catenin N-terminus is bound to β-catenin on
  the cadherin tail — a constitutive, force-independent link).
* updater: a per-BATCH ``hoomd.custom.Action`` (mirroring the H.3
  ``XlinkBondUpdater`` and the H.4 ``IntegrinBondUpdater``) that, each batch
  tick, (a) breaks engaged couple bonds with a force-dependent CATCH off-rate
  (α-catenin unfolds and recruits vinculin → the bond is STABILISED by force
  up to a peak, then slips), and (b) binds free ``junc_actin`` heads to a
  nearby same-cell cortex bead found by a cKDTree broad-phase query.

The catch behaviour is the documented α-catenin/F-actin two-state catch bond
(Buckley 2014 Science): under load the α-catenin–actin interface switches to
a strongly-bound state, so the off-rate FALLS with force up to a peak force
F* and rises beyond it (a two-pathway catch-slip, same functional family as
the H.3 filamin / H.4 integrin Pereverzev law). Vinculin recruitment
(Yonemura 2010; le Duc 2010) is the molecular cause: α-catenin's M-domain
unfurls under tension, exposing the vinculin-binding site, and bound vinculin
adds a second actin-anchoring arm that deepens the catch.

STUB status — why the enabled path is NotImplementedError
---------------------------------------------------------
The TOPOLOGY (the explicit particles + bonds + updater scaffold), the
resolver, the Sanity Gate, and the Performance Contract are all real and
authored here. But the QUANTITATIVE catch-bond constants for the
cadherin→cortex coupling clutch — the coupling stiffness ``k_couple``, the
zero-force catch/slip rates, and the catch/slip force scales — are NOT
anchored to a quantitative single-molecule measurement at this ×40
mesoscale. Buckley 2014 measured the α-catenin/F-actin catch bond
qualitatively (lifetime-vs-force shape, two-state model) but the published
constants are for the isolated α-catenin–actin interface under a specific
assay geometry, not a calibrated cadherin→cortex *coupling* stiffness for a
coarse-grained junction bead; le Duc 2010 / Yonemura 2010 establish the
force-dependent vinculin recruitment qualitatively. Per CLAUDE.md
no-magic-number, the coupling constants are therefore set to ``None`` with
``PI_DECISIONS`` entries, and any attempt to BUILD the enabled coupling
raises ``NotImplementedError`` until a PI sign-off anchors them. The
disabled path is a clean no-op; the enabled path is intentionally inert
until anchored. Honesty over completeness.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_junctional_actin.py``.*

1. **Dimensional analysis**
   - ``k_couple`` [N/m]; coupling-bond extension ``(|Δr| − r0)`` [m] →
     ``F = k_couple · Δr`` [N]. ✓ (constant is ``None`` in this stub; the
     dimensional contract is stated so the future anchored value plugs in.)
   - catch off-rate ``k_off(F) = k_catch0·exp(−F·x_catch/kT)
     + k_slip0·exp(+F·x_slip/kT)``: ``x_catch``/``x_slip`` [m], ``F`` [N],
     ``kT`` [J] → exponent dimensionless → ``k_off`` [1/s]. ✓
   - batch break/bind probability ``1 − exp(−k_off·Δt_batch)`` dimensionless,
     ``Δt_batch = batch_steps · dt`` [s]. ✓
   - All positions metres, forces Newtons, SI throughout.

2. **Boundary cases**
   - ``enabled=False`` (or a missing/false config): the resolver returns a
     zeroed ``ResolvedJunctionalActin(enabled=False, …)`` and the topology
     builder early-returns its input UNCHANGED (zero particles, zero bonds).
     The OFF path is a bit-for-bit no-op (RUNTIME + STATIC test).
   - ``enabled=True`` with the coupling constants still ``None``: the
     resolver succeeds (topology may be planned) but the BUILDER raises
     ``NotImplementedError`` — the un-anchored physics must not run.
   - No cadherin particles at the interface (n_cad = 0): nothing to couple;
     the builder adds zero heads even when enabled+anchored (trivially
     correct).
   - ``max_couple_dist ≤ 0`` / non-positive geometry: resolver ``ValueError``.

3. **Conservation invariants**
   - **Newton's 3rd law**: both the ``junc_actin_anchor`` (head↔cadherin) and
     the dynamic ``junc_actin_couple_b{i}`` (head↔cortex) bonds are
     ``md.bond.Harmonic`` symmetric pairs — they inject ZERO net momentum.
     The coupling is an INTERNAL load path between two existing particle
     populations, not an external field (unlike ``cortex/erm.py``). A
     break/bind event mutates topology but never injects momentum.
   - Particle count after a (future, anchored) extend is
     ``n_in + n_coupling_heads`` (cadherin + cortex populations unchanged;
     only ``junc_actin`` heads appended). The disabled builder adds zero.
   - The coupling bonds form a NON-CORTICAL load path → every bond type is
     prefixed ``junc_actin`` so the cortical-tension estimator denylists them
     (no γ contamination — Hard Rule 5; STATIC test).

4. **Numerical sanity (CFL)**
   - When anchored, the coupling-bond CFL is ``dt ≤ cfl_safety_factor ·
     γ_cortex / k_couple`` (the same bond-relaxation gate as
     ``cortex/crosslinkers`` / ``cortex/erm``), and the batch off-rate CFL is
     ``batch_steps · dt · k_off_max ≤ 1e-3`` (the shared D2 batched-updater
     contract). Both are stated here; the gates are dormant while
     ``k_couple`` / the rates are ``None`` (the build raises first).

5. **Sign / sense**
   - Coupling clutch is a TENSILE-only harmonic
     ``F = k_couple · max(0, |Δr| − r0)``: a STRETCHED head↔cortex bond pulls
     the cortex bead toward the junction and the junction head toward the
     cortex (mutual attraction; no compression push). Sign of the pure scalar
     law is tested analytically (``coupling_force_magnitude`` ≥ 0, and the
     vector pulls the two endpoints together).
   - CATCH signature: ``k_off(F)`` FALLS with force up to F* then RISES — the
     opposite sign-of-slope at low force to a slip bond (force STABILISES the
     α-catenin/actin link; Buckley 2014). The scalar ``catch_off_rate`` law
     is exposed so the catch-slip signature can be asserted once anchored;
     while the rates are ``None`` it raises ``NotImplementedError``.

6. **Measurement-protocol consistency**
   - At zero junction load (F ≈ 0) the engaged-coupling fraction approaches
     ``k_on/(k_on + k_catch0 + k_slip0)`` (two-state kinetics) — the
     production gate to assert once the rates are anchored. Stated, not yet
     asserted (stub).
   - The coupling transmits cadherin tension into cortical tension; because
     every coupling bond is on the ``junc_actin`` denylist it must NOT appear
     in the cortical-tension γ accounting (it is a junction load path, not a
     cortex-internal one). STATIC denylist-prefix test.

Compartment Performance Contract
--------------------------------
* **Particle types added**: ``junc_actin`` (α-catenin/vinculin coupling head,
  one per engaged cadherin particle at the interface).
* **Particle count** — one coupling head per cadherin particle that engages
  the cortex. At the mesoscale (n_fil = 1000) a cell-cell interface carries
  ``~n_cad`` cadherin particles (KU-4.17 ``N_cad_per_contact ≈ 100`` per
  mature contact); so ``~100`` ``junc_actin`` heads per interface (≪ the
  ~1000-filament cortex). At native ~38000-filament scale the cadherin count
  scales with interface area (×40 areal coarse-graining inverse) →
  ``~100·√40 ≈ 630`` heads per interface — still a thin belt, NOT a per-bead
  cost over the whole cortex.
* **Bond / angle count**: per engaged head, 1 permanent ``junc_actin_anchor``
  (head↔cadherin) + ≤ 1 dynamic ``junc_actin_couple_b{i}`` (head↔cortex);
  so ≤ ``2·n_head`` bonds (``~200`` per interface at mesoscale). 0 angles.
* **Per-step force**: NO. The coupling is carried by ``md.bond.Harmonic``
  (HOOMD's built-in bonded force), not a per-step ``md.force.Custom``.
* **Per-batch updater**: YES — a ``hoomd.custom.Action`` runs every
  ``batch_steps`` integration steps to break (catch off-rate) / bind
  coupling bonds (the H.3 ``XlinkBondUpdater`` / H.4 ``IntegrinBondUpdater``
  cadence; gated by the shared D2 batch CFL).
* **Uses cpu_local_snapshot**: YES (the batch updater reads positions via
  ``sim.state.get_snapshot()`` / ``cpu_local_snapshot`` like the H.3/H.4
  binders — a GPU→host sync per batch tick, not per step).
* **Uses cKDTree / broad-phase**: YES — a same-cell cortex-acceptor query
  (``scipy.spatial.cKDTree``) per batch tick, identical to the crosslinker
  binder. Restricted to the small junction belt, not the whole cortex.
* **Hot-path priority**: **P2** (per the compartment registry). The belt is
  thin and the updater is batched — meaningful but not the dominant cost.
* **GPU path now**: **CPU** (the batch updater is host-resident, sharing the
  GPU-MAIN port debt with the other binders). ``native ForceCompute
  candidate``: NO — the harmonic coupling already runs as a builtin HOOMD
  bonded force; only the binder needs the eventual cupy/native port.
* **Bottleneck risk**: LOW — small belt, batched updater, reuses the proven
  clutch pattern; the dominant junction cost is the cadherin cross-cell
  partner search, not this same-cell coupling.

References
----------
- Yonemura, Wada, Miyake, Suzuki, Kinoshita & Shibata (2010)
  "α-Catenin as a tension transducer that induces adherens junction
  development", *Nat. Cell Biol.* 12:533–542. doi:10.1038/ncb2055. —
  force-dependent vinculin recruitment / α-catenin tension sensing at
  adherens junctions. *In corpus (Yonemura2010_NCB).*
- le Duc, Shi, Blonk, Sonnenberg, Wang, Leckband & de Rooij (2010)
  "Vinculin potentiates E-cadherin mechanosensing and is recruited to
  actin-anchored sites within adherens junctions in a myosin II–dependent
  manner", *J. Cell Biol.* 189:1107–1115. doi:10.1083/jcb.201001149. —
  the vinculin force-recruitment step. *In SE candidate set (leDuc2010_JCB).*
- Buckley, Tan, Palmer, Doyle, Dunn, Sawyer, Manibog, Pielak, Astrof,
  Cheung, Liu, Fierke, Schepartz, Selvin, Frydman & Dunn (2014)
  "The minimal cadherin-catenin complex binds to actin filaments under
  force", *Science* 346:1254211. doi:10.1126/science.1254211. — THE
  α-catenin/F-actin CATCH bond (two-state, force-strengthened); the catch
  mechanism cited for the coupling clutch.
- Yao, Qiu, Chen, Zhang, Wen, Mi, Xie, Ma, Liu, Lei, Chen & Le (2014)
  "Force-dependent conformational switch of α-catenin controls vinculin
  binding", *Nat. Commun.* 5:4525. doi:10.1038/ncomms5525. — the
  force-dependent α-catenin conformational switch exposing the vinculin
  site (the molecular basis of the catch reinforcement).
- ``ffn_sim/cortex/crosslinkers.py`` — the dynamic per-r0-binned attach-bond
  family + batched Bell-Evans / Pereverzev catch-slip Updater pattern this
  coupling clutch mirrors.
- ``ffn_sim/cell/cell.py`` ``SubstrateLigandPin`` — ``hoomd.custom.Action``
  updater pattern.
- ``ffn_sim/junction/cadherin.py`` — the required ``cadherin_junction``
  compartment that owns the ``cadherin`` particles this module couples.
- Compartment contract: ``ffn_sim/cell/compartment_registry.py`` →
  ``junctional_actin`` (STUB, requires ``cadherin_junction``, denylist
  ``junc_actin``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Module-level contracts (read by the cortical-tension estimator + the Lead)
# ---------------------------------------------------------------------------
#: Bond-type prefix for every bond this module creates. The cortical-tension
#: estimator denylists this prefix so the junction load path NEVER contaminates
#: the emergent cortical-tension γ (Hard Rule 5). Mirrors the compartment
#: registry's ``denylist_bond_types=("junc_actin",)`` for this compartment.
GAMMA_DENYLIST_PREFIX: str = "junc_actin"

#: Open PI decisions — every entry is a constant we refuse to invent.
#: The enabled BUILD path raises NotImplementedError while any of these stand.
PI_DECISIONS: list[str] = [
    "k_couple (cadherin->cortex coupling-clutch stiffness, N/m) is not anchored "
    "to a single-molecule measurement at the x40 mesoscale. Buckley 2014 gives "
    "the alpha-catenin/F-actin catch-bond SHAPE, not a calibrated coupling "
    "stiffness for a coarse-grained junction bead. PI must anchor it.",
    "catch off-rate constants (k_catch0 [1/s], x_catch [m], k_slip0 [1/s], "
    "x_slip [m]) for the alpha-catenin/F-actin catch bond are qualitative in "
    "Buckley 2014 / Yao 2014 (two-state, force-strengthened) but not given as "
    "a calibrated two-pathway Pereverzev parameter set for THIS coupling. PI "
    "must anchor them (or sanction a transfer from the H.3 filamin catch set).",
    "k_on (coupling-head -> cortex single-head binding rate, 1/s) and "
    "max_couple_dist (acceptor search radius, m) for the junctional belt are "
    "not anchored. The vinculin-recruitment ON step (Yonemura 2010 / le Duc "
    "2010) is qualitative. PI must anchor or sanction a transfer.",
]


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class ResolvedJunctionalActin:
    """Resolved junctional-actin coupling parameters (all SI).

    STUB: the quantitative coupling/catch constants are ``None`` until a PI
    sign-off anchors them (see :data:`PI_DECISIONS`). The geometry /
    bookkeeping fields are real; the un-anchored physics fields are ``None``
    and any build that needs them raises ``NotImplementedError``.

    Attributes:
        enabled: Master on/off. When ``False`` the builder is a no-op.
        n_couplings_max: Upper bound on coupling heads to seed (one per
            engaged cadherin particle at the interface). 0 when disabled.
        max_couple_dist: Cortex-acceptor search radius [m] for the coupling
            clutch (head -> same-cell cortex bead). ``None`` until anchored.
        anchor_r0: Rest length [m] of the permanent ``junc_actin_anchor``
            bond (junc_actin head <-> parent cadherin); the α-catenin
            N-terminus to β-catenin/cadherin-tail constitutive link.
        k_anchor: Stiffness [N/m] of that permanent anchor bond. ``None``
            until anchored (constitutive, force-independent, but its
            stiffness still needs a literature value).
        k_couple: Dynamic coupling-clutch stiffness [N/m]
            (head <-> cortex). ``None`` (Buckley 2014 gives shape not value).
        k_catch0: Catch-pathway zero-force off-rate [1/s]. ``None``.
        x_catch: Catch-pathway force scale [m]. ``None``.
        k_slip0: Slip-pathway zero-force off-rate [1/s]. ``None``.
        x_slip: Slip-pathway force scale [m]. ``None``.
        k_on: Single-head binding rate when an acceptor is present [1/s].
            ``None``.
        n_bins: Number of per-r0 bins for force-free dynamic-bond placement
            (mirrors the H.3 ``xlink_attach_b{i}`` binning).
        batch_steps: BAOAB steps between updater ticks (D2 cadence).
        dt: Host sim timestep [s] (for the batch CFL).
        kT: Thermal energy [J] (for the catch off-rate exponent).
        seed: Deterministic RNG seed for the binder.
        extras: Free-form diagnostics bag.
    """

    enabled: bool

    # Geometry / bookkeeping (real even in the stub).
    n_couplings_max: int
    n_bins: int
    batch_steps: int
    dt: float
    kT: float
    seed: int

    # Un-anchored physics — ``None`` until PI sign-off (see PI_DECISIONS).
    max_couple_dist: float | None = None     # m   cortex-acceptor radius
    anchor_r0: float | None = None           # m   head<->cadherin rest length
    k_anchor: float | None = None            # N/m head<->cadherin stiffness
    k_couple: float | None = None            # N/m head<->cortex clutch stiffness
    k_catch0: float | None = None            # 1/s catch-pathway zero-force rate
    x_catch: float | None = None             # m   catch-pathway force scale
    k_slip0: float | None = None             # 1/s slip-pathway zero-force rate
    x_slip: float | None = None              # m   slip-pathway force scale
    k_on: float | None = None                # 1/s single-head binding rate

    extras: dict[str, Any] = field(default_factory=dict)

    @property
    def is_anchored(self) -> bool:
        """True iff every un-anchored physics constant has been supplied.

        While this is ``False`` the enabled BUILD path raises
        ``NotImplementedError`` (the un-anchored physics must not run).
        """
        return all(
            v is not None
            for v in (
                self.max_couple_dist,
                self.anchor_r0,
                self.k_anchor,
                self.k_couple,
                self.k_catch0,
                self.x_catch,
                self.k_slip0,
                self.x_slip,
                self.k_on,
            )
        )

    @property
    def batch_dt(self) -> float:
        """Batch interval ``batch_steps · dt`` [s]."""
        return self.batch_steps * self.dt


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def _disabled(
    *, dt: float, kT: float, n_bins: int, batch_steps: int, seed: int
) -> ResolvedJunctionalActin:
    """Zeroed OFF-identity result (every physics field ``None``)."""
    return ResolvedJunctionalActin(
        enabled=False,
        n_couplings_max=0,
        n_bins=int(n_bins),
        batch_steps=int(batch_steps),
        dt=float(dt),
        kT=float(kT),
        seed=int(seed),
    )


def resolve_junctional_actin(
    cfg: dict,
    *,
    kT: float,
    dt: float,
    n_cadherin: int | None = None,
) -> ResolvedJunctionalActin:
    """Resolve the ``junctional_actin`` config block (DEFAULT-OFF; STUB).

    ``cfg`` may be the YAML root, the ``junction`` sub-dict, or the
    ``junctional_actin`` sub-dict. The compartment is DEFAULT-OFF: a missing
    block, or ``enabled: false``, returns a zeroed
    ``ResolvedJunctionalActin(enabled=False, …)`` whose builder is a no-op.

    When ``enabled: true`` the resolver still succeeds (so the topology can be
    planned and the contract surfaced), but it does NOT invent the coupling /
    catch constants — those stay ``None`` (see :data:`PI_DECISIONS`). The
    BUILD path then raises ``NotImplementedError`` until a PI sign-off anchors
    them (``ResolvedJunctionalActin.is_anchored``). The resolver only validates
    quantities the host genuinely supplies (kT, dt, bin/batch counts, and any
    anchored constants passed through ``cfg``).

    Args:
        cfg: Config mapping (root, ``junction``, or ``junctional_actin``).
        kT: Thermal energy [J] (> 0; for the catch off-rate exponent).
        dt: Host simulation timestep [s] (> 0; for the batch CFL).
        n_cadherin: Number of cadherin particles at the interface (one
            coupling head per cadherin). ``None`` or 0 → no heads.

    Returns:
        Resolved junctional-actin parameters (SI). ``.enabled`` reflects the
        config; physics constants are ``None`` unless explicitly anchored via
        ``cfg``.

    Raises:
        ValueError: on non-positive ``kT`` / ``dt`` / bin / batch counts, a
            negative ``n_cadherin``, or an explicitly-supplied non-positive
            physics constant.
    """
    if "junction" in cfg:
        cfg = cfg["junction"]
    if "junctional_actin" in cfg:
        cfg = cfg["junctional_actin"]

    n_bins = int(cfg.get("n_bins", 10))
    batch_steps = int(cfg.get("batch_steps", 100))
    seed = int(cfg.get("seed", 42))

    # Host-supplied scalars are validated whether enabled or not.
    _require_finite_positive("kT", kT)
    _require_finite_positive("dt", dt)
    if n_bins < 1:
        raise ValueError(f"n_bins must be ≥ 1; got {n_bins}")
    if batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {batch_steps}")
    if n_cadherin is not None and int(n_cadherin) < 0:
        raise ValueError(f"n_cadherin must be ≥ 0; got {n_cadherin}")

    enabled = bool(cfg.get("enabled", False))
    if not enabled:
        return _disabled(
            dt=dt, kT=kT, n_bins=n_bins, batch_steps=batch_steps, seed=seed
        )

    # Enabled: plan the topology bound, but do NOT invent un-anchored physics.
    # Any anchored constant a PI later supplies in ``cfg`` is range-checked.
    n_couplings_max = int(n_cadherin) if n_cadherin else 0

    def _opt_positive(key: str) -> float | None:
        if key not in cfg or cfg[key] is None:
            return None
        v = float(cfg[key])
        _require_finite_positive(key, v)
        return v

    p = ResolvedJunctionalActin(
        enabled=True,
        n_couplings_max=n_couplings_max,
        n_bins=n_bins,
        batch_steps=batch_steps,
        dt=float(dt),
        kT=float(kT),
        seed=seed,
        max_couple_dist=_opt_positive("max_couple_dist"),
        anchor_r0=_opt_positive("anchor_r0"),
        k_anchor=_opt_positive("k_anchor"),
        k_couple=_opt_positive("k_couple"),
        k_catch0=_opt_positive("k_catch0"),
        x_catch=_opt_positive("x_catch"),
        k_slip0=_opt_positive("k_slip0"),
        x_slip=_opt_positive("x_slip"),
        k_on=_opt_positive("k_on"),
    )
    return p


# ---------------------------------------------------------------------------
# Per-r0 binning helpers (force-free dynamic attachment; mirrors H.3 cortex)
# ---------------------------------------------------------------------------
def junc_actin_couple_bin_names(n_bins: int) -> list[str]:
    """HOOMD bond-type names for the per-r0-binned coupling clutch.

    Every name carries the :data:`GAMMA_DENYLIST_PREFIX` so the cortical-
    tension estimator denylists the whole junction load path.

    Args:
        n_bins: Number of rest-length bins (≥ 1).

    Returns:
        ``["junc_actin_couple_b0", …, "junc_actin_couple_b{n-1}"]``.
    """
    if int(n_bins) < 1:
        raise ValueError(f"n_bins must be ≥ 1; got {n_bins}")
    return [f"{GAMMA_DENYLIST_PREFIX}_couple_b{i}" for i in range(int(n_bins))]


def junc_actin_bond_type_names(n_bins: int) -> list[str]:
    """All bond-type names this module would register (anchor + clutch bins).

    The permanent head<->cadherin anchor plus the per-r0 coupling-clutch bins.
    Used by the OFF-identity / denylist Sanity-Gate tests and by the Lead's
    topology registration. Every name starts with
    :data:`GAMMA_DENYLIST_PREFIX`.

    Args:
        n_bins: Number of coupling-clutch rest-length bins (≥ 1).

    Returns:
        ``["junc_actin_anchor", "junc_actin_couple_b0", …]``.
    """
    return [f"{GAMMA_DENYLIST_PREFIX}_anchor", *junc_actin_couple_bin_names(n_bins)]


def junc_actin_couple_bin_rest_lengths(
    n_bins: int, max_couple_dist: float
) -> np.ndarray:
    """Bin-center rest lengths uniformly across ``(0, max_couple_dist]`` [m].

    Mirrors ``cortex.crosslinkers.xlink_attach_bin_rest_lengths`` so a coupling
    bond placed at distance ``d`` gets a per-bin ``r0 ≈ d`` (force-free at
    construction). The first bin center is > 0 (a length-0 r0 would be a
    zero-force degeneracy).

    Args:
        n_bins: Number of bins (≥ 1).
        max_couple_dist: Acceptor search radius [m] (> 0).

    Returns:
        ``(n_bins,)`` float64 array of bin-center rest lengths [m].
    """
    if int(n_bins) < 1:
        raise ValueError(f"n_bins must be ≥ 1; got {n_bins}")
    _require_finite_positive("max_couple_dist", max_couple_dist)
    edges = np.linspace(0.0, float(max_couple_dist), int(n_bins) + 1, dtype=np.float64)
    return 0.5 * (edges[:-1] + edges[1:])


# ---------------------------------------------------------------------------
# Pure scalar force / rate laws (testable analytically; raise while un-anchored)
# ---------------------------------------------------------------------------
def coupling_force_magnitude(p: ResolvedJunctionalActin, extension: float) -> float:
    """Tensile coupling-clutch force magnitude ``F = k_couple·max(0, Δr)`` [N].

    The dynamic head<->cortex coupling is a TENSILE-only harmonic: a stretched
    bond pulls the endpoints together; a compressed bond carries no force
    (``Δr = |r| − r0 < 0`` → 0). Sign-sense: ``F ≥ 0`` always, pulling the
    cortex bead and the junction head toward each other.

    Args:
        p: Resolved parameters.
        extension: Signed bond extension ``|r| − r0`` [m].

    Returns:
        Tensile force magnitude [N] (≥ 0).

    Raises:
        NotImplementedError: if ``k_couple`` is un-anchored (``None``) — the
            un-anchored physics must not run (see :data:`PI_DECISIONS`).
    """
    if p.k_couple is None:
        raise NotImplementedError(
            "junctional_actin coupling stiffness k_couple is not anchored "
            "(STUB). See PI_DECISIONS; surface to PI before enabling."
        )
    return float(p.k_couple) * max(0.0, float(extension))


def catch_off_rate(p: ResolvedJunctionalActin, F: float) -> float:
    """α-catenin/F-actin two-pathway CATCH-slip off-rate ``k_off(F)`` [1/s].

    ``k_off(F) = k_catch0·exp(−F·x_catch/kT) + k_slip0·exp(+F·x_slip/kT)``.

    The catch pathway DECREASES the off-rate with force (force STABILISES the
    α-catenin/actin link — Buckley 2014; the vinculin-reinforced strongly-
    bound state), the slip pathway increases it. With ``k_catch0 > k_slip0``
    and ``x_catch > x_slip`` the catch branch dominates at low force, giving
    the catch-slip signature (off-rate falls to a minimum at the peak force
    ``F* = (kT/(x_catch+x_slip))·ln((k_catch0·x_catch)/(k_slip0·x_slip))`` then
    rises). Same two-pathway functional structure as the H.3 filamin / H.4
    integrin Pereverzev law.

    Args:
        p: Resolved parameters.
        F: Bond tensile force [N] (≥ 0).

    Returns:
        Off-rate [1/s].

    Raises:
        NotImplementedError: if any catch/slip constant is un-anchored
            (``None``) — the un-anchored physics must not run.
    """
    if None in (p.k_catch0, p.x_catch, p.k_slip0, p.x_slip):
        raise NotImplementedError(
            "junctional_actin catch/slip constants are not anchored (STUB). "
            "Buckley 2014 / Yao 2014 give the catch-bond SHAPE, not a "
            "calibrated parameter set. See PI_DECISIONS; surface to PI."
        )
    kT = p.kT
    return (
        float(p.k_catch0) * math.exp(-float(F) * float(p.x_catch) / kT)
        + float(p.k_slip0) * math.exp(+float(F) * float(p.x_slip) / kT)
    )


# ---------------------------------------------------------------------------
# Topology extension builder (DEFAULT-OFF no-op; enabled path = NotImplemented)
# ---------------------------------------------------------------------------
def extend_snapshot_with_junctional_actin(
    snapshot: Any,
    p: ResolvedJunctionalActin,
    *,
    cadherin_tags: np.ndarray | None = None,
    cortex_positions: np.ndarray | None = None,
    cortex_tag_start: int | None = None,
) -> Any:
    """Append the ``junc_actin`` coupling belt to a HOOMD/gsd snapshot.

    DEFAULT-OFF contract (Hard Rule 2): when ``p.enabled`` is ``False`` this
    returns the snapshot UNCHANGED — zero ``junc_actin`` particles, zero
    bonds, no type registration — a bit-for-bit no-op the Lead can call
    unconditionally.

    When ``p.enabled`` is ``True`` but the coupling/catch constants are NOT
    yet anchored (``not p.is_anchored``), this raises ``NotImplementedError``
    — the explicit topology is designed (heads at the cadherin tails, dynamic
    catch clutch to same-cell cortex, per-r0 binning) but the un-anchored
    physics must not be built (STUB; see :data:`PI_DECISIONS`).

    Args:
        snapshot: A ``hoomd``/``gsd`` snapshot to extend (returned unchanged
            when disabled).
        p: Resolved junctional-actin parameters.
        cadherin_tags: Global tags of the cadherin particles to couple (one
            ``junc_actin`` head per tag). Required when enabled+anchored.
        cortex_positions: Same-cell cortex bead positions ``(n, 3)`` [m] for
            initial force-free head placement. Required when enabled+anchored.
        cortex_tag_start: Global tag offset of the same-cell cortex beads.
            Required when enabled+anchored.

    Returns:
        The snapshot (UNCHANGED when disabled; otherwise it would be the
        extended snapshot — unreachable while the stub is un-anchored).

    Raises:
        NotImplementedError: when enabled but the coupling/catch constants are
            un-anchored (STUB). This is the intended state until PI sign-off.
    """
    # --- DEFAULT-OFF: bit-for-bit no-op (Hard Rule 2) ---
    if not p.enabled:
        return snapshot

    # --- enabled but un-anchored: the explicit topology is designed, but the
    #     physics must not run until PI anchors the constants (STUB). ---
    if not p.is_anchored:
        raise NotImplementedError(
            "junctional_actin is enabled but its cadherin<->cortex coupling / "
            "alpha-catenin catch-bond constants are NOT anchored (STUB). The "
            "explicit topology (junc_actin coupling heads + junc_actin_anchor "
            "+ junc_actin_couple_b{i} dynamic catch clutch) is designed and "
            "the Sanity Gate + Performance Contract are authored, but per "
            "CLAUDE.md no-magic-number the coupling stiffness k_couple and the "
            "catch/slip rates are None pending a PI literature anchor "
            "(Buckley 2014 gives the catch-bond SHAPE, not a calibrated "
            "parameter set for a coarse-grained junction bead). See "
            "junctional_actin.PI_DECISIONS. Do not invent these to make the "
            "build pass; surface to PI."
        )

    # --- enabled AND anchored: the real explicit build would go here. It is
    #     deliberately unreachable until a PI sign-off supplies the constants
    #     (is_anchored == True). The build would:
    #       1. append one `junc_actin` head per cadherin tag at the cadherin
    #          position (anchor_r0 force-free along the cadherin->cortex line);
    #       2. register `junc_actin_anchor` + `junc_actin_couple_b{i}` bond
    #          types (all GAMMA_DENYLIST_PREFIX-prefixed);
    #       3. add the permanent head<->cadherin anchor bonds;
    #       4. seed force-free head<->cortex coupling bonds (per-r0-binned),
    #          for the JunctionalActinCouplingUpdater to maintain.
    raise NotImplementedError(  # pragma: no cover - unreachable in the stub
        "junctional_actin enabled+anchored explicit build path is reserved "
        "for the PI-sanctioned anchored release; not implemented in the STUB."
    )


# ---------------------------------------------------------------------------
# Per-batch coupling Updater (catch off-rate break + bind); STUB scaffold
# ---------------------------------------------------------------------------
class JunctionalActinCouplingUpdater:
    """Per-batch α-catenin/vinculin coupling-clutch updater (STUB scaffold).

    Mirrors ``cortex.crosslinkers.XlinkBondUpdater`` / the H.4
    ``IntegrinBondUpdater``: a ``hoomd.custom.Action`` that, every
    ``batch_steps`` integration steps, breaks engaged
    ``junc_actin_couple_b{i}`` bonds with the force-dependent CATCH off-rate
    (:func:`catch_off_rate`) and binds free ``junc_actin`` heads to a nearby
    SAME-CELL cortex bead found by a ``scipy.spatial.cKDTree`` broad-phase
    query.

    STUB: the constructor accepts a ``ResolvedJunctionalActin`` but refuses to
    operate until the constants are anchored — instantiating it on an
    un-anchored ``p`` raises ``NotImplementedError`` so a caller cannot
    silently run the un-anchored physics. (It does NOT subclass
    ``hoomd.custom.Action`` at import so the resolver-level tests stay
    import-light; the anchored release wires it onto the Action base.)

    Args:
        p: Resolved parameters (must satisfy ``p.is_anchored``).
        cadherin_tags: Global tags of the coupled cadherin particles.
        n_same_cell_cortex: Count of same-cell cortex acceptor beads.

    Raises:
        NotImplementedError: if ``p`` is not anchored (STUB) — the un-anchored
            physics must not run.
        ValueError: on an inconsistent tag set.
    """

    def __init__(
        self,
        *,
        p: ResolvedJunctionalActin,
        cadherin_tags: np.ndarray,
        n_same_cell_cortex: int,
    ) -> None:
        if not p.enabled:
            raise ValueError(
                "JunctionalActinCouplingUpdater requires an enabled "
                "ResolvedJunctionalActin (got enabled=False)."
            )
        if not p.is_anchored:
            raise NotImplementedError(
                "JunctionalActinCouplingUpdater cannot run: the "
                "cadherin<->cortex coupling / alpha-catenin catch-bond "
                "constants are not anchored (STUB). See "
                "junctional_actin.PI_DECISIONS; surface to PI."
            )
        if int(n_same_cell_cortex) < 0:
            raise ValueError(
                f"n_same_cell_cortex must be ≥ 0; got {n_same_cell_cortex}"
            )
        self.p = p
        self.cadherin_tags = np.asarray(cadherin_tags, dtype=np.int64)
        self.n_same_cell_cortex = int(n_same_cell_cortex)
