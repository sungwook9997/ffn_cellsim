"""Explicit E-cadherin trans-dimer CATCH-BOND cell-cell junction (KU-4.2).

One-line purpose
----------------
A fine-grained, fully-explicit E-cadherin trans-dimer adherens junction
between TWO cells: every cadherin is a HOOMD particle on a cell's surface,
and every engaged cell-cell adhesion is an EXPLICIT dynamic ``cadherin_trans``
harmonic bond across the interface whose stochastic lifetime follows the
faithful Rakshit-2012 sliding-rebinding catch-slip off-rate. No lumped line
tension, no slip-only shortcut, no effective-force adiabatic average — the
trans-dimers are resolved individually as stochastic bonds (the architectural
rule's full-fidelity option, KU-4.2).

Explicit mechanism (what particles / bonds)
-------------------------------------------
* **Particles** — one ``cadherin`` particle type. Each of the two cells
  carries ``n_cad_per_cell`` cadherin ectodomain-tip beads, seeded on (or just
  outside) that cell's surface. They are EXPLICIT degrees of freedom: the
  BAOAB heat bath integrates them; they diffuse with the cell surface they are
  anchored to (the anchoring spring to the parent cortex is a SEPARATE concern
  owned by the host cell build — this module supplies the trans-cell adhesion
  load path only).
* **Bonds** — ``cadherin_trans`` (dynamic, harmonic). A bond joins a cadherin
  on cell A to a cadherin on cell B across the interface. It is created /
  destroyed by :class:`CadherinTransJunctionUpdater` per batch tick:
    - *bind*: an unbound cadherin on A that has an unbound cadherin on B within
      ``r_bind`` forms a trans-dimer with probability
      ``1 − exp(−k_on · Δt_batch)`` (one trans-dimer per cadherin — a single
      EC1 strand-swap interface).
    - *break*: an engaged trans-dimer carrying elastic force
      ``F = k_trans·max(0, L − r0_trans)`` ruptures with probability
      ``1 − exp(−k_off(F)·Δt_batch)`` where ``k_off(F) = 1/τ(F)`` is the
      faithful sliding-rebinding effective off-rate
      (``validation.cadherin_sliding_rebinding.effective_k_off``).
  The harmonic ``cadherin_trans`` bond carries the trans-dimer's elastic
  tension (a symmetric Newton-3rd-law pair force pulling the two cells
  together), exactly as the v1 single-cadherin spring (Rakshit elastic
  contact) but resolved per dimer rather than per contact ensemble.

REUSE (cited)
-------------
* ``ffn_sim.validation.cadherin_sliding_rebinding`` — the faithful Rakshit-2012
  sliding-rebinding closed-form catch-slip ``k_off(F)`` / ``mean_lifetime(F)``
  (runtime-importable; the same module ``spheroid/cadherin_bonds.py`` samples).
  This junction samples the SAME ``effective_k_off`` per trans-dimer.
* ``ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds`` — the existing single-cell-aggregate
  catch-bond cohesion (``CadherinBondUpdater``-class adiabatic contact law,
  ``N_cad = F_detach/f0 ≈ 223`` scale bridge). This module is the EXPLICIT
  two-cell-interface counterpart: it replaces the adiabatic per-contact
  effective force with per-dimer stochastic bonds (the timescale gap the
  spheroid module averaged out is resolved here because a single-cell
  junction has ~hundreds of dimers, not ~10⁴ growing cells).
* ``ffn_sim.archive.hoomd_legacy.cortex.crosslinkers`` — the D2 dynamic-bond ``Action`` rebuild
  idiom (cKDTree partner search + per-batch break/bind + ``bonds.group``
  rewrite) is mirrored here for the trans-cell partner search.
* ``ffn_sim.archive.hoomd_legacy.cell.cell.SubstrateLigandPin`` — the ``hoomd.custom.Action`` +
  ``attach``/``act`` + ``_sim_ref`` lifecycle.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_cadherin_junction.py``.*

1. **Dimensional analysis**
   - ``k_trans`` [N/m]; extension ``(L − r0_trans)`` [m]; bond force
     ``F = k_trans·Δ`` [N]. ✓
   - ``F → k_off(F) = effective_k_off(F)`` [s⁻¹] (Rakshit closed form,
     argument in N). ✓
   - ``Δt_batch = batch_steps · dt`` [s]; ``p = 1 − exp(−k·Δt_batch)`` [–]. ✓
   - ``k_on`` [s⁻¹]; ``r_bind`` [m]; ``r0_trans`` [m]. ✓
   - Bell-Evans/Rakshit batch CFL: ``batch_steps·dt·k_off_max ≤ 1e−3``
     enforced at resolve time (mirrors the D2 xlink CFL).

2. **Boundary cases**
   - ``enabled=False`` (or cfg missing / ``enabled: false``) → resolver
     returns ``ResolvedCadherinJunction(enabled=False, …)`` with zeros and the
     builder is a strict no-op (snapshot/sim unchanged, 0 particles, 0 bonds).
   - Negative / zero ``k_trans``, ``r0_trans``, ``r_bind``, ``k_on``,
     ``n_cad_per_cell``, ``batch_steps`` → ``ValueError`` at resolve time.
   - No B-cadherin within ``r_bind`` of any A-cadherin → 0 binds (the
     updater forms no bonds; force-free).
   - Single cell present (no partner cell tags) → trivially 0 trans bonds.

3. **Conservation / no-net-force**
   - ``cadherin_trans`` is an ``md.bond.Harmonic`` (symmetric pair): each
     engaged dimer applies equal-and-opposite forces to its two cadherins →
     the junction injects ZERO net momentum into the two-cell system
     (internal adhesion, not an external field). ✓
   - A cadherin engages at most ONE trans-dimer at a time (single EC1
     interface) → trans-bond count ≤ ``n_cad_per_cell`` (≤ min over the two
     cells' free cadherins).

4. **Sign / sense**
   - Stretched trans-dimer (``L > r0_trans``): harmonic force pulls the two
     cadherins (hence the two cells) TOGETHER (inward, cohesive). STATIC test
     on a 2-particle snapshot. ``L ≤ r0_trans`` → the binder samples no
     tensile rupture (compression is held by excluded volume elsewhere; here
     ``F = k_trans·max(0, L−r0)`` clamps the off-rate argument to ≥ 0, matching
     ``single_pair_off_rate``'s tensile-only domain).

5. **CFL / stiffness note**
   - ``cadherin_trans`` is a stiff harmonic; its overdamped CFL
     ``dt ≤ cfl_safety·γ_cad/k_trans`` is gated by
     :func:`attach_cadherin_junction` when a per-bead drag is supplied
     (mirrors ``cortex/erm.py`` CFL gate). The catch-bond batch CFL above
     gates the binder's rate accuracy.

6. **Measurement-protocol consistency**
   - The per-dimer force ``F`` that enters ``k_off(F)`` is the SAME elastic
     bond force HOOMD applies (``k_trans·max(0, L−r0)``), so the stochastic
     lifetime is consistent with the mechanical load — no second force model.
   - The ENSEMBLE rupture force of a fully-engaged junction ≈
     ``n_cad·f0`` matches the Iturri-2020 MCF7–MCF7 de-adhesion force
     (6.5 nN) when ``n_cad ≈ 223`` (the scale bridge; reported, not tuned).

Compartment Performance Contract
--------------------------------
* **Particle types added**: ``cadherin`` (1 new type).
* **Particle count**: ``2 · n_cad_per_cell`` for a two-cell junction. At the
  default ``n_cad_per_cell`` derived from the Iturri scale bridge
  (``N_cad ≈ 223``) → ~446 cadherin particles for a cell pair. This is a
  per-INTERFACE count and does NOT scale with ``n_fil``: at the ×40 meso
  scale (n_fil=1000) and at native ~38000 filaments the cadherin count is the
  SAME (~223/cell) — it is set by the measured adhesion force, not by the
  actin discretisation.
* **Bond / angle count**: ``cadherin_trans`` dynamic bonds, ≤ ``n_cad_per_cell``
  (≤ ~223 per interface) at any instant; typically a bound fraction
  ``φ ≈ k_on/(k_on+k_off)`` of that. 0 angles, 0 dihedrals.
* **Per-step force**: NO. The trans-dimer elastic force is a standard
  ``md.bond.Harmonic`` evaluated by HOOMD's built-in bond kernel each step —
  this module adds no Python per-step ``md.force.Custom``.
* **Per-batch updater**: YES — :class:`CadherinTransJunctionUpdater`
  (``hoomd.custom.Action`` in a ``CustomUpdater``) runs every
  ``batch_steps`` steps to break/bind trans-dimers (D2 batched-updater
  convention, gated by the catch-bond batch CFL).
* **Uses cpu_local_snapshot**: NO in the hot path — the per-batch updater
  uses ``sim.state.get_snapshot()`` (rank-0 gather) like the D2 xlink binder,
  NOT a per-step ``cpu_local_snapshot``. (The static sign test reads forces via
  ``cpu_local_snapshot`` once, off the hot path.)
* **Uses cKDTree / broad-phase**: YES — per batch tick the partner search
  builds one ``scipy.spatial.cKDTree`` over cell-B cadherins and queries
  cell-A cadherins within ``r_bind`` (same idiom as the D2 xlink binder).
* **Hot-path priority**: **P2** (per-prompt). It is a per-batch event, not
  per-step; the per-step cost is just the built-in harmonic bond kernel.
* **GPU path now**: built-in (the harmonic bond force is native HOOMD GPU); the
  Python binder runs on the host each batch (CPU). A future GPU-resident
  partner search (cupy/native) is the optimisation lever, not required for P2.
* **Native ForceCompute candidate**: NO for the elastic force (already a
  native ``md.bond.Harmonic``). The stochastic break/bind binder is a native
  CUDA candidate only at native-N with many interfaces (not the single-pair
  H.7 case).
* **Bottleneck risk**: LOW at single-cell scale (≤ ~223 dimers/interface,
  batched). The cKDTree rebuild each batch is O(n_cad log n_cad) — negligible
  vs the cortex. Risk rises only if many cell-cell interfaces are active
  simultaneously (multi-cell, Layer-2-fine) — that is the native-port trigger.

References
----------
- Rakshit S, Zhang Y, Manibog K, Shafraz O, Sivasankar S (2012) "Ideal,
  catch, and slip bonds in cadherin adhesion." PNAS 109(46):18815-18820.
  DOI 10.1073/pnas.1208349109. (SI Table S1: k-1⁰=30.4 s⁻¹, x=0.34 nm,
  k+1=5.3 s⁻¹, k+2=1985.9 s⁻¹, f0=29.2 pN, n=4.8.) — via the reused
  ``validation/cadherin_sliding_rebinding.py``.
- Lou J, Zhu C (2007) Biophys J 92(5):1471-1485, DOI
  10.1529/biophysj.106.097048. (sliding-rebinding model.)
- Iturri J, Toca-Herrera JL, et al. (2020) "Single-Cell Probe Force Studies to
  Identify Sox2 Overexpression-Promoted Cell Adhesion in MCF7 Breast Cancer
  Cells." Cells 9(4):935. DOI 10.3390/cells9040935 (PMC7227807). — MCF7–MCF7
  SCFS de-adhesion ~6–7 nN @120 s (6.5 nN midpoint); value is a Fig. 5 read-off,
  not a prose figure; E-cadherin not separately characterised in their cells;
  corroborated Omidvar 2014/2016. The ``N_cad = F_detach/f0 ≈ 223`` scale bridge
  (reused from ``spheroid/cadherin_bonds.py``).
  Provenance note: this anchor sets only a particle COUNT (n_cad ≈ 223) and the
  ensemble-force REPORT — both grid-invariant and overridable via
  ``cfg['n_cad_per_cell']``/``cfg['k_trans']``; the per-dimer mechanics run off
  the fully-registered Rakshit f0/catch-bond. SourceEvidence registration in the
  Notion SoT is a Lead/PI registry action (see PI_DECISIONS); confirm its
  verdict is OK before citing this anchor in a deliverable (CLAUDE.md).
- KU-4.2 (E-cadherin full catch-bond, Phase-1 contract).
- Reused repo modules (cited above): ``validation/cadherin_sliding_rebinding.py``,
  ``spheroid/cadherin_bonds.py``, ``cortex/crosslinkers.py``,
  ``cell/cell.py`` (SubstrateLigandPin), ``cortex/erm.py`` (CFL gate idiom).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd

from ffn_sim.validation.cadherin_sliding_rebinding import (
    RAKSHIT_W2A,
    CadherinCatchParams,
    effective_k_off,
    mean_lifetime,
)

__all__ = [
    "GAMMA_DENYLIST_PREFIX",
    "PI_DECISIONS",
    "CADHERIN_TRANS_BOND",
    "CADHERIN_PARTICLE_TYPE",
    "ResolvedCadherinJunction",
    "resolve_cadherin_junction",
    "extend_snapshot_with_cadherins",
    "CadherinTransJunctionUpdater",
    "attach_cadherin_junction",
]

# ---------------------------------------------------------------------------
# γ-contamination guard: every bond this module creates carries a NON-cortical
# (cell-cell adhesion) load path and must be denylisted by the cortical-tension
# estimator (cortex/cortical_tension.py). All bond type names start with this.
# ---------------------------------------------------------------------------
GAMMA_DENYLIST_PREFIX: str = "cadherin_"

# The single dynamic trans-dimer bond type (note the prefix).
CADHERIN_TRANS_BOND: str = "cadherin_trans"
# The cadherin particle type added by this compartment.
CADHERIN_PARTICLE_TYPE: str = "cadherin"

# Per-batch rupture-rate accuracy budget [–]: the maximum k_off·Δt for a single
# first-order Poisson break draw (mirrors the D2 xlink Bell-Evans batch CFL).
# Used by BOTH the resolver (to size batch_steps for the low-force regime) and
# the runtime binder (to sub-step the break Bernoulli when a STRETCHED dimer's
# k_off·Δt_batch exceeds it, so the per-batch rupture probability stays a valid
# first-order rate at any reached force). Not a physical magic number — it is a
# numerical integration-accuracy tolerance (a discretisation budget), grid-/
# parameter-invariant, identical to the dynamic-bond CFL convention in cortex/.
_BATCH_CFL_BUDGET: float = 1.0e-3

# Open decisions surfaced to PI (empty list ⇒ none).
PI_DECISIONS: list[str] = [
    # k_trans (trans-dimer ectodomain elastic stiffness): no single clean
    # MEASURED single-molecule axial stiffness for the E-cadherin EC1-5
    # ectodomain at this coarse scale is in the repo KB. We therefore DERIVE
    # it from the same measured anchors the reused spheroid module uses
    # (F_detach / f0 = N_cad cadherins; the per-dimer elastic force reaches f0
    # over one contact-zone width) — see resolve_cadherin_junction. If PI wants
    # a directly-measured ectodomain spring constant instead, supply it via
    # cfg['k_trans'] (N/m) and it overrides the derivation. Flagged as a
    # derivation-vs-direct-measurement choice, NOT an invented number.
    "k_trans: derived from Iturri F_detach + Rakshit f0 over one contact-zone "
    "(grid-invariant); pass cfg['k_trans'] [N/m] to override with a direct "
    "single-molecule ectodomain stiffness if PI has one.",
    # r_bind (trans-dimer capture radius): the strand-swap engagement distance.
    # Set to the EC1 trans-interface reach; no tight single literature value at
    # ×40 scale, so it is a REQUIRED-with-documented-default parameter
    # (default = contact_zone_width, the same engagement length the spheroid
    # module uses). Flagged; override via cfg['r_bind'].
    "r_bind: defaults to contact_zone_width (engagement length reused from the "
    "spheroid scale bridge); override via cfg['r_bind'] [m] if a measured EC1 "
    "capture radius is available.",
    # SourceEvidence registration gap (Lead/PI registry action, NOT a code/number
    # change). The Iturri-2020 de-adhesion anchor (Cells 9(4):935, DOI
    # 10.3390/cells9040935, PMC7227807; MCF7–MCF7 ~6.5 nN, a Fig.5 read-off,
    # Omidvar 2014/2016-corroborated) is cited in this module/tests/configs and in
    # the paper-ready Layer2_Report.docx, but has NO SourceEvidence row in the
    # Notion SoT and so no anchor_status/verdict. Per CLAUDE.md ("confirm verdict
    # is OK before citing in a deliverable") this must be registered: author an
    # SE_REGISTRATION_CANDIDATES_<date>.md entry, PI-review, add the SourceEvidence
    # row in Notion (SoT), then `bash outputs/tag_kb/refresh.sh` + verify_sources.py
    # so it gets a verdict. Real source, correct number — registration hygiene, not
    # a fabrication. This module owns no registry file, so the action is the Lead's.
    "Iturri-2020 de-adhesion anchor (Cells 9(4):935, DOI 10.3390/cells9040935) "
    "is UNREGISTERED in the Notion SoT (no SourceEvidence row / no verdict) yet "
    "cited in deliverables: Lead/PI must stage SE_REGISTRATION_CANDIDATES + add "
    "the SourceEvidence row + refresh so the CLAUDE.md verdict-OK gate passes.",
]


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class ResolvedCadherinJunction:
    """Resolved explicit E-cadherin trans-dimer junction parameters (SI units).

    Attributes:
        enabled: master switch. When False this is an all-zero OFF record and
            every builder/attacher is a strict no-op.
        catch: faithful single-cadherin sliding-rebinding kinetics
            (Rakshit 2012 Table S1). The trans-dimer lifetime / off-rate is
            sampled from this (``effective_k_off``).
        n_cad_per_cell: cadherin particles seeded per cell [count]. From the
            Iturri/Rakshit scale bridge ``N_cad = F_detach / f0`` (~223) by
            default.
        k_trans: trans-dimer ectodomain elastic stiffness [N/m] (harmonic
            ``cadherin_trans`` bond). Derived (see resolver) or overridden.
        r0_trans: trans-dimer rest length [m] (force-free ectodomain bridge
            span across the interface).
        r_bind: trans-dimer capture radius [m] — an unbound A↔B cadherin pair
            within this distance can form a dimer.
        k_on: trans-dimer on-rate [s⁻¹] (reformation rate; rest-symmetric
            ``= k_off(0)`` by default, the same convention the spheroid module
            uses for the contact reformation rate).
        batch_steps: integration steps between binder ticks [count].
        n_cad: scale-bridge cadherins-per-interface used to report the ensemble
            rupture force ``n_cad·f0`` [count] (= n_cad_per_cell here).
    """

    enabled: bool
    catch: CadherinCatchParams
    n_cad_per_cell: int
    k_trans: float            # N/m
    r0_trans: float           # m
    r_bind: float             # m
    k_on: float               # s⁻¹
    batch_steps: int
    n_cad: float              # count (scale-bridge report)


def resolve_cadherin_junction(
    cfg: dict | None,
    *,
    kT: float,
    dt: float,
    contact_zone_width: float,
    deadhesion_force: float | None = None,
    catch: CadherinCatchParams = RAKSHIT_W2A,
) -> ResolvedCadherinJunction:
    """Resolve the cadherin-junction config block (default-OFF).

    Args:
        cfg: the ``junction.cadherin`` sub-dict (or a parent dict containing
            ``junction``/``cadherin``). ``None``, missing, or
            ``{'enabled': False}`` → an OFF record with zeros.
        kT: thermal energy [J] (for the binder; passed through for parity with
            the catch params' kT and reserved for future thermal terms).
        dt: integrator timestep [s] (for the batch CFL gate).
        contact_zone_width: engagement length [m] of the trans-interface — the
            extension over which the per-dimer elastic force reaches ``f0``.
            From the host cell's measured contact-zone (the spheroid scale
            bridge). Used for the derived ``k_trans`` and the default
            ``r_bind`` / ``r0_trans``.
        deadhesion_force: measured single-cell de-adhesion force [N]
            (Iturri 2020 ≈ 6.5e-9 N for MCF7). Sets the default
            ``n_cad_per_cell = round(F_detach/f0)`` and the ensemble report.
            If ``None`` and ``cfg`` does not give ``n_cad_per_cell``, the
            resolver raises (no invented count).
        catch: single-cadherin sliding-rebinding kinetics (Rakshit default).

    Returns:
        ResolvedCadherinJunction. OFF when disabled.

    Raises:
        ValueError: on negative/zero physical parameters, or when neither a
            de-adhesion force nor an explicit ``n_cad_per_cell`` is available
            to set the cadherin count (no invented number).
    """
    # ---- unwrap nesting ----
    if cfg is None:
        cfg = {}
    if "junction" in cfg:
        cfg = cfg["junction"]
    if "cadherin" in cfg:
        cfg = cfg["cadherin"]

    enabled = bool(cfg.get("enabled", False))

    if not enabled:
        # OFF-identity record: zeros, but a valid CadherinCatchParams so the
        # dataclass type is uniform. No builder will touch a sim with this.
        return ResolvedCadherinJunction(
            enabled=False,
            catch=catch,
            n_cad_per_cell=0,
            k_trans=0.0,
            r0_trans=0.0,
            r_bind=0.0,
            k_on=0.0,
            batch_steps=0,
            n_cad=0.0,
        )

    # ---- context validation ----
    if not (math.isfinite(kT) and kT > 0.0):
        raise ValueError(f"kT must be finite and > 0; got {kT}")
    if not (math.isfinite(dt) and dt > 0.0):
        raise ValueError(f"dt must be finite and > 0; got {dt}")
    if not (math.isfinite(contact_zone_width) and contact_zone_width > 0.0):
        raise ValueError(
            f"contact_zone_width must be finite and > 0; got {contact_zone_width}"
        )

    # ---- cadherin count per cell (measured scale bridge, no invented number) ----
    if "n_cad_per_cell" in cfg:
        n_cad_per_cell = int(cfg["n_cad_per_cell"])
    elif deadhesion_force is not None:
        if not (math.isfinite(deadhesion_force) and deadhesion_force > 0.0):
            raise ValueError(
                f"deadhesion_force must be finite and > 0; got {deadhesion_force}"
            )
        # N_cad = F_detach / f0  (Iturri 2020 / Rakshit f0; the reused spheroid
        # scale bridge). round() → an integer particle count.
        n_cad_per_cell = int(round(deadhesion_force / catch.f0))
    else:
        raise ValueError(
            "cadherin_junction needs a cadherin count: pass cfg['n_cad_per_cell'] "
            "or a measured deadhesion_force (Iturri 2020 ≈ 6.5e-9 N). No default "
            "count is invented (CLAUDE.md no-magic-number rule)."
        )
    if n_cad_per_cell <= 0:
        raise ValueError(
            f"n_cad_per_cell must be > 0; got {n_cad_per_cell}"
        )

    # ---- rest length & capture radius (engagement length; PI-flagged) ----
    r0_trans = float(cfg.get("r0_trans", contact_zone_width))
    r_bind = float(cfg.get("r_bind", contact_zone_width))
    if not (math.isfinite(r0_trans) and r0_trans > 0.0):
        raise ValueError(f"r0_trans must be finite and > 0; got {r0_trans}")
    if not (math.isfinite(r_bind) and r_bind > 0.0):
        raise ValueError(f"r_bind must be finite and > 0; got {r_bind}")

    # ---- trans-dimer elastic stiffness (derived, grid-invariant; overridable) ----
    # Per-dimer elastic force f reaches the catch transition f0 at an extension
    # of one contact-zone width: k_trans = f0 / contact_zone_width. This is the
    # SAME force-reaches-f0-over-one-zone construction the reused spheroid
    # module uses for its contact stiffness (there k_bond = N_cad·f0/zone is the
    # ENSEMBLE spring; here we resolve a SINGLE dimer so the per-dimer spring is
    # f0/zone). It is grid-invariant (no n_fil, no n_cad dependence) and NOT
    # tuned to a gate. cfg['k_trans'] overrides with a direct measurement.
    if "k_trans" in cfg:
        k_trans = float(cfg["k_trans"])
    else:
        k_trans = float(catch.f0 / contact_zone_width)
    if not (math.isfinite(k_trans) and k_trans > 0.0):
        raise ValueError(f"k_trans must be finite and > 0; got {k_trans}")

    # ---- on-rate (rest-symmetric reformation = k_off(0); spheroid convention) ----
    if "k_on" in cfg:
        k_on = float(cfg["k_on"])
    else:
        k_on = 1.0 / mean_lifetime(0.0, catch)   # = effective_k_off(0)
    if not (math.isfinite(k_on) and k_on > 0.0):
        raise ValueError(f"k_on must be finite and > 0; got {k_on}")

    # ---- batch CFL gate (mirror the D2 xlink Bell-Evans batch CFL) ----
    # batch_steps · dt · k_off_max ≤ 1e-3 sizes batch_steps for the *low-force*
    # regime the bond DWELLS in. It does NOT (and cannot) bound the slip tail,
    # because a bound trans-dimer is stretched by the dynamics with NO ceiling at
    # r_bind (r_bind only governs binding CAPTURE, never the bound-state stretch),
    # so the per-dimer force the runtime samples is unbounded above. The previous
    # `f_envelope = k_trans·(r_bind − r0_trans)` was the WRONG bound on two counts:
    # (1) it is identically 0 at the default r_bind == r0_trans, collapsing the
    # envelope to k_off(0); (2) even when r_bind > r0_trans it bounds the capture
    # force, not the runtime stretch. The slip tail is instead handled at RUNTIME
    # by the binder's adaptive sub-stepping (see CadherinTransJunctionUpdater.act:
    # any dimer whose k_off(F)·Δt_batch exceeds the budget has its break Bernoulli
    # sub-stepped into first-order draws), which keeps the per-batch rupture
    # probability a valid first-order rate at ANY reached force WITHOUT an invented
    # max-force magic number. Here we therefore size batch_steps on the largest
    # k_off the bond's LIFETIME landscape exposes at low force — the rest rate
    # k_off(0) and the catch peak f0 — and let the runtime sub-stepper cover the
    # tail. (k_off(0) > k_off(f0) because this is a CATCH bond: lifetime is LONGER
    # at f0, so the rest rate is the relevant low-force bound.) Shrink batch_steps
    # to satisfy the bound (never grow it past the requested value).
    requested_batch = int(cfg.get("batch_steps", 100))
    if requested_batch <= 0:
        raise ValueError(f"batch_steps must be > 0; got {requested_batch}")
    k_off_max = max(
        effective_k_off(0.0, catch),
        effective_k_off(float(catch.f0), catch),
    )
    cfl_budget = _BATCH_CFL_BUDGET
    max_batch = int(cfl_budget / (dt * k_off_max))
    if max_batch < 1:
        # Even a single-step batch (the minimum) overshoots the rate-accuracy
        # budget: the timestep itself is too coarse for the catch-bond kinetics.
        # Mirror the D2 xlink resolver — raise rather than silently return a
        # CFL-violating value (no gate-loosening).
        raise ValueError(
            f"cadherin_junction batch CFL impossible to satisfy: dt={dt:.3e} s "
            f"> 1e-3 / k_off_max ({cfl_budget / k_off_max:.3e} s). Reduce dt "
            "(the catch-bond kinetics need a finer timestep) or surface to PI."
        )
    batch_steps = min(requested_batch, max_batch)

    return ResolvedCadherinJunction(
        enabled=True,
        catch=catch,
        n_cad_per_cell=int(n_cad_per_cell),
        k_trans=k_trans,
        r0_trans=r0_trans,
        r_bind=r_bind,
        k_on=k_on,
        batch_steps=int(batch_steps),
        n_cad=float(n_cad_per_cell),
    )


# ---------------------------------------------------------------------------
# Snapshot extension: append cadherin particles for two cells' surfaces
# ---------------------------------------------------------------------------
def extend_snapshot_with_cadherins(
    snap,
    p: ResolvedCadherinJunction,
    *,
    cell_a_surface_points: np.ndarray | None = None,
    cell_b_surface_points: np.ndarray | None = None,
):
    """Append ``cadherin`` particles for cell A and cell B to ``snap`` (in place).

    Mirrors ``cell/nucleus.py:_extend_snapshot_with_nucleus`` (append particles,
    register the new type, return tag bookkeeping). A strict NO-OP when
    ``p.enabled`` is False or no surface points are given.

    Args:
        snap: a ``gsd.hoomd.Frame`` (or HOOMD snapshot) being assembled. Must
            already have ``particles.N``, ``.types``, ``.typeid``, ``.position``
            set. Modified in place.
        p: resolved junction params.
        cell_a_surface_points: shape (n_cad_per_cell, 3) seed positions for
            cell A's cadherin tips [m]. If None, no A cadherins are added.
        cell_b_surface_points: shape (n_cad_per_cell, 3) seed positions for
            cell B's cadherin tips [m]. If None, no B cadherins are added.

    Returns:
        dict with keys:
            'cadherin_tag_start', 'n_cadherin', 'cell_a_tags', 'cell_b_tags',
            'cadherin_typeid'.
        On OFF / no-op: n_cadherin=0 and empty tag arrays; ``snap`` unchanged.
    """
    n_before = int(snap.particles.N)
    empty = {
        "cadherin_tag_start": n_before,
        "n_cadherin": 0,
        "cell_a_tags": np.empty((0,), dtype=np.int64),
        "cell_b_tags": np.empty((0,), dtype=np.int64),
        "cadherin_typeid": -1,
    }
    if not p.enabled:
        return empty
    if cell_a_surface_points is None and cell_b_surface_points is None:
        return empty

    pts_a = (
        np.empty((0, 3), dtype=np.float64)
        if cell_a_surface_points is None
        else np.asarray(cell_a_surface_points, dtype=np.float64)
    )
    pts_b = (
        np.empty((0, 3), dtype=np.float64)
        if cell_b_surface_points is None
        else np.asarray(cell_b_surface_points, dtype=np.float64)
    )
    for nm, pts in (("cell_a_surface_points", pts_a), ("cell_b_surface_points", pts_b)):
        if pts.ndim != 2 or (pts.size and pts.shape[1] != 3):
            raise ValueError(f"{nm} must have shape (n, 3); got {pts.shape}")
    new_pts = np.vstack([pts_a, pts_b]) if (pts_a.size or pts_b.size) else None
    if new_pts is None:
        return empty
    n_new = int(new_pts.shape[0])

    # Register the cadherin particle type (append if absent).
    types = list(snap.particles.types)
    if CADHERIN_PARTICLE_TYPE not in types:
        types.append(CADHERIN_PARTICLE_TYPE)
    cad_typeid = types.index(CADHERIN_PARTICLE_TYPE)
    snap.particles.types = types

    old_pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    old_typeid = np.asarray(snap.particles.typeid).reshape(-1)
    new_typeid_block = np.full(n_new, cad_typeid, dtype=old_typeid.dtype)

    snap.particles.N = n_before + n_new
    snap.particles.position = np.vstack([old_pos, new_pts])
    snap.particles.typeid = np.concatenate([old_typeid, new_typeid_block])

    # Extend any other per-particle arrays that are already present, so the
    # frame stays internally consistent (mass/velocity/etc.).
    for attr, fill in (("mass", 1.0), ("charge", 0.0), ("diameter", 0.0)):
        arr = getattr(snap.particles, attr, None)
        if arr is not None and np.asarray(arr).shape[0] == n_before:
            setattr(
                snap.particles,
                attr,
                np.concatenate(
                    [np.asarray(arr, dtype=np.float64), np.full(n_new, fill)]
                ),
            )
    vel = getattr(snap.particles, "velocity", None)
    if vel is not None and np.asarray(vel).reshape(-1, 3).shape[0] == n_before:
        snap.particles.velocity = np.vstack(
            [np.asarray(vel, dtype=np.float64).reshape(-1, 3), np.zeros((n_new, 3))]
        )
    # image flags (int32, (n,3)): extend too, or a host frame that carries an
    # image array desyncs from N and crashes Snapshot.from_gsd_frame at state
    # creation (broadcast (n_before,3) into (N,3)). Mirrors the velocity carry +
    # the nucleus/MT/IF extenders which all carry image.
    img = getattr(snap.particles, "image", None)
    if img is not None and np.asarray(img).reshape(-1, 3).shape[0] == n_before:
        snap.particles.image = np.vstack(
            [np.asarray(img, dtype=np.int32).reshape(-1, 3),
             np.zeros((n_new, 3), dtype=np.int32)]
        )

    a_tags = np.arange(n_before, n_before + pts_a.shape[0], dtype=np.int64)
    b_tags = np.arange(
        n_before + pts_a.shape[0], n_before + n_new, dtype=np.int64
    )
    return {
        "cadherin_tag_start": n_before,
        "n_cadherin": n_new,
        "cell_a_tags": a_tags,
        "cell_b_tags": b_tags,
        "cadherin_typeid": int(cad_typeid),
    }


# ---------------------------------------------------------------------------
# Dynamic trans-dimer break/bind updater (per-batch Action)
# ---------------------------------------------------------------------------
class CadherinTransJunctionUpdater(hoomd.custom.Action):
    """Per-batch explicit E-cadherin trans-dimer break/bind across two cells.

    Mirrors the D2 ``cortex/crosslinkers.py`` dynamic-bond Action: every
    ``batch_steps`` it (1) ruptures engaged ``cadherin_trans`` bonds with the
    faithful sliding-rebinding ``k_off(F)`` and (2) forms new trans-dimers
    between unbound A↔B cadherins within ``r_bind`` via a cKDTree partner
    search — then rewrites ``bonds.group`` / ``bonds.typeid``.

    The cell-A and cell-B cadherin tag sets define the interface; a trans-dimer
    only ever joins an A cadherin to a B cadherin (no intra-cell cadherin
    bonds), so this is a genuine TWO-CELL junction.

    Parameters
    ----------
    p : ResolvedCadherinJunction
    cell_a_tags, cell_b_tags : np.ndarray
        Global HOOMD tags of the two cells' cadherin particles.
    seed : int
        RNG seed for the stochastic break/bind draws.
    """

    def __init__(
        self,
        p: ResolvedCadherinJunction,
        *,
        cell_a_tags: np.ndarray,
        cell_b_tags: np.ndarray,
        seed: int = 0,
    ) -> None:
        super().__init__()
        if not p.enabled:
            raise ValueError(
                "CadherinTransJunctionUpdater requires an ENABLED "
                "ResolvedCadherinJunction (default-off contract)."
            )
        self.p = p
        self._a_tags = np.asarray(cell_a_tags, dtype=np.int64)
        self._b_tags = np.asarray(cell_b_tags, dtype=np.int64)
        self._rng = np.random.default_rng(seed)
        self._sim_ref: hoomd.Simulation | None = None
        self.batch_dt = float(p.batch_steps) * 0.0  # set in attach() from dt
        # break/bind counters (diagnostics).
        self.n_break_total = 0
        self.n_bind_total = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation
        ig = simulation.operations.integrator
        dt = float(ig.dt) if ig is not None else 0.0
        self.batch_dt = float(self.p.batch_steps) * dt

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        if sim is None:
            return
        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32).copy()
        bond_type_names = list(read_snap.bonds.types)

        if CADHERIN_TRANS_BOND not in bond_type_names:
            raise RuntimeError(
                "snap.bonds.types must include "
                f"'{CADHERIN_TRANS_BOND}'; got {bond_type_names}. "
                "Call attach_cadherin_junction() (it registers the bond type) "
                "before running the updater."
            )
        trans_typeid = bond_type_names.index(CADHERIN_TRANS_BOND)

        is_trans = bt == trans_typeid
        trans_bonds = bg[is_trans]
        other_bg = bg[~is_trans]
        other_bt = bt[~is_trans]

        # Track which cadherins are currently engaged (one dimer per cadherin).
        a_set = set(self._a_tags.tolist())
        b_set = set(self._b_tags.tolist())
        bound = set()
        for row in trans_bonds:
            bound.add(int(row[0]))
            bound.add(int(row[1]))

        # ---- Step 1: break engaged trans-dimers via faithful k_off(F) ----
        # ADAPTIVE batch-CFL (slip-tail-correct): a bound trans-dimer is stretched
        # by the dynamics with NO ceiling at r_bind, so F_mag — hence k_off(F) —
        # can exceed the value the resolver sized batch_steps for. If we drew a
        # single per-batch Bernoulli p = 1 − exp(−k_off·Δt_batch) at a slip-tail
        # k_off, that draw would violate the rate-accuracy budget the CFL is meant
        # to enforce (k_off·Δt_batch ≫ _BATCH_CFL_BUDGET). We therefore sub-step
        # the break decision PER DIMER: split Δt_batch into n_sub equal slices so
        # each slice's k_off·Δt_sub ≤ _BATCH_CFL_BUDGET, and rupture if ANY slice
        # fires. This keeps every Bernoulli a valid first-order rate at any reached
        # force, with no invented max-force constant (the bound is derived from the
        # actual current force each tick). The composite survival over n_sub equal
        # constant-rate slices equals the single exact-survival 1 − exp(−k·Δt) (the
        # rate is constant within a batch because positions are frozen during act),
        # so sub-stepping changes only the per-draw accuracy, not the expected rate.
        kept_rows = []
        if trans_bonds.shape[0] > 0:
            ta = trans_bonds[:, 0]
            tb = trans_bonds[:, 1]
            r = np.linalg.norm(pos[ta] - pos[tb], axis=1)
            ext = np.clip(r - self.p.r0_trans, 0.0, None)
            F_mag = self.p.k_trans * ext
            # Faithful sliding-rebinding effective off-rate per dimer (catch-slip).
            k_off = np.array(
                [effective_k_off(float(f), self.p.catch) for f in F_mag],
                dtype=np.float64,
            )
            # Per-dimer sub-step count so each slice obeys k_off·Δt_sub ≤ budget.
            kdt = k_off * self.batch_dt
            n_sub = np.maximum(
                1, np.ceil(kdt / _BATCH_CFL_BUDGET).astype(np.int64)
            )
            # Defensive cap: an EXTREME-slip dimer (force-runaway under an
            # un-equilibrated dynamic run) drives k_off → ∞, so the exact n_sub
            # would overflow the RNG allocation. Such a dimer ruptures with
            # probability ≈ 1 regardless, so cap the sub-step count and let the
            # per-slice draw certify the (certain) rupture rather than crash.
            _NSUB_MAX = 100_000
            for i, row in enumerate(trans_bonds):
                ns = min(int(n_sub[i]), _NSUB_MAX)
                # First-order per-slice rupture probability; rupture if any slice
                # fires. n_sub draws keep each slice's k·Δt_sub within budget.
                p_slice = 1.0 - math.exp(-k_off[i] * self.batch_dt / ns)
                u = self._rng.uniform(0.0, 1.0, size=ns)
                if np.any(u < p_slice):
                    self.n_break_total += 1
                    bound.discard(int(row[0]))
                    bound.discard(int(row[1]))
                else:
                    kept_rows.append((int(row[0]), int(row[1])))

        # ---- Step 2: bind unbound A↔B cadherins within r_bind (cKDTree) ----
        free_a = np.array(
            [t for t in self._a_tags.tolist() if t not in bound], dtype=np.int64
        )
        free_b = np.array(
            [t for t in self._b_tags.tolist() if t not in bound], dtype=np.int64
        )
        new_rows = []
        if free_a.size > 0 and free_b.size > 0:
            try:
                from scipy.spatial import cKDTree
            except ImportError as exc:  # pragma: no cover - env guard
                raise RuntimeError(
                    "scipy.spatial.cKDTree required for the cadherin trans-dimer "
                    "partner search."
                ) from exc
            tree_b = cKDTree(pos[free_b])
            nbr_lists = tree_b.query_ball_point(pos[free_a], r=self.p.r_bind)
            p_bind = 1.0 - math.exp(-self.p.k_on * self.batch_dt)
            u2 = self._rng.uniform(0.0, 1.0, size=free_a.size)
            b_used: set[int] = set()
            for k, nbrs in enumerate(nbr_lists):
                if len(nbrs) == 0 or u2[k] >= p_bind:
                    continue
                nbrs_arr = np.asarray(nbrs, dtype=np.int64)
                # nearest free B cadherin not already claimed this tick.
                d = np.linalg.norm(pos[free_b[nbrs_arr]] - pos[free_a[k]], axis=1)
                order = np.argsort(d)
                chosen = -1
                for j in order:
                    b_tag = int(free_b[nbrs_arr[j]])
                    if b_tag not in b_used:
                        chosen = b_tag
                        break
                if chosen < 0:
                    continue
                a_tag = int(free_a[k])
                if a_tag == chosen:
                    continue  # defensive (A and B sets are disjoint by construction)
                b_used.add(chosen)
                new_rows.append((a_tag, chosen))
                self.n_bind_total += 1

        # ---- Step 3: rebuild the bond topology ----
        all_trans = kept_rows + new_rows
        if all_trans:
            trans_arr = np.array(all_trans, dtype=np.int64)
            trans_tids = np.full(trans_arr.shape[0], trans_typeid, dtype=np.uint32)
        else:
            trans_arr = np.empty((0, 2), dtype=np.int64)
            trans_tids = np.empty((0,), dtype=np.uint32)

        if other_bg.shape[0] > 0:
            new_bg = np.concatenate(
                [trans_arr.astype(np.int64), other_bg.astype(np.int64)], axis=0
            )
            new_bt = np.concatenate([trans_tids, other_bt.astype(np.uint32)])
        else:
            new_bg = trans_arr.astype(np.int64)
            new_bt = trans_tids

        write_snap = hoomd.Snapshot()
        N_part = int(read_snap.particles.N)
        write_snap.particles.N = N_part
        write_snap.particles.types = list(read_snap.particles.types)
        write_snap.particles.typeid[:] = np.asarray(read_snap.particles.typeid)
        write_snap.particles.position[:] = pos
        write_snap.particles.velocity[:] = np.asarray(read_snap.particles.velocity)
        write_snap.particles.mass[:] = np.asarray(read_snap.particles.mass)
        write_snap.particles.image[:] = np.asarray(read_snap.particles.image)
        write_snap.configuration.box = list(read_snap.configuration.box)

        write_snap.bonds.N = int(new_bg.shape[0])
        write_snap.bonds.types = list(bond_type_names)
        if new_bg.shape[0] > 0:
            write_snap.bonds.group[:] = new_bg
            write_snap.bonds.typeid[:] = new_bt

        for grp_name in ("angles", "dihedrals", "impropers"):
            src = getattr(read_snap, grp_name)
            dst = getattr(write_snap, grp_name)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)

        sim.state.set_snapshot(write_snap)


# ---------------------------------------------------------------------------
# Public helper: register the trans bond force + attach the binder
# ---------------------------------------------------------------------------
def attach_cadherin_junction(
    sim: hoomd.Simulation,
    p: ResolvedCadherinJunction,
    *,
    cell_a_tags: np.ndarray,
    cell_b_tags: np.ndarray,
    seed: int = 0,
    gamma_cad: float | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
):
    """Wire the explicit trans-dimer junction into a built simulation.

    Adds (a) an ``md.bond.Harmonic`` carrying the ``cadherin_trans`` elastic
    force and (b) a periodic ``CustomUpdater`` running
    :class:`CadherinTransJunctionUpdater` every ``p.batch_steps`` steps.
    A strict NO-OP returning ``(None, None)`` when ``p.enabled`` is False.

    The caller is responsible for having appended ``cadherin_trans`` to the
    snapshot's ``bonds.types`` at build time (so HOOMD knows the type); this
    helper sets that bond type's harmonic params. Because a HOOMD 7.0.1
    ``md.bond.Harmonic`` demands params for EVERY bond type registered in the
    state (not only the one it logically owns — an absent type raises
    ``IncompleteSpecificationError`` at ``sim.run``), this helper zero-stiffness
    fills (``k=0, r0=0``) every OTHER registered bond type on its standalone
    Harmonic before setting the real ``cadherin_trans`` stiffness, so it is safe
    to append to an integrator that already carries cortex / xlink / FA / etc.
    bonds.

    Args:
        sim: a built simulation with an Integrator and cadherin particles.
        p: resolved (enabled) junction params.
        cell_a_tags, cell_b_tags: the two cells' cadherin global tags.
        seed: binder RNG seed.
        gamma_cad: per-cadherin Stokes drag [N·s/m]; if given, gates the
            overdamped CFL ``dt ≤ cfl_safety·γ_cad/k_trans`` (cortex/erm.py
            idiom). If None the elastic-bond CFL gate is skipped.
        cfl_safety_factor: CFL safety (default 0.1, the D3 convention).
        cfl_strict: raise on CFL violation when True.

    Returns:
        (harmonic_force, updater) — the attached ``md.bond.Harmonic`` and the
        ``CustomUpdater`` wrapping the binder. ``(None, None)`` when disabled.

    Raises:
        RuntimeError: no Integrator, or the elastic CFL is violated (strict).
    """
    if not p.enabled:
        return None, None

    import hoomd.md as md

    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching the "
            "cadherin junction."
        )

    # Overdamped elastic CFL on the stiff trans-dimer harmonic (erm.py idiom).
    if gamma_cad is not None:
        dt = float(ig.dt)
        tau_cad = gamma_cad / p.k_trans
        dt_cfl = cfl_safety_factor * tau_cad
        if dt > dt_cfl and cfl_strict:
            raise RuntimeError(
                f"cadherin_trans CFL violated: dt = {dt:.3e} s > "
                f"{cfl_safety_factor:.2f}·τ = {dt_cfl:.3e} s "
                f"(τ = γ_cad/k_trans = {tau_cad:.3e} s). Reduce dt or soften "
                "k_trans (the latter changes physics — PI sign-off). Pass "
                "cfl_strict=False for diagnostic runs."
            )

    # Harmonic force for the dynamic trans-dimer bond. The bond type must
    # already exist in the state; we set its params here. CRITICAL HOOMD 7.0.1
    # contract: an md.bond.Harmonic ForceCompute demands params for EVERY bond
    # type registered in sim.state — NOT just the ones it "owns" — or sim.run()
    # raises IncompleteSpecificationError. A real cell build carries many bond
    # types (cortex-bond, arp_branch, xlink_intra/_attach, FA integrin/ligand,
    # lamellipodium, …), so this standalone Harmonic must zero-stiffness EVERY
    # other registered type (k=0, r0=0 → no force, no rest length) and then set
    # the real cadherin_trans stiffness. (Mirrors cell.py:1099-1160 and
    # connected_mesh.py:162-175, which register every present bond type on one
    # Harmonic; cell.py:1111 carries the same "else md.bond.Harmonic would demand
    # params for an absent type" note.) The cleanest long-run wiring is for the
    # host cell build to register cadherin_trans on its SINGLE Harmonic; this
    # zero-fill path keeps the standalone helper safe when that has not happened.
    harmonic = md.bond.Harmonic()
    for bt in sim.state.bond_types:
        if bt != CADHERIN_TRANS_BOND:
            harmonic.params[bt] = dict(k=0.0, r0=0.0)
    harmonic.params[CADHERIN_TRANS_BOND] = dict(k=p.k_trans, r0=p.r0_trans)
    ig.forces.append(harmonic)

    binder = CadherinTransJunctionUpdater(
        p, cell_a_tags=cell_a_tags, cell_b_tags=cell_b_tags, seed=seed
    )
    updater = hoomd.update.CustomUpdater(
        action=binder, trigger=hoomd.trigger.Periodic(p.batch_steps)
    )
    sim.operations.updaters.append(updater)
    sim._cadherin_binder = binder  # noqa: SLF001 lifetime anchor
    return harmonic, updater
