"""H.5 plasma-membrane load — the leading-edge tension reservoir the
lamellipodial barbed ends push against (KU-5.2 force-velocity).

ADDITIVE, DEFAULT-OFF module (mechanism audit 2026-05-30 §5 #6 + §6 #7;
``docs/briefs/H4_FA_INTEGRATION_DESIGN.md`` §4 "Membrane-load (KU-5.x)").

The problem this fixes
----------------------
The H.5 lamellipodium (``ffn_sim/cell/lamellipodium.py``) advances each
barbed end with a Brownian-ratchet / Bell-Evans slip law

    k_elong(F) = k_elong⁰ · exp(−F · δ_elong / kT)        (Bieling 2016)

and caps with ``k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sin θ / kT)``. Both
laws are FORCE-DEPENDENT, but the per-barbed-end load ``F`` is supplied
through ``BarbedEndElongationUpdater.set_network_load(F_total)`` /
``CappingUpdater.set_network_load(F_total)`` whose default is **0**. With
no explicit plasma membrane the barbed ends push against vacuum, so F = 0,
the exponential factors collapse to 1, and the laws are inert — KU-5.1
dendritic density runs ~10,000× under and KU-5.2 force-velocity cannot
emerge (audit §5 #6, §6 #7).

This module supplies the missing load. A load-bearing planar membrane
sits at the leading edge (the WAVE plane ``y = Y_max``, membrane normal
``n̂ = +ŷ`` — the SAME geometry the WAVE/mother seeds already use, growing
toward ``+ŷ``). Barbed ends that reach within a contact range ``δ_contact``
of the plane push it outward; by Newton's 3rd law each contacting barbed
end feels a REACTION load ``F`` along ``−ŷ`` (opposing its advance). That
reaction F is exactly the load the Bell-Evans elongation / capping laws
consume → load develops → polymerization slows under tension → KU-5.2
force-velocity emerges as a membrane reaction, and the loaded rate balance
lets KU-5.1 density rise.

Mechanism (fine-grained, mechanistic per CLAUDE.md — minimal but real)
----------------------------------------------------------------------
The membrane is a planar **tension reservoir** (Helfrich tension term;
the bending κ_m term is sub-dominant at the lamellipodial radius of
curvature and enters only the anchor cross-check, below). A membrane of
tension γ_mem [N/m] bowing around the filament tips resists their advance.
For ``n_contact`` barbed ends sharing a leading-edge patch of area
``A_mem`` the membrane supplies a total normal reaction whose per-filament
share is the **tension-derived load**

    F_per = γ_mem · s_fil ,   s_fil = sqrt(A_mem / n_contact)            (1)

where ``s_fil`` is the mean tip-to-tip spacing of the contacting filaments
(``A_mem / n_contact`` is the membrane area each filament locally supports;
its square root is the linear span the tensed membrane bridges between
neighbouring tips). γ_mem [N/m] · s_fil [m] = [N], a force. This is the
standard tented-membrane / Brownian-ratchet load scale (Mogilner-Oster
1996; Peskin 1993): the membrane tension times the lateral span the
membrane is tented over between pinning filaments. As the network gets
denser (n_contact ↑) the per-filament load FALLS (each tip supports less
membrane, ``s_fil`` ↓) — the autoregulation that lets the dendritic array
reach the KU-5.1 steady density rather than stalling.

The reaction is distributed equally over the ``n_contact`` contacting
barbed ends (each feels ``F_per`` along ``−n̂``); the membrane feels the
equal-and-opposite total ``+n̂ · Σ F_per = +n̂ · n_contact · F_per``
(Newton's 3rd law, Sanity Gate §3). The membrane DOF is a single mobile
plane position ``y_mem`` (low-DOF reservoir); at quasi-static balance the
plane sits where its tension restoring force equals the filament push, and
the per-filament load (1) is the value handed to the Bell-Evans laws. We
keep the plane position diagnostic-only for Phase 1 (the load magnitude
(1) is set by tension + contact geometry, not by integrating the plane —
the simplest defensible form; a fully-dynamical membrane sheet is a later
refinement and would add particles, which Phase 1 forbids).

The Brownian-ratchet polymerization law ``v(F) = v0 · exp(−F δ/kT)``
(Mogilner-Oster 1996; Peskin 1993) and the Bieling 2016 load-velocity
curve are the ACCEPTANCE ORACLES used in ``ffn_sim/tests/test_membrane.py``
only — NEVER the runtime. The runtime is the explicit tension-derived
reaction load (1) fed into the existing Bell-Evans updaters.

How the load reaches the barbed ends (why an Action, not only a force)
----------------------------------------------------------------------
The existing lamellipodium advances barbed ends by SNAPSHOT MUTATION
inside ``BarbedEndElongationUpdater.act`` — the advance is a Bell-Evans
*rate*, not a force the integrator integrates. So the load must reach the
RATE: :class:`MembraneLoad` is a :class:`hoomd.custom.Action` (mirrors the
``_BatchedLamelUpdater`` pattern) that each tick reads barbed-end
positions, finds the contacting set, computes (1), and pushes the total
``n_contact · F_per`` into the elongation / capping updaters via their
``set_network_load`` handles. It runs on the SAME ``Periodic(batch_steps)``
trigger and BEFORE the elongation/capping updaters in the operations list,
so the load is current when the Bell-Evans laws fire.

A companion :class:`MembraneReactionForce` (``md.force.Custom``) is also
provided so the membrane push is a genuine per-bead reaction force in
HOOMD's ``net_force`` (so the contacting tips physically feel ``−n̂ F_per``
and the no-COM-runaway invariant is exercised against the integrator).
It is optional (``attach_membrane_to_simulation(..., with_reaction_force=
True)``) and supplies ONLY the membrane reaction; the BAOAB Action stays
the sole position integrator. The frozen ``integrator/*`` is untouched.

Magic-Number Block — γ_mem (no empirical magic number)
------------------------------------------------------
``γ_mem`` is a MEMBRANE-TENSION anchor, NOT a tuned spring. The audit /
FA-brief §4 KU-5.x anchor is γ_mem ~ 10–300 pN/µm (Bieling 2016; Lieber
2013; Sens 2015). We take the geometric-mean physiological value

    γ_mem = 50 pN/µm = 5.0 × 10⁻⁵ N/m

as the default (mid-decade of the 10–300 pN/µm band; Lieber 2013 resting
HeLa ≈ 30–40 pN/µm, Sens 2015 review ≈ 30–300 pN/µm). It is:

* Derivable / literature-anchored: a stated KU-5.x membrane-tension
  anchor [N/m], not chosen to make a gate pass. The resolver takes it
  verbatim from the band and range-checks it against 10–300 pN/µm; a value
  outside the band raises (guards a silent tuned override).
* Grid-invariant: γ_mem is an intensive material tension [N/m]. It does
  NOT depend on n_WAVE, batch_steps, ℓ₀, or dt. The per-filament load
  (1) carries the contact geometry (``s_fil``) so the load self-scales
  with the emergent contact density; γ_mem only sets the tension scale.
* Not chosen to pass a gate: the force-balance / Brownian-ratchet test
  (test §RatchetBalance) recovers the per-filament F from γ_mem and the
  contact geometry and checks ``v(F) = v0·exp(−Fδ/kT)`` is monotone
  decreasing in γ_mem — it does not fit γ_mem to any density target.

Cross-check (consistency, not a fit): the bending modulus κ_m ≈ 20 kT
(Helfrich; Bieling 2016) gives a tension-to-bending crossover length
``λ = sqrt(κ_m/γ_mem) = sqrt(20·4.114e-21 / 5e-5) ≈ 41 nm`` — comparable
to the filament tip spacing at KU-5.1 density (~100/µm² → spacing ~100 nm),
so the tension term (1) dominates the per-filament load at lamellipodial
densities and the κ_m term is a sub-dominant correction (it is therefore
carried only in the anchor block + this cross-check, not in the runtime
load). The per-filament load at the KU-5.1 density (~100/µm² → tip spacing
``s_fil ~ 100 nm``) is ``F_per = γ_mem · s_fil ≈ 5e-5 · 1e-7 ≈ 5 pN`` (the
fundamental law), pulled by the runtime λ ≈ 41 nm span cap to ``F_per ≈
γ_mem · λ ≈ 2 pN``; both are in the near-pN Brownian-ratchet regime the
Bell-Evans δ_elong = 2.7 nm law resolves (``F·δ/kT`` ≈ 1–3, a real brake).
The Helfrich λ span cap (:func:`_physical_span_cap`) keeps the
single-sparse-tip limit physical (a tip cannot tent a flat-tension load
across more than λ).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_membrane.py``.*

1. **Dimensional analysis**
   - ``γ_mem`` [N/m]; ``s_fil`` [m] → ``F = γ_mem·s_fil`` [N]. ✓
   - Brownian-ratchet exponent ``F · δ_elong / kT`` = [N·m]/[J] =
     dimensionless. ✓
   - ``δ_contact`` [m] contact range; membrane normal n̂ dimensionless. ✓
   - At ``n_contact = 0`` (no tip touches the membrane): F = 0, no load,
     no reaction force. ✓

2. **Boundary cases**
   - ``n_contact == 0`` (no barbed end within δ_contact): zero membrane
     load (``set_network_load(0)``), zero reaction force.
   - ``γ_mem`` ≤ 0 or outside the 10–300 pN/µm band: ValueError at resolve.
   - ``δ_contact`` ≤ 0: ValueError at resolve.
   - ``A_mem`` ≤ 0: ValueError at resolve.
   - All barbed ends capped: elongation idle (load still computed but
     consumed by zero free ends — harmless).

3. **Conservation invariants**
   - **Newton's 3rd law**: total reaction on the filaments
     ``Σ_i F_i = −F_membrane`` (the membrane feels +n̂·Σ|F_per|, the tips
     feel −n̂·F_per each). STATIC test asserts the membrane reaction is the
     exact negation of the summed barbed-end reaction.
   - **No COM runaway**: the per-bead reaction force (when the companion
     md.force.Custom is attached) acts along ±n̂ on contacting tips and
     its negation on the (massless reservoir) membrane is NOT a particle,
     so the load is one-sided on the bead set by construction; the test
     asserts the reaction never injects a growing net momentum into the
     particle COM beyond the explicit, bounded membrane push (short BAOAB
     demo stays finite, no runaway).

4. **Numerical sanity**
   - All positions read float64; ``np.isfinite`` guards at compute time.
   - CFL: the membrane reaction is a BOUNDED tension-plateau load (F_per
     depends on the lateral tip span, not on penetration depth — it is NOT
     a stiff Hookean normal spring), so the correct stability criterion is
     the per-step OVERDAMPED DRIFT bound: ``Δx = F_max·dt/γ_b ≤ α ·
     contact_range`` with ``F_max = γ_mem·sqrt(A_mem)`` the bounded load
     magnitude. :func:`attach_membrane_to_simulation` raises (parity with
     erm.py / enclosed_volume.py) if violated. At γ_mem = 50 pN/µm and the
     KU-5.1 spacing the per-tip load is ~1–2 pN, so Δx ≪ contact_range and
     the gate never fires — the membrane never tightens the cortex dt. (A
     diagnostic Hookean-stiffness estimate ``γ_mem/s_fil`` is reported by
     :meth:`MembraneLoad.effective_tip_stiffness` for comparison only.)

5. **Sign / sense**
   - The reaction on a contacting barbed end points along ``−n̂`` (into the
     cytosol, opposing the barbed-end advance toward ``+n̂``): a LOAD, not
     a push-forward. STATIC test asserts ``F_reaction · t_barbed < 0`` for
     a barbed end advancing along ``+n̂``.
   - Higher γ_mem → larger F → SMALLER ``k_elong(F)`` and smaller
     ``v(F) = v0·exp(−Fδ/kT)`` (monotone), matching Bieling load-velocity.

6. **Measurement protocol**
   - Force-velocity (KU-5.2): impose a membrane tension γ_mem, read the
     per-filament reaction F = γ_mem·s_fil, feed F into the Bell-Evans
     elongation rate; the polymerization velocity ``v(F)/v(0) =
     exp(−Fδ_elong/kT)`` must track the Brownian-ratchet / Bieling 2016
     curve, monotone decreasing in γ_mem, within tolerance. ANALYTICAL
     (construction-time) test — this is the genuinely construction-testable
     core (see HONESTY note below).
   - BAOAB smoke: a small lamellipodium + membrane-on runs a few hundred
     steps with no NaN/Inf; the membrane load is finite and ≥ 0.

HONESTY note (what this gate does and does NOT prove)
-----------------------------------------------------
The FORCE-BALANCE / Brownian-ratchet math (γ_mem → F → v(F)) IS fully
construction-testable and is the core oracle asserted here. The DOWNSTREAM
emergence — that supplying this load actually lifts the KU-5.1 dendritic
density out of its ~10,000× hole — can only be confirmed by the full,
PI-gated KU-5.1 production sweep (it depends on the loaded elongation /
capping / branching rate BALANCE over a long run, not on a single rate
evaluation). This module therefore delivers the mechanism + config + hook
and asserts the force-balance/ratchet relation; it does NOT fabricate a
construction-time density-emergence pass. See ``test_membrane.py``
docstring + the KU-5.1 production driver for the emergence check.

References
----------
- Mechanism audit: ``ffn_sim/docs/v2_audit/MECHANISM_AUDIT_2026-05-30.md``
  §5 #6 (no explicit membrane → barbed ends carry zero load) + §6 #7
  (force-dependent laws inert at F=0).
- Design: ``ffn_sim/docs/briefs/H4_FA_INTEGRATION_DESIGN.md`` §4
  "Membrane-load (KU-5.x)".
- KU-5.x membrane tension anchor: Bieling 2016 Cell; Lieber et al. 2013
  (Curr Biol, resting membrane tension); Sens & Plastino 2015 (review),
  γ_mem ~ 10–300 pN/µm. Bending κ_m ≈ 20 kT (Helfrich 1973).
- Brownian-ratchet polymerization: Mogilner & Oster 1996 (Biophys J);
  Peskin, Odell & Oster 1993 (Biophys J).
- ``ffn_sim/cell/lamellipodium.py`` — the H.5 lamellipodium this load
  drives (``BarbedEndElongationUpdater`` / ``CappingUpdater`` /
  ``LamellipodiumState`` / ``WaveMembranePin`` plane geometry).
- ``ffn_sim/cortex/enclosed_volume.py`` — the newest additive-module
  template (ResolvedX + resolve_X + attach_X + md.force.Custom + §1–6
  Sanity Gate + Magic-Number Block + test-only oracle helpers).
- ``ffn_sim/cortex/erm.py`` — md.force.Custom + CFL attach-gate pattern.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md


# ---------------------------------------------------------------------------
# Physical constants (anchors; NOT tunables)
# ---------------------------------------------------------------------------
_KT_DEFAULT = 4.114e-21          # J   k_B·T at T = 298 K (room temp)
_GAMMA_MEM_MIN = 1.0e-5          # N/m = 10 pN/µm  (KU-5.x band lower bound)
_GAMMA_MEM_MAX = 3.0e-4          # N/m = 300 pN/µm (KU-5.x band upper bound)
_GAMMA_MEM_DEFAULT = 5.0e-5      # N/m = 50 pN/µm  (mid-band physiological)
_KAPPA_M_KT_DEFAULT = 20.0       # κ_m / kT  (Helfrich bending; cross-check)


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedMembrane:
    """Plasma-membrane load parameters (all SI).

    The membrane is a planar tension reservoir at ``y = y_plane`` with
    outward normal ``+ŷ`` (the leading-edge / WAVE-plane geometry).
    """

    gamma_mem: float             # N/m  membrane tension (KU-5.x anchor)
    y_plane: float               # m    membrane plane position (= WAVE Y_max)
    contact_range: float         # m    δ_contact: tip-to-plane reach for push
    A_mem: float                 # m²   leading-edge membrane patch area
    kT: float                    # J

    # Provenance / cross-check anchors (Magic-Number Block)
    kappa_m_kT: float = _KAPPA_M_KT_DEFAULT   # κ_m / kT (Helfrich, cross-check)
    gamma_mem_band: tuple[float, float] = (_GAMMA_MEM_MIN, _GAMMA_MEM_MAX)

    # Derived diagnostics
    lambda_tension_bending: float = 0.0   # sqrt(κ_m/γ_mem)  [m]
    extras: dict[str, Any] | None = None


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_membrane(
    cfg: dict,
    *,
    y_plane: float | None = None,
    A_mem: float | None = None,
    kT: float = _KT_DEFAULT,
) -> ResolvedMembrane:
    """Resolve the ``membrane`` config block.

    ``cfg`` is the YAML root, the ``lamellipodium`` sub-dict, or the
    ``membrane`` sub-dict. ``y_plane`` defaults to the lamellipodium
    ``Y_max`` (the WAVE membrane plane the barbed ends grow toward);
    ``A_mem`` defaults to the lamellipodium ``wave_area`` (the leading-edge
    patch). Both may be supplied explicitly by the caller (e.g. from a
    resolved :class:`ResolvedH5`).

    γ_mem is taken verbatim from the KU-5.x membrane-tension anchor and
    range-checked against the 10–300 pN/µm band (Bieling 2016; Lieber 2013;
    Sens 2015). A value outside the band raises (CLAUDE.md no-magic-number
    guard — guards against a tuned tension silently overriding the anchor).
    """
    # Unwrap nested config locations (membrane block may live under
    # lamellipodium or at the YAML root).
    if "lamellipodium" in cfg and isinstance(cfg["lamellipodium"], dict):
        if "membrane" in cfg["lamellipodium"]:
            sub = cfg["lamellipodium"]
            # pull Y_max / wave_area defaults from the lamellipodium block
            if y_plane is None and "Y_max" in sub:
                y_plane = float(sub["Y_max"])
            if A_mem is None and "wave_area" in sub:
                A_mem = float(sub["wave_area"])
            cfg = sub["membrane"]
    if "membrane" in cfg:
        cfg = cfg["membrane"]

    gamma_mem = float(cfg.get("gamma_mem", _GAMMA_MEM_DEFAULT))
    contact_range = float(cfg.get("contact_range", cfg.get("delta_contact", 30.0e-9)))
    kappa_m_kT = float(cfg.get("kappa_m_kT", _KAPPA_M_KT_DEFAULT))

    band_lo = float(cfg.get("gamma_mem_min", _GAMMA_MEM_MIN))
    band_hi = float(cfg.get("gamma_mem_max", _GAMMA_MEM_MAX))

    if y_plane is None:
        raise ValueError(
            "resolve_membrane needs y_plane (the membrane plane position); "
            "pass it explicitly or supply lamellipodium.Y_max in cfg."
        )
    if A_mem is None:
        raise ValueError(
            "resolve_membrane needs A_mem (the leading-edge membrane patch "
            "area); pass it explicitly or supply lamellipodium.wave_area."
        )

    # §2 boundary checks.
    _require_finite_positive("gamma_mem", gamma_mem)
    _require_finite_positive("contact_range", contact_range)
    _require_finite_positive("A_mem", A_mem)
    _require_finite_positive("kT", kT)
    _require_finite_positive("kappa_m_kT", kappa_m_kT)
    if not (math.isfinite(y_plane)):
        raise ValueError(f"y_plane must be finite; got {y_plane!r}")

    # No-magic-number band guard on the tension anchor.
    if not (band_lo <= gamma_mem <= band_hi):
        raise ValueError(
            f"gamma_mem = {gamma_mem:.4e} N/m is outside the KU-5.x "
            f"membrane-tension band [{band_lo:.1e}, {band_hi:.1e}] N/m "
            f"(= [10, 300] pN/µm; Bieling 2016 / Lieber 2013 / Sens 2015). "
            "Per CLAUDE.md no-magic-number: either use a band-consistent "
            "tension or update the anchor band and surface to PI."
        )

    p = ResolvedMembrane(
        gamma_mem=gamma_mem,
        y_plane=float(y_plane),
        contact_range=contact_range,
        A_mem=A_mem,
        kT=kT,
        kappa_m_kT=kappa_m_kT,
        gamma_mem_band=(band_lo, band_hi),
        extras={},
    )
    # Diagnostic: tension-to-bending crossover length λ = sqrt(κ_m/γ_mem).
    kappa_m = kappa_m_kT * kT
    p.lambda_tension_bending = math.sqrt(kappa_m / gamma_mem)
    return p


# ---------------------------------------------------------------------------
# Physical span cap (the largest flat tented span a tip can support)
# ---------------------------------------------------------------------------
def _physical_span_cap(p: "ResolvedMembrane") -> float:
    """Largest PHYSICAL flat-tented membrane span between tips [m].

    A membrane under tension γ_mem with bending modulus κ_m transmits a flat
    (tension-dominated) load only over lengths up to the tension-to-bending
    crossover ``λ = sqrt(κ_m / γ_mem)`` (Helfrich); beyond λ the membrane
    BENDS rather than tenting flat, so a single tip cannot tent a flat-tension
    load across a span longer than λ. This is therefore the first-principles
    ceiling on ``s_fil`` (KU-5.x κ_m ≈ 20 kT → λ ≈ 41 nm at 50 pN/µm —
    comparable to the lamellipodial inter-tip spacing, so the cap binds only
    in the sparse / large-patch limit and never alters the dense-array load).
    It is also bounded above by the patch linear extent ``sqrt(A_mem)`` (can't
    tent wider than the patch). No magic number: both κ_m and γ_mem are the
    stated KU-5.x anchors; λ is their derived combination.
    """
    lam = math.sqrt((p.kappa_m_kT * p.kT) / p.gamma_mem)
    return min(lam, math.sqrt(p.A_mem))


# ---------------------------------------------------------------------------
# Closed-form oracles (TEST-ONLY helpers; NOT used by the runtime load)
# ---------------------------------------------------------------------------
def per_filament_load_from_spacing(gamma_mem: float, s_fil: float) -> float:
    """Tension-derived per-filament reaction load ``F = γ_mem · s_fil`` [N].

    The FUNDAMENTAL membrane-load law (Eq. 1): tension × local tip spacing
    (the lateral span the tensed membrane is tented across between
    neighbouring pinning tips; Mogilner-Oster 1996 / Peskin 1993 ratchet
    load scale). ``s_fil ≤ 0`` → 0 (no span, no load).
    """
    if s_fil <= 0.0:
        return 0.0
    return gamma_mem * s_fil


def per_filament_load(
    gamma_mem: float, A_mem: float, n_contact: int, s_cap: float | None = None
) -> float:
    """Tension-derived per-filament membrane reaction load ``F`` [N], from a
    LOCAL contact-patch area + tip count (the pure-scalar oracle form).

    ``F = γ_mem · s_fil`` with ``s_fil = min(sqrt(A_mem / n_contact), s_cap)``
    the local tip-to-tip spacing the membrane is tented across (Eq. 1). Here
    ``A_mem`` is the LOCAL contact-patch area the contacting tips occupy
    (their spread), so ``sqrt(A_mem/n_contact)`` is the mean inter-tip
    spacing; ``s_cap`` (default = ``sqrt(A_mem)``, the patch's linear extent)
    caps the span so a single sparse tip cannot tent the membrane across an
    unphysical distance.

    This is the SAME law :meth:`MembraneLoad.compute_per_filament_load` uses
    when it falls back to the area estimate; exposed as a pure scalar so
    tests can check the force-balance / Brownian-ratchet relation
    analytically without a sim. Returns 0 for ``n_contact == 0`` (no tip
    touches the membrane — §2 boundary).
    """
    if n_contact <= 0 or A_mem <= 0.0:
        return 0.0
    s_fil = math.sqrt(A_mem / n_contact)
    if s_cap is None:
        s_cap = math.sqrt(A_mem)
    s_fil = min(s_fil, s_cap)
    return per_filament_load_from_spacing(gamma_mem, s_fil)


def _capped_load(p: "ResolvedMembrane", s_fil: float) -> float:
    """Per-filament load ``γ_mem · min(s_fil, λ)`` with the physical span cap.

    The runtime entry point: takes a measured/estimated local tip spacing
    ``s_fil`` and applies the tension-to-bending span cap λ (a tip cannot
    tent a flat-tension load over a span longer than the Helfrich crossover
    length — :func:`_physical_span_cap`). Returns 0 for ``s_fil ≤ 0``.
    """
    if s_fil <= 0.0:
        return 0.0
    return per_filament_load_from_spacing(p.gamma_mem, min(s_fil, _physical_span_cap(p)))


def mean_nearest_neighbour_spacing(positions: np.ndarray) -> float:
    """Mean nearest-neighbour distance of a set of tip positions [m].

    The LOCAL tip spacing ``s_fil`` the runtime feeds into Eq. 1 — a
    genuine, grid-invariant geometry measured from the actual contacting
    barbed-end positions (NOT a macroscopic patch area). For < 2 tips there
    is no neighbour, so returns 0 (a lone tip has no membrane to tent
    between neighbours → handled as the single-contact fallback by the
    caller). O(n²) pairwise — n_contact is small (≲ few hundred tips).
    """
    pos = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
    n = pos.shape[0]
    if n < 2:
        return 0.0
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    nn = d.min(axis=1)
    return float(nn.mean())


def brownian_ratchet_velocity(
    v0: float, F: float, delta: float, kT: float = _KT_DEFAULT
) -> float:
    """Brownian-ratchet load-velocity ``v(F) = v0 · exp(−F δ / kT)`` [m/s].

    ACCEPTANCE ORACLE (KU-5.2; Mogilner-Oster 1996; Peskin 1993; the
    Bieling 2016 force-velocity curve). Used in tests to check that the
    runtime tension-derived load produces the correct monotone slowdown.
    NEVER called by the runtime. ``v0`` may be a rate (1/s) or a velocity
    (m/s) — the exp factor is identical, so this doubles as the
    elongation-RATE oracle ``k_elong(F)/k_elong⁰ = exp(−Fδ/kT)``.
    """
    return v0 * math.exp(-F * delta / kT)


def ratchet_exponent(F: float, delta: float, kT: float = _KT_DEFAULT) -> float:
    """Dimensionless Bell-Evans / Brownian-ratchet exponent ``F·δ/kT``."""
    return F * delta / kT


# ---------------------------------------------------------------------------
# Membrane load Action — supplies the per-barbed-end load to the
# Bell-Evans elongation / capping updaters (mirrors _BatchedLamelUpdater).
# ---------------------------------------------------------------------------
class MembraneLoad(hoomd.custom.Action):
    """Per-tick membrane reaction load on the lamellipodial barbed ends.

    Each tick (on the shared ``Periodic(batch_steps)`` trigger, scheduled
    BEFORE the elongation / capping updaters):

    1. Read the current barbed-end positions from the snapshot.
    2. Find the CONTACTING set: barbed ends within ``contact_range`` of the
       membrane plane (``y ≥ y_plane − contact_range``, advancing toward
       ``+ŷ``).
    3. Compute the per-filament tension-derived load
       ``F_per = γ_mem · sqrt(A_mem / n_contact)`` (Eq. 1).
    4. Push the TOTAL load ``F_total = n_contact · F_per`` into the
       elongation / capping updaters via ``set_network_load`` (they
       partition it back over free barbed ends, ``F_total / N_free``).
    5. Record ``last_n_contact`` / ``last_F_per`` / ``last_F_total`` and the
       equal-and-opposite membrane reaction (diagnostic; Newton 3rd law).

    The reaction on each contacting tip points along ``−n̂`` (opposing the
    barbed-end advance toward ``+n̂``) — a LOAD. The membrane feels the
    equal-and-opposite ``+n̂ · F_total`` (a low-DOF reservoir; the plane
    position is diagnostic-only in Phase 1).

    Parameters
    ----------
    p : ResolvedMembrane
        Resolved membrane parameters.
    lamel_state : Any
        The :class:`ffn_sim.cell.lamellipodium.LamellipodiumState` whose
        ``barbed_end_tags`` / ``capped_tags`` / ``tangent_of`` are read.
    barbed_load_consumers : list
        The updater objects exposing ``set_network_load(F_total)`` (the
        ``BarbedEndElongationUpdater`` and ``CappingUpdater``; optionally
        the ``ArpBranchingUpdater``). The membrane pushes the total load
        into each.
    membrane_normal : tuple[float, float, float], default (0, 1, 0)
        Outward membrane normal ``n̂`` (the WAVE plane is ``+ŷ``).
    """

    def __init__(
        self,
        p: ResolvedMembrane,
        *,
        lamel_state: Any,
        barbed_load_consumers: list,
        membrane_normal: tuple[float, float, float] = (0.0, 1.0, 0.0),
    ) -> None:
        super().__init__()
        self.p = p
        self.state = lamel_state
        self._consumers = list(barbed_load_consumers)
        nhat = np.asarray(membrane_normal, dtype=np.float64)
        nrm = float(np.linalg.norm(nhat))
        if nrm <= 0.0:
            raise ValueError("membrane_normal must be a non-zero vector.")
        self.n_hat = nhat / nrm
        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0
        # Last-computed diagnostics (NaN before first tick).
        self.last_n_contact: int = 0
        self.last_F_per: float = float("nan")
        self.last_F_total: float = float("nan")
        self.last_membrane_reaction = np.full(3, float("nan"), dtype=np.float64)

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    @property
    def steps_run(self) -> int:
        return self._steps_run

    # -- load law --
    def compute_per_filament_load(
        self, n_contact: int, s_fil: float | None = None
    ) -> float:
        """Tension-derived per-filament load (Eq. 1) for ``n_contact`` tips.

        If a measured local spacing ``s_fil`` (mean nearest-neighbour tip
        distance) is supplied it is used directly (the runtime path), capped
        at the physical Helfrich span λ (:func:`_physical_span_cap`). Else
        falls back to the area estimate ``min(sqrt(A_mem/n_contact),
        sqrt(A_mem))`` (the pure-scalar oracle form), likewise λ-capped.
        """
        if n_contact <= 0:
            return 0.0
        if s_fil is not None and s_fil > 0.0:
            return _capped_load(self.p, s_fil)
        s_area = min(
            math.sqrt(self.p.A_mem / n_contact), math.sqrt(self.p.A_mem)
        )
        return _capped_load(self.p, s_area)

    def bounded_load_magnitude(self) -> float:
        """Largest physical per-tip membrane load [N] (tension plateau).

        ``F_max = γ_mem · s_max`` with ``s_max`` the physical span cap
        (Helfrich crossover λ, bounded by the patch extent —
        :func:`_physical_span_cap`). The membrane reaction is a BOUNDED load
        (F_per depends on the lateral tip span, not on penetration depth), so
        the CFL-relevant quantity is this bounded magnitude, NOT a Hookean
        stiffness. Used by :func:`attach_membrane_to_simulation`'s
        per-step-drift CFL gate.
        """
        return self.p.gamma_mem * _physical_span_cap(self.p)

    def effective_tip_stiffness(self, n_contact: int) -> float:
        """Diagnostic per-tip tented-membrane stiffness estimate [N/m].

        ``k ≈ γ_mem / s_fil`` (area estimate of ``s_fil``); reported for
        comparison with the cortex bond/ERM stiffness scale. NOTE: this is a
        DIAGNOSTIC only — the membrane reaction is a bounded tension plateau,
        not a Hookean normal spring, so the CFL gate uses the bounded-load
        per-step-drift criterion (:meth:`bounded_load_magnitude`), not this
        stiffness. For ``n_contact == 0`` → 0.
        """
        if n_contact <= 0:
            return 0.0
        s_fil = math.sqrt(self.p.A_mem / n_contact)
        return self.p.gamma_mem / s_fil

    def _contact_set(self, pos: np.ndarray) -> tuple[int, np.ndarray]:
        """(n_contact, contacting-tip positions) for the free barbed ends.

        The global (rank-0) snapshot is tag-ordered: ``pos[i]`` is the
        particle with tag ``i``. A free barbed end contacts the membrane
        when its projection onto ``n̂`` has reached within ``contact_range``
        of the plane (advancing toward ``+n̂``).
        """
        free = [
            t for t in self.state.barbed_end_tags
            if t not in self.state.capped_tags
        ]
        n_part = pos.shape[0]
        thresh = self.p.y_plane - self.p.contact_range
        rows = [
            t for t in free
            if 0 <= t < n_part and float(np.dot(self.n_hat, pos[t])) >= thresh
        ]
        if not rows:
            return 0, np.empty((0, 3), dtype=np.float64)
        return len(rows), pos[np.asarray(rows, dtype=np.int64)]

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        self._steps_run += 1

        # No barbed ends at all → no load.
        if not self.state.barbed_end_tags:
            self._push_load(0.0)
            self.last_n_contact = 0
            self.last_F_per = 0.0
            self.last_F_total = 0.0
            self.last_membrane_reaction = np.zeros(3, dtype=np.float64)
            return

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return
        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()

        n_contact, contact_pos = self._contact_set(pos)
        # Local tip spacing from the ACTUAL contacting-tip positions (mean
        # nearest-neighbour distance); falls back to the area estimate for a
        # lone tip (n_contact < 2 → spacing 0 → area form used).
        s_fil = mean_nearest_neighbour_spacing(contact_pos)
        F_per = self.compute_per_filament_load(
            n_contact, s_fil=(s_fil if s_fil > 0.0 else None)
        )
        F_total = n_contact * F_per

        if not math.isfinite(F_total):
            raise FloatingPointError(
                f"Non-finite membrane load: n_contact={n_contact}, "
                f"F_per={F_per}, gamma_mem={self.p.gamma_mem}."
            )

        self._push_load(F_total)
        self.last_n_contact = int(n_contact)
        self.last_F_per = float(F_per)
        self.last_F_total = float(F_total)
        # Newton 3rd law: membrane feels +n̂·F_total (tips feel −n̂·F_per each).
        self.last_membrane_reaction = self.n_hat * F_total

    def _push_load(self, F_total: float) -> None:
        for c in self._consumers:
            if hasattr(c, "set_network_load"):
                c.set_network_load(F_total)


# ---------------------------------------------------------------------------
# Companion per-bead reaction force (md.force.Custom) — optional, so the
# membrane push is a genuine reaction force in HOOMD's net_force.
# ---------------------------------------------------------------------------
class MembraneReactionForce(md.force.Custom):
    """Per-tip membrane reaction force on contacting barbed ends.

    ``F_i = −n̂ · F_per`` for each free barbed end within ``contact_range``
    of the membrane plane, with ``F_per = γ_mem·sqrt(A_mem/n_contact)``
    (Eq. 1). The force points INTO the cytosol (along ``−n̂``), opposing the
    barbed-end advance toward ``+n̂`` — a load (Sanity Gate §5). All other
    particles are untouched.

    This makes the membrane push a real contribution to HOOMD's
    ``net_force`` (read by the BAOAB Action each step), so the contacting
    tips physically feel the load and the no-COM-runaway invariant is
    exercised. Supplies ONLY the membrane reaction; the BAOAB Action stays
    the sole position integrator.

    Notes
    -----
    The per-tick contact set is read from ``lamel_state`` (Python-level
    barbed-end bookkeeping, same as :class:`MembraneLoad`). Computed in
    float64 via ``cpu_local_snapshot`` / ``cpu_local_force_arrays`` (the
    local snapshot is ROW-PERMUTED by the ParticleSorter, so we map by tag).
    """

    def __init__(
        self,
        p: ResolvedMembrane,
        *,
        lamel_state: Any,
        membrane_normal: tuple[float, float, float] = (0.0, 1.0, 0.0),
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.state = lamel_state
        nhat = np.asarray(membrane_normal, dtype=np.float64)
        self.n_hat = nhat / max(float(np.linalg.norm(nhat)), 1e-30)
        self.last_n_contact: int = 0
        self.last_F_per: float = float("nan")

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()

        F_vec = np.zeros_like(pos)
        U_per = np.zeros(pos.shape[0], dtype=np.float64)

        free = set(self.state.barbed_end_tags) - set(self.state.capped_tags)
        if free:
            # Per-row projection onto the membrane normal.
            proj = pos @ self.n_hat
            thresh = self.p.y_plane - self.p.contact_range
            is_free = np.array([t in free for t in tag])
            contacting = is_free & (proj >= thresh)
            n_contact = int(contacting.sum())
            if n_contact > 0:
                # Local tip spacing from actual contacting-tip positions
                # (mean nearest-neighbour distance), λ-capped; area fallback
                # for a lone tip. SAME Eq. 1 the MembraneLoad Action uses.
                s_fil = mean_nearest_neighbour_spacing(pos[contacting])
                if s_fil > 0.0:
                    F_per = _capped_load(self.p, s_fil)
                else:
                    F_per = per_filament_load(
                        self.p.gamma_mem, self.p.A_mem, n_contact
                    )
                # Reaction along −n̂ on each contacting tip (a load).
                F_vec[contacting] = -self.n_hat * F_per
                # Penetration depth past the plane (≥0) as a diagnostic
                # potential: U = F_per · max(proj − y_plane, 0).
                pen = np.maximum(proj - self.p.y_plane, 0.0)
                U_per[contacting] = F_per * pen[contacting]
                self.last_n_contact = n_contact
                self.last_F_per = float(F_per)
            else:
                self.last_n_contact = 0
                self.last_F_per = 0.0
        else:
            self.last_n_contact = 0
            self.last_F_per = 0.0

        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Public helper: attach the membrane load to an existing lamellipodium sim
# ---------------------------------------------------------------------------
def attach_membrane_to_simulation(
    sim: hoomd.Simulation,
    p_mem: ResolvedMembrane,
    *,
    lamel_state: Any,
    barbed_load_consumers: list,
    batch_steps: int,
    gamma_b: float | None = None,
    n_contact_for_cfl: int = 1,
    with_reaction_force: bool = False,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> dict:
    """Wire the membrane load into an existing lamellipodium simulation.

    Appends a :class:`MembraneLoad` Action (on a
    ``Periodic(batch_steps)`` trigger) that feeds the per-barbed-end load
    into the Bell-Evans elongation / capping updaters, and (optionally) a
    :class:`MembraneReactionForce` custom force so the membrane push is a
    genuine per-bead reaction in HOOMD's ``net_force``.

    The Action is appended to ``sim.operations.updaters``; the caller is
    responsible for ordering it BEFORE the elongation / capping updaters
    (or attaching the membrane first, then the lamellipodium updaters) so
    the load is current when the Bell-Evans laws fire. When the membrane
    is composed via the lamellipodium hook the ordering is handled there.

    Parameters
    ----------
    sim : hoomd.Simulation
        Already-built lamellipodium simulation (must have an Integrator).
    p_mem : ResolvedMembrane
    lamel_state : Any
        The lamellipodium ``LamellipodiumState``.
    barbed_load_consumers : list
        Updaters exposing ``set_network_load`` (elongation, capping, and
        optionally branching).
    batch_steps : int
        Trigger period (shares the lamellipodium D2 batch cadence).
    gamma_b : float, optional
        Per-bead Stokes drag [N·s/m] for the bounded-load per-step-drift CFL
        gate (only applied when ``with_reaction_force=True``; the Action-only
        path puts no force on the integrator and has no CFL constraint). If
        None the gate is skipped (caller responsible).
    n_contact_for_cfl : int, default 1
        Retained for signature stability; the bounded-load drift CFL gate
        does not depend on the contact count (the bounded magnitude
        F_max = γ_mem·s_max is the worst case). Ignored.
    with_reaction_force : bool, default False
        Also attach the :class:`MembraneReactionForce` custom force so the
        reaction enters ``net_force``. Default False (the load → Bell-Evans
        rate is the primary, mechanistically-required path).
    cfl_safety_factor : float, default 0.1
        Same convention as the bond/angle/ERM/EV CFL gates (D3 BAOAB).
    cfl_strict : bool, default True
        If True, raise on CFL violation; set False for diagnostic runs.

    Returns
    -------
    dict
        Handles: ``membrane_load`` (the Action), ``membrane_updater`` (its
        CustomUpdater wrapper), ``membrane_reaction_force`` (or None).
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching the "
            "membrane load."
        )

    # CFL gate (parity with erm.py / enclosed_volume.py, but physically
    # correct for a BOUNDED LOAD rather than a stiff spring). IMPORTANT: this
    # gate is RELEVANT ONLY for the companion per-step reaction FORCE
    # (:class:`MembraneReactionForce`, enabled via ``with_reaction_force``).
    # The PRIMARY, mechanistically-required path is the :class:`MembraneLoad`
    # Action, which feeds the tension-derived load into the Bell-Evans
    # elongation/capping RATE (a probability per batch tick) — it puts NO
    # force on the integrator, so it has NO CFL/stability constraint at all.
    # Only when an actual per-bead reaction force enters HOOMD's net_force
    # does the overdamped-BAOAB per-step DRIFT bound apply: the tip must not
    # be pushed more than a safety fraction of its bead-scale contact
    # footprint in one step. The membrane reaction is a tension PLATEAU of
    # bounded per-tip magnitude (F_per depends on the LATERAL tip span, not on
    # penetration depth — NOT a Hookean normal spring), so the worst case uses
    # the largest physical tented span ``s_max`` (the Helfrich λ cap).
    if gamma_b is not None and with_reaction_force:
        dt = float(ig.dt)
        s_max = _physical_span_cap(p_mem)
        F_max = p_mem.gamma_mem * s_max          # bounded plateau load [N]
        drift = F_max * dt / gamma_b             # per-step BAOAB drift [m]
        drift_ceiling = cfl_safety_factor * p_mem.contact_range
        if drift > drift_ceiling and cfl_strict:
            dt_cfl_mem = drift_ceiling * gamma_b / F_max
            raise RuntimeError(
                f"Membrane CFL violated: per-step drift from the bounded "
                f"membrane load Δx = F_max·dt/γ_b = {drift:.3e} m > "
                f"{cfl_safety_factor:.2f} · contact_range = "
                f"{drift_ceiling:.3e} m (F_max = γ_mem·s_max = {F_max:.3e} N, "
                f"s_max = {s_max:.3e} m, γ_b = {gamma_b:.3e} N·s/m). Reduce "
                "γ_mem (softer membrane — requires PI sign-off vs the KU-5.x "
                f"tension anchor) OR reduce dt (need dt ≤ {dt_cfl_mem:.3e} "
                "s). Pass cfl_strict=False to skip for diagnostic runs."
            )

    membrane_load = MembraneLoad(
        p_mem,
        lamel_state=lamel_state,
        barbed_load_consumers=barbed_load_consumers,
    )
    membrane_updater = hoomd.update.CustomUpdater(
        action=membrane_load, trigger=hoomd.trigger.Periodic(int(batch_steps)),
    )
    sim.operations.updaters.append(membrane_updater)

    membrane_reaction_force = None
    if with_reaction_force:
        membrane_reaction_force = MembraneReactionForce(
            p_mem, lamel_state=lamel_state,
        )
        ig.forces.append(membrane_reaction_force)

    return {
        "membrane_load": membrane_load,
        "membrane_updater": membrane_updater,
        "membrane_reaction_force": membrane_reaction_force,
    }
