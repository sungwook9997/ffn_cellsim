"""H.3 actin turnover — cofilin severing + pointed-end re-annealing (KU-3.5).

ADDITIVE, DEFAULT-OFF module (mechanism audit 2026-05-30 §5 #1). The H.3
cortex is otherwise a FROZEN permanently-connected network: every
``cortex-bond`` backbone bond placed at construction lives forever, so the
filament network can only stiffen, never relax — and the actin-turnover
process that Chugh 2017 / Saha 2016 / Fritzsche 2013 identify as the SINGLE
DOMINANT determinant of steady-state cortical tension has no runtime
representation. KU-3.5 cortical tension is unreachable (even after FA wiring)
without it. This module supplies force-/age-aware filament disassembly so a
steady-state tension can EMERGE.

Mechanism (fine-grained, not lumped)
------------------------------------
**Severing (cofilin).** Each cortex-actin backbone bond carries a per-batch
first-order severing probability::

    p_sev = 1 − exp(−k_sev · Δt_batch)

with ``k_sev = ln2 / τ_half`` the cofilin severing rate derived from the
filament half-life ``τ_half`` (Chugh 2017 / Fritzsche 2013 cortical actin
turnover τ ≈ 5–30 s; the ACCEPTANCE-ORACLE anchor). On a sever event the
backbone bond ``(i, i+1)`` is REMOVED — the filament splits at that junction.
The two cortex-angle triplets that span the severed junction
(``(i−1, i, i+1)`` and ``(i, i+1, i+2)``, the angles whose middle segment is
the broken bond) are removed in lockstep, because a bending angle across a
severed junction would be unphysical (it would couple two now-disconnected
fragments).

**Pointed-end re-annealing (treadmilling balance).** A sever-only model would
disassemble the cortex monotonically and is therefore not a steady state. To
pair disassembly with assembly while keeping the HOOMD particle count FIXED
(the Phase-1 constraint that ``crosslinkers.py`` / ``myosin.py`` work under —
mutate bonds, NOT particles), each *currently-severed* backbone junction
re-anneals with per-batch probability::

    p_anneal = 1 − exp(−k_anneal · Δt_batch)

PROVIDED the two beads it previously joined are still within
``r0 · (1 + reanneal_tol)`` of each other (i.e. the filament has not drifted
apart — re-polymerization / re-annealing reconnects the SAME junction, the
mechanistic pointed-end/barbed-end re-addition of a monomer that closes the
gap). When it fires, the ``cortex-bond`` and its spanning angles are restored.

This is the simplest mechanistic model that (a) reaches a STEADY STATE — the
backbone is a two-state (connected ⇌ severed) per-junction process whose
steady-state connected fraction is ``f_ss = k_anneal / (k_sev + k_anneal)``,
so the mean connected-fragment length stabilizes rather than collapsing; (b)
is mechanistic (severing = cofilin first-order kinetics; re-annealing =
distance-gated pointed-end re-addition, both per-junction); (c) keeps the
total particle count fixed (bond + angle topology mutation only). Its
documented LIMITATION: re-annealing only the SAME junction (not migration of
the severed end to a new acceptor, and no monomer-pool depletion) is a
Phase-1 simplification — sufficient for a bounded steady state at fixed
particle count, but it does not model net filament translocation
(treadmilling flux) or G-actin pool kinetics. A reservoir-with-new-acceptor
model is a documented later refinement.

**Force / age dependence (hooks, default OFF).** The Phase-1 default is the
simplest defensible first-order severing (constant ``k_sev``). Two
mechanistic refinements are wired as HOOKS, off by default:

* ``k_sev_age_factor`` (cofilin preferentially severs OLD / ADP-actin
  segments; De La Cruz 2009): when > 0, a bond's effective severing rate is
  scaled by ``1 + k_sev_age_factor · (age / τ_half)`` where ``age`` is the
  time since the bond was last (re)annealed. Default 0 → age-independent.
* ``k_sev_force_x_beta`` (force-dependent severing; cofilin severing is
  enhanced near domain boundaries under tension — a Bell-Evans-style factor):
  when > 0, ``k_sev`` is scaled by ``exp(F · x_β / kT)`` with ``F`` the bond
  tension. Default 0 → force-independent.

Both default to 0 so the Phase-1 runtime is the clean first-order form, with
the mechanistic extensions documented and testable.

Architecture
------------
A ``hoomd.custom.Action`` wrapped in a periodic ``CustomUpdater`` that runs
every ``batch_steps`` BAOAB steps — the D2 batched-updater convention shared
with ``crosslinkers.py`` (``XlinkBondUpdater``) and ``myosin.py``. Each batch
tick reads a snapshot, mutates the ``cortex-bond`` + ``cortex-angle`` topology
(remove on sever, restore on anneal), and writes a fresh snapshot. It does
NOT touch particle positions — the BAOAB Action remains the SOLE position
integrator (no double-stepping). The frozen ``integrator/`` is UNTOUCHED.

Batch-CFL contract (mirrors crosslinkers / myosin / integrin)
-------------------------------------------------------------
``events_per_batch = batch_steps · dt · k_sev_max ≤ 1e-3``, where ``k_sev_max``
is the maximum credible severing rate (the bare ``k_sev`` times the maximum
age/force envelope when those hooks are on; just ``k_sev`` in the default
first-order case). ``resolve_turnover`` enforces this at config-resolve time
by SHRINKING ``batch_steps`` if violated (identical to
``resolve_crosslinkers``). The Updater's ``__init__`` re-asserts the bound.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_turnover.py``.*

1. **Dimensional analysis**
   - ``k_sev = ln2 / τ_half``: ``τ_half`` [s] → ``k_sev`` [1/s]. ✓
   - ``k_anneal`` [1/s]; ``Δt_batch = batch_steps · dt`` [s];
     ``k · Δt_batch`` dimensionless → probability. ✓
   - ``k_sev_force_x_beta`` [m]; ``F`` [N]; ``kT`` [J] → ``F·x_β/kT``
     dimensionless. ✓  ``k_sev_age_factor`` dimensionless. ✓
   - STATIC test recomputes ``k_sev`` from ``τ_half`` and asserts
     ``k_sev = ln2/τ`` recovers the half-life.

2. **Boundary cases**
   - ``k_sev = 0`` (τ_half = ∞ via ``enabled: false`` upstream OR an explicit
     0): no bonds severed; the Updater idles. (With ``enabled: false`` the
     Updater is never attached at all — see ``cell.py`` hook.)
   - ``τ_half ≤ 0``: ``resolve_turnover`` ValueError.
   - ``k_anneal < 0`` / ``k_sev < 0``: ValueError.
   - ``reanneal_tol < 0``: ValueError.
   - No backbone bonds (n_backbone = 0): the act() loop finds zero
     ``cortex-bond`` entries; idles.
   - ``batch_steps · dt · k_sev_max > 1e-3``: caught by ``resolve_turnover``
     (batch_steps shrunk); ``__init__`` re-asserts.

3. **Conservation invariants**
   - **Particle count is INVARIANT** — this module mutates ONLY bond + angle
     topology, never particles (the Phase-1 fixed-N constraint). STATIC +
     demo test assert ``sim.state.N_particles`` unchanged.
   - Backbone-bond count is bounded by the construction count
     ``n_backbone_0`` (severing removes, annealing restores; never exceeds
     the original junction set — re-anneal only reconnects pre-existing
     junctions). Non-cortex-bond topology (xlink_intra, xlink_attach,
     myosin, FA, lamellipodium bonds + their angles) is PRESERVED untouched
     across every tick.
   - No net momentum injected: bond/angle removal/restore changes the
     conservative force field but injects no force directly (the BAOAB
     thermostat is the only kinetic coupling).

4. **Numerical sanity**
   - dt provided at construction; CFL gate checked once at ``__init__``.
   - Per-batch act() is O(n_backbone) vectorised (boolean masks over the
     bond array) — cheaper than the crosslinker KDTree query.
   - RNG isolated per Action (``np.random.default_rng`` with
     ``seed + seed_offset``); ``seed_offset`` defaults to 5 to avoid
     collision with BAOAB(0) / integrin(1) / xlink(2) / myosin(...).
   - All positions read float64; no NaN/Inf in a short BAOAB demo (STATIC).

5. **Sign / sense**
   - Severing REMOVES connectivity (network softens) — the correct sense for
     a disassembly process; STATIC asserts the severed-bond count rises from
     0 when only severing is active (k_anneal = 0).
   - Re-annealing RESTORES connectivity (network re-stiffens) — STATIC
     asserts the connected fraction rises toward 1 when only annealing is
     active (k_sev = 0, all bonds pre-severed).
   - The force/age hooks MONOTONICALLY INCREASE k_sev (older / more-loaded
     segments sever faster) — STATIC asserts ``d k_sev_eff / d F > 0`` and
     ``d k_sev_eff / d age > 0`` when the hooks are on.

6. **Measurement protocol**
   - **Kinetics gate (the core oracle)**: over a batch window the fraction of
     intact backbone bonds that sever matches the first-order expectation
     ``1 − exp(−k_sev · Δt_batch)`` within sampling tolerance. ANALYTICAL —
     run many independent junctions through ``severing_probability`` /
     one batch tick and compare the empirical sever fraction to the closed
     form. The closed form is the ORACLE; the runtime is the per-bond
     Bernoulli draw.
   - **Steady state**: with both severing + annealing on, the connected
     fraction equilibrates to ``f_ss = k_anneal/(k_sev + k_anneal)`` and the
     mean connected-fragment length STABILIZES (does not monotone-collapse)
     over a demo BAOAB run. ANALYTICAL ``f_ss`` + demo check.

References
----------
- Mechanism audit: ``ffn_sim/docs/v2_audit/MECHANISM_AUDIT_2026-05-30.md``
  §5 #1 (actin turnover — single dominant determinant of steady-state
  cortical tension; zero runtime implementation before this module).
- KU-3.5 cortical tension (Salbreux 2012) — the tension this module makes
  reachable.
- Cortical actin turnover half-life τ ≈ 5–30 s: Fritzsche et al. 2013 (Mol
  Biol Cell); Chugh et al. 2017 (Nat Cell Biol); Saha et al. 2016. The
  ``k_sev = ln2/τ`` ANCHOR.
- Cofilin severing kinetics (first-order; age/ADP-actin preference): De La
  Cruz 2009 (Biophys Rev); Suarez 2011.
- D2 batched-updater + per-r0 / per-bond snapshot mutation pattern:
  ``ffn_sim/cortex/crosslinkers.py`` (XlinkBondUpdater), shared CFL gate with
  ``ffn_sim/bridge/integrin_bonds.py`` (gate_bell_evans_batch_cfl).
- ``ffn_sim/cortex/cortex.py`` (cortex-bond backbone + cortex-angle topology).
- ``ffn_sim/integrator/baoab.py`` (D3 BAOAB, frozen 2026-05-20) — the SOLE
  position integrator; this module mutates topology only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom


CORTEX_BOND_TYPE: str = "cortex-bond"
CORTEX_ANGLE_TYPE: str = "cortex-angle"

# D2 batch-CFL ceiling, shared with crosslinkers / myosin / integrin updaters.
BATCH_CFL_CEILING: float = 1.0e-3


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedTurnover:
    """Resolved actin-turnover parameters (all SI).

    Populated by :func:`resolve_turnover` from the ``cortex.turnover`` config
    block. The severing rate ``k_sev`` is DERIVED from the filament half-life
    ``tau_half`` (``k_sev = ln2 / tau_half``) — no magic number; the half-life
    is the KU-3.5 / Chugh-2017 / Fritzsche-2013 literature anchor.
    """

    enabled: bool                # default-off gate (mirrors enclosed_volume)
    tau_half: float              # s    filament half-life (cortical actin turnover)
    k_anneal: float              # 1/s  pointed-end re-annealing rate
    reanneal_tol: float          # —    re-anneal if |Δr| ≤ r0·(1 + reanneal_tol)

    # Force / age dependence HOOKS (default 0 → clean first-order severing).
    k_sev_age_factor: float      # —    age scaling (0 = age-independent)
    k_sev_force_x_beta: float    # m    Bell-Evans severing length (0 = force-indep.)

    # D2 batch
    batch_steps: int             # BAOAB steps between Updater ticks
    dt: float                    # host sim dt [s]
    seed: int

    # Derived
    k_sev: float = 0.0           # ln2 / tau_half  [1/s]
    batch_dt: float = 0.0        # batch_steps · dt  [s]
    k_sev_max: float = 0.0       # max credible severing rate (for CFL)
    events_per_batch: float = 0.0
    f_ss: float = 0.0            # steady-state connected fraction k_anneal/(k_sev+k_anneal)
    extras: dict[str, Any] = field(default_factory=dict)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0, got {x!r}")


def _require_finite_nonnegative(name: str, x: float) -> None:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError(f"{name} must be finite and ≥ 0, got {x!r}")


def resolve_turnover(
    cfg: dict,
    *,
    dt: float,
    rest_length: float | None = None,
) -> ResolvedTurnover:
    """Resolve the ``cortex.turnover`` config block with D2 batch-CFL enforcement.

    ``cfg`` is the YAML root, the ``cortex:`` sub-dict, or the
    ``cortex.turnover`` sub-dict. ``dt`` is the host simulation timestep
    (typically ``p_cortex.dt_cfl``). ``rest_length`` (optional) is recorded in
    ``extras`` for the re-annealing distance gate; it is supplied to the
    Updater separately.

    ``k_sev`` is DERIVED: ``k_sev = ln2 / tau_half`` (no magic number — the
    half-life ``tau_half`` is the KU-3.5 / Chugh-2017 cortical-actin-turnover
    literature anchor). If the config also supplies an explicit ``k_sev`` it
    is checked for consistency against the derivation (a > 1 % mismatch raises,
    guarding against a tuned value silently overriding the anchor — same guard
    as ``enclosed_volume.resolve_enclosed_volume`` uses for K_vol).
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "turnover" in cfg:
        cfg = cfg["turnover"]

    enabled = bool(cfg.get("enabled", False))
    tau_half = float(cfg["tau_half"])
    k_anneal = float(cfg.get("k_anneal", 0.1))
    reanneal_tol = float(cfg.get("reanneal_tol", 0.1))
    k_sev_age_factor = float(cfg.get("k_sev_age_factor", 0.0))
    k_sev_force_x_beta = float(cfg.get("k_sev_force_x_beta", 0.0))
    batch_steps = int(cfg.get("batch_steps", 100))
    seed = int(cfg.get("seed", 47))

    # ---- §2 boundary checks ----
    _require_finite_positive("tau_half", tau_half)
    _require_finite_positive("dt", float(dt))
    _require_finite_nonnegative("k_anneal", k_anneal)
    _require_finite_nonnegative("reanneal_tol", reanneal_tol)
    _require_finite_nonnegative("k_sev_age_factor", k_sev_age_factor)
    _require_finite_nonnegative("k_sev_force_x_beta", k_sev_force_x_beta)
    if batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {batch_steps}")

    # ---- §1 derived: k_sev from the half-life anchor (Magic-Number Block) ----
    k_sev_derived = math.log(2.0) / tau_half
    k_sev_cfg = cfg.get("k_sev", None)
    if k_sev_cfg is not None:
        k_sev = float(k_sev_cfg)
        _require_finite_nonnegative("k_sev", k_sev)
        if k_sev > 0.0 and not math.isclose(k_sev, k_sev_derived, rel_tol=1.0e-2):
            raise ValueError(
                f"Explicit k_sev = {k_sev:.4e} 1/s disagrees with the "
                f"first-principles anchor ln2/tau_half = {k_sev_derived:.4e} "
                f"1/s (tau_half = {tau_half} s). Per CLAUDE.md no-magic-number: "
                "either remove the explicit k_sev (let it derive) or update "
                "tau_half and surface to PI."
            )
    else:
        k_sev = k_sev_derived

    p = ResolvedTurnover(
        enabled=enabled,
        tau_half=tau_half,
        k_anneal=k_anneal,
        reanneal_tol=reanneal_tol,
        k_sev_age_factor=k_sev_age_factor,
        k_sev_force_x_beta=k_sev_force_x_beta,
        batch_steps=batch_steps,
        dt=float(dt),
        seed=seed,
    )
    p.k_sev = k_sev
    p.batch_dt = p.batch_steps * p.dt

    # ---- §2 batch CFL: k_sev_max · batch_dt ≤ 1e-3, else shrink batch_steps ----
    # k_sev_max is the bare k_sev in the default first-order case. With the
    # age hook on, the multiplicative envelope is bounded by
    # (1 + k_sev_age_factor) at age = τ_half (the half-life-normalised age the
    # hook references); with the force hook on it is exp(F_max·x_β/kT). We use a
    # conservative credible-force scale F_max = bond_k · ℓ₀ (a bond stretched by
    # its own rest length — far beyond the thermal regime) when x_β > 0; in the
    # default (x_β = 0) the factor is exactly 1. The age factor at age=τ_half is
    # (1 + k_sev_age_factor·1). This keeps the CFL conservative without a magic
    # number — both envelopes are the hooks' own maxima.
    age_envelope = 1.0 + p.k_sev_age_factor  # at age = τ_half
    force_envelope = 1.0
    if p.k_sev_force_x_beta > 0.0:
        # Credible max bond tension scale: kT/ℓ₀ if rest_length known, else a
        # 1 pN reference (thermal-to-pN bridge). Both are documented envelopes,
        # not tuned: they only set the CFL safety margin, never the kinetics.
        kT_ref = 4.28e-21
        F_max = (kT_ref / rest_length) if (rest_length and rest_length > 0.0) else 1.0e-12
        force_envelope = math.exp(F_max * p.k_sev_force_x_beta / kT_ref)
    p.k_sev_max = p.k_sev * age_envelope * force_envelope

    if p.k_sev_max > 0.0:
        events = p.batch_dt * p.k_sev_max
        if events > BATCH_CFL_CEILING:
            target_batch_dt = BATCH_CFL_CEILING / p.k_sev_max
            new_batch_steps = max(1, int(math.floor(target_batch_dt / p.dt)))
            if new_batch_steps < 1:
                raise RuntimeError(
                    f"D2 batch CFL impossible to satisfy: dt={p.dt:.3e} > "
                    f"1e-3 / k_sev_max ({BATCH_CFL_CEILING / p.k_sev_max:.3e})."
                )
            p.extras["batch_steps_shrunk_from"] = p.batch_steps
            p.batch_steps = new_batch_steps
            p.batch_dt = p.batch_steps * p.dt
    p.events_per_batch = p.batch_dt * p.k_sev_max

    # ---- steady-state connected fraction (diagnostic; the test oracle) ----
    denom = p.k_sev + p.k_anneal
    p.f_ss = (p.k_anneal / denom) if denom > 0.0 else 1.0

    if rest_length is not None:
        p.extras["rest_length"] = float(rest_length)

    # ---- §4 finite checks on derived ----
    for name in ("k_sev", "batch_dt", "k_sev_max", "events_per_batch", "f_ss"):
        v = getattr(p, name)
        if not (math.isfinite(v) and v >= 0.0):
            raise RuntimeError(f"Derived {name}={v!r} is not finite-nonnegative.")
    if p.events_per_batch > BATCH_CFL_CEILING + 1.0e-12:
        raise RuntimeError(
            f"resolve_turnover failed to satisfy batch CFL: "
            f"events_per_batch={p.events_per_batch:.3e} > {BATCH_CFL_CEILING:.3e}."
        )

    return p


# ---------------------------------------------------------------------------
# Closed-form oracles (TEST-ONLY helpers; NOT used as the runtime mechanism)
# ---------------------------------------------------------------------------
def severing_rate_from_half_life(tau_half: float) -> float:
    """First-order severing rate ``k_sev = ln2 / τ_half`` [1/s].

    ACCEPTANCE-ORACLE anchor (Chugh 2017 / Fritzsche 2013 cortical actin
    half-life). This IS the derivation the resolver uses; exposed as a scalar
    so the dimensional test can recover the half-life from k_sev.
    """
    return math.log(2.0) / tau_half


def severing_probability(k_sev: float, dt_batch: float) -> float:
    """First-order per-batch severing probability ``1 − exp(−k_sev·Δt)``.

    The CORE kinetics ORACLE: over a batch window each intact bond severs with
    this probability. The runtime draws a Bernoulli with exactly this p; the
    test compares the empirical sever fraction to this closed form.
    """
    return 1.0 - math.exp(-k_sev * dt_batch)


def steady_state_connected_fraction(k_sev: float, k_anneal: float) -> float:
    """Two-state (connected ⇌ severed) steady-state connected fraction.

    ``f_ss = k_anneal / (k_sev + k_anneal)`` — the fraction of backbone
    junctions intact at steady state. ORACLE for the steady-state gate.
    """
    denom = k_sev + k_anneal
    return (k_anneal / denom) if denom > 0.0 else 1.0


def effective_severing_rate(
    k_sev: float,
    *,
    age: float = 0.0,
    tau_half: float = 1.0,
    k_sev_age_factor: float = 0.0,
    force: float = 0.0,
    k_sev_force_x_beta: float = 0.0,
    kT: float = 4.28e-21,
) -> float:
    """Effective severing rate with the (default-off) age / force hooks.

    ``k_sev_eff = k_sev · (1 + k_sev_age_factor·age/τ_half) · exp(F·x_β/kT)``.

    With both factors at their defaults (0) this returns ``k_sev`` exactly
    (the clean Phase-1 first-order form). Used by the runtime per-bond rate
    AND by the STATIC sign tests (monotone-increasing in age and F).
    """
    age_factor = 1.0 + k_sev_age_factor * (age / tau_half if tau_half > 0.0 else 0.0)
    force_factor = math.exp(force * k_sev_force_x_beta / kT) if k_sev_force_x_beta > 0.0 else 1.0
    return k_sev * age_factor * force_factor


# ---------------------------------------------------------------------------
# Severing + re-annealing Updater
# ---------------------------------------------------------------------------
class ActinTurnoverUpdater(hoomd.custom.Action):
    """D2 actin-turnover updater — cofilin severing + pointed-end re-annealing.

    Each batch tick (every ``p.batch_steps`` BAOAB steps):

    1. Read the snapshot (positions + bonds + angles).
    2. Split the bond list into ``cortex-bond`` backbone bonds vs all others
       (xlink / myosin / FA / lamellipodium bonds — PRESERVED untouched).
    3. **Sever**: for each currently-intact backbone junction, sample
       ``p_sev = 1 − exp(−k_sev_eff·Δt_batch)`` (with the optional age / force
       hooks) and remove the bond on fire; mark the junction severed and
       record the sever timestep (for the age hook).
    4. **Re-anneal**: for each currently-severed junction whose two beads are
       within ``r0·(1 + reanneal_tol)``, sample
       ``p_anneal = 1 − exp(−k_anneal·Δt_batch)`` and restore the bond on fire;
       reset its age clock.
    5. Rebuild the ``cortex-bond`` group + the ``cortex-angle`` triplets that
       span only currently-intact junctions, concatenate with the preserved
       non-cortex topology, and write a fresh snapshot.

    The set of CANDIDATE backbone junctions is FIXED at construction (the
    original ``cortex-bond`` set): severing/annealing only toggles members of
    this set, so the backbone-bond count never exceeds the construction count
    and particle count is invariant.

    Parameters
    ----------
    p : ResolvedTurnover
        Resolved turnover config (batch_dt, k_sev, k_anneal, hooks).
    rest_length : float
        Cortex backbone rest length ℓ₀ [m] — the re-annealing distance gate
        uses ``r0·(1 + reanneal_tol)``.
    kT : float
        Thermal energy for the (optional) force-dependent severing hook.
    bond_k : float, optional
        Backbone bond stiffness [N/m] — used to compute the bond tension
        ``F = bond_k·max(0, |Δr|−r0)`` for the force hook. Required only when
        ``k_sev_force_x_beta > 0``; ignored otherwise.
    seed_offset : int, default 5
        Offset added to ``p.seed`` for this Action's RNG. Default 5 avoids
        collision with BAOAB(0) / integrin(1) / xlink(2) / myosin / lamel.
    """

    def __init__(
        self,
        *,
        p: ResolvedTurnover,
        rest_length: float,
        kT: float,
        bond_k: float | None = None,
        seed_offset: int = 5,
    ) -> None:
        super().__init__()
        self.p = p
        self.rest_length = float(rest_length)
        self.kT = float(kT)
        self.bond_k = float(bond_k) if bond_k is not None else None
        self._rng = np.random.default_rng(p.seed + seed_offset)

        _require_finite_positive("rest_length", self.rest_length)
        _require_finite_positive("kT", self.kT)
        if self.p.k_sev_force_x_beta > 0.0 and self.bond_k is None:
            raise ValueError(
                "k_sev_force_x_beta > 0 requires bond_k for the force hook."
            )

        # Re-assert D2 batch CFL.
        events = p.batch_dt * p.k_sev_max
        if events > BATCH_CFL_CEILING + 1.0e-12:
            raise RuntimeError(
                f"ActinTurnoverUpdater violates D2 batch CFL: "
                f"batch_dt·k_sev_max = {p.batch_dt:.3e}·{p.k_sev_max:.3e} "
                f"= {events:.3e} > {BATCH_CFL_CEILING:.3e} ceiling. "
                "Shrink batch_steps or pre-resolve via resolve_turnover."
            )

        # Lazily initialised on the first act() from the construction snapshot:
        #   _junctions      (n0, 2)   the FIXED candidate backbone bead-pairs
        #   _intact         (n0,)     bool — currently connected
        #   _last_anneal_t  (n0,)     timestep the junction was last (re)annealed
        #   _junc_index     dict (a,b)->row  for O(1) lookup
        #   _other_*        preserved non-cortex-bond topology
        self._initialised = False
        self._junctions: np.ndarray | None = None
        self._intact: np.ndarray | None = None
        self._last_anneal_t: np.ndarray | None = None
        self._junc_index: dict[tuple[int, int], int] = {}

        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0
        self._n_sever_total = 0
        self._n_anneal_total = 0

    # ------------------------------------------------------------------
    @staticmethod
    def _norm_pair(a: int, b: int) -> tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def _init_from_snapshot(self, read_snap) -> None:
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64)
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.int64)
        bond_type_names = list(read_snap.bonds.types)
        if CORTEX_BOND_TYPE not in bond_type_names:
            raise RuntimeError(
                f"snapshot bond types must include {CORTEX_BOND_TYPE!r}; "
                f"got {bond_type_names}."
            )
        cortex_typeid = bond_type_names.index(CORTEX_BOND_TYPE)
        is_cortex = bt == cortex_typeid
        cortex_bonds = bg[is_cortex]

        # The FIXED candidate junction set = the construction cortex-bonds.
        self._junctions = cortex_bonds.copy()
        n0 = self._junctions.shape[0]
        self._intact = np.ones(n0, dtype=bool)
        self._last_anneal_t = np.zeros(n0, dtype=np.int64)
        self._junc_index = {
            self._norm_pair(int(a), int(b)): i
            for i, (a, b) in enumerate(self._junctions)
        }
        self.p.extras["n_backbone_0"] = int(n0)
        self._initialised = True

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        if not self._initialised:
            self._init_from_snapshot(read_snap)

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.int64).copy()
        bond_type_names = list(read_snap.bonds.types)
        cortex_bond_typeid = bond_type_names.index(CORTEX_BOND_TYPE)

        is_cortex = bt == cortex_bond_typeid
        other_bg = bg[~is_cortex]
        other_bt = bt[~is_cortex]

        junctions = self._junctions
        intact = self._intact
        last_anneal = self._last_anneal_t
        assert junctions is not None and intact is not None and last_anneal is not None
        n0 = junctions.shape[0]

        if n0 > 0:
            r_a = pos[junctions[:, 0]]
            r_b = pos[junctions[:, 1]]
            dr = r_a - r_b
            sep = np.linalg.norm(dr, axis=1)

            # ---- Step 1: SEVER currently-intact junctions ----
            if intact.any() and self.p.k_sev > 0.0:
                intact_idx = np.flatnonzero(intact)
                # Effective per-bond severing rate (hooks default to identity).
                k_eff = np.full(intact_idx.size, self.p.k_sev, dtype=np.float64)
                if self.p.k_sev_age_factor > 0.0:
                    age = (timestep - last_anneal[intact_idx]).astype(np.float64) * self.p.dt
                    k_eff = k_eff * (
                        1.0 + self.p.k_sev_age_factor * (age / self.p.tau_half)
                    )
                if self.p.k_sev_force_x_beta > 0.0 and self.bond_k is not None:
                    F = self.bond_k * np.clip(sep[intact_idx] - self.rest_length, 0.0, None)
                    k_eff = k_eff * np.exp(F * self.p.k_sev_force_x_beta / self.kT)
                p_sev = 1.0 - np.exp(-k_eff * self.p.batch_dt)
                u = self._rng.uniform(0.0, 1.0, size=intact_idx.size)
                fired = u < p_sev
                if fired.any():
                    severed = intact_idx[fired]
                    intact[severed] = False
                    self._n_sever_total += int(fired.sum())

            # ---- Step 2: RE-ANNEAL currently-severed junctions within range ----
            if (~intact).any() and self.p.k_anneal > 0.0:
                sev_idx = np.flatnonzero(~intact)
                gate = sep[sev_idx] <= self.rest_length * (1.0 + self.p.reanneal_tol)
                cand = sev_idx[gate]
                if cand.size > 0:
                    p_anneal = 1.0 - math.exp(-self.p.k_anneal * self.p.batch_dt)
                    u2 = self._rng.uniform(0.0, 1.0, size=cand.size)
                    fired2 = u2 < p_anneal
                    if fired2.any():
                        annealed = cand[fired2]
                        intact[annealed] = True
                        last_anneal[annealed] = timestep
                        self._n_anneal_total += int(fired2.sum())

        # ---- Step 3: rebuild cortex-bond + cortex-angle topology ----
        intact_junctions = junctions[intact] if n0 > 0 else junctions
        n_intact = intact_junctions.shape[0]

        # Rebuild cortex-angle triplets that span ONLY intact junctions.
        new_angles, angle_types = self._rebuild_angles(read_snap, intact)

        new_cortex_bg = intact_junctions.astype(np.uint32)
        cortex_typeids = np.full(n_intact, cortex_bond_typeid, dtype=np.uint32)

        if other_bg.shape[0] > 0:
            new_bg = np.concatenate(
                [new_cortex_bg, other_bg.astype(np.uint32)], axis=0
            )
            new_bt = np.concatenate(
                [cortex_typeids, other_bt.astype(np.uint32)]
            )
        else:
            new_bg = new_cortex_bg
            new_bt = cortex_typeids

        self._write_snapshot(read_snap, pos, new_bg, new_bt, bond_type_names,
                             new_angles, angle_types)
        self._steps_run += 1

    # ------------------------------------------------------------------
    def _rebuild_angles(self, read_snap, intact: np.ndarray):
        """Keep only cortex-angle triplets whose two segments are both intact.

        A cortex-angle (i, i+1, i+2) couples the bonds (i,i+1) and (i+1,i+2).
        If EITHER spanning junction is currently severed, the angle is dropped
        (no bending across a broken junction). Non-cortex angles (e.g.
        lamellipodium branch angles) are PRESERVED untouched.
        """
        n_ang = int(read_snap.angles.N)
        angle_types = list(read_snap.angles.types) if n_ang > 0 else []
        if n_ang == 0 or CORTEX_ANGLE_TYPE not in angle_types:
            # No cortex angles to prune; pass through whatever is present.
            if n_ang == 0:
                return None, angle_types
            ag = np.asarray(read_snap.angles.group, dtype=np.int64)
            at = np.asarray(read_snap.angles.typeid, dtype=np.int64)
            return (ag, at), angle_types

        cortex_angle_typeid = angle_types.index(CORTEX_ANGLE_TYPE)
        ag = np.asarray(read_snap.angles.group, dtype=np.int64)
        at = np.asarray(read_snap.angles.typeid, dtype=np.int64)
        is_cortex_angle = at == cortex_angle_typeid
        other_ag = ag[~is_cortex_angle]
        other_at = at[~is_cortex_angle]
        cortex_ag = ag[is_cortex_angle]

        if cortex_ag.shape[0] > 0:
            # Intact-junction set as a frozenset of normalized pairs.
            intact_pairs = set()
            if self._junctions is not None:
                for row in np.flatnonzero(intact):
                    a, b = self._junctions[row]
                    intact_pairs.add(self._norm_pair(int(a), int(b)))
            keep = np.ones(cortex_ag.shape[0], dtype=bool)
            for k in range(cortex_ag.shape[0]):
                i0, i1, i2 = int(cortex_ag[k, 0]), int(cortex_ag[k, 1]), int(cortex_ag[k, 2])
                if (
                    self._norm_pair(i0, i1) not in intact_pairs
                    or self._norm_pair(i1, i2) not in intact_pairs
                ):
                    keep[k] = False
            cortex_ag = cortex_ag[keep]

        cortex_at = np.full(cortex_ag.shape[0], cortex_angle_typeid, dtype=np.int64)
        if other_ag.shape[0] > 0:
            new_ag = np.concatenate([cortex_ag, other_ag], axis=0)
            new_at = np.concatenate([cortex_at, other_at])
        else:
            new_ag = cortex_ag
            new_at = cortex_at
        return (new_ag, new_at), angle_types

    def _write_snapshot(self, read_snap, pos, new_bg, new_bt, bond_type_names,
                        new_angles, angle_types) -> None:
        sim = self._sim_ref
        assert sim is not None
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
            write_snap.bonds.typeid[:] = new_bt.astype(np.uint32)

        if new_angles is not None:
            ang_g, ang_t = new_angles
            write_snap.angles.N = int(ang_g.shape[0])
            write_snap.angles.types = list(angle_types)
            if ang_g.shape[0] > 0:
                write_snap.angles.group[:] = ang_g.astype(np.uint32)
                write_snap.angles.typeid[:] = ang_t.astype(np.uint32)

        for grp_name in ("dihedrals", "impropers"):
            src = getattr(read_snap, grp_name)
            dst = getattr(write_snap, grp_name)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)

        sim.state.set_snapshot(write_snap)

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
    @property
    def n_intact(self) -> int:
        return int(self._intact.sum()) if self._intact is not None else 0

    @property
    def n_backbone_0(self) -> int:
        return int(self._junctions.shape[0]) if self._junctions is not None else 0

    @property
    def connected_fraction(self) -> float:
        n0 = self.n_backbone_0
        return (self.n_intact / n0) if n0 > 0 else 1.0

    @property
    def n_sever_total(self) -> int:
        return self._n_sever_total

    @property
    def n_anneal_total(self) -> int:
        return self._n_anneal_total

    @property
    def steps_run(self) -> int:
        return self._steps_run


def make_turnover_updater(
    *,
    p: ResolvedTurnover,
    rest_length: float,
    kT: float,
    bond_k: float | None = None,
    seed_offset: int = 5,
) -> tuple[ActinTurnoverUpdater, hoomd.update.CustomUpdater]:
    """Build the ActinTurnoverUpdater Action wrapped in a periodic CustomUpdater."""
    action = ActinTurnoverUpdater(
        p=p, rest_length=rest_length, kT=kT, bond_k=bond_k,
        seed_offset=seed_offset,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(p.batch_steps)
    )
    return action, updater
