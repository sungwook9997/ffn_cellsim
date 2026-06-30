"""Hand kinetic layer — KMC binding / stepping / unbinding (NF2007 §10.1, grounded).

A *Hand* (Cytosim's term) is the binding element of a motor or crosslinker. Following Nédélec &
Foethke 2007 §10.1 (p20–21) a Hand has two procedures:

  * **attach(m)** — an unbound Hand binds the closest fiber site within a capture radius ``ε`` with
    per-tick probability ``1 − exp(−τ·k_on)`` (Poisson; ≈ τ·k_on for small τ·k_on).
  * **step(f)** — a bound Hand, given its load ``f`` from the mechanics, (i) advances its abscissa
    by the active step ``δa = τ·v_max·(1 − f/f_stall)`` (motors; crosslinkers have v_max=0), and
    (ii) detaches with a **force-dependent** off-rate ``p_off`` → per-tick detach probability
    ``1 − exp(−τ·p_off)``.

Two grounded off-rate laws:

  * **Bell slip** ``p_off = p₀ · exp(|f|/f₀)`` — NF2007 §10.1 (kinesin), also non-muscle myosin II
    head–actin and α-actinin (Bell length x_β ⇒ f₀ = k_BT/x_β).
  * **Pereverzev catch–slip** ``p_off = k_catch0·exp(−f·x_catch/k_BT) + k_slip0·exp(+f·x_slip/k_BT)``
    — filamin/F-actin (the project's faithful filamin catch bond, KU-3.19).

This module is the **regime-B kinetic layer** of the FF engine (ENGINE.md §1): it runs
event-driven *on top of* the mechanical equilibrium solve — each tick (a) solve the mechanics to
get loads, (b) run these KMC trials. It carries the kinetics + state only; the geometric site
search (nearest fiber point) and the link forces are supplied by the caller (cortex assembly /
γ-floor harness, Stage 6d). All parameters are FF units (pN, µm, s).

Parameter presets are LITERATURE-ANCHORED (none invented):
  * KINESIN — NF2007 Fig 8 worked example (the validation anchor): f_stall=5pN, v0=0.4µm/s,
    k_on=10/s, p₀=0.5/s, f₀=2.5pN, ε=10nm, link k=200pN/µm.
  * NMIIA_MYOSIN — v0=0.2µm/s (Kovács 2003), F_stall=0.5pN/head, k_off0=0.35/s (Stam-Hocky 2015 /
    Tam 2021), Bell x_β=0.6nm (Veigel 2002), k_on=50/s (configs/phase1_h3.yaml myosin block).
  * ALPHA_ACTININ — k_off0=0.066/s (Ferrer 2008 PNAS), x_β=0.4nm, k_on=10/s, ε=60nm (no stepping).
  * FILAMIN — Pereverzev catch–slip k_catch0=0.1/s (Furuike 2001), x_catch=0.8nm, k_slip0=0.02/s,
    x_slip=0.3nm (Pereverzev 2005 form), k_on=10/s, ε=60nm (no stepping).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ffn_sim.ff import units as U


# ── Off-rate laws (FF units: force pN, rate 1/s) ──────────────────────────────────────────────────
def bell_off_rate(f: np.ndarray | float, p0: float, f0: float) -> np.ndarray | float:
    """Bell slip off-rate ``p_off = p₀·exp(|f|/f₀)`` (NF2007 §10.1). ``f`` load [pN], ``f₀`` [pN]."""
    return p0 * np.exp(np.abs(f) / f0)


def pereverzev_off_rate(f: np.ndarray | float, k_catch0: float, x_catch_um: float,
                        k_slip0: float, x_slip_um: float, kT: float = U.KBT) -> np.ndarray | float:
    """Pereverzev two-pathway catch–slip off-rate (filamin/F-actin; KU-3.19).

    ``p_off = k_catch0·exp(−f·x_catch/k_BT) + k_slip0·exp(+f·x_slip/k_BT)`` — off-rate FALLS with
    load (catch) up to a peak F*, then RISES (slip). ``f`` [pN], ``x_*`` [µm], ``kT`` [pN·µm].
    """
    f = np.asarray(f, dtype=np.float64) if not np.isscalar(f) else f
    return k_catch0 * np.exp(-f * x_catch_um / kT) + k_slip0 * np.exp(+f * x_slip_um / kT)


def bell_f0_from_x_beta(x_beta_um: float, kT: float = U.KBT) -> float:
    """Bell characteristic force f₀ = k_BT/x_β [pN] from a Bell length x_β [µm]."""
    return kT / x_beta_um


# ── Hand type parameters ──────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True, slots=True)
class HandParams:
    """Kinetic parameters of one Hand type (FF units; pN, µm, s). Grounded inputs — see module doc.

    Attributes:
        name: label.
        k_on: per-tick binding rate when an acceptor is in range [1/s].
        capture_radius_um: binding capture radius ε [µm].
        p0: zero-force off-rate prefactor [1/s] (Bell) — used iff ``catch_slip`` is False.
        f0: Bell characteristic force [pN] — used iff ``catch_slip`` is False.
        v_max_um_s: unloaded stepping speed [µm/s] (0 for a pure crosslinker).
        f_stall: stall force [pN] (motors). Ignored if v_max=0.
        catch_slip: if True use the Pereverzev law with the ``cs_*`` fields instead of Bell.
        cs_k_catch0, cs_x_catch_um, cs_k_slip0, cs_x_slip_um: Pereverzev catch–slip parameters.
        link_k: interaction-link stiffness [pN/µm] (the bound Hand's spring to its anchor).
    """

    name: str
    k_on: float
    capture_radius_um: float
    p0: float = 0.0
    f0: float = 1.0
    v_max_um_s: float = 0.0
    f_stall: float = 1.0
    catch_slip: bool = False
    cs_k_catch0: float = 0.0
    cs_x_catch_um: float = 0.0
    cs_k_slip0: float = 0.0
    cs_x_slip_um: float = 0.0
    link_k: float = 0.0

    def off_rate(self, f):
        """Force-dependent off-rate p_off(f) [1/s] for this Hand type (Bell or catch–slip)."""
        if self.catch_slip:
            return pereverzev_off_rate(f, self.cs_k_catch0, self.cs_x_catch_um,
                                       self.cs_k_slip0, self.cs_x_slip_um)
        return bell_off_rate(f, self.p0, self.f0)

    def active_step(self, f, tau: float, *, clamp: bool = False):
        """Active abscissa advance δa = τ·v_max·(1 − f/f_stall) [µm] (NF2007 §10.1).

        ``clamp=True`` floors the velocity factor at 0 (no super-stall back-stepping); the faithful
        default (False) keeps the paper's linear form (back-stepping above stall is then possible,
        but such a Hand has a large p_off and detaches).
        """
        fac = 1.0 - np.asarray(f, dtype=np.float64) / self.f_stall
        if clamp:
            fac = np.maximum(fac, 0.0)
        return tau * self.v_max_um_s * fac


# ── Grounded presets (FF units) ───────────────────────────────────────────────────────────────────
KINESIN = HandParams(                 # NF2007 Fig 8 (validation anchor)
    name="kinesin", k_on=10.0, capture_radius_um=0.010, p0=0.5, f0=2.5,
    v_max_um_s=0.4, f_stall=5.0, link_k=200.0)

NMIIA_MYOSIN = HandParams(            # non-muscle myosin IIA head (configs/phase1_h3.yaml myosin)
    name="nmiia", k_on=50.0, capture_radius_um=0.210,   # head_actin_capture_perp 210 nm
    p0=0.35, f0=bell_f0_from_x_beta(0.6e-3),            # k_off0=0.35/s, x_β=0.6nm (→ f0≈7.1pN)
    v_max_um_s=0.2, f_stall=0.5, link_k=1.0)            # v0=0.2µm/s, F_stall=0.5pN, k_head=1pN/µm

ALPHA_ACTININ = HandParams(          # α-actinin crosslinker Hand (Ferrer 2008)
    name="alpha_actinin", k_on=10.0, capture_radius_um=0.060,
    p0=0.066, f0=bell_f0_from_x_beta(0.4e-3),           # k_off0=0.066/s, x_β=0.4nm (→ f0≈10.7pN)
    v_max_um_s=0.0, link_k=4.6e5)        # k_xl=455 pN/nm (Ferrer 2008 PNAS AFM; PI-approved 2026-06-30)
# CROSSLINK STIFFNESS (link_k) — lit-anchored: Ferrer 2008 PNAS AFM α-actinin 455 pN/nm = 4.6e5 pN/µm,
# filamin 820 pN/nm = 8.2e5 (SAME paper as the α-actinin k_off0 above). The old 0.1 was the broken
# AFINES soft-surrogate (~4.5e6× too soft; FF_KIM §c). The stiff value is now SAFE in production via the
# crosslink-turnover resting baseline (Stage 6N, gamma_floor.equilibrate(crosslink_turnover=True)):
# settle the shell bending-only (stiff crosslinks excluded → no CFL/overshoot) then bind crosslinks
# force-free (Hand §10.1 attach). The ACTIVE γ (myosin) is identical to the soft value — the γ-floor is
# force-magnitude-limited, robust to crosslink stiffness (last root-cause CLOSED, 6M). SE rows pending
# (Ferrer stiffness, SE_REGISTRATION_CANDIDATES §5). ⚠️ explicit/implicit relaxers WITHOUT
# crosslink_turnover still over-stretch at this stiffness — use crosslink_turnover (or kim_network's
# own solver). FF_STAGE6M/6N docs.

FILAMIN = HandParams(                 # filamin crosslinker Hand (Pereverzev catch–slip)
    name="filamin", k_on=10.0, capture_radius_um=0.060,
    catch_slip=True, cs_k_catch0=0.1, cs_x_catch_um=0.8e-3,
    cs_k_slip0=0.02, cs_x_slip_um=0.3e-3, v_max_um_s=0.0,
    link_k=8.2e5)                        # k_xl=820 pN/nm (Ferrer 2008 PNAS AFM, companion; PI-approved)

INTEGRIN_A5B1 = HandParams(           # integrin α5β1–fibronectin FA clutch (Kong 2009, KB-2.5) — PI-gated
    name="integrin_a5b1", k_on=1.0, capture_radius_um=0.300,   # h_c=300 nm FA gap (Kim 2012; PI-gated k_on)
    catch_slip=True,
    cs_k_catch0=0.4, cs_x_catch_um=U.KBT / 7.0,    # k_catch=0.4/s, Fc=7 pN  → catch term k_catch·e^(−F/Fc)
    cs_k_slip0=0.5, cs_x_slip_um=U.KBT / 30.0,     # k_slip=0.5/s, Fs=30 pN → slip term k_slip·e^(+F/Fs)
    v_max_um_s=0.0, link_k=1.0e3)                  # k_int=1 pN/nm = 1e3 pN/µm FA clutch spring (PI-gated)
# ⚠️ FA clutch = INTEGRIN-ECM (Kong 2009 α5β1-fibronectin catch-slip), NOT Rakshit 2012 / KU-4.2 (that
# is cell-cell E-CADHERIN). Same two-pathway pereverzev_off_rate FORM, KB-2.5 numbers (k_slip=0.5/s,
# Fs=30pN, k_catch=0.4/s, Fc=7pN). ⚠️ KB-CONSISTENCY FLAG (→PI): these recorded params give an interior
# lifetime peak F* = ln[(k_catch·Fs)/(k_slip·Fc)]/(1/Fc+1/Fs) ≈ 7 pN, NOT the "F*≈30 pN" the KB-2.5
# claim text states — the catch-bond (catch→slip) IS present, but the peak-force claim is inconsistent
# with its own params; used AS-RECORDED (no tuning). k_on, k_int = Phase-default (KB-2.4/2.18) → PI-gated.

PRESETS = {h.name: h for h in (KINESIN, NMIIA_MYOSIN, ALPHA_ACTININ, FILAMIN, INTEGRIN_A5B1)}


# ── Per-tick probabilities (Poisson) ──────────────────────────────────────────────────────────────
def attach_probability(tau: float, k_on: float) -> float:
    """Per-tick attach probability 1 − exp(−τ·k_on) (Poisson; ≈ τ·k_on for small τ·k_on)."""
    return float(-np.expm1(-tau * k_on))


def detach_probability(tau: float, p_off):
    """Per-tick detach probability 1 − exp(−τ·p_off) for off-rate ``p_off`` [1/s]."""
    return -np.expm1(-tau * np.asarray(p_off, dtype=np.float64))


# ── Vectorized Hand population (the regime-B kinetic state) ────────────────────────────────────────
@dataclass(slots=True)
class HandPopulation:
    """A population of Hands of one type — the kinetic state evolved by the KMC ticks.

    Geometry-agnostic: the caller supplies, each tick, the nearest-acceptor distance (for attach)
    and the per-Hand load (for step/detach). State arrays (length = n_hands):

        bound:    (n,) bool   — bound vs free.
        anchor:   (n,) int    — index of the bound fiber node/site (−1 if free). Caller-defined.
        abscissa: (n,) float  — curvilinear position advanced by active stepping [µm].
    """

    params: HandParams
    bound: np.ndarray
    anchor: np.ndarray
    abscissa: np.ndarray

    @classmethod
    def empty(cls, n: int, params: HandParams) -> "HandPopulation":
        """All-free population of ``n`` Hands."""
        return cls(params=params, bound=np.zeros(n, bool),
                   anchor=np.full(n, -1, np.int64), abscissa=np.zeros(n, np.float64))

    @property
    def n(self) -> int:
        return int(self.bound.shape[0])

    @property
    def bound_fraction(self) -> float:
        return float(self.bound.mean()) if self.n else 0.0

    def attach_step(self, nearest_dist: np.ndarray, nearest_idx: np.ndarray, tau: float,
                    rng: np.random.Generator) -> int:
        """Run one attachment sub-step (NF2007 attach(m)). Free Hands whose nearest acceptor is
        within ε bind it with prob 1−exp(−τk_on).

        Args:
            nearest_dist: (n,) distance to the closest acceptor site [µm] (∞ if none).
            nearest_idx: (n,) index of that acceptor site (caller-defined; used to set ``anchor``).
            tau: tick interval [s].
            rng: random generator.

        Returns:
            number of Hands that newly bound.
        """
        p = attach_probability(tau, self.params.k_on)
        eligible = (~self.bound) & (nearest_dist <= self.params.capture_radius_um)
        roll = rng.random(self.n) < p
        new = eligible & roll
        self.bound[new] = True
        self.anchor[new] = nearest_idx[new]
        self.abscissa[new] = 0.0
        return int(new.sum())

    def step_detach(self, loads: np.ndarray, tau: float, rng: np.random.Generator,
                    *, clamp_step: bool = False) -> int:
        """Run one step+detach sub-step (NF2007 step(f)) for the bound Hands.

        Advances each bound Hand's abscissa by the active step and detaches it with the
        force-dependent probability. ``loads`` is the per-Hand load magnitude [pN].

        Returns:
            number of Hands that detached.
        """
        bmask = self.bound
        if not bmask.any():
            return 0
        f = np.asarray(loads, dtype=np.float64)
        # active stepping (motors only; v_max=0 ⇒ no advance)
        if self.params.v_max_um_s != 0.0:
            self.abscissa[bmask] += self.params.active_step(f[bmask], tau, clamp=clamp_step)
        # force-dependent detachment
        p_off = self.params.off_rate(f[bmask])
        pdet = detach_probability(tau, p_off)
        roll = rng.random(int(bmask.sum())) < pdet
        idx = np.where(bmask)[0][roll]
        self.bound[idx] = False
        self.anchor[idx] = -1
        self.abscissa[idx] = 0.0
        return int(roll.sum())
