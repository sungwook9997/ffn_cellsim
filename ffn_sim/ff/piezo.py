r"""Piezo1 mechanosensitive channel — tension-gated open-probability reporter (FF, Warp, µm·pN·s).

The FF cell has no biochemical/mechanosensing readout layer; Piezo1 is the canonical membrane-tension→signal
transducer (the force-FELT gate, distinct from the force-BEARING compartments). It bears no mechanical load —
it READS the membrane tension the model already computes and reports an open probability. This is a pure
REPORTER: the mechano-chemical FEEDBACK (P_open → Ca²⁺ → contractility) is a DEFAULT-OFF stub (gain = 0) until
the channel density / conductance / Ca-gain are KB-registered (surfaced, not defaulted).

KB-3.10: ``P_open(γ) = 1 / (1 + exp(−(γ − γ_half)/γ_s))`` with half-activation tension γ_half ≈ 5 pN/nm and
sensitivity γ_s ≈ 1–2 pN/nm (Cox 2016 / Lewis-Grandl patch-tension gating). The FF membrane tension γ is carried
in **pN/µm**; KB-3.10 is in **pN/nm** — the unit conversion (÷1000) is the whole ballgame (the "unit trap").
The guard gate: at the resting bilayer tension γ_mem ≈ 10 pN/µm (= 0.01 pN/nm) Piezo is essentially CLOSED,
``P_open(rest) ≈ 0.035`` — if the conversion were dropped (10 treated as 10 pN/nm) it would read ~0.97 (open),
which is the bug this gate catches.

No magic numbers: γ_half, γ_s are the KB-3.10 values; the feedback gain is 0 (registration-gated).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

PN_UM_PER_PN_NM = 1000.0        # 1 pN/nm = 1000 pN/µm  (unit-trap conversion)

# KB-3.10 (kept in pN/nm as published; converted to FF pN/µm at resolve time)
GAMMA_HALF_PN_NM = 5.0          # half-activation membrane tension
GAMMA_S_PN_NM = 1.5            # activation sensitivity (band 1–2)
CA_FEEDBACK_GAIN = 0.0         # ⚠ mechano-chemical feedback DEFAULT-OFF (density/conductance/Ca-gain unregistered)


@dataclass
class PiezoParams:
    gamma_half_pn_um: float     # half-activation tension in FF units [pN/µm]
    gamma_s_pn_um: float        # sensitivity [pN/µm]
    ca_gain: float = CA_FEEDBACK_GAIN


def resolve_piezo(gamma_half_pn_nm: float = GAMMA_HALF_PN_NM, gamma_s_pn_nm: float = GAMMA_S_PN_NM) -> PiezoParams:
    """KB-3.10 gating tensions converted pN/nm → FF pN/µm (the unit-trap guard lives here)."""
    return PiezoParams(gamma_half_pn_um=gamma_half_pn_nm * PN_UM_PER_PN_NM,
                       gamma_s_pn_um=gamma_s_pn_nm * PN_UM_PER_PN_NM)


def p_open(gamma_pn_um, p: PiezoParams | None = None):
    """Piezo1 open probability at membrane tension ``gamma_pn_um`` [pN/µm] (KB-3.10 logistic)."""
    p = p or resolve_piezo()
    g = np.asarray(gamma_pn_um, float)
    return 1.0 / (1.0 + np.exp(-(g - p.gamma_half_pn_um) / p.gamma_s_pn_um))


@wp.kernel
def piezo_popen_kernel(tension: wp.array(dtype=wp.float64), area: wp.array(dtype=wp.float64),
                       gamma_half: wp.float64, gamma_s: wp.float64,
                       popen: wp.array(dtype=wp.float64), open_area: wp.array(dtype=wp.float64)):
    """Per-membrane-face open probability from the local tension (pN/µm) and area-weighted open area sum.
    ``open_area[0]`` accumulates Σ P_open·A → the whole-cell open fraction = open_area / total area."""
    i = wp.tid()
    g = tension[i]
    po = wp.float64(1.0) / (wp.float64(1.0) + wp.exp(-(g - gamma_half) / gamma_s))
    popen[i] = po
    wp.atomic_add(open_area, 0, po * area[i])


__all__ = ["PiezoParams", "resolve_piezo", "p_open", "piezo_popen_kernel",
           "GAMMA_HALF_PN_NM", "GAMMA_S_PN_NM", "PN_UM_PER_PN_NM", "CA_FEEDBACK_GAIN"]
