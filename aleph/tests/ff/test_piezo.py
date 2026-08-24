r"""Piezo1 tension-gated open-probability reporter (ff/piezo) — the KB-3.10 logistic + the unit-trap guard.

Piezo1 is a pure REPORTER: it reads the membrane tension the model already carries and returns an open
probability ``P_open(γ) = 1/(1 + exp(−(γ − γ_half)/γ_s))`` (Cox 2016 / Lewis-Grandl). It bears no mechanical
load, and the mechano-chemical feedback (P_open → Ca²⁺ → contractility) is DEFAULT-OFF (gain 0) until the
channel density / conductance / Ca-gain are KB-registered.

THE gate is the UNIT TRAP: KB-3.10 is published in pN/nm; the FF membrane tension is carried in pN/µm. The
conversion (×1000) is the whole ballgame — at the resting bilayer tension γ_mem ≈ 10 pN/µm (= 0.01 pN/nm)
Piezo must read essentially CLOSED (P_open ≈ 0.035). Drop the conversion (treat 10 as pN/nm vs γ_half=5 pN/nm)
and it reads ≈0.97 (wide open) — the bug this test catches. Checked against the ACTUAL resting membrane, not
just the docstring number.
"""

import numpy as np
import pytest
import warp as wp

from aleph.laws.compartments import resolve_membrane
from aleph.laws.piezo import (
    CA_FEEDBACK_GAIN,
    GAMMA_HALF_PN_NM,
    GAMMA_S_PN_NM,
    PN_UM_PER_PN_NM,
    p_open,
    piezo_popen_kernel,
    resolve_piezo,
)

D = "cpu"


def test_resolve_converts_pn_nm_to_pn_um():
    """The pN/nm → FF pN/µm conversion (×1000) happens once, at resolve time (the unit-trap guard lives here)."""
    pz = resolve_piezo()
    assert PN_UM_PER_PN_NM == 1000.0
    assert pz.gamma_half_pn_um == pytest.approx(GAMMA_HALF_PN_NM * 1000.0)   # 5 pN/nm → 5000 pN/µm
    assert pz.gamma_s_pn_um == pytest.approx(GAMMA_S_PN_NM * 1000.0)         # 1.5 pN/nm → 1500 pN/µm


def test_logistic_midpoint_is_half():
    """At γ = γ_half the channel is exactly half-open (logistic centre)."""
    pz = resolve_piezo()
    assert float(p_open(pz.gamma_half_pn_um, pz)) == pytest.approx(0.5, abs=1e-12)


def test_unit_trap_guard_resting_is_closed():
    """UNIT-TRAP GUARD: at the resting membrane tension Piezo reads closed (P_open < 0.1); the dropped-conversion
    bug would read wide-open (> 0.9). Tested against the ACTUAL resting bilayer tension, not only the ~10 pN/µm
    docstring example."""
    pz = resolve_piezo()
    # docstring example: 10 pN/µm resting bilayer
    assert float(p_open(10.0, pz)) == pytest.approx(0.0347, abs=5e-3)
    assert float(p_open(10.0, pz)) < 0.1
    # the real model membrane (resolve_membrane default = 10 pN/µm, KB-3.B1.1 band [3,40]) is also closed
    g_rest = resolve_membrane().gamma_mem
    assert float(p_open(g_rest, pz)) < 0.1
    # the BUG: if the ×1000 conversion were dropped, 10 would be compared to γ_half=5 pN/nm → wide open
    p_bug = 1.0 / (1.0 + np.exp(-(10.0 - GAMMA_HALF_PN_NM) / GAMMA_S_PN_NM))
    assert p_bug > 0.9                                         # ~0.966 — the failure mode the conversion prevents


def test_p_open_monotone_and_saturates():
    """P_open is monotone increasing in tension and saturates 0→1 across the gating range."""
    pz = resolve_piezo()
    gs = np.array([0.0, 1000.0, 3000.0, 5000.0, 7000.0, 9000.0, 20000.0])
    po = np.array([float(p_open(g, pz)) for g in gs])
    assert np.all(np.diff(po) > 0)
    assert po[0] < 0.05 and po[-1] > 0.99


def test_feedback_is_default_off():
    """The mechano-chemical feedback gain is 0 (registration-gated) — Piezo is a reporter, not yet a driver."""
    assert CA_FEEDBACK_GAIN == 0.0
    assert resolve_piezo().ca_gain == 0.0


def test_kernel_matches_closed_form_and_area_weights():
    """piezo_popen_kernel per-face P_open matches the closed form; open_area accumulates Σ P_open·A so the
    whole-cell open fraction = open_area / total_area."""
    pz = resolve_piezo()
    tens = np.array([10.0, 5000.0, 9000.0, 20000.0])          # closed, half, mostly-open, saturated
    areas = np.array([0.4, 0.1, 0.25, 0.05])
    tension = wp.array(tens, dtype=wp.float64, device=D)
    area = wp.array(areas, dtype=wp.float64, device=D)
    popen = wp.zeros(len(tens), dtype=wp.float64, device=D)
    open_area = wp.zeros(1, dtype=wp.float64, device=D)
    wp.launch(piezo_popen_kernel, dim=len(tens),
              inputs=[tension, area, wp.float64(pz.gamma_half_pn_um), wp.float64(pz.gamma_s_pn_um),
                      popen, open_area], device=D)
    want = np.array([float(p_open(g, pz)) for g in tens])
    assert np.allclose(popen.numpy(), want, atol=1e-12)
    assert float(open_area.numpy()[0]) == pytest.approx(float((want * areas).sum()), rel=1e-12)
    open_fraction = float(open_area.numpy()[0]) / areas.sum()
    assert 0.0 < open_fraction < 1.0
