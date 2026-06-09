"""Tests for the H.7 fine-grained traction → Layer-2 CBM scale-bridge (h7_traction_seam).

The bridge is a PURE arithmetic/IO resolver (no sim): it turns H.7's mechanistic single-SF
coherent contractile traction (+131 pN) into the CBM per-cell f_active via the FA-pairing SF
count, and supplies the physiological (Gil-Redondo) anchor. These tests pin the dimensional
contract, the anchors (no magic numbers), the seam-record precedence, and the provisional flag.
"""

from __future__ import annotations

import json

import pytest

from ffn_sim.spheroid import h7_traction_seam as seam
from ffn_sim.spheroid.h7_traction_seam import (
    FA_COUNT_PER_CELL,
    GIL_REDONDO_CONTACT_AREA_M2,
    GIL_REDONDO_TOTAL_CONTRACTILITY_N,
    N_SF_PER_CELL,
    resolve_h7_traction,
)


def test_f_cell_test_is_per_sf_times_n_sf():
    """Dimensional contract: f_cell_test = per_SF · N_SF."""
    r = resolve_h7_traction()
    assert r.f_cell_test == pytest.approx(r.per_sf_traction * r.n_sf, rel=1e-12)


def test_per_sf_is_the_decisive_sarcomeric_value():
    """per_SF = +131 pN (decisive single sarcomeric-SF), CONTRACTILE (positive)."""
    r = resolve_h7_traction()
    assert r.per_sf_traction == pytest.approx(131.08e-12, rel=0.01)
    assert r.per_sf_traction > 0.0          # contractile (the mixed-polarity − value is NOT consumed)


def test_n_sf_is_fa_pairing_anchor():
    """N_SF = FA_count/2 (each ventral SF spans 2 FAs) — derived, not tuned."""
    assert N_SF_PER_CELL == pytest.approx(FA_COUNT_PER_CELL / 2.0, rel=1e-12)
    assert N_SF_PER_CELL == pytest.approx(21.5, rel=1e-9)


def test_test_density_lands_in_the_few_nN_regime():
    """f_cell_test ≈ 2.8 nN — near the B1 ceiling, between the old heuristics (1.6/9.4 nN)."""
    r = resolve_h7_traction()
    assert r.f_cell_test == pytest.approx(2.82e-9, rel=0.05)
    assert 1.6e-9 < r.f_cell_test < 9.4e-9


def test_physiological_anchor_is_gil_redondo():
    """f_cell_physio = Gil-Redondo 102 nN; stress ≈ 56 Pa (≈ their 63 Pa)."""
    r = resolve_h7_traction()
    assert r.f_cell_physio == pytest.approx(GIL_REDONDO_TOTAL_CONTRACTILITY_N, rel=1e-12)
    assert r.stress_physio_Pa == pytest.approx(
        GIL_REDONDO_TOTAL_CONTRACTILITY_N / GIL_REDONDO_CONTACT_AREA_M2, rel=1e-12
    )
    assert r.stress_physio_Pa == pytest.approx(56.0, rel=0.05)


def test_density_gap_factor():
    """Physiological ÷ test-density ≈ ×36 (≈ sanctioned ×40 mesoscale; the honest residual)."""
    r = resolve_h7_traction()
    assert r.density_gap_factor == pytest.approx(r.f_cell_physio / r.f_cell_test, rel=1e-12)
    assert 25.0 < r.density_gap_factor < 45.0


def test_provisional_without_seam_record():
    """No h7_traction_seam.json present ⇒ provisional from the single-SF value."""
    r = resolve_h7_traction()
    assert r.source == "provisional_single_sf"
    assert r.provisional is True


def test_boundary_zero_per_sf_via_n_sf_scaling():
    """N_SF must be > 0; a positive n_sf scales f_cell_test linearly (boundary/monotone)."""
    r1 = resolve_h7_traction(n_sf=10.0)
    r2 = resolve_h7_traction(n_sf=20.0)
    assert r2.f_cell_test == pytest.approx(2.0 * r1.f_cell_test, rel=1e-12)
    with pytest.raises(ValueError):
        resolve_h7_traction(n_sf=0.0)


def test_seam_record_takes_precedence(tmp_path, monkeypatch):
    """When the H.7 seam record is present, the consumer reads its measured array aggregate."""
    rec = {
        "per_SF_traction_N": 150.0e-12,
        "per_SF_traction_sem_N": 9.0e-12,
        "single_cell_total_traction_N": 5.0e-9,
        "n_SF": 20,
        "provenance": "H.7 array aggregate (test)",
    }
    d = tmp_path / "seam_inputs"
    d.mkdir()
    (d / "h7_traction_seam.json").write_text(json.dumps(rec))
    monkeypatch.setattr(seam, "_SEAM_INPUTS", d)
    monkeypatch.setattr(seam, "_H7_PROD", tmp_path / "nonexistent")
    r = resolve_h7_traction()
    assert r.source == "seam_record"
    assert r.provisional is False
    assert r.per_sf_traction == pytest.approx(150.0e-12, rel=1e-12)
    assert r.f_cell_test == pytest.approx(5.0e-9, rel=1e-12)   # measured aggregate, not per_SF×N_SF
    assert r.n_sf == pytest.approx(20.0, rel=1e-12)
