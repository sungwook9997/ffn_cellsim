"""Tests for the lamellipodium → CBM active-motility scale-bridge (motility_bridge).

The bridge is a PURE arithmetic resolver (no sim): it maps an anchored single-cell speed to the
CBM active force via the overdamped relation f = v0·γ_cell, for two candidate anchors. These
tests pin the bridge contract + the load-bearing discriminating numbers (protrusion ~9.4 nN
exceeds the 3 nN ceiling; whole-cell ~1.6 nN is in-band; ratio ~6× ≈ the magnitude gap).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.motility_bridge import (
    V0_MCF7_WHOLECELL,
    resolve_active_traction,
)
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def _resolved():
    return resolve_layer2(yaml.safe_load(_CFG.read_text()))


def test_f_active_is_v0_times_gamma():
    """Dimensional/arithmetic contract: f = v0·γ for both anchors."""
    m = resolve_active_traction(_resolved())
    assert m.f_protrusion == pytest.approx(m.v0_protrusion * m.gamma_cell, rel=1e-12)
    assert m.f_wholecell == pytest.approx(m.v0_wholecell * m.gamma_cell, rel=1e-12)


def test_protrusion_anchor_value_and_ceiling():
    """Bieling v0 = k_elong0·δ_elong ⇒ ~31 nm/s ⇒ ~9.4 nN, which EXCEEDS the 3 nN ceiling."""
    m = resolve_active_traction(_resolved())
    assert m.v0_protrusion == pytest.approx(31.3e-9, rel=0.02)   # 11.6/s × 2.7 nm
    assert m.f_protrusion == pytest.approx(9.4e-9, rel=0.05)
    assert m.protrusion_exceeds_ceiling is True
    assert m.f_protrusion > m.ceiling


def test_wholecell_anchor_value_in_band():
    """MCF7 whole-cell 19 µm/h ⇒ ~1.6 nN, inside the B1 stable band."""
    m = resolve_active_traction(_resolved())
    assert m.v0_wholecell == pytest.approx(19.0e-6 / 3600.0, rel=1e-9)  # 19 µm/h in m/s
    assert m.f_wholecell == pytest.approx(1.6e-9, rel=0.06)
    assert m.f_wholecell < m.ceiling


def test_active_ratio_matches_the_magnitude_gap():
    """The protrusion/whole-cell ratio ~6× ≈ the observed ~5–9× magnitude gap (the active-side fork)."""
    m = resolve_active_traction(_resolved())
    assert m.active_ratio == pytest.approx(5.9, rel=0.1)
    assert 5.0 <= m.active_ratio <= 9.0


def test_bieling_constant_read_from_h5_sot():
    """The protrusion v0 tracks the H.5 config (single source of truth), not a duplicated number."""
    h5 = yaml.safe_load((Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml").read_text())

    def _find(key, node):
        if isinstance(node, dict):
            if key in node:
                return float(node[key])
            for v in node.values():
                r = _find(key, v)
                if r is not None:
                    return r
        return None

    ke, de = _find("k_elong_0", h5), _find("delta_elong", h5)
    m = resolve_active_traction(_resolved())
    assert m.v0_protrusion == pytest.approx(ke * de, rel=1e-12)


def test_overlay_v0_override():
    """An alternate measured whole-cell speed flows through (override path)."""
    m = resolve_active_traction(_resolved(), v0_wholecell=2.0 * V0_MCF7_WHOLECELL)
    assert m.f_wholecell == pytest.approx(2.0 * V0_MCF7_WHOLECELL * m.gamma_cell, rel=1e-12)
