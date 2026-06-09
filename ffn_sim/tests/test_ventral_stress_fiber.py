"""Sanity gates for the ventral stress-fiber layout (ADHERENT pivot, Stage 2a).

Geometry only — explicit bead-spring actin bundle, mixed polarity, FA-anchored ends on the
flat ventral plane. No mechanics yet (myosin/α-actinin in Stage 2b).
"""
from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.cortex.ventral_stress_fiber import generate_ventral_sf_layout


def _layout(n_fil=8, L=5.0e-6, ell0=5.0e-7, z_basal=-3.75e-6, br=2.0e-7, seed=0):
    return generate_ventral_sf_layout(
        n_fil=n_fil, fiber_length=L, ell0=ell0, z_basal=z_basal,
        bundle_radius=br, rng=np.random.default_rng(seed),
    )


def test_dimensions_and_topology():
    lay = _layout()
    nb = lay.n_beads_per_fil
    assert nb == int(round(5.0e-6 / 5.0e-7)) + 1 == 11
    assert lay.positions.shape == (lay.n_fil * nb, 3)
    assert lay.backbone_bonds.shape == (lay.n_fil * (nb - 1), 2)
    # bond length ≈ ell0 along the fiber
    d = lay.positions[lay.backbone_bonds[:, 0]] - lay.positions[lay.backbone_bonds[:, 1]]
    seg = np.linalg.norm(d, axis=1)
    assert np.allclose(seg, 5.0e-7, atol=1e-12)


def test_ventral_above_substrate():
    """Every bead sits ON/ABOVE the ventral plane (adherent, never below the substrate)."""
    lay = _layout()
    assert np.all(lay.positions[:, 2] >= lay.z_basal)


def test_aligned_along_axis():
    """Filament tangents are along the fiber axis ±x̂ (aligned bundle, not isotropic)."""
    lay = _layout()
    for f in range(lay.n_fil):
        beads = np.where(lay.filament_idx == f)[0]
        t = lay.positions[beads[-1]] - lay.positions[beads[0]]
        t /= np.linalg.norm(t)
        assert abs(abs(np.dot(t, lay.axis)) - 1.0) < 1e-9


def test_mixed_polarity_present():
    """Mixed polarity (both +x̂ and -x̂ filaments) so bipolar myosin has antiparallel overlap."""
    lay = _layout()
    assert set(np.unique(lay.polarity)) == {-1, 1}


def test_anchor_beads_are_fiber_ends():
    """FA-anchor beads are the extreme-x (end) beads — the traction reaction set."""
    lay = _layout()
    assert lay.anchor_beads.size == 2 * lay.n_fil
    xs = lay.positions[:, 0]
    xmin, xmax = xs.min(), xs.max()
    for b in lay.anchor_beads:
        assert np.isclose(xs[b], xmin) or np.isclose(xs[b], xmax)


def test_boundary_cases():
    one = _layout(n_fil=1)          # single chain: valid, no bundle
    assert one.n_fil == 1 and one.polarity.tolist() == [1]
    with pytest.raises(ValueError):
        _layout(L=1.0e-7, ell0=5.0e-7)   # fiber shorter than one segment
    with pytest.raises(ValueError):
        _layout(n_fil=0)
