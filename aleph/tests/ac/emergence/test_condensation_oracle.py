"""Self-test of the local condensation fields + bundle localization (pure NumPy, no Warp/CUDA).

The core I5 deliverable: a real replacement for the vacuous ``bundle_count`` that LOCALIZES a condensed
bundle out of an isotropic background, is rotation-equivariant + permutation-invariant, holds across the
neighborhood-scale plateau (grid-invariance), and rejects crowding-without-order.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.emergence import synthetic as syn
from aleph.components.emergence.condensation import (
    NEIGHBORHOOD_SCALE,
    compute_local_fields,
    localize_bundles,
    neighborhood_radius,
)

_NB, _NBUN = 400, 40  # background / bundle fiber counts for the planted fixtures


def _score(members: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """(precision, recall) of a recovered member set against the planted ground-truth mask."""
    rec = set(members.tolist())
    planted = set(np.flatnonzero(mask).tolist())
    tp = len(rec & planted)
    return tp / max(len(rec), 1), tp / max(len(planted), 1)


def test_planted_bundle_is_localized() -> None:
    """A compact aligned patch is recovered as ONE bundle: high recall/precision, centroid within a spacing."""
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3)
    fields = compute_local_fields(pos, off)
    bundles = localize_bundles(fields)
    assert len(bundles) >= 1
    top = max(bundles, key=lambda b: b.size)
    precision, recall = _score(top.members, mask)
    assert recall >= 0.8
    assert precision >= 0.7
    assert float(np.linalg.norm(top.center - center)) < fields.radius


def test_recovered_axis_matches_planted_director() -> None:
    """The recovered bundle's nematic axis matches the planted director (up to nematic sign)."""
    d = np.array([0.0, 0.0, 1.0])
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3, director=d)
    top = max(localize_bundles(compute_local_fields(pos, off)), key=lambda b: b.size)
    assert abs(abs(float(top.axis @ d)) - 1.0) < 0.05
    assert top.s_bundle > 0.5


def test_localization_is_rotation_equivariant() -> None:
    """Rotate the config -> the recovered centroid rotates by the SAME R (equivariance), same member count."""
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3)
    top = max(localize_bundles(compute_local_fields(pos, off)), key=lambda b: b.size)
    rng = np.random.default_rng(0)
    r = syn.random_rotation(rng)
    pos_r = syn.rotate_config(pos, r, about=center)
    top_r = max(localize_bundles(compute_local_fields(pos_r, off)), key=lambda b: b.size)
    expected = center + (top.center - center) @ r.T
    assert float(np.linalg.norm(top_r.center - expected)) < 1e-6
    assert top_r.size == top.size


def test_localization_is_permutation_invariant() -> None:
    """Relabel the fibers -> the SAME recovered member set (mapped through the permutation)."""
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3)
    top = set(max(localize_bundles(compute_local_fields(pos, off)), key=lambda b: b.size).members.tolist())
    rng = np.random.default_rng(1)
    perm = rng.permutation(off.shape[0] - 1)
    pos_p, off_p = syn.permute_fibers(pos, off, perm)
    top_p = max(localize_bundles(compute_local_fields(pos_p, off_p)), key=lambda b: b.size)
    remapped = set(perm[top_p.members].tolist())
    assert remapped == top


@pytest.mark.parametrize("scale", [2.5, 3.0, 3.5, 4.0])
def test_recovery_holds_across_neighborhood_scale_plateau(scale: float) -> None:
    """The planted bundle is recovered across the NEIGHBORHOOD_SCALE plateau (grid-invariance, not a
    single tuned k)."""
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3)
    fields = compute_local_fields(pos, off, scale=scale)
    bundles = localize_bundles(fields)
    assert len(bundles) >= 1
    top = max(bundles, key=lambda b: b.size)
    _, recall = _score(top.members, mask)
    assert recall >= 0.75
    assert float(np.linalg.norm(top.center - center)) < fields.radius


def test_crowding_without_order_is_not_a_bundle() -> None:
    """A dense but ISOTROPIC patch is NOT flagged: condensation requires ORDER, not mere crowding."""
    pos, off, mask, center = syn.dense_isotropic_patch(_NB, 60, seed=17)
    bundles = localize_bundles(compute_local_fields(pos, off))
    # no bundle should sit on the dense isotropic patch; if anything survives it must be tiny + off-center
    for b in bundles:
        _, recall = _score(b.members, mask)
        assert recall < 0.3


def test_pure_isotropic_recovers_no_bundle() -> None:
    """A homogeneous isotropic network yields no localized bundle (specificity of localization)."""
    pos, off = syn.isotropic_network(440, seed=23)
    assert localize_bundles(compute_local_fields(pos, off)) == []


def test_neighborhood_radius_scales_with_spacing() -> None:
    """The derived radius scales linearly with the config's own spacing (grid-invariant construction):
    uniformly scaling every length by 2x doubles the radius, and r = scale * median-NN-spacing."""
    from aleph.components.emergence.nematic import fiber_centroids

    pos, off = syn.isotropic_network(400, seed=24, box=10.0)
    cent = fiber_centroids(pos, off)
    r1 = neighborhood_radius(cent)
    r2 = neighborhood_radius(cent * 2.0)  # double every length -> spacing doubles -> radius doubles
    assert r2 == pytest.approx(2.0 * r1, rel=1e-9)
    # r == NEIGHBORHOOD_SCALE * median nearest-neighbor spacing, computed independently here
    d = np.linalg.norm(cent[None] - cent[:, None], axis=-1)
    np.fill_diagonal(d, np.inf)
    s_med = np.median(d.min(axis=1))
    assert r1 == pytest.approx(NEIGHBORHOOD_SCALE * s_med, rel=1e-9)
    assert r1 > 0.0


def test_local_field_shapes_and_density() -> None:
    """Local fields have the right shapes; the dense patch has higher local density than the background."""
    pos, off, mask, center = syn.planted_bundle(_NB, _NBUN, seed=3)
    fields = compute_local_fields(pos, off)
    f = off.shape[0] - 1
    assert fields.density.shape == (f,) and fields.s_local.shape == (f,) and fields.align.shape == (f,)
    assert fields.density[mask].mean() > fields.density[~mask].mean()  # bundle patch is denser
