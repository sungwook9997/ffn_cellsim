"""KU-3.5-active GATE-B: bipolar-completion antiparallel-availability metric.

Locks in the pure-geometry root-cause metric
(``h3_ku35_completion_diag.antiparallel_partner_availability``) on KNOWN
synthetic lattices. This is the structural ceiling on
``frac_complete_pairs``: a minifilament whose + side grips filament A can only
COMPLETE a bipolar contractile pair if its − side finds an ANTIPARALLEL partner
filament (opposite minus-end-ward tangent projection on the rod axis) within the
head reach. The metric must satisfy:

  * all-parallel actin (no antiparallel anywhere) → availability == 0 even
    though actin IS densely in reach (the limiter is POLARITY, not reach).
  * checkerboard polarity, neighbours in reach → availability == 1.
  * checkerboard polarity but spacing ≫ reach → availability == 0 with the reach
    EMPTY (the limiter is OVERLAP-SCARCITY: antiparallel actin exists but is out
    of reach in a sparse mesoscale cortex).
  * growing the reach on the sparse-but-antiparallel lattice RECOVERS
    availability (reach monotonicity) — distinguishing a reach limit from a
    polarity limit.

MEASUREMENT-ONLY; changes no physics; HOOMD-free.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.scripts.h3_ku35_completion_diag import (
    CompletionGeometry,
    antiparallel_partner_availability,
    _lattice,
    _minus_tangents,
    self_test,
)


# Production-derived geometry constants (phase1_h3.yaml myosin block).
HEAD_OFF = 200e-9          # head_rest_length
CAPTURE_PERP = 210e-9      # head_actin_capture_perp
BACKBONE_HALF = 350e-9     # backbone_length 700 nm / 2
REACH = HEAD_OFF + CAPTURE_PERP + BACKBONE_HALF   # ~760 nm − side reach


def _geo(F_side, spacing, polarity, length=300e-9, seed=0):
    pos, tang = _lattice(F_side, spacing, length, polarity,
                         np.random.default_rng(seed))
    return CompletionGeometry(pos, tang, HEAD_OFF, CAPTURE_PERP, BACKBONE_HALF)


def test_minus_tangent_is_negated_plus_tangent():
    """m̂ = −(plus-end-ward tangent), unit-normalised (cortex Option-A polarity)."""
    g = _geo(3, 0.3 * REACH, "all_plus")
    m = _minus_tangents(g)
    assert np.allclose(m, -g.tangents)
    assert np.allclose(np.linalg.norm(m, axis=1), 1.0)


def test_all_parallel_lattice_zero_availability_despite_dense_actin():
    """POLARITY limit: dense actin in reach but ALL same polarity → avail == 0."""
    g = _geo(6, 0.4 * REACH, "all_plus")
    r = antiparallel_partner_availability(g, rng=np.random.default_rng(1))
    assert r["availability"] == 0.0
    assert r["frac_reach_empty"] == 0.0       # actin IS densely present
    assert r["n_probed"] > 0


def test_checkerboard_in_reach_full_availability():
    """Antiparallel neighbours within reach → availability == 1."""
    g = _geo(6, 0.4 * REACH, "checker")
    r = antiparallel_partner_availability(g, rng=np.random.default_rng(2))
    assert r["availability"] == 1.0
    assert r["mean_n_partners"] > 0.0


def test_sparse_checkerboard_zero_availability_overlap_scarcity():
    """OVERLAP-SCARCITY: antiparallel actin EXISTS but spacing ≫ reach → 0."""
    g = _geo(6, 5.0 * REACH, "checker")
    r = antiparallel_partner_availability(g, rng=np.random.default_rng(3))
    assert r["availability"] == 0.0
    assert r["frac_reach_empty"] > 0.5        # nothing near the rod centre


def test_reach_scale_monotonic_recovery():
    """Reach monotonicity: on a sparse-but-antiparallel lattice, growing the −
    side reach RECOVERS availability (so reach, not polarity, was the limiter
    there) — distinguishing a reach limit from a structural polarity limit."""
    g = _geo(6, 5.0 * REACH, "checker")
    base = antiparallel_partner_availability(
        g, rng=np.random.default_rng(3), reach_scale=1.0)["availability"]
    grown = antiparallel_partner_availability(
        g, rng=np.random.default_rng(3), reach_scale=6.0)["availability"]
    assert grown > base
    assert grown == pytest.approx(1.0)


def test_availability_is_a_fraction():
    """Metric is a well-defined fraction in [0, 1] for an isotropic random mix."""
    rng = np.random.default_rng(11)
    F, N = 200, 4
    pos = rng.uniform(-2e-6, 2e-6, size=(F, N, 3))
    tang = rng.normal(size=(F, 3))
    tang /= np.linalg.norm(tang, axis=1, keepdims=True)
    g = CompletionGeometry(pos, tang, HEAD_OFF, CAPTURE_PERP, BACKBONE_HALF)
    r = antiparallel_partner_availability(g, rng=np.random.default_rng(12))
    assert 0.0 <= r["availability"] <= 1.0
    assert r["R_reach"] == pytest.approx(REACH)


def test_self_test_passes():
    """The module's bundled self-test (synthetic lattices) passes end-to-end."""
    assert self_test(verbose=False) is True
