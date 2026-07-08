r"""Crawl-driver knobs (scripts/ff_crawl_on_substrate) — pure-helper coverage for the 2026-07-09 overnight additions.

The crawl driver's full build()/run() loop is an integration entry point (heavy: cortex build + Warp), but the
discrete-FA site selector ``_fps_subsample`` is a pure, deterministic function and IS unit-testable. It implements
the physical `--n-fa` knob (place N discrete focal-adhesion sites instead of one clutch per cortex node), so its
spread + determinism are contracts worth pinning.
"""

import numpy as np

from ffn_sim.scripts.ff_crawl_on_substrate import _fps_subsample


def test_fps_returns_all_when_k_ge_n():
    """k ≥ n ⇒ every point is a site (identity), no resampling."""
    pts = np.random.default_rng(0).standard_normal((7, 3))
    assert np.array_equal(_fps_subsample(pts, 7), np.arange(7))
    assert np.array_equal(_fps_subsample(pts, 99), np.arange(7))


def test_fps_selects_k_unique_indices():
    pts = np.random.default_rng(1).standard_normal((200, 3))
    sel = _fps_subsample(pts, 30)
    assert sel.shape == (30,)
    assert len(np.unique(sel)) == 30                          # no duplicates
    assert sel.min() >= 0 and sel.max() < 200


def test_fps_is_deterministic():
    """No RNG inside — identical input ⇒ identical sites (reproducible FA placement)."""
    pts = np.random.default_rng(2).standard_normal((150, 3))
    assert np.array_equal(_fps_subsample(pts, 20), _fps_subsample(pts, 20))


def test_fps_is_spread_not_clustered():
    """Farthest-point picks well-separated sites: on a 1-D line it grabs BOTH ends before the middle, and the
    minimum inter-site spacing beats a naive contiguous slice of the same count."""
    line = np.stack([np.linspace(0.0, 1.0, 100), np.zeros(100), np.zeros(100)], axis=1)
    sel = _fps_subsample(line, 5)
    xs = np.sort(line[sel, 0])
    assert xs[0] < 0.02 and xs[-1] > 0.98                     # both extremes chosen
    fps_min_gap = np.diff(xs).min()
    contiguous_min_gap = np.diff(np.sort(line[:5, 0])).min()  # first-5 slice = tightly clustered
    assert fps_min_gap > contiguous_min_gap                   # FPS is genuinely more spread
