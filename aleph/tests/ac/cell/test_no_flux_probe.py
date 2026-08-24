"""Non-vacuity gate for the nucleus relative-no-flux probe (``_nucleus_suppressed_flux_kernel``).

The previous probe computed ``mobility·(p−p)/dx`` — a self-minus-self difference that is IDENTICALLY 0 for
ANY field, so ``leaked_flux()`` could never detect a leak (vacuous). The fix reads the REAL nucleus-neighbour
cell, so the probe is a genuine measurement: 0 iff the field is uniform across the envelope, > 0 under a
transmembrane gradient. I0-A requires these kernel tests to run on CUDA; they skip on the dev Mac.
"""
from __future__ import annotations

import numpy as np
import pytest

wp = pytest.importorskip("warp")

from aleph.components.fluid.boundary import _nucleus_suppressed_flux_kernel  # noqa: E402

wp.init()
_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)
pytestmark = pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
DEV = str(_CUDA_DEVICE) if _CUDA_DEVICE is not None else ""
FLUID, OUTSIDE, NUCLEUS = 1, 0, 2


def _probe(p: np.ndarray, mask: np.ndarray, mobility: float = 5e-3, dx: float = 0.5) -> float:
    pd = wp.array(np.ascontiguousarray(p, np.float64), dtype=wp.float64, device=DEV)
    md = wp.array(np.ascontiguousarray(mask, np.int32), dtype=wp.int32, device=DEV)
    out = wp.zeros(1, dtype=wp.float64, device=DEV)
    wp.launch(_nucleus_suppressed_flux_kernel, dim=p.shape,
              inputs=[pd, md, wp.float64(mobility), wp.float64(dx)], outputs=[out], device=DEV)
    wp.synchronize_device(DEV)
    return float(out.numpy()[0])


def _sphere_nucleus_mask(n: int = 12) -> np.ndarray:
    """A grid with a central NUCLEUS blob surrounded by FLUID (like the real classification)."""
    mask = np.full((n, n, n), FLUID, np.int32)
    c = (n - 1) / 2.0
    idx = np.indices((n, n, n)).transpose(1, 2, 3, 0) - c
    r = np.linalg.norm(idx, axis=-1)
    mask[r <= n * 0.22] = NUCLEUS
    mask[r >= n * 0.48] = OUTSIDE
    return mask


def test_uniform_field_gives_zero_flux():
    """A uniform pore-pressure field ⇒ no gradient across the envelope ⇒ suppressed flux exactly 0."""
    mask = _sphere_nucleus_mask()
    p = np.full(mask.shape, 40.0)                         # Π₀ everywhere
    assert _probe(p, mask) == 0.0


def test_gradient_field_gives_nonzero_flux():
    """A field with a transmembrane gradient ⇒ the mask holds a NONZERO would-be flux (the probe has teeth).

    This is exactly what the old p−p probe could not detect: it returned 0 here too. The fix makes it > 0.
    """
    mask = _sphere_nucleus_mask()
    xr = np.linspace(-1.0, 1.0, mask.shape[0])
    p = 40.0 + 10.0 * xr[:, None, None] * np.ones(mask.shape)   # linear ramp along x
    assert _probe(p, mask) > 1e-6


def test_no_nucleus_faces_gives_zero():
    """With no NUCLEUS cells adjacent to fluid, there are no envelope faces ⇒ 0 (nothing to suppress)."""
    mask = np.full((10, 10, 10), FLUID, np.int32)
    mask[0, :, :] = OUTSIDE
    p = np.random.default_rng(0).random(mask.shape)
    assert _probe(p, mask) == 0.0


def test_old_probe_would_be_vacuous():
    """Document the bug: the self-difference the OLD probe computed is 0 for the SAME gradient field the fixed
    probe flags as nonzero — proving the fix added falsifiability, not just a rename."""
    mask = _sphere_nucleus_mask()
    xr = np.linspace(-1.0, 1.0, mask.shape[0])
    p = 40.0 + 10.0 * xr[:, None, None] * np.ones(mask.shape)
    old_vacuous = float(np.sum(np.abs(p - p)))           # what mobility·(p−p)/dx summed to: identically 0
    assert old_vacuous == 0.0
    assert _probe(p, mask) > old_vacuous                 # fixed probe strictly separates the two
