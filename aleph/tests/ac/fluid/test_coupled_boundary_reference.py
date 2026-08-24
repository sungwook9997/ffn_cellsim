"""NumPy gate: coupled enclosed-volume fluid candidate (mass + no-flux + adjoint-work close together).

Pure-host acceptance oracle for the ``cytosol`` FLUID_VOLUME CONNECTED milestone (roadmap §3.1): the Biot
field plus BOTH moving boundaries (semipermeable membrane + impermeable nucleus) close the three ledger
channels simultaneously.  Reuses the ``fv_reference`` stencils the Warp-CUDA solver ports.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.coupled_boundary_reference import EnclosedVolumeCoupledReference
from aleph.components.fluid.fv_reference import FLUID, NUCLEUS, OUTSIDE


def _enclosed_mask(n=10):
    """A cubic domain: OUTSIDE shell (membrane), FLUID cytosol, NUCLEUS core."""
    mask = np.full((n, n, n), FLUID, dtype=np.int_)
    mask[0, :, :] = mask[-1, :, :] = OUTSIDE
    mask[:, 0, :] = mask[:, -1, :] = OUTSIDE
    mask[:, :, 0] = mask[:, :, -1] = OUTSIDE
    c = n // 2
    mask[c - 1:c + 1, c - 1:c + 1, c - 1:c + 1] = NUCLEUS
    return mask


def _ref(mask):
    return EnclosedVolumeCoupledReference(
        shape=mask.shape, dx=0.5, mobility=1.0, storage_S=1e-3, mask=mask, L_p=2.0, sigma_refl=1.0
    )


def test_coupled_mass_channel_content_equals_membrane_permeation() -> None:
    mask = _enclosed_mask()
    ref = _ref(mask)
    rng = np.random.default_rng(1)
    p = np.where(mask == FLUID, 30.0 + rng.random(mask.shape), 0.0)
    dt = 0.4 * ref._biot.cfl_dt()
    _, ledger = ref.coupled_step(p, d_pi_osm=40.0, dt=dt)
    # content change equals integrated membrane permeation; nucleus contributes no permeation term.
    assert ledger.mass_residual < 1e-9
    assert abs(ledger.membrane_outward_flux) > 0.0  # membrane is actually exchanging
    assert ledger.nucleus_applied_flux == 0.0


def test_no_flux_teeth_are_load_bearing() -> None:
    mask = _enclosed_mask()
    ref = _ref(mask)
    # a transmembrane gradient across the nuclear envelope -> teeth > 0 (mask holds the flux at zero)
    p = np.zeros(mask.shape)
    p[mask == FLUID] = np.arange(np.count_nonzero(mask == FLUID))  # non-uniform across the envelope
    assert ref.nucleus_teeth(p) > 0.0
    assert ref.nucleus_applied_flux() == 0.0
    # a field uniform ACROSS the envelope (nucleus cells at the same pressure) has nothing to suppress
    assert ref.nucleus_teeth(np.full(mask.shape, 5.0)) == pytest.approx(0.0)


def test_adjoint_work_conjugacy_is_a_transpose_pair() -> None:
    rng = np.random.default_rng(2)
    n_surface, n_grid = 12, 40
    # normalized Peskin-like spread: each surface node spreads to a few grid cells, columns sum to 1.
    spread = np.zeros((n_grid, n_surface))
    for s in range(n_surface):
        cells = rng.choice(n_grid, size=4, replace=False)
        w = rng.random(4)
        spread[cells, s] = w / w.sum()
    v = rng.random(n_surface)
    g = rng.random(n_grid)
    assert EnclosedVolumeCoupledReference.adjoint_work_residual(spread, v, g) < 1e-12


def test_moving_membrane_remap_conserves_content_with_swept_term() -> None:
    mask = _enclosed_mask()
    ref = _ref(mask)
    rng = np.random.default_rng(3)
    p = np.where(mask == FLUID, 20.0 + rng.random(mask.shape), 0.0)
    old_total = ref.total_content(p)
    # the membrane advances: one FLUID cell just inside the shell becomes OUTSIDE (envelope moved outward->in).
    new_mask = mask.copy()
    new_mask[1, 5, 5] = OUTSIDE
    p_new, swept = ref.remap_swept_content(p, new_mask)
    cell_vol = ref.dx**ref.dim
    new_total = float(ref.storage_S * np.sum(p_new[new_mask == FLUID]) * cell_vol)
    # total content over the new domain equals the old total plus the swept (crossed-face) content.
    assert abs(new_total - (old_total + swept)) < 1e-9


def test_mask_shape_must_match_grid() -> None:
    with pytest.raises(ValueError, match="mask shape"):
        EnclosedVolumeCoupledReference(
            shape=(4, 4, 4), dx=0.5, mobility=1.0, storage_S=1e-3,
            mask=np.full((3, 3, 3), FLUID, dtype=np.int_),
        )
