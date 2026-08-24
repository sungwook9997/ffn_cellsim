"""ac/ intermediate-filament cage — modulus presets (CPU) + force-free rest & nucleus coupling (CUDA).

The IF cage's job is to UNFREEZE the nucleus: a cortex deformation must reach the nucleus through the radial
spokes. These tests certify the ac/ wrapper (a) carries the keratin<vimentin EMT modulus split and (b) is
force-free at build yet transmits a cortex displacement to the nucleus beads.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.solid.intermediate_filament import IF_E_KERATIN_PA, IF_E_VIMENTIN_PA

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)


def _sphere_points(n: int, R: float) -> np.ndarray:
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    th = np.pi * (1.0 + 5.0 ** 0.5) * k
    return R * np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), np.cos(phi)], axis=1)


def test_emt_modulus_split() -> None:
    """Keratin (MCF7) is the softer preset; vimentin (MDA-231) the stiffer large-strain one; both ~MPa."""
    assert IF_E_KERATIN_PA < IF_E_VIMENTIN_PA
    assert 1e6 <= IF_E_KERATIN_PA <= 2e7 and 1e6 <= IF_E_VIMENTIN_PA <= 2e7


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_cage_is_force_free_at_rest_and_couples_nucleus() -> None:
    from aleph.components.solid.intermediate_filament import build_if_compartment
    dev = str(_CUDA_DEVICE)
    R_nuc, R_cortex = 5.0, 7.4
    nuc_pos = _sphere_points(200, R_nuc)
    cortex_pos = _sphere_points(400, R_cortex)
    n_nuc, n_cortex = nuc_pos.shape[0], cortex_pos.shape[0]
    cage = build_if_compartment(
        centre=(0.0, 0.0, 0.0), R_nuc_um=R_nuc, R_cortex_um=R_cortex, nuc_pos=nuc_pos, cortex_pos=cortex_pos,
        n_fil=40, cell_type="keratin", nuc_offset=0, cortex_offset=n_nuc, if_offset=n_nuc + n_cortex, device=dev)
    assert cage.n_backbone > 0 and cage.n_linc == 40 and cage.n_anchor == 40

    pos_all = np.concatenate([nuc_pos, cortex_pos, cage.pos0], axis=0)
    n_total = pos_all.shape[0]
    with wp.ScopedDevice(dev):
        pos = wp.array(pos_all, dtype=wp.vec3d, device=dev)
        f = wp.zeros(n_total, dtype=wp.vec3d, device=dev)
        cage.accumulate(pos, f)
        f_rest = np.linalg.norm(f.numpy(), axis=1).max()
        assert f_rest < 1e-6, f"cage must be force-free at build, got max|F|={f_rest:.3g} pN"

        # Single force-eval only stretches bonds whose endpoints moved (full transmission to the far side needs
        # relaxation = the outer solver). So test the two ENDS of the connected load path directly:
        # (i) cortex→cage: displacing cortex nodes outward loads the anchor bonds → the IF spoke beads feel it.
        if_slice = slice(n_nuc + n_cortex, n_total)
        p = pos_all.copy()
        cortex_slice = slice(n_nuc, n_nuc + n_cortex)
        r = np.linalg.norm(p[cortex_slice], axis=1, keepdims=True)
        p[cortex_slice] = p[cortex_slice] * (1.0 + 0.5 / r)
        fd = wp.zeros(n_total, dtype=wp.vec3d, device=dev)
        cage.accumulate(wp.array(p, dtype=wp.vec3d, device=dev), fd)
        f_if = np.linalg.norm(fd.numpy()[if_slice], axis=1).max()
        assert f_if > 1e-3, f"a cortex displacement must load the cage via the anchors, got {f_if:.3g} pN"

        # (ii) cage→nucleus: displacing the nucleus inward loads the LINC bonds → the nucleus beads feel it
        # (the un-freezing of the previously-decoupled nucleus).
        p2 = pos_all.copy()
        rn = np.linalg.norm(p2[:n_nuc], axis=1, keepdims=True)
        p2[:n_nuc] = p2[:n_nuc] * (1.0 - 0.4 / rn)
        fd2 = wp.zeros(n_total, dtype=wp.vec3d, device=dev)
        cage.accumulate(wp.array(p2, dtype=wp.vec3d, device=dev), fd2)
        f_nuc = np.linalg.norm(fd2.numpy()[:n_nuc], axis=1).max()
        assert f_nuc > 1e-3, f"a nucleus displacement must be resisted via LINC (nucleus coupled), got {f_nuc:.3g} pN"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_vimentin_preset_builds() -> None:
    from aleph.components.solid.intermediate_filament import build_if_compartment
    dev = str(_CUDA_DEVICE)
    nuc_pos, cortex_pos = _sphere_points(200, 5.0), _sphere_points(400, 7.4)
    cage = build_if_compartment(centre=(0.0, 0.0, 0.0), R_nuc_um=5.0, R_cortex_um=7.4, nuc_pos=nuc_pos,
                                cortex_pos=cortex_pos, n_fil=40, cell_type="vimentin",
                                nuc_offset=0, cortex_offset=200, if_offset=600, device=dev)
    assert cage.cell_type == "vimentin" and cage.E_if_Pa == IF_E_VIMENTIN_PA
