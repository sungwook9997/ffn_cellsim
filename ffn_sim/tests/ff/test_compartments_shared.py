"""Shared engine-agnostic compartments (common/compartments) — bit-parity + lit-anchored resolvers.

The radial-shell 3-law kernel (nucleus/membrane/turgor) is a VERBATIM port of the validated
dcm/radial_shell_warp; this guards that it stays bit-identical (so the FF and DCM engines share ONE
implementation, PI 2026-07-01) and that resolve_nucleus honors the KU-3.B2 bands.
"""

import numpy as np
import pytest

from ffn_sim.common.compartments import resolve_nucleus, run_radial_shell_warp


def _nucleus_cloud(n=60, R=2.5, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v * (R + 0.3 * rng.standard_normal((n, 1))), np.zeros(n, dtype=np.uint32)


def test_shared_kernel_bit_identical_to_dcm_original():
    """common/compartments.run_radial_shell_warp == dcm/radial_shell_warp (verbatim port → bit-identical)."""
    from ffn_sim.dcm.radial_shell_warp import run_radial_shell_warp as rs_dcm
    n = resolve_nucleus(R_nuc_um=2.5, n_beads=60)
    pos, tag = _nucleus_cloud()
    kw = dict(pos=pos, tag=tag, tag_range=(0, 1), law=0, R0=n.R_nuc_um,
              pa=n.k_chrom, pb=n.k_lamin, pc=n.d_knee_um, pd=n.F_knee_pN, device="cpu")
    fc = run_radial_shell_warp(**kw)["force"]
    fd = rs_dcm(**kw)["force"]
    assert np.array_equal(fc, fd)                     # exactly identical, not just close


def test_resolve_nucleus_bridge_and_bands():
    """KU-3.B2 dimensional bridge + band guards (E_nuc 1-10 kPa; ratio_lamin 1.4-5×, not isolated 10×)."""
    n = resolve_nucleus(R_nuc_um=2.5, n_beads=60, E_nuc_Pa=5.0e3, ratio_lamin=3.0, knee_strain=0.10)
    assert n.k_chrom == pytest.approx(4 * np.pi * 5.0e3 * 2.5 / 60)   # 1 Pa ≡ 1 pN/µm²
    assert n.k_lamin == pytest.approx(2.0 * n.k_chrom)               # (3−1)×
    assert n.d_knee_um == pytest.approx(0.25)                         # 0.10·2.5
    assert n.F_knee_pN == pytest.approx(n.k_chrom * n.d_knee_um)
    # grid-invariance: n_beads·k_chrom = 4π·E·R (intensive)
    n2 = resolve_nucleus(R_nuc_um=2.5, n_beads=120)
    assert n.n_beads * n.k_chrom == pytest.approx(n2.n_beads * n2.k_chrom)
    with pytest.raises(ValueError):
        resolve_nucleus(R_nuc_um=2.5, n_beads=60, ratio_lamin=10.0)   # isolated value rejected
    with pytest.raises(ValueError):
        resolve_nucleus(R_nuc_um=2.5, n_beads=60, E_nuc_Pa=1.0e5)     # out of KU band


def test_nucleus_strain_stiffens_past_knee():
    """The bilinear law stiffens past the lamin knee: force grows faster once |d|>d_knee (strain-stiffening)."""
    n = resolve_nucleus(R_nuc_um=2.5, n_beads=60)
    # one bead pushed inward by δ (compression): interior slope k_chrom, exterior slope k_chrom+k_lamin
    def fmag(delta):
        pos = np.array([[2.5 - delta, 0.0, 0.0], [-(2.5 - delta), 0.0, 0.0],
                        [0, 2.5, 0], [0, -2.5, 0], [0, 0, 2.5], [0, 0, -2.5]])
        tag = np.zeros(len(pos), dtype=np.uint32)
        f = run_radial_shell_warp(pos=pos, tag=tag, tag_range=(0, 1), law=0, R0=n.R_nuc_um,
                                  pa=n.k_chrom, pb=n.k_lamin, pc=n.d_knee_um, pd=n.F_knee_pN, device="cpu")["force"]
        return abs(f[0, 0])
    small = fmag(0.1)                    # inside knee (0.1 < 0.25)
    big = fmag(0.5)                      # past knee (0.5 > 0.25)
    assert big / 0.5 > small / 0.1       # stiffer per unit strain past the knee
