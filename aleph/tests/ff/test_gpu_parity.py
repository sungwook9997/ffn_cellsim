"""GPU↔CPU bit-parity of every FF on-device kernel + integrator — runs on the A5000, skips offline.

FF is GPU-native (Warp): the production path is `device="cuda:0"`. This pins that path as
bit-faithful to the `device="cpu"` reference — the same skip-if-absent pattern as the Cytosim parity
oracle (no CUDA on the dev Mac → skipped; on the gbook A5000 → runs). Validated 2026-07-01
(FF_STAGE6R): at native scale every force kernel agrees to ≤2.2e-11 (float64 atomic-add ordering) and
the full relax integrator to ~1e-14 after 400 steps. Tolerances here are loose headroom over those.
"""

import numpy as np
import pytest

wp = pytest.importorskip("warp")
if not wp.is_cuda_available():
    pytest.skip("no CUDA device (FF GPU-parity runs on the A5000, skipped offline)",
                allow_module_level=True)

from aleph.laws.gamma_floor import CortexParams, build_crosslinked_cortex
from aleph.laws.network_warp import (
    link_spring_force_np,
    myosin_force_np,
    relax_on_device,
    reshape_np,
)


def _maxabs(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


@pytest.fixture(scope="module")
def cortex():
    return build_crosslinked_cortex(CortexParams(), n_filaments=4000, n_xl=4000, n_myo=400,
                                    rng=np.random.default_rng(1))


def test_reshape_kernel_cpu_cuda_parity(cortex):
    pos, foff, srest = cortex.net.pos, cortex.net.fiber_offsets, cortex.net.seg_rest
    c = reshape_np(pos, foff, srest, n_iter=4, device="cpu")
    g = reshape_np(pos, foff, srest, n_iter=4, device="cuda:0")
    assert _maxabs(c, g) < 1e-12          # constraint projection is exact-arithmetic identical


def test_link_spring_kernel_cpu_cuda_parity(cortex):
    links = np.stack([cortex.xl_i, cortex.xl_j], 1).astype(np.int64)
    c = link_spring_force_np(cortex.net.pos, links, cortex.xl_k, cortex.xl_rest, device="cpu")
    g = link_spring_force_np(cortex.net.pos, links, cortex.xl_k, cortex.xl_rest, device="cuda:0")
    assert _maxabs(c, g) < 1e-8           # float64 atomic-add ordering only


def test_myosin_kernel_cpu_cuda_parity(cortex):
    mlinks = np.stack([cortex.myo_i, cortex.myo_j], 1).astype(np.int64)
    c = myosin_force_np(cortex.net.pos, mlinks, 5.0, device="cpu")
    g = myosin_force_np(cortex.net.pos, mlinks, 5.0, device="cuda:0")
    assert _maxabs(c, g) < 1e-12


def test_relax_integrator_cpu_cuda_parity(cortex):
    """The FULL on-device integrator (bending + crosslink + reshape, 400 steps) stays bit-faithful —
    the per-step float64 atomic differences do NOT amplify over the trajectory."""
    import copy
    links = np.stack([cortex.xl_i, cortex.xl_j], 1).astype(np.int64)
    a, b = copy.deepcopy(cortex.net), copy.deepcopy(cortex.net)
    pc = relax_on_device(a, links=links, k_xl=cortex.xl_k, xl_rest=cortex.xl_rest,
                         n_steps=400, device="cpu")
    pg = relax_on_device(b, links=links, k_xl=cortex.xl_k, xl_rest=cortex.xl_rest,
                         n_steps=400, device="cuda:0")
    assert _maxabs(pc, pg) < 1e-10        # ~1e-14 in practice; loose headroom
