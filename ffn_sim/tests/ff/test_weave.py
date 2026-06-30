"""Unified actin-architecture builder (ff/weave) — task (a) increment 1 validation.

Proves the unification: ONE ``weave()`` produces both the isotropic cortical shell (CORTEX) and a
parallel bundle (FILOPODIUM) from a spec, with NO structure-specific branch beyond manifold + bind
mode. CORTEX must reproduce the H.3 γ-floor cortex; FILOPODIUM must hit the bundle architectural metrics.
"""

import numpy as np
import pytest

from ffn_sim.ff.architecture_metrics import (
    bundle_count,
    cortex_metrics,
    inter_filament_spacing_nm,
    parallel_order_parameter,
)
from ffn_sim.ff.architecture_spec import CORTEX, FILOPODIUM
from ffn_sim.ff.weave import weave


def test_weave_cortex_reproduces_gamma_floor():
    """weave(CORTEX) reproduces gamma_floor.build_crosslinked_cortex (same RNG order ⇒ bit-exact γ) —
    so routing the cortex through the unified builder does NOT regress the γ-floor production number."""
    from ffn_sim.ff.gamma_floor import (
        NMIIA_MINIFIL_STALL_PN,
        CortexParams,
        build_crosslinked_cortex,
        equilibrate,
        measure_gamma,
    )

    cw = weave(CORTEX, rng=np.random.default_rng(0))
    cb = build_crosslinked_cortex(CortexParams(), n_filaments=1000, n_xl=1000, n_myo=100,
                                  rng=np.random.default_rng(0))
    # structure parity
    assert cw.net.n_nodes == cb.net.n_nodes
    assert cw.xl_i.size == cb.xl_i.size and cw.myo_i.size == cb.myo_i.size
    # γ parity (bit-exact — RNG call order preserved)
    equilibrate(cw, 0.0, n_steps=300, method="device", device="cpu")
    equilibrate(cb, 0.0, n_steps=300, method="device", device="cpu")
    gw = measure_gamma(cw, NMIIA_MINIFIL_STALL_PN, turgor=False)["gamma_active"]
    gb = measure_gamma(cb, NMIIA_MINIFIL_STALL_PN, turgor=False)["gamma_active"]
    assert abs(gw - gb) / gb < 1e-9

    m = cortex_metrics(cw)
    assert m["n_filaments"] == 1000
    assert 0.7 < m["areal_density_um2"] < 0.9              # ≈ 0.8/µm² (H.3 cortex)
    assert m["connectivity_z"] > 1.0


def test_weave_filopodium_bundle_metrics():
    """weave(FILOPODIUM) — the SAME builder, manifold='bundle' — produces a parallel bundle with the
    literature architectural metrics (bundle count 10–30, nematic order ≈1, ~7–8 nm spacing)."""
    cf = weave(FILOPODIUM, rng=np.random.default_rng(0))
    assert 10 <= bundle_count(cf.net) <= 30                # filopodium bundle count
    assert parallel_order_parameter(cf.net) > 0.95         # tightly parallel (uniform polarity)
    assert 6.0 <= inter_filament_spacing_nm(cf.net) <= 10.0 # fascin cross-bridge spacing ~7–8 nm
    assert cf.myo_i.size == 0                               # filopodial core has no motor


def test_weave_one_builder_two_architectures():
    """The unification claim: the same weave() gives DISTINCT architectures (isotropic vs parallel)."""
    cortex = weave(CORTEX, rng=np.random.default_rng(0))
    filo = weave(FILOPODIUM, rng=np.random.default_rng(0))
    assert parallel_order_parameter(cortex.net) < 0.3      # cortex ≈ isotropic
    assert parallel_order_parameter(filo.net) > 0.95       # filopodium ≈ parallel


def test_weave_lamellipodium_dendritic_metrics():
    """weave(LAMELLIPODIUM) — the SAME builder, manifold='patch' — produces an Arp2/3 dendritic array
    with the literature architecture: branch junctions at θ₀=70° (Fäßler 2020 68±9°) AND a ±35°
    two-mode filament orientation about the protrusion axis (Mueller 2017)."""
    from ffn_sim.ff.architecture_metrics import branch_angle_distribution, two_mode_orientation
    from ffn_sim.ff.architecture_spec import LAMELLIPODIUM

    cx = weave(LAMELLIPODIUM, rng=np.random.default_rng(0))
    assert cx.branch_triples.shape[0] > 0                   # Arp2/3 branch junctions emitted
    ba = branch_angle_distribution(cx)
    assert abs(ba["mean_deg"] - 70.0) < 5.0                 # branch junctions at ~70° (Fäßler 2020)
    tm = two_mode_orientation(cx)
    assert abs(tm["plus_mode_deg"] - 35.0) < 8.0           # +35° mode
    assert abs(tm["minus_mode_deg"] + 35.0) < 8.0          # −35° mode (Mueller 2017 ±35°)
    assert tm["two_mode_frac"] > 0.6                        # dendritic two-mode signature


def test_weave_lamellipodium_relax_stable_with_branch_kernel():
    """The lamellipodium relaxes STABLY on-device with the Arp2/3 angle-harmonic branch kernel wired in
    (bending + crosslink/anchor springs + branch kernel + reshape): finite, and the branch junctions
    stay near θ₀=70° (within the soft k_angle spread)."""
    from ffn_sim.ff.architecture_metrics import branch_angle_distribution
    from ffn_sim.ff.architecture_spec import (
        ARP23_BRANCH_ANGLE_RAD,
        ARP23_BRANCH_K,
        LAMELLIPODIUM,
    )
    from ffn_sim.ff.network_warp import relax_on_device

    cx = weave(LAMELLIPODIUM, rng=np.random.default_rng(0))
    links = np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64) if cx.xl_i.size else None
    cx.net.pos = relax_on_device(cx.net, links=links, k_xl=cx.xl_k, xl_rest=cx.xl_rest,
                                 branch_triples=cx.branch_triples, branch_theta0=ARP23_BRANCH_ANGLE_RAD,
                                 branch_k=ARP23_BRANCH_K, n_steps=800, device="cpu")
    assert np.isfinite(cx.net.pos).all()                    # stable (CFL incl branch k_eff)
    ba = branch_angle_distribution(cx)
    assert 55.0 < ba["mean_deg"] < 85.0                    # branches held near 70° by the kernel
