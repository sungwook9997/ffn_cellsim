"""Unified actin-architecture builder (ff/weave) — task (a) increment 1 validation.

Proves the unification: ONE ``weave()`` produces both the isotropic cortical shell (CORTEX) and a
parallel bundle (FILOPODIUM) from a spec, with NO structure-specific branch beyond manifold + bind
mode. CORTEX must reproduce the H.3 γ-floor cortex; FILOPODIUM must hit the bundle architectural metrics.
"""

import dataclasses

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


def _small(spec, n):
    """A small-N copy of an ArchitectureSpec for fast CPU unit tests (production CORTEX is now native
    ~38000; these tests validate the unification MECHANISM, not scale)."""
    return dataclasses.replace(spec, filament=dataclasses.replace(spec.filament, n_filaments=n))


CORTEX_SMALL = _small(CORTEX, 1000)


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

    # weave maps spec densities → counts: n_xl=round(N·xl.density_per_fil), n_myo=round(N·motor.density_per_fil).
    # At the lit-faithful CORTEX (actin 100/µm², α-actinin 1:1, NMIIA Nie 0.625/µm²) the small-N (1000)
    # cortex has n_xl=1000, n_myo=round(1000·0.00625)=6 — match build_crosslinked_cortex to those for parity.
    n_xl_s = int(round(1000 * CORTEX.crosslinker.density_per_fil))
    n_myo_s = int(round(1000 * CORTEX.motor.density_per_fil))
    cw = weave(CORTEX_SMALL, rng=np.random.default_rng(0))
    cb = build_crosslinked_cortex(CortexParams(), n_filaments=1000, n_xl=n_xl_s, n_myo=n_myo_s,
                                  rng=np.random.default_rng(0))
    # structure parity
    assert cw.net.n_nodes == cb.net.n_nodes
    assert cw.xl_i.size == cb.xl_i.size and cw.myo_i.size == cb.myo_i.size
    # γ parity (bit-exact — RNG call order preserved)
    equilibrate(cw, 0.0, n_steps=300, method="device", device="cpu")
    equilibrate(cb, 0.0, n_steps=300, method="device", device="cpu")
    gw = measure_gamma(cw, NMIIA_MINIFIL_STALL_PN, turgor=False)["gamma_active"]
    gb = measure_gamma(cb, NMIIA_MINIFIL_STALL_PN, turgor=False)["gamma_active"]
    assert abs(gw - gb) <= 1e-12 + 1e-9 * abs(gb)   # bit-exact (robust to the sparse-Nie small-N gb≈0)

    m = cortex_metrics(cw)
    assert m["n_filaments"] == 1000
    assert 1.2 < m["areal_density_um2"] < 1.6              # 1000/(4π·7.5²) ≈ 1.41/µm² (small-N test at MCF7 R)
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
    cortex = weave(CORTEX_SMALL, rng=np.random.default_rng(0))
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


def test_weave_stress_fiber_sarcomeric():
    """weave(STRESS_FIBER) — manifold='bundle' — produces a SARCOMERIC FA–FA bundle: α-actinin Z-bodies
    periodic at ~1 µm (Hotulainen 2006) with NMIIA bands ANTI-registered (Murrell 2015), from the SAME
    builder (periodic node mask, no new kernel)."""
    from ffn_sim.ff.architecture_metrics import sarcomeric_period_um
    from ffn_sim.ff.architecture_spec import STRESS_FIBER

    sf = weave(STRESS_FIBER, rng=np.random.default_rng(0))
    assert sf.myo_i.size > 0 and sf.xl_i.size > 0
    sp = sarcomeric_period_um(sf)
    assert abs(sp["period_um"] - 1.0) < 0.25                # sarcomere period ~1 µm (band 0.5–1.4)
    assert sp["anti_registered"]                            # NMIIA bands anti-registered to Z-bodies


def test_weave_microvillus_bundle_metrics():
    """weave(MICROVILLUS) — manifold='bundle' — a tight parallel finger bundle: count 20–30, ~12 nm
    packing, nematic S≈1 (no core motor)."""
    from ffn_sim.ff.architecture_metrics import bundle_count, inter_filament_spacing_nm
    from ffn_sim.ff.architecture_spec import MICROVILLUS

    mv = weave(MICROVILLUS, rng=np.random.default_rng(0))
    assert 20 <= bundle_count(mv.net) <= 30
    assert 9.0 <= inter_filament_spacing_nm(mv.net) <= 15.0  # ~12 nm lateral c2c (PI-gated)
    assert parallel_order_parameter(mv.net) > 0.95
    assert mv.myo_i.size == 0                                # microvillus core has no contractile motor


def test_weave_all_five_architectures_distinct():
    """The full unified table: ONE weave() spans 5 distinct architectures by orientation order S +
    branch + sarcomere — cortex(isotropic) / filopodium+microvillus(parallel) / lamellipodium(branched)
    / stress-fiber(sarcomeric)."""
    from ffn_sim.ff.architecture_metrics import sarcomeric_period_um
    from ffn_sim.ff.architecture_spec import LAMELLIPODIUM, MICROVILLUS, STRESS_FIBER

    assert parallel_order_parameter(weave(CORTEX_SMALL, rng=np.random.default_rng(0)).net) < 0.3
    assert parallel_order_parameter(weave(MICROVILLUS, rng=np.random.default_rng(0)).net) > 0.95
    assert weave(LAMELLIPODIUM, rng=np.random.default_rng(0)).branch_triples.shape[0] > 0
    assert np.isfinite(sarcomeric_period_um(weave(STRESS_FIBER, rng=np.random.default_rng(0)))["period_um"])
