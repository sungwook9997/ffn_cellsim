r"""CPU tests for the MIXED formin + Arp2/3 cortex sub-population (FLAGGED prototype, 2026-07-24).

Native cortex is ~2/3 formin (long) + ~1/3 Arp2/3-nucleated (short, branched) BY MASS (Bovellan 2014,
KB-3.18); the default cortex is formin-only. This asserts the OPT-IN mixed-cortex split — WITHOUT a CUDA
``build_cell`` (guarded out; region/weave construction is pure host NumPy):

  1. **Default = bit-parity.** ``cortex_arp23_fraction == 0.0`` ⇒ ``_cortex_regions`` returns the SINGLE
     historical formin cortex region, unchanged.
  2. **Mass→count derivation.** ``derive_cortex_arp23_split`` turns the Bovellan ~1/3 MASS fraction into a
     count split that PRESERVES the total budget (areal density) and matches the closed-form arithmetic —
     including the task's illustrative "5× more Arp2/3 by count" case.
  3. **The Arp2/3 sub-population carries branches.** The built region emits ``branch_triples`` (a valid
     mother→branch→daughter tree) whose junction angles cluster at the Arp2/3 θ₀ (70°), and it is REUSING the
     shared ``branch_angle`` path. A small ``weave_cell`` integration confirms the two regions concat into one
     network with the derived counts + branches.

CUDA is guarded: ``build_cell`` is never called. The math + host region builders run on the dev Mac.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.assemble import CellConfig, R_CORTEX_UM, _cortex_region, _cortex_regions
from aleph.components.weave.branch_angle import ARP23_THETA0_RAD, branch_angle
from aleph.components.weave.regions import (
    build_region,
    cortex_arp23_region,
    derive_cortex_arp23_split,
)
from aleph.components.weave.woven_cell import weave_cell

N_FIL = 70686  # FULL native cortical F-actin budget (the production baseline)


# ── 1. default (fraction 0) is bit-identical to today's single formin cortex ───────────────────────


def test_default_config_has_arp23_off():
    cfg = CellConfig()
    assert cfg.cortex_arp23_fraction == 0.0  # OFF by default → pure-formin cortex


def test_fraction_zero_returns_single_legacy_region():
    """cortex_arp23_fraction=0.0 ⇒ exactly one region == the historical _cortex_region call (bit-parity)."""
    cfg = CellConfig()
    regions = _cortex_regions(cfg)
    assert len(regions) == 1
    legacy = _cortex_region(
        cfg.n_filaments, seg_um=cfg.cortex_seg_um, length_um=cfg.cortex_length_um,
        density_per_fil=cfg.cortex_density_per_fil)
    assert regions[0].arch == legacy.arch  # frozen dataclass value-equality ⇒ identical build
    assert regions[0].arch.filament.n_filaments == N_FIL
    assert regions[0].arch.filament.nucleator == "formin"


# ── 2. mass → count derivation (preserves areal density) ───────────────────────────────────────────


def test_split_preserves_total_budget():
    """The split preserves the cortex filament BUDGET exactly ⇒ areal density (~100/µm²) unchanged."""
    n_formin, n_arp23 = derive_cortex_arp23_split(N_FIL, 0.33, formin_length_um=3.0, arp23_length_um=0.15)
    assert n_formin + n_arp23 == N_FIL
    assert n_formin >= 1 and n_arp23 >= 1


def test_split_recovers_the_mass_fraction():
    """Round-trip: the derived counts reproduce the requested Arp2/3 MASS fraction (mass ∝ n·length)."""
    L_f, L_a, f = 3.0, 0.15, 0.33
    n_formin, n_arp23 = derive_cortex_arp23_split(N_FIL, f, L_f, L_a)
    arp23_mass = n_arp23 * L_a
    formin_mass = n_formin * L_f
    assert arp23_mass / (arp23_mass + formin_mass) == pytest.approx(f, abs=1e-3)


def test_illustrative_five_x_count_case():
    """The task's illustrative numbers: L_a=0.1, L_f=1.0, f=1/3 ⇒ ~5× more Arp2/3 filaments BY COUNT."""
    n_formin, n_arp23 = derive_cortex_arp23_split(N_FIL, 1.0 / 3.0, formin_length_um=1.0, arp23_length_um=0.1)
    assert n_arp23 / n_formin == pytest.approx(5.0, rel=1e-3)
    # count fraction (≈0.83) is FAR above the 1/3 mass fraction — the whole point of the mass→count map
    assert n_arp23 / N_FIL == pytest.approx(5.0 / 6.0, rel=1e-3)


def test_baseline_default_lengths_split():
    """The CellConfig defaults (L_f=3.0, L_a=0.15, f=0.33) give the documented ~9.85× count ratio."""
    cfg = CellConfig(cortex_arp23_fraction=0.33)
    n_formin, n_arp23 = derive_cortex_arp23_split(
        cfg.n_filaments, cfg.cortex_arp23_fraction, cfg.cortex_length_um, cfg.cortex_arp23_length_um)
    assert (n_formin, n_arp23) == (6514, 64172)
    assert n_arp23 / n_formin == pytest.approx(9.85, rel=1e-2)


def test_split_rejects_bad_inputs():
    with pytest.raises(ValueError):
        derive_cortex_arp23_split(N_FIL, 0.0, 3.0, 0.15)   # fraction must be in (0,1)
    with pytest.raises(ValueError):
        derive_cortex_arp23_split(N_FIL, 1.0, 3.0, 0.15)
    with pytest.raises(ValueError):
        derive_cortex_arp23_split(N_FIL, 0.33, 3.0, 0.0)   # length must be > 0


# ── 3. _cortex_regions produces two regions with the derived split ─────────────────────────────────


def test_mixed_config_returns_two_regions_with_derived_counts():
    cfg = CellConfig(cortex_arp23_fraction=0.33)
    regions = _cortex_regions(cfg)
    assert len(regions) == 2
    formin, arp23 = regions
    assert formin.arch.filament.nucleator == "formin"
    assert arp23.arch.filament.nucleator == "arp23"
    assert arp23.arch.manifold == "sphere_dendritic"
    # total Arp2/3 count = mothers (arch.n_filaments) + daughters (n_daughter_pool); total budget preserved
    n_arp23 = arp23.arch.filament.n_filaments + arp23.n_daughter_pool
    n_formin = formin.arch.filament.n_filaments
    assert n_formin + n_arp23 == N_FIL
    _, n_arp23_expected = derive_cortex_arp23_split(N_FIL, 0.33, cfg.cortex_length_um, cfg.cortex_arp23_length_um)
    assert n_arp23 == n_arp23_expected
    # the Arp2/3 shell sits just inside the membrane (same shell as the formin cortex)
    assert arp23.arch.R_um == R_CORTEX_UM == formin.arch.R_um


# ── 4. the built Arp2/3 region carries a valid branched topology at θ₀ ─────────────────────────────


def _build_small_arp23_leaf(n_arp23=240, mother_fraction=0.2, length_um=0.15, seg_um=0.05, seed=0):
    spec = cortex_arp23_region(
        n_arp23, length_um=length_um, seg_um=seg_um, mother_fraction=mother_fraction, R_um=R_CORTEX_UM)
    return spec, build_region(spec, np.random.default_rng(seed))


def test_arp23_region_carries_branches():
    """The Arp2/3 sub-population is genuinely BRANCHED (not just short formin): one branch per daughter."""
    spec, leaf = _build_small_arp23_leaf()
    n_mothers = spec.arch.filament.n_filaments
    n_daughters = spec.n_daughter_pool
    n_fib = leaf.fiber_offsets.shape[0] - 1
    assert n_fib == n_mothers + n_daughters
    # every daughter carries exactly one branch junction (a mother→branch→daughter tree)
    assert leaf.branch_triples.shape[0] == n_daughters
    assert leaf.branch_anchors.shape[0] == n_daughters
    assert bool(leaf.branch_active.all())
    # triples reference valid nodes
    assert leaf.branch_triples.min() >= 0
    assert leaf.branch_triples.max() < leaf.pos.shape[0]


def test_arp23_branch_angles_cluster_at_theta0():
    """Measured junction angles cluster at the Arp2/3 θ₀ (70°) — the REUSED angle-harmonic branch path."""
    _, leaf = _build_small_arp23_leaf(n_arp23=400)
    angles = np.array([
        branch_angle(leaf.pos[i], leaf.pos[j], leaf.pos[k]) for (i, j, k) in leaf.branch_triples
    ])
    # mean within a few degrees of θ₀ (the solid-angle measure shifts the peak slightly above θ₀)
    assert abs(np.degrees(angles.mean()) - np.degrees(ARP23_THETA0_RAD)) < 8.0


def test_arp23_region_rejects_unbranchable_discretization():
    """A too-coarse seg (< 3 nodes/filament) cannot carry a branch → explicit error, not a silent rod."""
    spec = cortex_arp23_region(100, length_um=0.1, seg_um=0.1, mother_fraction=0.2, R_um=R_CORTEX_UM)
    with pytest.raises(ValueError):
        build_region(spec, np.random.default_rng(0))


# ── 5. weave_cell integration: two regions concat into ONE network with the branches ──────────────


def test_weave_cell_integrates_formin_and_arp23():
    """Small end-to-end weave: formin sphere cortex + Arp2/3 branched region → one woven network."""
    n_total = 600
    n_formin, n_arp23 = derive_cortex_arp23_split(n_total, 0.33, 3.0, 0.15)
    formin = _cortex_region(n_formin, seg_um=0.5, length_um=3.0, density_per_fil=20.0)
    arp23 = cortex_arp23_region(
        n_arp23, length_um=0.15, seg_um=0.05, mother_fraction=0.2, R_um=R_CORTEX_UM)
    wc = weave_cell([formin, arp23], rng=np.random.default_rng(0))  # overlap_free default False (fast)

    assert wc.n_fibers == n_formin + n_arp23
    # branches come ONLY from the Arp2/3 region and are present in the unified network
    assert wc.branch_triples.shape[0] > 0
    assert wc.branch_triples.max() < wc.n_nodes
    ledger = wc.population_ledger()
    assert ledger["N_unique_active_fibers"] == n_formin + n_arp23
    assert ledger["per_region"]["cortex"]["total_fibers"] == n_formin
    assert ledger["per_region"]["cortex_arp23"]["total_fibers"] == n_arp23
    assert ledger["N_branch_junctions_active"] == wc.branch_triples.shape[0]


def test_weave_cell_fraction_zero_is_single_region():
    """weave_cell of the fraction-0 cortex list has NO branches (pure formin) — the parity baseline."""
    formin = _cortex_region(300, seg_um=0.5, length_um=3.0, density_per_fil=20.0)
    wc = weave_cell([formin], rng=np.random.default_rng(0))
    assert wc.branch_triples.shape[0] == 0
    assert wc.n_fibers == 300
