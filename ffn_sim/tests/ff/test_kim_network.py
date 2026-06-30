"""FF actin network vs Taeyoon Kim's BD model (ff/kim_network) — Stage 6i peer-model validation.

The single-cell FF network's correct oracle is a peer fine-grained model (Kim 2007 MIT thesis), NOT
the multi-cell spheroid law. Anchors: the semiflexible mesh ξ(151µM)≈64nm (Kim/actin lit ~50-100nm);
the FF box network reproduces the ξ ∝ C_A^(−1/2) scaling (Schmidt 1989 / MacKintosh / Kim); and the
network connectivity scales as z = 2R with the crosslinker ratio.
"""

import numpy as np
import pytest

from ffn_sim.ff.kim_network import (
    build_box_network,
    connectivity_z,
    length_density_per_um2,
    measure_mesh_size_um,
    mesh_size_theory_um,
)


def test_mesh_theory_matches_actin_literature():
    xi = mesh_size_theory_um(151.0) * 1e3      # nm
    assert 50.0 < xi < 100.0                    # Kim cortex mesh ~50-100 nm; analytic ~64 nm
    # ρ_L scales linearly in C_A → ξ = ρ_L^(-1/2) ∝ C_A^(-1/2)
    assert length_density_per_um2(302.0) == pytest.approx(2 * length_density_per_um2(151.0))


def test_ff_box_reproduces_mesh_scaling():
    """The FF 3D box actin network reproduces ξ ∝ C_A^(−1/2) (Kim / Schmidt-MacKintosh)."""
    Cs = [75.0, 150.0, 300.0, 600.0, 1200.0]
    xi = []
    for C in Cs:
        net, xl, meta = build_box_network(C, box_um=1.0, rng=np.random.default_rng(0))
        xi.append(measure_mesh_size_um(net, box_um=1.0))
    slope = np.polyfit(np.log(Cs), np.log(xi), 1)[0]
    assert slope == pytest.approx(-0.5, abs=0.06)   # measured exponent ≈ −0.5


def test_connectivity_scales_with_crosslinker_ratio():
    """Mean crosslinks/filament z = 2R (connectivity linear in the ACP ratio R)."""
    for R in (0.1, 0.3, 0.5, 1.0):
        net, xl, meta = build_box_network(151.0, box_um=1.0, R_acp=R, rng=np.random.default_rng(0))
        assert connectivity_z(net, xl) == pytest.approx(2 * R, abs=0.05)
