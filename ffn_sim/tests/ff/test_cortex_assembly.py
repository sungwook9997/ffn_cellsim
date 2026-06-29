"""FF cortex assembly (ff/cortex_assembly) — Stage 6c-b.

Anchors: parameters match the H.3 production cortex (R=10µm, 1000 filaments, 7 beads, ℓ₀=0.5µm,
ℓ_p=17µm); fibers lie on the shell (radius residual ≪ ℓ₀); segment rest lengths ≈ ℓ₀; areal
density ≈ 0.8 µm⁻² (the cortex.py coverage check); the network topology is consistent
(nodes = F·beads, bending triples per fiber = beads−2); and the bending force on the assembled
cortex is finite and (since fibers are curved by the sphere) nonzero.
"""

import numpy as np
import pytest

from ffn_sim.ff import units as U
from ffn_sim.ff.cortex_assembly import CortexParams, build_cortex_network, equatorial_circumference_um
from ffn_sim.ff.forces_warp import bending_energy, bending_force


def test_default_params_match_h3_config():
    p = CortexParams()
    assert p.R_um == 10.0
    assert p.n_filaments == 1000
    assert p.beads_per_filament == 7
    assert p.seg_um == 0.5
    assert p.L_filament_um == pytest.approx(3.0)
    assert p.persistence_length_um == 17.0
    assert p.kappa == pytest.approx(U.KAPPA_ACTIN)
    assert p.areal_density_um2 == pytest.approx(1000 / (4 * np.pi * 100), rel=1e-9)
    assert 0.7 < p.areal_density_um2 < 0.85       # cortex.py coverage ≈ 0.8 µm⁻²


def test_small_cortex_geometry():
    """A small prototype cortex: fibers on the shell, seg ≈ ℓ₀, topology consistent."""
    net, meta = build_cortex_network(n_filaments=50, rng=np.random.default_rng(1))
    nb = 7
    assert net.n_fibers == 50
    assert net.n_nodes == 50 * nb
    assert net.bend_triples.shape[0] == 50 * (nb - 2)       # beads−2 triples/fiber
    # all model-points lie on the shell (radius residual tiny vs ℓ₀=0.5µm)
    assert meta["on_shell_residual_um"] < 1e-6
    # segment lengths ≈ ℓ₀ (great-circle chord of a 0.5µm arc on R=10µm)
    assert meta["seg_len_mean_um"] == pytest.approx(0.5, rel=2e-3)
    assert meta["seg_len_std_um"] < 1e-6


def test_bending_force_finite_and_curved():
    """The cortex fibers are curved by the sphere ⇒ nonzero bending energy + finite forces."""
    net, _ = build_cortex_network(n_filaments=20, rng=np.random.default_rng(2))
    F = bending_force(net)
    assert np.all(np.isfinite(F))
    assert bending_energy(net) > 0.0                        # curved on the shell → bent
    # whole-fiber net force ≈ 0 (each bending triple is internal torque only)
    off = net.fiber_offsets
    for f in range(net.n_fibers):
        sl = slice(int(off[f]), int(off[f + 1]))
        assert np.allclose(F[sl].sum(axis=0), 0.0, atol=1e-9)


def test_equatorial_circumference():
    assert equatorial_circumference_um(CortexParams()) == pytest.approx(2 * np.pi * 10.0)
