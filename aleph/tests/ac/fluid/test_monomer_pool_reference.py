"""NumPy gate: multi-consumer G-actin conserved pool (A_total invariant, separable, positive).

Pure-host acceptance oracle for the CHEMICAL_FLUX conserved-pool convention (spec §4).  Verifies the discrete
invariant the Warp-CUDA ``monomer_flux_commit_kernel`` ports to the device: ``A_total = integral phi c dV +
sum_c bound_c`` conserved to round-off across N consumers, with per-consumer separability (no shared pool).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.monomer_pool_reference import MonomerPoolConservationReference
from aleph.components.fluid.transport_reference import rad_cfl_dt


def _pool(consumers=("cortex", "sf_arc", "lamellipodium", "filopodium")):
    return MonomerPoolConservationReference(
        shape=(8, 8, 8), dx=0.5, phi=0.7, d_c=3.0, consumers=consumers
    )


def test_total_actin_conserved_across_consumers_and_transport() -> None:
    ref = _pool()
    rng = np.random.default_rng(0)
    u = np.abs(rng.random((8, 8, 8))) * 10.0  # free-monomer density
    a0 = ref.total_actin(u)
    dt = rad_cfl_dt(ref.dx, ref.d_c, v_max=0.2, dim=3) * 0.5
    v_f = np.zeros((8, 8, 8, 3))
    v_f[..., 0] = 0.1
    for _ in range(20):
        # each consumer draws (negative reaction) or returns (positive) a distinct amount
        reactions = {
            "cortex": -0.02 * u,
            "sf_arc": -0.01 * u,
            "lamellipodium": +0.005 * u,
            "filopodium": -0.003 * u,
        }
        u = ref.step(u, dt, per_consumer_reactions=reactions, v_f=v_f)
        ref.assert_conserved(u, a0, atol=1e-8)


def test_per_consumer_separability_no_shared_pool() -> None:
    ref = _pool(("cortex", "sf_arc"))
    u = np.full((8, 8, 8), 5.0)
    dt = 1e-3
    # only cortex draws; sf_arc's pool must be untouched
    ref.step(u, dt, per_consumer_reactions={"cortex": -np.ones_like(u)})
    assert ref.bound["cortex"] > 0.0
    assert ref.bound["sf_arc"] == 0.0


def test_balanced_treadmill_leaves_pools_unchanged() -> None:
    ref = _pool(("cortex",))
    u = np.full((8, 8, 8), 5.0)
    dt = 1e-3
    r = np.zeros((8, 8, 8))
    r[2, 2, 2] = 3.0      # pointed release
    r[5, 5, 5] = -3.0     # barbed consumption (balanced)
    u2 = ref.step(u, dt, per_consumer_reactions={"cortex": r})
    assert abs(ref.bound["cortex"]) < 1e-12          # net zero exchange
    assert abs(float(u2.sum()) - float(u.sum())) < 1e-9


def test_draw_against_unknown_consumer_is_rejected() -> None:
    ref = _pool(("cortex",))
    u = np.full((8, 8, 8), 5.0)
    with pytest.raises(KeyError, match="unknown monomer consumer"):
        ref.step(u, 1e-3, per_consumer_reactions={"filopodium": -np.ones_like(u)})


def test_duplicate_consumer_is_rejected() -> None:
    with pytest.raises(ValueError, match="unique"):
        MonomerPoolConservationReference(
            shape=(4, 4, 4), dx=0.5, phi=0.7, d_c=3.0, consumers=("cortex", "cortex")
        )
