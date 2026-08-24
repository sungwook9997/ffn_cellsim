"""Structural gates for the G-actin CHEMICAL_FLUX conserved-pool connector (whole-cell slice-3).

CUDA-free gates for ``ac/engine/monomer_flux.py``:

* the pure-host :class:`MonomerFluxAccount` conserved-pool discipline (draw/release, separability, finite pool),
* the ``ConnectorFamily.CHEMICAL_FLUX`` capability + the connector's declared endpoints/kinetics shape,
* the CUDA-residency guards on the device-lane owner/connector (rejected with metadata doubles on the Mac),
* the connector's full transaction + event API surface (the actor's kinetic-connector binding contract).

The device physics (the exchange kernels) run on the native lane; here the host account IS the accounting the
kernel mirrors, so the conservation invariant is genuinely exercised without CUDA.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from aleph.engine.contracts import ConnectorFamily
from aleph.engine.monomer_flux import (
    MonomerFluxAccount,
    MonomerFluxConnector,
    MonomerPoolOwner,
)


# --------------------------------------------------------------------------------------------------
# 1. MonomerFluxAccount — pure-host conserved pool
# --------------------------------------------------------------------------------------------------
def test_account_draw_release_conserves_a_total() -> None:
    acc = MonomerFluxAccount(field_monomer=100.0, consumers=("cortex", "sf_arc", "lamellipodium"))
    a0 = acc.a_total
    acc.draw("cortex", 10.0)
    acc.draw("sf_arc", 5.0)
    acc.release("cortex", 3.0)
    assert acc.bound["cortex"] == pytest.approx(7.0)
    assert acc.bound["sf_arc"] == pytest.approx(5.0)
    assert acc.field_monomer == pytest.approx(100.0 - 7.0 - 5.0)
    assert acc.a_total == pytest.approx(a0)
    acc.assert_conserved()


def test_account_consumer_pools_are_disjoint() -> None:
    acc = MonomerFluxAccount(field_monomer=50.0, consumers=("cortex", "filopodium"))
    acc.draw("cortex", 20.0)
    assert acc.bound["filopodium"] == 0.0  # cortex draw never touches filopodium's pool


def test_account_cannot_draw_more_than_free_pool() -> None:
    acc = MonomerFluxAccount(field_monomer=5.0, consumers=("cortex",))
    with pytest.raises(ValueError, match="exceeds the free-monomer pool"):
        acc.draw("cortex", 6.0)  # finite G-actin, never resampled


def test_account_cannot_release_more_than_bound() -> None:
    acc = MonomerFluxAccount(field_monomer=5.0, consumers=("cortex",))
    acc.draw("cortex", 2.0)
    with pytest.raises(ValueError, match="exceeds its bound pool"):
        acc.release("cortex", 3.0)


def test_account_unknown_consumer_is_rejected() -> None:
    acc = MonomerFluxAccount(field_monomer=5.0, consumers=("cortex",))
    with pytest.raises(KeyError, match="unknown monomer consumer"):
        acc.draw("sf_arc", 1.0)


def test_account_requires_at_least_one_consumer() -> None:
    with pytest.raises(ValueError, match="at least one"):
        MonomerFluxAccount(field_monomer=1.0, consumers=())


# --------------------------------------------------------------------------------------------------
# 2. CHEMICAL_FLUX family + device-lane guards (metadata doubles; no CUDA)
# --------------------------------------------------------------------------------------------------
def test_chemical_flux_family_exists() -> None:
    assert ConnectorFamily.CHEMICAL_FLUX.value == "chemical_flux"


@dataclass
class _Dev:
    is_cuda: bool


@dataclass
class _ArrDouble:
    dtype: object
    shape: tuple
    device: _Dev
    ptr: int


class _FieldDouble:
    """Minimal MonomerField stand-in exposing the attributes the owner probes."""

    def __init__(self) -> None:
        class _Grid:
            dx = 0.5
            dim = 3
            shape = (4, 4, 4)
            mask = object()
        self.grid = _Grid()
        self.u = object()


def _cuda_scalar(ptr: int, shape=(1,)):
    import warp as wp

    return _ArrDouble(dtype=wp.float64, shape=shape, device=_Dev(is_cuda=True), ptr=ptr)


def test_owner_rejects_non_cuda_pool() -> None:
    import warp as wp

    host = _ArrDouble(dtype=wp.float64, shape=(1,), device=_Dev(is_cuda=False), ptr=1)
    with pytest.raises(ValueError, match="CUDA device array"):
        MonomerPoolOwner(_FieldDouble(), field_pool_d=host, a_total_d=_cuda_scalar(2))


def test_owner_rejects_aliased_scalars() -> None:
    with pytest.raises(ValueError, match="distinct storage"):
        MonomerPoolOwner(_FieldDouble(), field_pool_d=_cuda_scalar(7), a_total_d=_cuda_scalar(7))


def _owner_double() -> MonomerPoolOwner:
    return MonomerPoolOwner(_FieldDouble(), field_pool_d=_cuda_scalar(10), a_total_d=_cuda_scalar(11))


def test_connector_rejects_bound_pool_aliasing_field_pool() -> None:
    owner = _owner_double()
    with pytest.raises(ValueError, match="aliases the Cytosol field pool"):
        MonomerFluxConnector(
            name="cortex_gactin_flux",
            component_a="cortex",
            pool=owner,
            consumer_bound_d=_cuda_scalar(10),   # SAME ptr as field_pool_d -> shared array
            drawn_candidate_d=_cuda_scalar(20),
            net_drawn_d=_cuda_scalar(21),
            rate_d=_cuda_scalar(22),
            snapshot_d=_cuda_scalar(23, shape=(3,)),
        )


def test_connector_declares_chemical_flux_family_and_cytosol_endpoint() -> None:
    owner = _owner_double()
    conn = MonomerFluxConnector(
        name="cortex_gactin_flux",
        component_a="cortex",
        pool=owner,
        consumer_bound_d=_cuda_scalar(30),
        drawn_candidate_d=_cuda_scalar(31),
        net_drawn_d=_cuda_scalar(32),
        rate_d=_cuda_scalar(33),
        snapshot_d=_cuda_scalar(34, shape=(3,)),
    )
    assert conn.family is ConnectorFamily.CHEMICAL_FLUX
    assert conn.component_b == "cytosol"
    # the connector must expose the full transaction + event + mechanics API the actor binds.
    for hook in (
        "propose_events", "snapshot_candidate", "rollback",
        "commit_irreversible", "accumulate_ledger", "accumulate",
    ):
        assert callable(getattr(conn, hook))


def test_connector_consumer_must_not_be_cytosol() -> None:
    owner = _owner_double()
    with pytest.raises(ValueError, match="distinct component"):
        MonomerFluxConnector(
            name="bad", component_a="cytosol", pool=owner,
            consumer_bound_d=_cuda_scalar(40), drawn_candidate_d=_cuda_scalar(41),
            net_drawn_d=_cuda_scalar(42), rate_d=_cuda_scalar(43), snapshot_d=_cuda_scalar(44, shape=(3,)),
        )
