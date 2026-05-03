from __future__ import annotations

import pytest

from acs.v2.junction import JunctionState


def _junction(**overrides) -> JunctionState:
    base = dict(
        junction_id="j-1",
        cell_id_a="cell-1",
        cell_id_b="cell-2",
        contact_length_um=2.0,
        age_s=10.0,
        maturity=0.5,
    )
    base.update(overrides)
    return JunctionState(**base)


def test_minimal_junction_validates():
    junction = _junction()
    junction.validate()
    assert junction.state == "contacting"


def test_self_junction_rejected():
    junction = _junction(cell_id_b="cell-1")
    with pytest.raises(ValueError, match="itself"):
        junction.validate()


def test_non_canonical_pair_rejected():
    junction = _junction(cell_id_a="cell-2", cell_id_b="cell-1")
    with pytest.raises(ValueError, match="lexicographically"):
        junction.validate()


def test_from_cell_pair_canonicalizes():
    junction = JunctionState.from_cell_pair(
        junction_id="j-2",
        cell_id_pair=("zebra", "alpha"),
        contact_length_um=1.0,
        age_s=0.0,
        maturity=0.1,
    )
    junction.validate()
    assert junction.cell_id_a == "alpha"
    assert junction.cell_id_b == "zebra"


def test_negative_contact_length_rejected():
    junction = _junction(contact_length_um=-1.0)
    with pytest.raises(ValueError, match="contact_length_um"):
        junction.validate()


def test_negative_age_rejected():
    junction = _junction(age_s=-0.1)
    with pytest.raises(ValueError, match="age_s"):
        junction.validate()


def test_maturity_out_of_range_rejected():
    junction = _junction(maturity=1.5)
    with pytest.raises(ValueError, match="maturity"):
        junction.validate()


def test_invalid_state_rejected():
    junction = _junction(state="floating")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="state"):
        junction.validate()


def test_negative_edge_index_rejected():
    junction = _junction(contact_edge_indices_a=(-1,))
    with pytest.raises(ValueError, match="contact_edge_indices_a"):
        junction.validate()


def test_duplicate_edge_index_rejected():
    junction = _junction(contact_edge_indices_b=(2, 2))
    with pytest.raises(ValueError, match="contact_edge_indices_b"):
        junction.validate()


def test_cadherin_proxy_out_of_range_rejected():
    junction = _junction(cadherin_proxy=1.5)
    with pytest.raises(ValueError, match="cadherin_proxy"):
        junction.validate()


def test_contact_inhibition_out_of_range_rejected():
    junction = _junction(contact_inhibition_signal=-0.01)
    with pytest.raises(ValueError, match="contact_inhibition_signal"):
        junction.validate()


def test_tension_proxy_finite_required():
    junction = _junction(tension_proxy_nN=float("inf"))
    with pytest.raises(ValueError, match="tension_proxy_nN"):
        junction.validate()


def test_full_optional_fields_validate():
    junction = _junction(
        state="cortex_coupled",
        contact_edge_indices_a=(0, 1, 2),
        contact_edge_indices_b=(5, 6),
        cadherin_proxy=0.7,
        tension_proxy_nN=0.5,
        contact_inhibition_signal=0.3,
    )
    junction.validate()
