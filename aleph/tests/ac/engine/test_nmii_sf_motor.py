"""Head exclusivity for the two NMII motor edges — derived from geometry, and what it refuses.

One myosin head engages ONE actin filament. With a cortex motor and an SF motor over the same head
array, nothing structural stops a head from binding both and applying its stall force twice. These
tests pin that the partition cannot produce that state, and that it refuses rather than guesses when
the inputs cannot support a partition at all.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.nmii_sf_motor import HeadPartition, partition_minifilaments_by_proximity


class _Heads:
    """Stand-in for the actuator's head→node map; only its length is read."""

    def __init__(self, n):
        self.shape = (n,)


def _two_cluster_geometry(n_mini=6, heads_per=4):
    """Half the minifilaments sit near x=0 (cortex), half near x=10 (sf) — an unambiguous split."""
    centres = np.array([[0.0, 0.0, 0.0]] * (n_mini // 2) + [[10.0, 0.0, 0.0]] * (n_mini // 2))
    heads = np.repeat(centres, heads_per, axis=0).astype(np.float64)
    heads[:, 1] += np.tile(np.arange(heads_per) * 1e-3, n_mini)   # distinct, but not enough to reorder
    return heads, np.array([[0.0, 0.0, 0.0]]), np.array([[10.0, 0.0, 0.0]])


def test_each_head_serves_exactly_one_target() -> None:
    heads, cortex, sf = _two_cluster_geometry()
    part = partition_minifilaments_by_proximity(
        head_node_d=_Heads(len(heads)), head_positions_um=heads,
        cortex_positions_um=cortex, sf_positions_um=sf, heads_per_minifilament=4,
    )
    assert part.cortex_heads.size + part.sf_heads.size == len(heads)
    assert not np.intersect1d(part.cortex_heads, part.sf_heads).size
    assert part.sf_fraction == pytest.approx(0.5)


def test_a_head_in_both_sets_is_REFUSED() -> None:
    """The constructor's own guard: a double-served head applies its stall force twice."""
    with pytest.raises(ValueError, match="BOTH targets"):
        HeadPartition(
            cortex_minifilaments=np.array([0]), sf_minifilaments=np.array([0]),
            cortex_heads=np.array([0, 1]), sf_heads=np.array([1, 2]),
            n_heads=4, tie_broken_to_cortex=0,
        )


def test_a_dropped_head_is_REFUSED() -> None:
    """A head in neither set is a motor that stops existing with nothing saying so."""
    with pytest.raises(ValueError, match="silently dropped"):
        HeadPartition(
            cortex_minifilaments=np.array([0]), sf_minifilaments=np.array([], dtype=np.int64),
            cortex_heads=np.array([0, 1]), sf_heads=np.array([], dtype=np.int64),
            n_heads=4, tie_broken_to_cortex=0,
        )


def test_the_partition_is_a_distance_ORDERING_so_rescaling_cannot_change_it() -> None:
    """There is no radius to choose here. Scale the whole geometry and the split must be identical."""
    heads, cortex, sf = _two_cluster_geometry()
    base = partition_minifilaments_by_proximity(
        head_node_d=_Heads(len(heads)), head_positions_um=heads,
        cortex_positions_um=cortex, sf_positions_um=sf, heads_per_minifilament=4)
    scaled = partition_minifilaments_by_proximity(
        head_node_d=_Heads(len(heads)), head_positions_um=heads * 1000.0,
        cortex_positions_um=cortex * 1000.0, sf_positions_um=sf * 1000.0, heads_per_minifilament=4)
    np.testing.assert_array_equal(base.sf_heads, scaled.sf_heads)


def test_a_tie_goes_to_the_cortex_and_is_COUNTED() -> None:
    """Equidistant is a real case; the rule is stated and the count is reported, not hidden.

    A large tie count would mean the two targets are not separated in space, so the partition is not
    carrying the physical meaning it claims — which the caller can only notice if it is reported.
    """
    heads = np.zeros((8, 3))
    part = partition_minifilaments_by_proximity(
        head_node_d=_Heads(8), head_positions_um=heads,
        cortex_positions_um=np.array([[1.0, 0.0, 0.0]]),
        sf_positions_um=np.array([[-1.0, 0.0, 0.0]]), heads_per_minifilament=4)
    assert part.tie_broken_to_cortex == 2
    assert part.sf_heads.size == 0 and part.cortex_heads.size == 8


def test_an_empty_target_is_refused_rather_than_handed_every_head() -> None:
    heads, cortex, _ = _two_cluster_geometry()
    with pytest.raises(ValueError, match="empty target"):
        partition_minifilaments_by_proximity(
            head_node_d=_Heads(len(heads)), head_positions_um=heads,
            cortex_positions_um=cortex, sf_positions_um=np.zeros((0, 3)), heads_per_minifilament=4)


def test_a_head_count_that_would_split_a_minifilament_is_refused() -> None:
    heads = np.zeros((10, 3))
    with pytest.raises(ValueError, match="whole number"):
        partition_minifilaments_by_proximity(
            head_node_d=_Heads(10), head_positions_um=heads,
            cortex_positions_um=np.array([[0.0, 0.0, 0.0]]),
            sf_positions_um=np.array([[1.0, 0.0, 0.0]]), heads_per_minifilament=4)


def test_heads_per_minifilament_must_be_read_not_assumed() -> None:
    heads = np.zeros((8, 3))
    with pytest.raises(ValueError, match="must be positive"):
        partition_minifilaments_by_proximity(
            head_node_d=_Heads(8), head_positions_um=heads,
            cortex_positions_um=np.array([[0.0, 0.0, 0.0]]),
            sf_positions_um=np.array([[1.0, 0.0, 0.0]]), heads_per_minifilament=0)


def test_the_chunked_nearest_search_agrees_with_a_dense_one() -> None:
    """The chunking exists so a native cortex does not materialise an n_mini x n_nodes matrix.

    It must not change the answer — a blocking bug here would silently mis-assign heads near the
    boundary between the two targets, which is exactly where the assignment matters.
    """
    rng = np.random.default_rng(3)
    heads = rng.normal(0.0, 5.0, (12 * 4, 3))
    cortex = rng.normal(-3.0, 1.0, (9000, 3))       # > the 4096 chunk, so several blocks are taken
    sf = rng.normal(+3.0, 1.0, (9000, 3))
    part = partition_minifilaments_by_proximity(
        head_node_d=_Heads(len(heads)), head_positions_um=heads,
        cortex_positions_um=cortex, sf_positions_um=sf, heads_per_minifilament=4)

    centres = heads.reshape(12, 4, 3).mean(axis=1)
    dense_sf = (
        ((centres[:, None, :] - sf[None, :, :]) ** 2).sum(2).min(1)
        < ((centres[:, None, :] - cortex[None, :, :]) ** 2).sum(2).min(1)
    )
    np.testing.assert_array_equal(np.flatnonzero(dense_sf), part.sf_minifilaments)


class _FakeSegId:
    def __init__(self, values):
        self._v = list(values)
        self.shape = (len(values),)
        self.device = "none"


class _FakeState:
    def __init__(self, seg_id):
        self.query_seg_id_d = seg_id


class _FakeInner:
    """Records that the real query ran, and exposes the scratch the wrapper vetoes."""

    def __init__(self, seg_id):
        self._state = _FakeState(seg_id)
        self.filled = 0
        self.some_other_attribute = "forwarded"

    def fill_attach_scratch(self, actuator, port) -> None:
        self.filled += 1


def test_the_mask_vetoes_through_the_gate_the_attach_kernel_ALREADY_honours() -> None:
    """`q_seg_id < 0` is the sentinel `segment_query` writes for "nothing in range".

    Reusing it is why exclusivity needs no change to the connector, the attach kernel or the
    crossbridge: a foreign head looks, to this edge, exactly like a head with no target.
    """
    from aleph.engine.nmii_sf_motor import HeadMaskedAttachQuery

    launched = {}
    inner = _FakeInner(_FakeSegId([7, 7, 7, 7]))

    class _Recording(HeadMaskedAttachQuery):
        pass

    query = _Recording(inner=inner, serves_d=object())
    import aleph.engine.nmii_sf_motor as mod
    real_launch = mod.wp.launch
    mod.wp.launch = lambda kernel, dim, inputs, device: launched.update(
        kernel=kernel, dim=dim, inputs=inputs)
    try:
        query.fill_attach_scratch(object(), object())
    finally:
        mod.wp.launch = real_launch

    assert inner.filled == 1, "the real query must still run; the wrapper only adds a veto"
    assert launched["kernel"] is mod._veto_foreign_heads_kernel
    assert launched["dim"] == 4


def test_the_wrapper_forwards_everything_else_untouched() -> None:
    from aleph.engine.nmii_sf_motor import HeadMaskedAttachQuery

    query = HeadMaskedAttachQuery(inner=_FakeInner(_FakeSegId([0])), serves_d=object())
    assert query.some_other_attribute == "forwarded"
