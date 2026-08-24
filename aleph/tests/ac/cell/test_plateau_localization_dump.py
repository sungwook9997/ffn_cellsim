"""CPU test for the GATE-A fine-mesh plateau residual-localization diagnostic (``--dump-plateau``).

Exercises the pure-numpy classifier + the npz dump writer of ``scripts/ac_gate_a_fq_coarse_test.py`` on a tiny
SYNTHETIC cortex (no CUDA, no ``build_cell``): a handful of fibers with ONE engineered node per structural
feature (steric contact / crosslink / ERM anchor / Arp2/3 branch / fiber boundary / interior). Asserts the
top-K are actually the K largest ``|PF|`` (self-consistency), that each engineered node lands in its expected
primary category by the documented priority, and that the writer round-trips a self-consistent npz + summary.

The production run is default-off: with no ``--dump-plateau`` flag the descent is byte-identical, so this test
only imports the diagnostic helpers — it never builds the CUDA cell.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

# The script is written for the gbook run convention (``ac``/``scripts`` top-level, PYTHONPATH=…/ffn_sim). Put
# ffn_sim on the path and load it by file so the CPU test can import the pure-numpy classifier without CUDA.
FFN_SIM = Path(__file__).resolve().parents[3]
if str(FFN_SIM) not in sys.path:
    sys.path.insert(0, str(FFN_SIM))
_spec = importlib.util.spec_from_file_location(
    "ac_gate_a_fq_coarse_test", FFN_SIM / "scripts" / "ac_gate_a_fq_coarse_test.py")
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


class _Arr:
    """Minimal warp-array stand-in: ``.numpy()`` returns the held ndarray (CPU-only fake cell)."""

    def __init__(self, a):
        self._a = np.asarray(a)

    def numpy(self):
        return self._a


def _synthetic_cortex():
    """6 fibers × 10 nodes with one engineered node per feature. Returns (arrays, expected primary per top node)."""
    n_fib, per = 6, 10
    n = n_fib * per
    fiber_offsets = np.arange(0, n + 1, per, dtype=np.int64)
    node_fiber = np.repeat(np.arange(n_fib), per).astype(np.int64)

    # positions: each fiber a separate line, fibers >> WCA cutoff apart so there are NO accidental steric pairs
    pos = np.zeros((n, 3), np.float64)
    for f in range(n_fib):
        for kk in range(per):
            pos[f * per + kk] = (0.1 * kk, 2.0 * f, 0.0)

    # engineered feature nodes (all INTERIOR of their fiber except the boundary one)
    S1, XL, ERM, BR, BND, INT = 3, 24, 34, 44, 50, 5   # 50 = fiber5 start ⇒ a fiber-boundary node
    S2 = 13                                             # fiber1 partner placed inside the WCA cutoff of S1
    pos[S2] = pos[S1] + np.array([0.003, 0.0, 0.0])     # 3 nm < cutoff ≈ 7.85 nm ⇒ a cross-fiber steric contact

    xl_pairs = np.array([[XL, 26]], np.int64)           # crosslink touches node XL
    erm_cortex_idx = np.array([ERM], np.int64)
    branch_triples = np.array([[43, BR, 45]], np.int64)  # Arp2/3 junction touches node BR

    pf_mag = np.full(n, 0.05, np.float64)
    pf_mag[S1] = 5.0
    pf_mag[XL] = 4.0
    pf_mag[ERM] = 3.0
    pf_mag[BR] = 2.0
    pf_mag[BND] = 1.5
    pf_mag[INT] = 1.2                                    # a plain interior node, the 6th-worst
    expected = {S1: "steric", XL: "crosslink", ERM: "erm", BR: "branch", BND: "boundary", INT: "interior"}
    return dict(pos=pos, pf_mag=pf_mag, node_fiber=node_fiber, xl_pairs=xl_pairs,
                erm_cortex_idx=erm_cortex_idx, branch_triples=branch_triples,
                fiber_offsets=fiber_offsets), expected


def test_classify_topk_self_consistent_and_feature_labels():
    arrs, expected = _synthetic_cortex()
    sigma = 0.007
    k = 6
    res = mod.classify_plateau_nodes(
        arrs["pos"], arrs["pf_mag"], arrs["node_fiber"],
        xl_pairs=arrs["xl_pairs"], erm_cortex_idx=arrs["erm_cortex_idx"],
        branch_triples=arrs["branch_triples"], fiber_offsets=arrs["fiber_offsets"],
        sigma_ev_um=sigma, k=k)

    # (1) top-K are ACTUALLY the K largest |PF|, in descending order
    ref = np.argsort(arrs["pf_mag"])[::-1][:k]
    assert np.array_equal(res["topk_idx"], ref), (res["topk_idx"], ref)
    assert np.all(np.diff(res["topk_pf"]) <= 0.0)                 # descending |PF|
    assert float(res["topk_pf"][0]) == float(arrs["pf_mag"].max())

    # (2) each engineered node lands in its expected PRIMARY category (documented priority)
    primary = {int(i): str(p) for i, p in zip(res["topk_idx"], res["topk_primary"], strict=True)}
    for node, cat in expected.items():
        assert primary[node] == cat, (node, primary[node], cat)

    # (3) breakdown is a partition of the top-K
    assert sum(res["breakdown"].values()) == k
    for cat in ("steric", "crosslink", "erm", "branch", "boundary", "interior"):
        assert res["breakdown"][cat] == 1

    # (4) deciles are a well-formed, monotone 11-point summary of the |PF| distribution
    assert res["deciles"].shape == (11,)
    assert np.all(np.diff(res["deciles"]) >= -1e-12)
    assert float(res["deciles"][-1]) == float(arrs["pf_mag"].max())

    # (5) spatial-spread + summary fields present and finite
    assert np.isfinite(res["r_gyration"]) and res["r_gyration"] >= 0.0
    assert np.isfinite(res["bbox_diag"]) and res["bbox_diag"] >= 0.0
    assert 0.0 <= res["frac_within_0p5um_of_worst"] <= 1.0
    assert "top-6" in res["summary_text"] and "feature breakdown" in res["summary_text"]


def _synthetic_term_forces(n=40, seed=0):
    """Build per-term (N,3) force arrays whose vector sum defines the total F_raw, plus a PF that differs from
    F_raw by a nonzero constraint (Jᵀλ = F_raw − PF). Returns (order, term_forces, raw_vec, pf_vec)."""
    rng = np.random.default_rng(seed)
    order = ["bending", "crosslink", "erm", "steric", "pressure"]
    term_forces = {nm: rng.standard_normal((n, 3)) for nm in order}
    raw_vec = np.zeros((n, 3), np.float64)
    for nm in order:
        raw_vec += term_forces[nm]
    # a constraint force removed by the projection: make one node dominated by a big pull the projector cancels
    constraint = rng.standard_normal((n, 3)) * 0.1
    constraint[0] = np.array([20.0, 0.0, 0.0])           # engineered worst node: a constraint force (Jᵀλ) that
    # dominates every term (~1-2 pN) and the ~4 pN raw total, so this node is unambiguously the worst by |PF|
    pf_vec = raw_vec - constraint
    return order, term_forces, raw_vec, pf_vec


def test_force_decomposition_sums_to_total_and_reports_dominant():
    """Σ_terms F_term == F_raw to round-off (complete + faithful), constraint = F_raw − PF, magnitudes exact."""
    n = 40
    order, term_forces, raw_vec, pf_vec = _synthetic_term_forces(n=n)
    # rank by |PF| so topk_idx matches how the real dump ranks the plateau nodes
    pf_mag = np.linalg.norm(pf_vec, axis=1)
    topk_idx = np.argsort(pf_mag)[::-1][:10].astype(np.int64)

    res = mod.summarize_force_decomposition(order, term_forces, raw_vec, pf_vec, topk_idx)

    # (1) completeness: the isolated terms sum back to the total assembled force to round-off
    assert res["completeness"] < 1e-9, res["completeness"]

    # (2) if a term is DROPPED, the identity must break (proves the check actually detects incompleteness)
    dropped = {nm: term_forces[nm] for nm in order if nm != "steric"}
    res_bad = mod.summarize_force_decomposition(
        [nm for nm in order if nm != "steric"], dropped, raw_vec, pf_vec, topk_idx)
    assert res_bad["completeness"] > 1e-6

    # (3) per-node magnitudes are exactly the vector norms of each isolated term
    col = list(res["col_names"])
    for nm in order:
        exp = np.linalg.norm(term_forces[nm], axis=1)
        got = res["topk_table"][:, col.index(nm)]
        assert np.allclose(got, exp[topk_idx]), nm

    # (4) the constraint column equals |F_raw − PF|
    con = res["topk_table"][:, col.index("JTlambda")]
    assert np.allclose(con, np.linalg.norm(raw_vec - pf_vec, axis=1)[topk_idx])

    # (5) totals columns are correct
    assert np.allclose(res["topk_table"][:, col.index("raw_total")], np.linalg.norm(raw_vec, axis=1)[topk_idx])
    assert np.allclose(res["topk_table"][:, col.index("PF_total")], np.linalg.norm(pf_vec, axis=1)[topk_idx])

    # (6) the engineered worst node (id 0, |PF| dominated by the constraint pull) is ranked first and its
    #     dominant physical term is correctly reported as the constraint (JTlambda)
    assert int(topk_idx[0]) == 0
    assert res["worst_dominant"] == "JTlambda"
    assert "DOMINANT term at worst node = JTlambda" in res["summary_text"]


def test_dump_writer_roundtrip(tmp_path):
    """The npz writer emits every array + the top-K stay the K largest |PF| after a reload (self-consistent)."""
    arrs, _ = _synthetic_cortex()
    n = arrs["pos"].shape[0]
    from types import SimpleNamespace
    fake = SimpleNamespace(
        n_actin=n,
        pos_d=_Arr(arrs["pos"]),
        node_fiber_d=_Arr(arrs["node_fiber"].astype(np.int32)),
        xl_d=_Arr(arrs["xl_pairs"].astype(np.int32)),
        n_xl=int(arrs["xl_pairs"].shape[0]),
        branch_triples_d=_Arr(arrs["branch_triples"].astype(np.int32)),
        foff_d=_Arr(arrs["fiber_offsets"].astype(np.int32)),
        membrane=None,
    )
    projected = _Arr(np.concatenate([arrs["pf_mag"][:, None] * np.array([[1.0, 0.0, 0.0]]), ], axis=0))
    out = tmp_path / "plateau.npz"
    mod._dump_plateau_from_cell(fake, projected, str(out), label="TEST", k=6)

    assert out.exists()
    d = np.load(out, allow_pickle=True)
    for key in ("pos", "pf_mag", "node_fiber", "topk_idx", "topk_pf", "topk_pos",
                "topk_primary", "feature_names", "deciles", "summary_text"):
        assert key in d.files, key
    # reload self-consistency: the dumped |PF| reproduces the top-K ranking
    pf = d["pf_mag"]
    ref = np.argsort(pf)[::-1][:6]
    assert np.array_equal(d["topk_idx"], ref)
    assert pf.shape[0] == n
    assert str(d["label"]) == "TEST"
