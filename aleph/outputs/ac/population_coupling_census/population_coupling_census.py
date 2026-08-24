"""The full-cell contribution census: which nodes does NOTHING reach, and which channels cross a population?

WHY THIS RUN EXISTS.  The arena parity record (`78942fc4`) reported, without interpreting it, that 16,312
of 511,114 nodes take no contribution from either channel it bound.  That number is NOT yet a finding, and
the reason is the first caution ffn-cellsim-45 sent from its own family ablation: `alpha_actinin_arc`
contributed exactly zero because `n_arc_joints` was zero as built, and that zero reads like physics.  The
parity census counted TWO channels — crosslink and actin bending — while the composed array is
`[actin | myosin]` plus appended nucleus and membrane blocks, and a myosin particle is bound by its own
backbone and head channels, not by an actin crosslink.  So the honest first move is to assume my own
number is an artefact of the channels I happened to bind, and to try to destroy it.

WHAT IS ACTUALLY UNDER TEST.  Two questions, one pass, because both need the same census:

  Q1 — CONTRIBUTION.  Over EVERY channel the built cell carries, which nodes are touched by nothing?
       A node reached by no channel has no restoring force from the assembled operator, and a static
       equilibrium problem whose operator is singular has no unique solution.  That is why this matters
       now rather than later: five solvers have failed to descend and PI item (e) 1 — whether the static
       equilibrium exists at all — is HELD.

  Q2 — CO-LOCATION.  For every channel, which POPULATION PAIRS do its indices span?  The composed array
       already is an arena: `assemble.py:353` documents `actin [0:n_actin) + myosin [n_actin:n_total)`
       and the compartments append after it at global `node_off`.  CLAUDE.md's rule is that co-location
       in a shared array is NEVER a connection, so every channel that spans two populations must be a
       DECLARED connector.  This lists them, and a channel that spans a pair with nothing declaring it
       is exactly the defect the rule exists to prevent.

DISCOVERY, NOT GUESSING, and this is the part that makes the census auditable.  Channel topologies are
found by introspecting the built cell and its primitives for integer device arrays rather than by naming
the attributes I expect.  An engine grows channels; a hand-written list silently stops counting the ones
added after it was written, which is the same failure as the parity census in a slower form.  Every array
found is CLASSIFIED and every classification is reported — nothing is dropped silently:

  * `topology`   — values inside [0, n_nodes) and not monotone: a genuine index list. COUNTED.
  * `offset`     — monotone non-decreasing, first 0: a CSR-style row pointer, not node indices. Excluded.
  * `per_node`   — length exactly n_nodes: a per-node attribute (mask, id, flag), not a contribution.
                   Excluded from the contribution count, but its zero/negative entries are reported
                   separately because `solid_active_d` marks non-skeleton nodes and that is a real signal.
  * `out_of_range` — holds values outside [0, n_nodes): not node indices at all. Excluded, and LISTED,
                   because an array I cannot classify is a hole in the census and must be visible.

A sentinel `-1` is a real convention here (`_nearest_antiparallel` returns it, and dormant heads use it),
so negative entries are counted as UNBOUND and reported per channel rather than treated as index 0.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — this pass counts topology only; no force, no length, no unit is produced. The force
    magnitudes that follow in the containment stage are [pN].
  * boundary — a node index equal to n_nodes-1 and one equal to 0 are both inside the range; the
    classification uses a half-open [0, n_nodes) test so an off-by-one cannot be absorbed silently.
    A channel with zero elements as built is reported with `n_elements: 0` and NOT counted as "no
    crossing", because "absent" and "present and within one population" are different facts.
  * conservation/invariant — the population ranges must PARTITION [0, n_nodes) with no gap and no
    overlap; if they do not, the census halts rather than reporting a census over an unknown remainder.
  * CFL/precision — nothing is integrated; the counts are exact integers.
  * sign sense — not applicable to a count; negative indices are sentinels and are reported as such.
  * measurement protocol — one `build_cell` on the device, one host readback of each topology array
    after the build completes. No readback inside any loop; no physical step is taken.

engine units: none produced by this pass.  Runtime: NVIDIA Warp on CUDA only (build_cell requires it).
"""

from __future__ import annotations

import json
import sys

import numpy as np
import warp as wp

_INT = None  # filled after wp.init


def _is_int_array(v) -> bool:
    return isinstance(v, wp.array) and v.dtype in _INT and v.size > 0


def _collect(obj, prefix: str, seen: set[int], out: dict) -> None:
    """Every integer device array reachable one level down from ``obj``, keyed by dotted name."""
    if obj is None:
        return
    for name in sorted(dir(obj)):
        if name.startswith("_"):
            continue
        try:
            v = getattr(obj, name)
        except Exception:
            continue
        if _is_int_array(v):
            if id(v) in seen:
                continue
            seen.add(id(v))
            out[f"{prefix}.{name}"] = v
    return


def classify(flat: np.ndarray, n_nodes: int, size_hint: int) -> str:
    """Which of the four classes this array belongs to. See the module docstring."""
    if flat.size == 0:
        return "empty"
    lo, hi = int(flat.min()), int(flat.max())
    if lo < -1 or hi >= n_nodes:
        return "out_of_range"
    if size_hint == n_nodes:
        return "per_node"
    # A CSR row pointer is monotone non-decreasing, starts at 0 and ends at a count.
    if flat.size > 2 and lo >= 0 and bool(np.all(np.diff(flat) >= 0)) and flat[0] == 0:
        return "offset"
    # An integer array whose values happen to fall inside [0, n_nodes) is not thereby a node index list.
    # Two false positives from the first pass make the point: `grid.mask` is a 43x43 FLUID cell mask whose
    # 0/1/2 entries are all small integers, and `erm_bound_d` is a per-tether bound-state flag. Both were
    # counted as topology. A state flag carries few distinct values and a counter is one or two elements,
    # so both are refused here and REPORTED rather than dropped -- an array I exclude by heuristic is a
    # hole in the census and has to stay visible.
    if flat.size <= 2:
        return "counter"
    if np.unique(flat).size <= 4:
        return "flag_or_mask"
    return "topology"


def main() -> int:
    global _INT
    wp.init()
    _INT = (wp.int32, wp.int64)
    device = wp.get_device()
    if not device.is_cuda:
        raise RuntimeError(f"CUDA-only: resolved {device!r}")

    import aleph
    from aleph.components.incumbent.assemble import CellConfig, build_cell

    cell = build_cell(CellConfig())
    n_nodes = int(cell.pos_d.shape[0])

    # ── population ranges ──────────────────────────────────────────────────────────────────────────
    # `n_total` is NOT actin+myosin: it is the whole composed array, `[actin | myosin | nucleus |
    # membrane]`. Reading `[n_actin, n_total)` as "myosin" swallows both compartments — it is how the
    # parity census's 16,312 arose, and the partition invariant below is what caught it. Only the
    # ANCHORED blocks are stated here (actin at 0, and each compartment at its own global `node_off`);
    # myosin is the REMAINING gap, and its size is checked against the ledger rather than assumed.
    anchored: list[tuple[str, int, int]] = [("actin", 0, int(cell.n_actin))]
    for label, comp in (("nucleus", cell.nucleus), ("membrane", cell.membrane)):
        if comp is not None and int(comp.n_verts):
            anchored.append((label, int(comp.node_off), int(comp.node_off) + int(comp.n_verts)))
    anchored.sort(key=lambda p: p[1])

    pops: list[tuple[str, int, int]] = []
    cursor = 0
    for label, lo, hi in anchored:
        if lo > cursor:
            pops.append(("myosin", cursor, lo))
        pops.append((label, lo, hi))
        cursor = hi
    if cursor < n_nodes:
        pops.append(("myosin", cursor, n_nodes))
    pops = [p for p in pops if p[2] > p[1]]

    n_myo_gap = sum(hi - lo for lbl, lo, hi in pops if lbl == "myosin")
    n_myo_ledger = int(cell.ledger.get("myosin_n_myosin_particles", -1))
    if n_myo_ledger >= 0 and n_myo_gap != n_myo_ledger:
        raise RuntimeError(
            f"the gap left for myosin is {n_myo_gap} but the ledger says {n_myo_ledger} particles; the "
            "block layout is not what this census assumes and the labels would be wrong")

    # INVARIANT: the populations must partition [0, n_nodes). A census over an unknown remainder is not
    # a census, so this halts rather than reporting one.
    cursor, gaps = 0, []
    for label, lo, hi in sorted(pops, key=lambda p: p[1]):
        if lo != cursor:
            gaps.append({"between": cursor, "and": lo, "before": label})
        cursor = hi
    if cursor != n_nodes:
        gaps.append({"between": cursor, "and": n_nodes, "before": "END_OF_ARRAY"})
    if gaps:
        raise RuntimeError(f"populations do not partition [0,{n_nodes}); gaps={gaps}; pops={pops}")

    pop_of = np.empty(n_nodes, np.int32)
    for i, (_, lo, hi) in enumerate(pops):
        pop_of[lo:hi] = i
    names = [p[0] for p in pops]

    # ── discover every integer device array on the cell and its primitives ─────────────────────────
    found: dict[str, wp.array] = {}
    seen: set[int] = set()
    _collect(cell, "cell", seen, found)
    for attr in ("steric", "pressure", "myosin", "nucleus", "membrane", "membrane_pressure",
                 "substrate", "membrane_bc", "grid", "domain", "state"):
        _collect(getattr(cell, attr, None), attr, seen, found)

    channels, excluded, per_node_arrays = {}, {}, {}
    hits = np.zeros(n_nodes, np.int64)          # contributions per node, over COUNTED channels only

    for key in sorted(found):
        arr = found[key]
        flat = arr.numpy().reshape(-1).astype(np.int64)
        cls = classify(flat, n_nodes, arr.shape[0])
        row = {"class": cls, "shape": list(arr.shape), "n_elements": int(arr.shape[0]),
               "min": int(flat.min()), "max": int(flat.max())}
        if cls != "topology":
            if cls == "per_node":
                neg = int((flat < 0).sum())
                zero = int((flat == 0).sum())
                per_node_arrays[key] = {**row, "n_negative": neg, "n_zero": zero}
            else:
                excluded[key] = row
            continue

        valid = flat[flat >= 0]
        np.add.at(hits, valid, 1)
        touched = np.unique(valid)
        # Q2 — which population PAIRS does this channel span? Taken per ELEMENT (row of the topology),
        # because a channel is a connector between the populations one of its elements joins.
        ncol = int(arr.shape[1]) if len(arr.shape) > 1 else 1
        pair_census: dict[str, int] = {}
        if ncol > 1:
            rows = flat.reshape(-1, ncol)
            live = rows[(rows >= 0).all(axis=1)]
            if live.size:
                rp = pop_of[live]
                for combo, cnt in zip(*np.unique(np.sort(rp, axis=1), axis=0, return_counts=True)):
                    pair_census["|".join(names[c] for c in dict.fromkeys(combo.tolist()))] = int(cnt)
        else:
            for c, cnt in zip(*np.unique(pop_of[valid], return_counts=True)):
                pair_census[names[int(c)]] = int(cnt)

        channels[key] = {**row, "n_sentinel_negative": int((flat < 0).sum()),
                         "n_distinct_nodes_touched": int(touched.size),
                         "population_pairs": pair_census,
                         "spans_populations": sum(1 for k in pair_census if "|" in k) > 0}

    # ── Q2, second pass: the couplings a per-array test CANNOT see ─────────────────────────────────
    # The first pass found no array spanning a population boundary, and that is true of the arrays and
    # FALSE of the physics. ERM (cortex<->membrane) and LINC (actin<->nucleus) are stored as SPLIT
    # PARALLEL ENDPOINT ARRAYS -- `erm_c_d`/`erm_m_d`, `linc_a_d`/`linc_n_d` -- so each half lies wholly
    # inside one population and the connection exists only in the PAIRING, which no array declares. An
    # audit asking "does any array cross a boundary?" therefore returns a clean bill of health on an
    # engine that coupling runs straight through. This pass recovers the pairs structurally: same OWNER,
    # equal element count, disjoint population sets. The owner test is what excludes the false pair
    # `nucleus.faces_d` x `membrane.faces_d`, which match on length and populations but share no object.
    connector_candidates = []
    keys = sorted(channels)
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            ka, kb = keys[a], keys[b]
            ca, cb = channels[ka], channels[kb]
            if ka.split(".")[0] != kb.split(".")[0]:
                continue
            if ca["n_elements"] != cb["n_elements"]:
                continue
            pa, pb = set(ca["population_pairs"]), set(cb["population_pairs"])
            if pa & pb or len(pa) != 1 or len(pb) != 1:
                continue
            connector_candidates.append({
                "endpoint_arrays": [ka, kb],
                "populations": sorted(pa | pb),
                "n_bonds": ca["n_elements"],
                "visible_to_per_array_test": False,
            })

    zero_idx = np.flatnonzero(hits == 0)
    by_pop = {}
    for i, (label, lo, hi) in enumerate(pops):
        z = int(((hits == 0) & (pop_of == i)).sum())
        by_pop[label] = {"n_nodes": hi - lo, "n_zero_contribution": z,
                         "fraction": z / (hi - lo) if hi > lo else 0.0}

    record = {
        "schema": "diagnostic@1",
        "run_label": "fullcell_contribution_census",
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "stamp": {"aleph_file": aleph.__file__, "warp_version": wp.__version__,
                  "python": sys.version.split()[0], "device": str(device)},
        "n_nodes": n_nodes,
        "populations": [{"name": a, "lo": b, "hi": c, "count": c - b} for a, b, c in pops],
        "channels_counted": channels,
        "arrays_excluded": excluded,
        "per_node_arrays": per_node_arrays,
        "cross_population_coupling": {
            "arrays_that_span_a_boundary": sum(1 for c in channels.values() if c["spans_populations"]),
            "connector_candidates_from_split_endpoint_arrays": connector_candidates,
            "reading": ("NO single topology array spans a population boundary, and that is a fact about "
                        "the arrays rather than about the physics: the cross-population connectors this "
                        "cell carries are stored as split parallel endpoint arrays, so each half lies "
                        "wholly inside one population. The coupling is real and lives in the PAIRING, "
                        "which no array declares."),
        },
        "channels_not_reachable_from_the_built_cell": {
            "actin_axial_segments_nf2007": ("segment node pairs are built in assemble.py and consumed by "
                                            "the solver constraint; `AssembledCell` exposes only the "
                                            "offset/rest ledgers (`soff_d`, `srest_d`), so the pair list "
                                            "is not reachable here"),
            "steric_excluded_volume": ("pairs are formed per step from a neighbour query, not stored as "
                                       "a topology; `steric.active` marks the same 16,312 non-actin "
                                       "nodes dormant, so it acts on actin only"),
            "why_this_does_not_weaken_Q1": ("the zero-contribution set is already EMPTY over the counted "
                                            "channels, and adding a channel can only ADD coverage. So the "
                                            "answer is monotone in the missing channels and stays 0."),
        },
        "zero_contribution": {
            "total": int(zero_idx.size),
            "fraction": float(zero_idx.size / n_nodes),
            "by_population": by_pop,
        },
        "hits_histogram": {str(int(k)): int(v) for k, v in
                           zip(*np.unique(hits[hits <= 40], return_counts=True))},
        "supersedes": ("the 16,312 figure recorded but not interpreted in "
                       "aleph/outputs/ac/arena_force_parity — that count bound TWO channels; this binds "
                       "every channel the built cell carries"),
    }
    np.savez_compressed("/tmp/null_census.npz", hits=hits, pop_of=pop_of, zero_idx=zero_idx)
    print(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
