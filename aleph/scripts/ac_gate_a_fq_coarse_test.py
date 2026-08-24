"""GATE-A native probe: fiber-quotient inter-fiber coarse ON vs OFF on the seeded resting cortex.

Modes (--mode): ``pathB`` (default, RECOMMENDED) assembles the crosslink-weighted BSR ``A_c`` and solves it
with a DEEP block-Jacobi CG (cheap ``bsr_mv`` SpMV) — the strong solve that closes the native plateau;
``pathA`` is the matrix-free ``P^T A P`` variant (correctness-first, too-weak inner budget).

Builds the EXACT seeded cell used by ``ac_gate_a_projected_residual.py`` (fraction 0.5 / f 1.5 / capture 0.6 /
subdiv 6, ``overlap_free_cortex`` + ``erm_radial_pairing``), then runs a faithful projected quasi-static
descent — accumulate raw F -> project to ``P F`` -> ``ProjectedAnalyticCG.solve`` (the exact
``(aI + P K P) dx = P F`` inner solve) -> geometric backtracking line search against the exact projected
residual — with the fiber-quotient coarse correction DISABLED then ENABLED, and reports the ``max|PF|``
trajectory + final for each.

The plateau to beat is ``max|PF| ~ 1.16 pN`` over the actin cortex, below which node-Jacobi /
per-fiber-block / global-``l<=2`` preconditioning cannot drive the resting residual (design:
``aleph/docs/v2_audit/FIBER_QUOTIENT_COARSE_PLAN_2026-07-24.md``, Path A milestone P1). This is a
PRECONDITIONER A/B: the fixed point and the residual gate are IDENTICAL ON vs OFF (CLAUDE.md no
gate-loosening); only the CG convergence path changes. Per-trial ``max|PF|`` host readback drives the line
search — this is a diagnostic probe, not the production hot loop.

CUDA-gated -> gbook A5000. Run:  python aleph/scripts/ac_gate_a_fq_coarse_test.py --native
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all
from aleph.components.incumbent.implicit_mechanics import (
    ProjectedAnalyticCG,
    compute_regularization_kernel,
    omitted_regularization_base,
)
from aleph.components.incumbent.inner_mechanics import (
    max_abs_pressure_kernel,
    max_force_kernel,
    project_constraint_forces_kernel,
)


@wp.kernel
def _axpy_kernel(
    base: wp.array(dtype=wp.vec3d), scale: wp.float64,
    direction: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    out[i] = base[i] + scale * direction[i]


def _build_regularization(cell, omitted_base, reg_d, max_pressure_d, d) -> None:
    """Compose the live regularizer ``a`` exactly as the runtime does (omitted-base vs live pressure jump)."""
    max_pressure_d.zero_()
    if cell.grid is not None and cell.membrane_pressure is not None:
        wp.launch(max_abs_pressure_kernel, dim=cell.grid.shape,
                  inputs=[cell.grid.p, cell.grid.mask, wp.float64(cell.membrane_pressure.p_ext)],
                  outputs=[max_pressure_d], device=d)
    wp.launch(compute_regularization_kernel, dim=1,
              inputs=[wp.float64(omitted_base), max_pressure_d, wp.float64(cell.pressure_edge_um), reg_d],
              device=d)


def _projected_force(cell, pos, f, projected, diag, rhs, finite, d) -> None:
    """Accumulate raw F at ``pos`` and write the exact per-fiber projected force ``P F`` into ``projected``."""
    _accumulate_all(cell, pos, f)
    wp.copy(projected, f)
    if cell.n_fibers:
        wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
                  inputs=[pos, f, cell.foff_d, cell.soff_d, projected, diag, rhs, finite], device=d)


def _max_pf(cell, projected, maxbuf, d) -> float:
    """Whole-cell ``max|PF|`` via the device reduction (single scalar readback)."""
    maxbuf.zero_()
    wp.launch(max_force_kernel, dim=cell.n_total, inputs=[projected], outputs=[maxbuf], device=d)
    return float(maxbuf.numpy()[0])


def _actin_max_pf(projected, n_actin) -> float:
    """Actin-cortex ``max|PF|`` (the ~1.16 pN plateau metric) via a host norm over the cortex block."""
    pf = projected.numpy()[:n_actin]
    return float(np.linalg.norm(pf, axis=1).max())


# ── per-node residual-localization diagnostic (--dump-plateau; default-off, out-of-hot-loop) ─────────────
# The plateau metric ``max|PF|`` is a SINGLE-worst-node scalar, so it cannot say whether the ~3.42 pN fine-mesh
# floor is one geometric hotspot (a fixable steric/branch/anchor locus) or a distributed fine-mesh floor. This
# classifier reads the exact per-actin-node projected force ``|PF|`` the plateau is taken over, ranks the top-K,
# and tags each of them with the nearest structural feature it sits on — so the Lead can SEE and CLASSIFY the
# residual. Pure host numpy, run ONCE at the end of a descent (never in the inner solve). PRIMARY-label priority
# (a node can touch several features): steric > branch > crosslink > erm > boundary > interior — steric contact
# (k_EV=1e3, the stiffest local coupling) is the most likely hotspot, interior backbone the genuine-floor case.
_FEATURE_NAMES = ("interior", "boundary", "crosslink", "erm", "branch", "steric")


def classify_plateau_nodes(
    pos: np.ndarray,
    pf_mag: np.ndarray,
    node_fiber: np.ndarray,
    *,
    xl_pairs: np.ndarray | None,
    erm_cortex_idx: np.ndarray | None,
    branch_triples: np.ndarray | None,
    fiber_offsets: np.ndarray | None,
    sigma_ev_um: float,
    k: int = 200,
) -> dict:
    r"""Rank the top-``k`` highest-``|PF|`` actin nodes and tag each with its nearest structural feature.

    Pure numpy so it runs on the dev Mac / in a CPU test. ``pos`` (N,3), ``pf_mag`` (N,) and ``node_fiber`` (N,)
    are the per-actin-node position / projected-force magnitude / fiber id. The structural arrays are the same
    ones the operator couples through: ``xl_pairs`` (M,2) crosslink node pairs, ``erm_cortex_idx`` the cortex-side
    ERM anchor nodes, ``branch_triples`` (B,3) Arp2/3 junctions, ``fiber_offsets`` (F+1,) for the free fiber-end
    nodes. A node is a STERIC-contact node if a DIFFERENT-fiber node sits within the WCA cutoff
    ``2^{1/6} sigma_ev`` (the exact repulsive support the operator uses). Returns the ranked indices, per-node
    feature booleans + primary category, the ``|PF|`` deciles, the top-k feature breakdown, and the top-k spatial
    spread (bbox / radius-of-gyration / fraction within 0.5 um of the single worst node) — localized-vs-distributed.
    """
    from scipy.spatial import cKDTree

    n = int(pf_mag.shape[0])
    k = int(min(k, n))
    order = np.argsort(pf_mag)                      # ascending
    topk_idx = order[::-1][:k].astype(np.int64)     # k largest |PF|, descending
    topk_pos = pos[topk_idx]
    topk_pf = pf_mag[topk_idx].astype(np.float64)

    # feature membership sets (node indices into the actin block)
    xl_nodes = (np.unique(np.asarray(xl_pairs, np.int64).reshape(-1))
                if xl_pairs is not None and np.asarray(xl_pairs).size else np.zeros(0, np.int64))
    erm_nodes = (np.unique(np.asarray(erm_cortex_idx, np.int64).reshape(-1))
                 if erm_cortex_idx is not None and np.asarray(erm_cortex_idx).size else np.zeros(0, np.int64))
    branch_nodes = (np.unique(np.asarray(branch_triples, np.int64).reshape(-1))
                    if branch_triples is not None and np.asarray(branch_triples).size else np.zeros(0, np.int64))
    if fiber_offsets is not None and np.asarray(fiber_offsets).size >= 2:
        off = np.asarray(fiber_offsets, np.int64)
        boundary_nodes = np.unique(np.concatenate([off[:-1], off[1:] - 1]))
    else:
        boundary_nodes = np.zeros(0, np.int64)

    # steric: a same-position tree, queried ONLY at the top-k nodes for a cross-fiber neighbour < WCA cutoff.
    cutoff = (2.0 ** (1.0 / 6.0)) * float(sigma_ev_um)
    tree = cKDTree(pos)
    steric_hit = np.zeros(k, bool)
    for r, idx in enumerate(topk_idx):
        fi = node_fiber[idx]
        for j in tree.query_ball_point(pos[idx], cutoff):
            if j != idx and node_fiber[j] != fi:
                steric_hit[r] = True
                break

    feat = np.zeros((k, len(_FEATURE_NAMES)), bool)
    feat[:, 0] = True                                              # interior baseline (always "present")
    feat[:, 1] = np.isin(topk_idx, boundary_nodes)
    feat[:, 2] = np.isin(topk_idx, xl_nodes)
    feat[:, 3] = np.isin(topk_idx, erm_nodes)
    feat[:, 4] = np.isin(topk_idx, branch_nodes)
    feat[:, 5] = steric_hit
    # primary = highest-priority feature present (steric=5 down to interior=0)
    prio = np.arange(len(_FEATURE_NAMES))
    primary_code = (feat * prio[None, :]).max(axis=1).astype(np.int64)
    primary = np.array(_FEATURE_NAMES, dtype="<U16")[primary_code]

    breakdown = {name: int((primary_code == c).sum()) for c, name in enumerate(_FEATURE_NAMES)}
    deciles = np.percentile(pf_mag, np.arange(0, 101, 10)).astype(np.float64)

    # spatial spread of the top-k → localized (one locus) vs distributed (over the sphere)
    worst = topk_pos[0]
    d_worst = np.linalg.norm(topk_pos - worst[None, :], axis=1)
    centroid = topk_pos.mean(axis=0)
    r_gyr = float(np.sqrt(np.mean(np.sum((topk_pos - centroid[None, :]) ** 2, axis=1))))
    bbox_lo = topk_pos.min(axis=0).astype(np.float64)
    bbox_hi = topk_pos.max(axis=0).astype(np.float64)
    bbox_diag = float(np.linalg.norm(bbox_hi - bbox_lo))
    frac_within_half_um = float(np.mean(d_worst < 0.5))
    shell_r = float(np.linalg.norm(pos, axis=1).mean())            # cortex shell radius (for scale)

    lines = [
        f"top-{k} residual-localization on the actin cortex ({n:,} nodes)",
        f"  max|PF|            = {float(pf_mag.max()):.4f} pN   (worst node id {int(topk_idx[0])})",
        f"  top-{k} mean|PF|   = {float(topk_pf.mean()):.4f} pN   (min in top-k {float(topk_pf.min()):.4f})",
        "  |PF| deciles [pN]  = " + " ".join(f"{v:.3f}" for v in deciles),
        "  top-k feature breakdown (primary): "
        + ", ".join(f"{breakdown[nm]} {nm}" for nm in reversed(_FEATURE_NAMES)),
        f"  top-k spread: R_gyration={r_gyr:.3f} um  bbox_diag={bbox_diag:.3f} um  "
        f"(cortex shell R~{shell_r:.2f} um)",
        f"  clustered@worst: {frac_within_half_um*100:.0f}% of top-k within 0.5 um of the single worst node",
        "  => LOCALIZED (fixable geometric hotspot) if R_gyration << shell R and one feature dominates; "
        "DISTRIBUTED (fine-mesh floor) if the top-k scatter over the shell and are mostly interior.",
    ]
    summary_text = "\n".join(lines)
    return {
        "topk_idx": topk_idx, "topk_pf": topk_pf, "topk_pos": topk_pos.astype(np.float32),
        "topk_feature_bool": feat, "topk_primary_code": primary_code, "topk_primary": primary,
        "feature_names": np.array(_FEATURE_NAMES, dtype="<U16"),
        "breakdown": breakdown, "deciles": deciles,
        "bbox_lo": bbox_lo, "bbox_hi": bbox_hi, "bbox_diag": bbox_diag,
        "r_gyration": r_gyr, "frac_within_0p5um_of_worst": frac_within_half_um,
        "cortex_shell_radius_um": shell_r, "cutoff_um": cutoff,
        "summary_text": summary_text,
    }


# ── per-FORCE-TERM decomposition at the plateau (--dump-plateau; default-off, out-of-hot-loop) ────────────
# The residual-localization classifier above tags WHERE the worst nodes sit (steric/crosslink/... loci); it
# cannot say WHICH force term is the ~3.4 pN unbalanced contribution there. This decomposition MEASURES that:
# it evaluates EACH additive force family of the inner assembly (driver._accumulate_all) SEPARATELY into its
# own zeroed buffer at the plateau positions, plus the constraint force the projection removes
# (Jᵀλ = F_raw − PF), so the Lead can read off which single term is ≈|PF| at the worst node. Pure diagnostic:
# no launch here touches the solve; the union of the isolated launches is byte-identical to _accumulate_all,
# so Σ_terms F_term == F_raw to round-off (the completeness self-check + the CPU test both assert this).


def summarize_force_decomposition(
    term_names: list[str],
    term_forces: dict[str, np.ndarray],
    raw_vec: np.ndarray,
    pf_vec: np.ndarray,
    topk_idx: np.ndarray,
) -> dict:
    r"""Build the per-node per-term force magnitudes + the constraint term and a top-K table (pure numpy).

    ``term_forces[name]`` is the ``(N,3)`` per-actin-node force of ONE isolated force family (bending,
    crosslink, branch, erm, steric, myosin, pressure, ...). ``raw_vec`` ``(N,3)`` is the TOTAL assembled force
    ``F`` BEFORE projection (``driver._accumulate_all`` output over the actin block); ``pf_vec`` ``(N,3)`` is the
    projected force ``PF`` the plateau ``max|PF|`` is taken over. The inextensibility / Lagrange constraint force
    the projection removes is exactly ``Jᵀλ = F_raw − PF`` (tracked as its own column, NOT part of the term sum).

    Returns the completeness residual ``max|Σ_terms F_term − F_raw|`` (must be ~machine round-off — this proves
    the decomposition is complete and faithful), per-node magnitude arrays for every term + constraint + totals,
    the top-K table aligned to ``topk_idx`` (rows = worst nodes, columns = the per-term magnitudes), and a
    printable summary naming the DOMINANT term at the worst node and the top-K mean. Runs on the dev Mac / in a
    CPU test — it consumes numpy arrays and never launches a kernel.
    """
    raw_vec = np.asarray(raw_vec, np.float64)
    pf_vec = np.asarray(pf_vec, np.float64)
    n = int(raw_vec.shape[0])
    topk_idx = np.asarray(topk_idx, np.int64)

    stacked = np.zeros((n, 3), np.float64)
    mags: dict[str, np.ndarray] = {}
    for nm in term_names:
        fv = np.asarray(term_forces[nm], np.float64)
        stacked = stacked + fv
        mags[nm] = np.linalg.norm(fv, axis=1)
    completeness = float(np.linalg.norm(stacked - raw_vec, axis=1).max()) if n else 0.0

    constraint_vec = raw_vec - pf_vec                       # Jᵀλ — the force the constraint projection removes
    mags["JTlambda"] = np.linalg.norm(constraint_vec, axis=1)
    raw_total = np.linalg.norm(raw_vec, axis=1)
    pf_total = np.linalg.norm(pf_vec, axis=1)

    phys_cols = list(term_names) + ["JTlambda"]             # the physical terms compared to find the dominant one
    col_names = phys_cols + ["raw_total", "PF_total"]
    table = np.zeros((topk_idx.shape[0], len(col_names)), np.float64)
    for c, nm in enumerate(phys_cols):
        table[:, c] = mags[nm][topk_idx]
    table[:, len(phys_cols)] = raw_total[topk_idx]
    table[:, len(phys_cols) + 1] = pf_total[topk_idx]

    worst_row = table[0] if table.shape[0] else np.zeros(len(col_names))
    topk_mean = table.mean(axis=0) if table.shape[0] else np.zeros(len(col_names))
    worst_dom_c = int(np.argmax(worst_row[:len(phys_cols)])) if table.shape[0] else 0
    mean_dom_c = int(np.argmax(topk_mean[:len(phys_cols)])) if table.shape[0] else 0

    width = max(9, *(len(c) for c in col_names))
    header = "  ".join(f"{c:>{width}}" for c in col_names)

    def _row(vals):
        return "  ".join(f"{v:>{width}.4f}" for v in vals)

    worst_id = int(topk_idx[0]) if topk_idx.shape[0] else -1
    lines = [
        f"force-term decomposition at the plateau (top-{topk_idx.shape[0]} worst |PF| nodes)",
        f"  completeness: max|Σ_terms F − F_raw| = {completeness:.3e} pN  (should be ~machine round-off)",
        "  per-term |force| [pN]:",
        "    " + header,
        f"  worst node {worst_id:>8d}:",
        "    " + _row(worst_row),
        f"    => DOMINANT term at worst node = {phys_cols[worst_dom_c]}  "
        f"({worst_row[worst_dom_c]:.4f} pN of the {worst_row[len(phys_cols)+1]:.4f} pN |PF|)",
        f"  top-{topk_idx.shape[0]} mean:",
        "    " + _row(topk_mean),
        f"    => DOMINANT term over top-K (mean) = {phys_cols[mean_dom_c]}  ({topk_mean[mean_dom_c]:.4f} pN)",
        "  => the term whose column is ≈ the PF_total column at the worst node IS the plateau residual.",
    ]
    summary_text = "\n".join(lines)
    return {
        "completeness": completeness,
        "col_names": np.array(col_names, dtype="<U24"),
        "phys_cols": np.array(phys_cols, dtype="<U24"),
        "topk_table": table,
        "worst_row": worst_row.astype(np.float64),
        "topk_mean_row": topk_mean.astype(np.float64),
        "worst_dominant": phys_cols[worst_dom_c] if table.shape[0] else "",
        "topk_mean_dominant": phys_cols[mean_dom_c] if table.shape[0] else "",
        "per_node_mag": {nm: mags[nm] for nm in phys_cols},
        "summary_text": summary_text,
    }


def _accumulate_terms(cell, pos, d) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    r"""Evaluate each additive force family of ``driver._accumulate_all`` SEPARATELY into its own zeroed buffer.

    Returns ``(order, term_forces, raw_vec)`` over the ACTIN block only: ``order`` is the term-name list,
    ``term_forces[name]`` the ``(n_actin,3)`` per-node force of that ONE family, and ``raw_vec`` the ``(n_actin,3)``
    TOTAL assembled force from a full ``_accumulate_all`` (for the completeness self-check ``Σ_terms == raw``).

    Native-only (launches Warp-CUDA kernels; the dev Mac cannot run it). Each entry launches EXACTLY the same
    kernel / accumulate call, with the same args and units, that ``_accumulate_all`` composes — the only change
    is that each is run into a freshly-zeroed buffer instead of the shared accumulator, so the isolated arrays
    sum back to the total. Membrane bending/area act only on membrane nodes, so on the actin block the membrane
    term is exactly the cortex-side ERM tether force (labelled ``erm``).
    """
    from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD
    from aleph.components.weave.branch_angle_warp import branch_angle_kernel
    from aleph.laws.forces_warp import cytosim_bending_kernel
    from aleph.laws.network_warp import _zero, link_spring_kernel

    n = cell.n_total
    na = int(cell.n_actin)
    fb = wp.zeros(n, dtype=wp.vec3d, device=d)
    order: list[str] = []
    term_forces: dict[str, np.ndarray] = {}

    def _emit(name: str, launch) -> None:
        wp.launch(_zero, dim=n, inputs=[fb], device=d)
        launch(fb)
        term_forces[name] = fb.numpy()[:na].astype(np.float64).copy()
        order.append(name)

    if cell.n_tri:
        _emit("bending", lambda f: wp.launch(cytosim_bending_kernel, dim=cell.n_tri,
              inputs=[pos, cell.tri_d, cell.alpha_d, f], device=d))
    if cell.n_xl:
        _emit("crosslink", lambda f: wp.launch(link_spring_kernel, dim=cell.n_xl,
              inputs=[pos, cell.xl_d, cell.kxl_d, cell.r0xl_d, f], device=d))
    if cell.branch_triples_d is not None and cell.n_branch:
        _emit("branch", lambda f: wp.launch(branch_angle_kernel, dim=cell.n_branch,
              inputs=[pos, cell.branch_triples_d, cell.branch_active_d,
                      ARP23_THETA0_RAD, ARP23_K_THETA, f], device=d))
    if cell.myosin is not None:
        _emit("myosin", lambda f: cell.myosin.accumulate(pos, f))
    if cell.nucleus is not None:
        _emit("nucleus", lambda f: cell.nucleus.accumulate(pos, f))
    if cell.membrane is not None:                            # on the actin block == cortex-side ERM tether force
        _emit("erm", lambda f: cell.membrane.accumulate(pos, f))
    if cell.membrane_pressure is not None:
        _emit("membrane_pressure", lambda f: cell.membrane_pressure.accumulate(pos, f))
    if cell.steric is not None:
        _emit("steric", lambda f: cell.steric.accumulate(cell.state, f))
    if cell.pressure is not None:
        _emit("pressure", lambda f: cell.pressure.accumulate(cell.state, f))

    raw = wp.zeros(n, dtype=wp.vec3d, device=d)
    _accumulate_all(cell, pos, raw)                          # byte-identical total for the completeness check
    raw_vec = raw.numpy()[:na].astype(np.float64).copy()
    return order, term_forces, raw_vec


def _dump_plateau_from_cell(cell, projected, path: str, *, label: str, k: int = 200,
                            term_forces: dict[str, np.ndarray] | None = None,
                            term_order: list[str] | None = None,
                            raw_vec: np.ndarray | None = None) -> None:
    """Extract the actin per-node ``|PF|`` + structural arrays from the composed cell and write the dump npz."""
    from aleph.components.incumbent.assemble import SIGMA_EV_UM

    n_actin = int(cell.n_actin)
    pf_vec = projected.numpy()[:n_actin]
    pf_mag = np.linalg.norm(pf_vec, axis=1).astype(np.float64)
    pos = cell.pos_d.numpy()[:n_actin].astype(np.float64)
    node_fiber = (cell.node_fiber_d.numpy().astype(np.int64) if cell.node_fiber_d is not None
                  else np.zeros(n_actin, np.int64))
    xl_pairs = cell.xl_d.numpy() if cell.n_xl else None
    branch_triples = (cell.branch_triples_d.numpy() if cell.branch_triples_d is not None else None)
    fiber_offsets = cell.foff_d.numpy()
    erm_cortex_idx = None
    if cell.membrane is not None and getattr(cell.membrane, "n_erm", 0):
        c = cell.membrane.erm_c_d.numpy().astype(np.int64)
        bound = cell.membrane.erm_bound_d.numpy().astype(np.int64)
        erm_cortex_idx = c[bound != 0]

    result = classify_plateau_nodes(
        pos, pf_mag, node_fiber, xl_pairs=xl_pairs, erm_cortex_idx=erm_cortex_idx,
        branch_triples=branch_triples, fiber_offsets=fiber_offsets, sigma_ev_um=SIGMA_EV_UM, k=k)

    print(f"\n  --- plateau dump [{label}] ---\n" + "\n".join("    " + ln for ln in result["summary_text"].split("\n")))

    extra: dict = {}
    if term_forces is not None and raw_vec is not None:
        order = term_order if term_order is not None else list(term_forces.keys())
        decomp = summarize_force_decomposition(order, term_forces, raw_vec, pf_vec, result["topk_idx"])
        print("\n" + "\n".join("    " + ln for ln in decomp["summary_text"].split("\n")))
        extra = {
            "decomp_col_names": decomp["col_names"],
            "decomp_phys_cols": decomp["phys_cols"],
            "decomp_topk_table": decomp["topk_table"],
            "decomp_worst_row": decomp["worst_row"],
            "decomp_topk_mean_row": decomp["topk_mean_row"],
            "decomp_worst_dominant": np.array(decomp["worst_dominant"], dtype=object),
            "decomp_topk_mean_dominant": np.array(decomp["topk_mean_dominant"], dtype=object),
            "decomp_completeness_pn": np.float64(decomp["completeness"]),
            "decomp_summary_text": np.array(decomp["summary_text"], dtype=object),
            **{f"term_mag_{nm}": decomp["per_node_mag"][nm].astype(np.float64)
               for nm in decomp["per_node_mag"]},
        }

    np.savez_compressed(
        path,
        pos=pos.astype(np.float32), pf_mag=pf_mag, node_fiber=node_fiber.astype(np.int32),
        topk_idx=result["topk_idx"], topk_pf=result["topk_pf"], topk_pos=result["topk_pos"],
        topk_feature_bool=result["topk_feature_bool"], topk_primary_code=result["topk_primary_code"],
        topk_primary=result["topk_primary"], feature_names=result["feature_names"],
        deciles=result["deciles"], bbox_lo=result["bbox_lo"], bbox_hi=result["bbox_hi"],
        bbox_diag=np.float64(result["bbox_diag"]), r_gyration=np.float64(result["r_gyration"]),
        frac_within_0p5um_of_worst=np.float64(result["frac_within_0p5um_of_worst"]),
        cortex_shell_radius_um=np.float64(result["cortex_shell_radius_um"]),
        cutoff_um=np.float64(result["cutoff_um"]),
        breakdown_json=np.array(str(result["breakdown"]), dtype=object),
        summary_text=np.array(result["summary_text"], dtype=object),
        label=np.array(label, dtype=object), k=np.int64(k),
        **extra,
    )
    print(f"    wrote {path}  (per-node pos/|PF|/fiber over {n_actin:,} actin nodes + top-{k} classification)")


def _descent(cell, ws, pos0, reg_d, finite_d, d, *, outer: int, label: str,
             dump_path: str | None = None, dump_k: int = 200) -> dict:
    """Projected quasi-static descent from ``pos0``; return the max|PF| trajectory (whole-cell + actin)."""
    wp.copy(cell.pos_d, pos0)
    pos = cell.pos_d
    n = cell.n_total
    f = wp.zeros(n, dtype=wp.vec3d, device=d)
    projected = wp.zeros(n, dtype=wp.vec3d, device=d)
    trial_pf = wp.zeros(n, dtype=wp.vec3d, device=d)
    pos_prev = wp.zeros(n, dtype=wp.vec3d, device=d)
    pos_trial = wp.zeros(n, dtype=wp.vec3d, device=d)
    pos_best = wp.zeros(n, dtype=wp.vec3d, device=d)
    diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
    rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
    maxbuf = wp.zeros(1, dtype=wp.float64, device=d)
    trials = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]

    whole_hist: list[float] = []
    actin_hist: list[float] = []
    t_start = time.time()
    for it in range(outer):
        # P F at the current position (the RHS for this linear solve and the residual we report/minimize).
        _projected_force(cell, pos, f, projected, diag, rhs, finite_d, d)
        whole = _max_pf(cell, projected, maxbuf, d)
        actin = _actin_max_pf(projected, cell.n_actin)
        whole_hist.append(whole)
        actin_hist.append(actin)

        # dx solves (aI + P K P) dx = P F at the frozen tangent; finite latch reset per solve.
        finite_d.fill_(1)
        dx = ws.solve(pos, projected, reg_d, finite_d)

        # Geometric backtracking against the EXACT projected residual; t=0 (no move) is the monotone baseline.
        wp.copy(pos_prev, pos)
        best_res = whole
        best_t = 0.0
        wp.copy(pos_best, pos_prev)
        for t in trials:
            wp.launch(_axpy_kernel, dim=n, inputs=[pos_prev, wp.float64(t), dx, pos_trial], device=d)
            _projected_force(cell, pos_trial, f, trial_pf, diag, rhs, finite_d, d)
            res = _max_pf(cell, trial_pf, maxbuf, d)
            if res < best_res:
                best_res = res
                best_t = t
                wp.copy(pos_best, pos_trial)
        wp.copy(cell.pos_d, pos_best)
        pos = cell.pos_d
        if it % 5 == 0 or it == outer - 1:
            print(f"    [{label}] it {it:3d}  max|PF| whole {whole:9.4f}  actin {actin:9.4f}  "
                  f"(step t={best_t:g}, {time.time() - t_start:6.1f}s)")

    # Final residual after the last accepted step.
    _projected_force(cell, cell.pos_d, f, projected, diag, rhs, finite_d, d)
    final_whole = _max_pf(cell, projected, maxbuf, d)
    final_actin = _actin_max_pf(projected, cell.n_actin)
    if dump_path is not None:
        # Per-force-term decomposition at the frozen plateau positions (native-only Warp launches). Each term
        # is isolated into its own zeroed buffer; the union is byte-identical to _accumulate_all so Σ_terms
        # == F_raw to round-off (the writer prints the completeness residual as a self-check).
        term_order, term_forces, raw_vec = _accumulate_terms(cell, cell.pos_d, d)
        _dump_plateau_from_cell(cell, projected, dump_path, label=label, k=dump_k,
                                term_forces=term_forces, term_order=term_order, raw_vec=raw_vec)
    return {
        "whole_hist": whole_hist, "actin_hist": actin_hist,
        "final_whole": final_whole, "final_actin": final_actin,
        "wall": time.time() - t_start,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native", action="store_true", help="full 70,686 cortical filaments (else 1500)")
    ap.add_argument("--mode", choices=["pathA", "pathB"], default="pathB",
                    help="pathB = strong BSR A_c + deep block-Jacobi CG (recommended); pathA = matrix-free")
    ap.add_argument("--outer", type=int, default=40, help="projected quasi-static descent iterations")
    ap.add_argument("--cg-iters", type=int, default=64,
                    help="ProjectedAnalyticCG.max_iterations (fine PCG needs ~60 with a strong coarse)")
    ap.add_argument("--fq-iters", type=int, default=0,
                    help="inner coarse CG budget (0 => mode default: pathB 1000, pathA 40)")
    ap.add_argument("--coarse-modes", type=int, default=0, help="global l<=2 coarse modes (0..12)")
    ap.add_argument("--coarse-iters", type=int, default=0, help="translation-coarse V-cycle iterations")
    ap.add_argument("--myosin-fraction", type=float, default=0.5,
                    help="resting bound-myosin fraction (0 => CLEAN turgor baseline, no discrete myosin)")
    ap.add_argument("--membrane-subdiv", type=int, default=6,
                    help="membrane subdivisions (8 = ratified grid-converged value)")
    ap.add_argument("--cortex-seg-um", type=float, default=0.5, dest="cortex_seg_um",
                    help="cortex segment rest length ℓ₀ [µm] → mesh/pore size (baseline 0.5; fine mesh 0.075)")
    ap.add_argument("--cortex-density", type=float, default=20.0, dest="cortex_density",
                    help="crosslinks / filament (baseline 20; fine mesh 40)")
    ap.add_argument("--disable-fiber-block", action="store_true", dest="disable_fiber_block",
                    help="drop the dense per-fiber block preconditioner (frees ~L_f²-scaling GPU mem; "
                         "essential at fine mesh where 41 nodes/fiber blows it to ~2.85 GB — pathB coarse "
                         "supplies the inter-fiber coupling)")
    ap.add_argument("--multigrid", action="store_true", dest="multigrid",
                    help="use the fiber-arclength geometric-multigrid V-cycle as M^-1 (design memo "
                         "FINE_MESH_MULTIGRID_DESIGN_2026-07-24): a symmetric SPD V(2,2) — line-smoother "
                         "pre/post + matrix-free Galerkin arclength coarse levels + additive fiber-quotient "
                         "coarsest. Spans the fine-mesh class-3 error (smooth-along-arclength × inter-fiber) "
                         "the per-fiber block + rigid coarse both miss; this is the wiring meant to let the "
                         "physiological ~75nm mesh converge. Implies --disable-fiber-block (the O(L²) Cholesky "
                         "is replaced by the O(L) line smoother).")
    ap.add_argument("--cg-check-every", type=int, default=0, dest="cg_check_every",
                    help="SPEED: outer-PCG host early-exit cadence. 0 = incumbent fixed --cg-iters budget "
                         "(byte-identical). >0 (e.g. 4) reads the device active/converged latch every K iters "
                         "and breaks once cleared — a strong MG converges in ~10-15 iters, so this skips the "
                         "wasted post-convergence V-cycles. dx is bit-identical to running to --cg-iters.")
    ap.add_argument("--mg-assembled-coarse", action="store_true", dest="mg_assembled_coarse",
                    help="SPEED: apply the CHEAP assembled per-fiber operator (ext + pentadiagonal bending, "
                         "O(n_nodes[l])) at MG coarse levels l>=1 instead of the full-fine R∘_operator∘P (a "
                         "2.9M-node apply PER LEVEL). Inter-fiber coupling stays with the additive fiber-quotient "
                         "coarsest level. Cuts full-fine applies/V-cycle from ~43 to ~5 (level-0 only). SPD + "
                         "h-independent on the CPU gate; native-validate the reduction factor with the Lead.")
    ap.add_argument("--no-overlap-free", action="store_true", dest="no_overlap_free",
                    help="DIAGNOSTIC: build overlap_free_cortex=False so cortex fibers stay SMOOTH great-circle "
                         "arcs (the WCA overlap relaxation, which pushes crossing nodes transversely off the arc "
                         "and creates ~5deg kinks at 75nm that drive the bending plateau, is skipped). Combine with "
                         "--no-steric to remove the ~64k t0 interpenetrations' force too: if the plateau then "
                         "collapses to ~coarse, the relaxation-kink IS the confirmed bending-plateau root cause. "
                         "NOT production (overlap-free intent must be preserved by a smoothness-safe fix).")
    ap.add_argument("--no-steric", action="store_true", dest="no_steric",
                    help="DIAGNOSTIC: build with_steric=False to isolate the node-based WCA excluded-volume "
                         "contribution to the fine-mesh max|PF| plateau. If the plateau collapses to ~coarse "
                         "levels, the fine-mesh floor is steric — and combined with the superlinear "
                         "max|PF|-vs-mesh divergence (500/200/75nm=0.14/0.21/3.42) confirms node-based steric is "
                         "grid-dependent (K_EV is a flagged PI discretization GAP). NOT production (EV must be on).")
    ap.add_argument("--cortex-length-um", type=float, default=3.0, dest="cortex_length_um",
                    help="representative cortical filament contour length L [µm] (baseline 3.0)")
    ap.add_argument("--cortex-arp23-fraction", type=float, default=0.0, dest="cortex_arp23_fraction",
                    help="Arp2/3 fraction of cortical actin BY MASS (~0.33 Bovellan; 0.0 = formin-only, "
                         "bit-identical baseline). >0 splits the fixed n_filaments budget into formin + short "
                         "branched Arp2/3 (areal density preserved). NOTE: the driver's _accumulate_all does "
                         "NOT launch branch_angle_kernel, so the 70° junctions carry only the anchor spring, "
                         "not the angle-harmonic restoring force (see report).")
    ap.add_argument("--cortex-overlap-mode", choices=["transverse", "radial_span"], default="transverse",
                    dest="cortex_overlap_mode",
                    help="HOW overlap_free_cortex resolves build interpenetrations. 'transverse' (default, "
                         "byte-identical) shoves crossing nodes IN the shell plane, kinking the arc (the fine-mesh "
                         "bending plateau). 'radial_span' separates crossing fibers OUT-of-plane with a smooth "
                         "±cortex-overlap-span-node radial bump so fine (75nm) fibers stay unkinked while still "
                         "starting overlap-free — the smoothness-preserving fix under test.")
    ap.add_argument("--cortex-overlap-span", type=int, default=2, dest="cortex_overlap_span",
                    help="radial_span bump half-window in arc-neighbours (default 2 => 5-node cosine bump).")
    ap.add_argument("--seed", type=int, default=0, help="cortex construction RNG seed (ensemble check)")
    ap.add_argument("--dump-plateau", type=str, default=None, dest="dump_plateau",
                    help="DIAGNOSTIC (default-off, no flag => byte-identical run): at the END of the chosen "
                         "descent write an npz with per-actin-node (x,y,z), per-node |PF| (the same actin "
                         "projected-force the max|PF| plateau is taken over), node->fiber id, and a nearest-"
                         "structural-feature classification of the top-K worst nodes (steric-contact / Arp2/3 "
                         "branch / crosslink / ERM anchor / fiber-boundary / interior backbone) + a text "
                         "summary (max, top-K mean, |PF| deciles, feature breakdown, spatial spread), PLUS a "
                         "per-FORCE-TERM decomposition at those worst nodes: each additive family (bending / "
                         "crosslink / branch / erm / steric / myosin / pressure) is evaluated in isolation and "
                         "the constraint force Jᵀλ = F_raw − PF is measured separately, so the Lead can read off "
                         "WHICH single term is ≈|PF| at the worst node. Lets the Lead SEE whether the fine-mesh "
                         "residual is a LOCALIZED hotspot or a DISTRIBUTED floor and MEASURE which term causes it.")
    ap.add_argument("--dump-plateau-which", choices=["ON", "OFF"], default="ON", dest="dump_plateau_which",
                    help="which descent to dump the plateau from (default ON = the coarse/MG-enabled run)")
    ap.add_argument("--dump-plateau-k", type=int, default=200, dest="dump_plateau_k",
                    help="number of highest-|PF| nodes to classify (default 200)")
    args = ap.parse_args()
    fq_iters = args.fq_iters or (1000 if args.mode == "pathB" else 40)

    n = 70686 if args.native else 1500
    print(f"=== GATE-A fiber-quotient coarse ({args.mode}) A/B — native={args.native}, n_filaments={n}, "
          f"fq_iters={fq_iters}, cg_iters={args.cg_iters} ===")
    print(f"    cortex mesh: seg_um={args.cortex_seg_um}  density_per_fil={args.cortex_density}  "
          f"length_um={args.cortex_length_um}  "
          f"({'FINE' if args.cortex_seg_um < 0.5 or args.cortex_density > 20.0 else 'COARSE baseline'})")
    print(f"    overlap resolution: mode={args.cortex_overlap_mode}  span={args.cortex_overlap_span}  "
          f"overlap_free={not args.no_overlap_free}")
    _seed = args.myosin_fraction > 0.0   # <=0 => CLEAN turgor baseline (no discrete resting myosin)
    cell = build_cell(CellConfig(
        n_filaments=n, overlap_free_cortex=(not args.no_overlap_free), erm_radial_pairing=True, seed=args.seed,
        cortex_overlap_mode=args.cortex_overlap_mode, cortex_overlap_span=args.cortex_overlap_span,
        with_steric=(not args.no_steric),   # diagnostic: --no-steric isolates the node-based WCA contribution to the plateau
        membrane_subdivisions=args.membrane_subdiv,
        cortex_seg_um=args.cortex_seg_um,
        cortex_density_per_fil=args.cortex_density,
        cortex_length_um=args.cortex_length_um,
        cortex_arp23_fraction=args.cortex_arp23_fraction,
        resting_bound_myosin_fraction=(args.myosin_fraction if _seed else None),
        resting_bound_myosin_force_pn=(1.5 if _seed else None),
        resting_bound_myosin_source=("fq_coarse_TEST" if _seed else ""),
        resting_bound_myosin_capture_um=(0.6 if _seed else None),
    ))
    d = cell.device
    pos0 = wp.zeros(cell.n_total, dtype=wp.vec3d, device=d)
    wp.copy(pos0, cell.pos_d)
    print(f"    n_total={cell.n_total}  n_actin={cell.n_actin}  n_fibers={cell.n_fibers}  device={d}")

    omitted_base = omitted_regularization_base(cell)
    reg_d = wp.zeros(1, dtype=wp.float64, device=d)
    max_pressure_d = wp.zeros(1, dtype=wp.float64, device=d)
    finite_d = wp.ones(1, dtype=wp.int32, device=d)
    _build_regularization(cell, omitted_base, reg_d, max_pressure_d, d)
    print(f"    regularizer a = {float(reg_d.numpy()[0]):.4e} pN/um   "
          f"(omitted_base {omitted_base:.4e})")

    results = {}
    for label, fq_on in (("OFF (baseline)", False), (f"ON  ({args.mode})", True)):
        wp.copy(cell.pos_d, pos0)
        print(f"\n  --- coarse {label} ---")
        ws = ProjectedAnalyticCG(
            cell, max_iterations=args.cg_iters, coarse_iterations=args.coarse_iters,
            coarse_modes=args.coarse_modes, fiber_quotient_coarse=(fq_on and not args.multigrid),
            fq_coarse_iterations=fq_iters, fq_coarse_mode=args.mode,
            disable_fiber_block=(args.disable_fiber_block or (args.multigrid and fq_on)),
            multigrid=(args.multigrid and fq_on),
            mg_assembled_coarse=args.mg_assembled_coarse, cg_check_every=args.cg_check_every)
        if args.multigrid and fq_on and ws.mg is not None:
            print(f"    MULTIGRID: {ws.mg.n_levels} arclength levels (lengths {ws.mg.level_lengths}), "
                  f"coarsest={'fiber-quotient' if ws.mg.fq is not None else 'line-smooth'}")
        if fq_on and ws.fq_coarse is not None:
            print(f"    fq[{args.mode}]: N_c={ws.fq_coarse.n_coarse}  "
                  f"block_rank={getattr(ws.fq_coarse, 'block_rank', ws.fq_coarse.max_rank)}  "
                  f"ranks={ws.fq_coarse.ranks_unique}  inner_iters={fq_iters}")
        _which_on = args.dump_plateau_which == "ON"
        _do_dump = args.dump_plateau is not None and (fq_on == _which_on)
        results[fq_on] = _descent(
            cell, ws, pos0, reg_d, finite_d, d, outer=args.outer, label=label,
            dump_path=(args.dump_plateau if _do_dump else None), dump_k=args.dump_plateau_k)

    off, on = results[False], results[True]
    print("\n" + "=" * 84)
    print("VERDICT — fiber-quotient coarse vs baseline (max|PF| over ACTIN cortex; plateau ~1.16 pN)")
    print("=" * 84)
    print(f"  baseline (OFF): start {off['actin_hist'][0]:9.4f} -> final {off['final_actin']:9.4f}   "
          f"({off['wall']:.1f}s)")
    print(f"  fiber-quotient  start {on['actin_hist'][0]:9.4f} -> final {on['final_actin']:9.4f}   "
          f"({on['wall']:.1f}s)")
    beat = on["final_actin"] < off["final_actin"] - 1e-9
    below = on["final_actin"] < 1.16
    print(f"  ON beats OFF: {beat}   ON below 1.16 plateau: {below}")
    if below:
        print(f"  => {args.mode} fiber-quotient coarse DRIVES the resting residual below the 1.16 pN plateau.")
    else:
        print("  => did NOT clear 1.16; inspect trajectory / raise --fq-iters or --cg-iters or --outer.")


if __name__ == "__main__":
    main()
