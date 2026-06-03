#!/usr/bin/env python
"""Track-1 myosin head->actin binding-fraction diagnostic (KU-3.5, READ-ONLY).

WHY: STAGE-2 eliminated percolation / turnover / measurement / coherence as the
cortical-tension (g_soft) wall and localized it UPSTREAM at force GENERATION.
The stage2 head-tension smoke run found bound heads ARE loaded to ~44% of stall
BUT only ~3% of heads are bound and s_grip ~ 0. This script DIAGNOSES *why* so
few heads bind, and whether s_grip ~ 0 is a downstream consequence of the low
bound fraction or an independent stepping bug.

It instruments the ACTUALLY-BUILT cortex+myosin cell (no physics change, no
parameter change). It does NOT raise k_on / stall / capture radius -- it only
measures the gates that the existing MyosinStepUpdater applies.

Measured panels (all from the real built state + the real updater internals):

  A. Geometry-at-construction:
       head -> nearest cortex-actin BEAD distance distribution
       head -> nearest cortex-actin SEGMENT (perpendicular) distance distribution
       vs head_actin_capture_perp (perp eligibility) and
          head_actin_max_bind_dist (KDTree bead-search ceiling).
       -> answers hypothesis (a) capture distance too small, and (d) geometry.

  B. Per-gate accounting of the binding Step-2 (one instrumented tick): of all
     unbound heads, how many
       (1) had >=1 candidate bead inside max_bind_dist (KDTree),
       (2) passed the Bernoulli p_bind = 1-exp(-k_on*batch_dt) draw,
       (3) found >=1 segment within capture_perp (perp-eligible),
       (4) passed the bipolar sidedness gate (grip_walk),
       (5) were blocked by the per-bead degree caps.
     -> isolates WHICH gate sheds the heads (hypotheses a/c/d/bipolar/caps).

  C. Kinetics: run the real updater for several ticks, record per-tick
     n_bind / n_break / bound-fraction trajectory -> steady-state bound frac and
     whether it matches the naive k_on/(k_on+k_off0) ceiling.
     -> answers hypothesis (b) rate*dt and (e) strip-faster-than-rebind.

  D. s_grip downstream test: among heads that ARE bound + loaded (F>0), does
     s_grip advance over ticks? -> tells us if s_grip~0 is just "few/young
     bonds" (downstream) or a stepping bug (independent).

Usage (small / fast, CPU):
    python ffn_sim/scripts/track1_binding_diag.py --n-fil 300 --ticks 8
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml
import hoomd
from scipy.spatial import cKDTree

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import (
    resolve_cortex_myosin,
    cortex_myosin_attach_bin_rest_lengths,
)
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cell.cell import build_cortex_full_simulation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
OUT = PKG / "outputs" / "h3" / "track1_binding_diag.json"


def _tagpos(snap) -> np.ndarray:
    """Return positions indexed by particle tag.

    ``sim.state.get_snapshot()`` on a single rank returns particles already in
    tag order (row i == tag i), and the updater's tag arithmetic indexes that
    same global snapshot, so position[row] == position of tag==row directly.
    """
    return np.asarray(snap.particles.position, dtype=np.float64)


def _percentiles(x: np.ndarray) -> dict:
    if x.size == 0:
        return {}
    return {
        "n": int(x.size),
        "min_nm": float(x.min() * 1e9),
        "p10_nm": float(np.percentile(x, 10) * 1e9),
        "median_nm": float(np.median(x) * 1e9),
        "p90_nm": float(np.percentile(x, 90) * 1e9),
        "max_nm": float(x.max() * 1e9),
    }


def geometry_panel(ma, snap) -> dict:
    """Panel A: head->nearest-bead and head->nearest-segment(perp) at build."""
    p = ma.p
    pos = _tagpos(snap)
    n_actin = ma.n_cortex_actin
    r_actin = pos[:n_actin]

    # All head global tags (every head, bound or not, at construction all free).
    n_heads = 2 * p.n_heads_per_side * p.n_motors_per_cell
    head_tags = np.array(
        [ma._head_global_tag(h) for h in range(n_heads)], dtype=np.int64
    )
    r_heads = pos[head_tags]

    # head -> nearest actin BEAD
    tree = cKDTree(r_actin)
    d_bead, _ = tree.query(r_heads, k=1)

    # head -> nearest actin SEGMENT (perpendicular distance to the segment line,
    # clamped to the segment) -- mirrors the updater's eligibility metric.
    bg = ma._cortex_bond_groups
    perp_min = np.full(r_heads.shape[0], np.inf)
    if bg is not None:
        A = r_actin[bg[:, 0]]
        B = r_actin[bg[:, 1]]
        seg = B - A
        L2 = np.einsum("ij,ij->i", seg, seg)
        valid = L2 > 0.0
        A = A[valid]; seg = seg[valid]; L2 = L2[valid]
        # For each head, loop only over segments whose endpoint is plausibly near
        # (use the same KDTree gate the updater uses to keep this O(heads*nbr)).
        nbr_lists = tree.query_ball_point(r_heads, r=p.head_actin_max_bind_dist)
        seg_of_bead = ma._bead_to_segs
        for k, nbrs in enumerate(nbr_lists):
            if not nbrs:
                continue
            cand = np.unique(np.concatenate(
                [seg_of_bead[int(b)] for b in nbrs] + [np.empty(0, np.int64)]
            ))
            if cand.size == 0:
                continue
            h = r_heads[k]
            a = r_actin[bg[cand, 0]]
            sv = r_actin[bg[cand, 1]] - a
            l2 = np.einsum("ij,ij->i", sv, sv)
            good = l2 > 0
            if not good.any():
                continue
            a = a[good]; sv = sv[good]; l2 = l2[good]
            t = np.clip(np.einsum("ij,ij->i", h - a, sv) / l2, 0.0, 1.0)
            closest = a + t[:, None] * sv
            perp_min[k] = float(np.linalg.norm(h - closest, axis=1).min())
    perp_min = perp_min[np.isfinite(perp_min)]

    cap = p.head_actin_capture_perp
    mbd = p.head_actin_max_bind_dist
    return {
        "head_actin_capture_perp_nm": cap * 1e9,
        "head_actin_max_bind_dist_nm": mbd * 1e9,
        "head_rest_length_nm": p.head_rest_length * 1e9,
        "nearest_bead_dist": _percentiles(d_bead),
        "nearest_segment_perp_dist": _percentiles(perp_min),
        "frac_heads_with_bead_in_max_bind_dist": float(
            (d_bead <= mbd).mean()
        ),
        "frac_heads_perp_eligible": float(
            (perp_min <= cap).sum() / max(d_bead.size, 1)
        ),
    }


def gate_accounting(ma, snap) -> dict:
    """Panel B: per-gate accounting of one binding Step-2 (no state mutation)."""
    p = ma.p
    pos = _tagpos(snap)
    n_actin = ma.n_cortex_actin
    r_actin = pos[:n_actin]
    bg = ma._cortex_bond_groups
    cap = p.head_actin_capture_perp
    mbd = p.head_actin_max_bind_dist

    unbound = np.flatnonzero(ma._head_bound_to_actin < 0)
    head_tags = np.array(
        [ma._head_global_tag(int(h)) for h in unbound], dtype=np.int64
    )
    r_heads = pos[head_tags]
    tree = cKDTree(r_actin)
    nbr_lists = tree.query_ball_point(r_heads, r=mbd)

    n_unbound = unbound.size
    n_has_candidate = 0      # gate 1: >=1 bead in max_bind_dist
    n_perp_eligible = 0      # gate 3: >=1 segment within capture_perp
    n_bipolar_pass = 0       # gate 4: bipolar sidedness accepts (grip_walk)
    n_degree_block = 0       # gate 5: best bead degree-capped
    p_bind = 1.0 - np.exp(-p.head_actin_k_on * p.batch_dt)

    seg_of_bead = ma._bead_to_segs
    # degree snapshot (myosin attach + everything) for the cap test
    bgrp = np.asarray(snap.bonds.group, dtype=np.int64).reshape(-1)
    deg = np.bincount(bgrp[bgrp < n_actin], minlength=n_actin)

    for k, nbrs in enumerate(nbr_lists):
        if not nbrs:
            continue
        n_has_candidate += 1
        if bg is None:
            continue
        cand = np.unique(np.concatenate(
            [seg_of_bead[int(b)] for b in nbrs] + [np.empty(0, np.int64)]
        ))
        h = r_heads[k]
        best_perp = np.inf
        best_bead = -1
        for s in cand:
            a_idx = int(bg[s, 0]); b_idx = int(bg[s, 1])
            A = r_actin[a_idx]; B = r_actin[b_idx]
            sv = B - A; l2 = float(sv @ sv)
            if l2 <= 0:
                continue
            t = max(0.0, min(1.0, float((h - A) @ sv) / l2))
            closest = A + t * sv
            perp = float(np.linalg.norm(h - closest))
            if perp <= cap and perp < best_perp:
                dA = float(np.linalg.norm(h - A))
                dB = float(np.linalg.norm(h - B))
                best_bead = a_idx if dA <= dB else b_idx
                best_perp = perp
        if best_bead < 0:
            continue
        n_perp_eligible += 1
        # degree cap (matches updater _MAX_CORTEX_BEAD_DEGREE = 6)
        if deg[best_bead] >= 6:
            n_degree_block += 1
            continue
        if ma.stepping_mode == "grip_walk":
            if ma._bipolar_accepts(int(unbound[k]), best_bead, pos):
                n_bipolar_pass += 1
        else:
            n_bipolar_pass += 1

    return {
        "n_unbound": int(n_unbound),
        "p_bind_per_tick": float(p_bind),
        "k_on_per_s": float(p.head_actin_k_on),
        "batch_dt_s": float(p.batch_dt),
        "gate1_has_bead_candidate": int(n_has_candidate),
        "gate3_perp_eligible": int(n_perp_eligible),
        "gate4_bipolar_pass": int(n_bipolar_pass),
        "gate5_degree_blocked": int(n_degree_block),
        "naive_equilib_bound_frac": float(
            p.head_actin_k_on / (p.head_actin_k_on + p.head_actin_k_off0)
        ),
        "effective_acceptors_per_unbound_head": (
            float(n_bipolar_pass / max(n_unbound, 1))
        ),
    }


def kinetics_panel(ma, sim, p_cortex, ticks: int) -> dict:
    """Panels C+D: run the REAL updater for `ticks` ticks; record kinetics +
    whether s_grip advances on bound+loaded heads."""
    p = ma.p
    batch_steps = p.batch_steps
    rows = []
    prev_bind = ma._n_bind_total
    prev_break = ma._n_break_total
    prev_adv = ma._n_step_advances_total
    n_heads = ma._head_bound_to_actin.size

    for t in range(ticks):
        sim.run(batch_steps)
        snap = sim.state.get_snapshot()
        pos = _tagpos(snap)
        bound = ma._head_bound_to_actin >= 0
        n_bound = int(bound.sum())
        # per-head bond stretch r and F for the bound set
        F_over_stall = []
        s_grip_bound = []
        if n_bound > 0:
            bidx = np.flatnonzero(bound)
            htags = np.array([ma._head_global_tag(int(h)) for h in bidx])
            atags = ma._head_bound_to_actin[bidx]
            r = np.linalg.norm(pos[htags] - pos[atags], axis=1)
            F = p.k_head_actin * np.clip(r, 0.0, None)  # r0~0 grip_walk
            F_over_stall = (F / p.F_stall_per_head).tolist()
            s_grip_bound = ma._head_grip_s[bidx].tolist()
        nb = ma._n_bind_total
        nk = ma._n_break_total
        na = ma._n_step_advances_total
        rows.append({
            "tick": t,
            "n_bound": n_bound,
            "bound_frac": float(n_bound / n_heads),
            "binds_this_tick": int(nb - prev_bind),
            "breaks_this_tick": int(nk - prev_break),
            "step_advances_this_tick": int(na - prev_adv),
            "mean_F_over_stall": (
                float(np.mean(F_over_stall)) if F_over_stall else 0.0
            ),
            "mean_s_grip_nm": (
                float(np.mean(s_grip_bound) * 1e9) if s_grip_bound else 0.0
            ),
            "max_s_grip_nm": (
                float(np.max(s_grip_bound) * 1e9) if s_grip_bound else 0.0
            ),
        })
        prev_bind, prev_break, prev_adv = nb, nk, na

    # s_grip downstream test: among heads bound at the end AND loaded, did
    # s_grip move off zero? (compare first vs last tick on the bound subset)
    s_grip_trend = [r["mean_s_grip_nm"] for r in rows]
    return {
        "ell0_cortex_nm": p_cortex.rest_length * 1e9,
        "rows": rows,
        "s_grip_mean_trend_nm": s_grip_trend,
        "s_grip_advanced": bool(
            len(s_grip_trend) >= 2 and s_grip_trend[-1] > s_grip_trend[0] + 1e-6
        ),
        "steady_bound_frac": float(np.mean([r["bound_frac"] for r in rows[-3:]]))
        if len(rows) >= 1 else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--ticks", type=int, default=8)
    ap.add_argument(
        "--stepping-mode", choices=["grip_walk", "binned_r0", "config"],
        default="grip_walk",
        help="override stepping_mode (default grip_walk = the production mode)",
    )
    args = ap.parse_args()

    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = args.n_fil
    cfg["cortex"]["demo_mode"] = True
    if args.stepping_mode != "config":
        cfg.setdefault("cortex", {}).setdefault("myosin", {})
        cfg["cortex"]["myosin"]["stepping_mode"] = args.stepping_mode

    p_cortex = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)
    # Crosslinkers are not needed for the head->actin binding diagnosis (myosin
    # heads bind cortex actin, not xlinks). Include them only if available so the
    # built cortex matches production topology; otherwise omit.
    try:
        p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    except Exception:
        p_xl = None

    print(
        f"[build] n_fil={args.n_fil} beads/fil={p_cortex.beads_per_filament} "
        f"n_motors={p_myo.n_motors_per_cell} heads={2*p_myo.n_heads_per_side*p_myo.n_motors_per_cell} "
        f"mode={p_myo.stepping_mode} k_on={p_myo.head_actin_k_on}/s "
        f"k_off0={p_myo.head_actin_k_off0}/s batch_dt={p_myo.batch_dt:.3e}s "
        f"capture_perp={p_myo.head_actin_capture_perp*1e9:.0f}nm "
        f"max_bind_dist={p_myo.head_actin_max_bind_dist*1e9:.0f}nm "
        f"hoomd={hoomd.version.version}",
        flush=True,
    )

    device = hoomd.device.CPU()
    hw = build_cortex_full_simulation(
        p_cortex, p_xlinks=p_xl, p_myosin=p_myo,
        device=device, with_baoab=True,
        rng=np.random.default_rng(p_cortex.seed),
    )
    # build_cortex_full_simulation returns a dict of handles.
    sim = hw["sim"]
    ma = hw["myosin_action"]
    if ma is None:
        raise RuntimeError("Could not locate MyosinStepUpdater on the sim.")

    # The builder already appended + attached the updater; run(0) initializes it.
    sim.run(0)
    snap0 = sim.state.get_snapshot()

    print("[panel A] geometry at construction ...", flush=True)
    panel_a = geometry_panel(ma, snap0)
    print(json.dumps(panel_a, indent=2), flush=True)

    print("[panel B] one-tick gate accounting ...", flush=True)
    panel_b = gate_accounting(ma, snap0)
    print(json.dumps(panel_b, indent=2), flush=True)

    print(f"[panel C/D] {args.ticks}-tick kinetics ...", flush=True)
    panel_cd = kinetics_panel(ma, sim, p_cortex, args.ticks)
    print(json.dumps({k: v for k, v in panel_cd.items() if k != "rows"},
                     indent=2), flush=True)
    for r in panel_cd["rows"]:
        print("   ", r, flush=True)

    result = {
        "config": {
            "n_fil": args.n_fil,
            "beads_per_filament": int(p_cortex.beads_per_filament),
            "n_motors": int(p_myo.n_motors_per_cell),
            "n_heads": int(2 * p_myo.n_heads_per_side * p_myo.n_motors_per_cell),
            "stepping_mode": p_myo.stepping_mode,
            "k_on_per_s": float(p_myo.head_actin_k_on),
            "k_off0_per_s": float(p_myo.head_actin_k_off0),
            "capture_perp_nm": p_myo.head_actin_capture_perp * 1e9,
            "max_bind_dist_nm": p_myo.head_actin_max_bind_dist * 1e9,
            "head_rest_length_nm": p_myo.head_rest_length * 1e9,
            "ell0_cortex_nm": p_cortex.rest_length * 1e9,
        },
        "panelA_geometry": panel_a,
        "panelB_gate_accounting": panel_b,
        "panelCD_kinetics": panel_cd,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    print("RESULT_FILE " + str(OUT), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
