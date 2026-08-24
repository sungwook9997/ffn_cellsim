#!/usr/bin/env python3
r"""Net-new state dumper for the composed Active Cell — geometry + composed per-node force, for rendering.

Builds the SAME assembled cell as :mod:`ac.cell.driver` (``build_cell`` → ``AssembledCell``), advances the ONE
``--from-resting`` outer physical step, and writes an ``.npz`` of the post-step geometry together with the
COMPOSED per-node force magnitude the mechanics actually sees (the sum of Cytosim bending + Hookean
link-spring + MyosinForce + StericForce + PressureCoupling, exactly what ``driver._accumulate_all`` sums).
The dev-Mac viewer ``aleph/scripts/ac_cell_assembled_viz.py`` consumes the ``.npz`` and colours the cortex
by |F| so the Lead/PI can SEE where the mechanics acts (in particular the excluded-volume hotspots at the
builder's interpenetrating nodes).

This module is a READ-ONLY consumer: it imports ``build_cell`` (assemble) and the driver's frozen
``_accumulate_all`` / ``_residual_host`` / ``make_inner_solve`` + the fluid ``PhysicalScheduler``, and edits no
track / assemble / ff module. It additionally isolates the StericForce-only per-node force (a fresh scratch
accumulator, no state mutation) so the viewer can highlight the interpenetration explicitly.

Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000 (``build_cell`` constructs device state); the dev Mac
cannot launch these kernels. Emits a plain ``.npz`` + sidecar ``.json`` that the Mac renders with numpy only.

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.dump_state \
        --from-resting --out aleph/outputs/ac/cell_assembled/assembled_state.npz
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import PI_0_PA, AssembledCell, CellConfig, build_cell
from aleph.components.incumbent.driver import _residual_host, make_inner_solve
from aleph.components.fluid.scheduler import PhysicalScheduler

__all__ = ["dump_state"]


def _steric_only_mag(cell: AssembledCell) -> np.ndarray:
    """Per-actin-node |F| from the StericForce primitive ALONE (fresh scratch; no state mutation) [pN].

    Isolates the excluded-volume contribution so the viewer can mark the ~interpenetrating nodes the cortex
    builder produced (it places filaments without the 7 nm EV shell). Reads current ``cell.state.pos``.
    """
    if cell.steric is None:
        return np.zeros(cell.n_actin, np.float64)
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    cell.steric.accumulate(cell.state, f)
    wp.synchronize_device(cell.device)
    return np.linalg.norm(f.numpy()[: cell.n_actin], axis=1)


def _radii(pos: np.ndarray):
    c = pos.mean(axis=0)
    r = np.linalg.norm(pos - c, axis=1)
    return float(r.mean()), float(r.min()), float(r.max()), c


def dump_state(cfg: CellConfig, out_npz: str, *, n_inner: int = 60, dt_phys: float = 0.05,
               reshape_every: int = 20, max_inner_retries: int = 0,
               biot_cfl_safety: float = 0.9) -> dict:
    """Build the composed cell, advance one outer step, and dump geometry + composed |F| to ``out_npz``."""
    t0 = time.time()
    cell = build_cell(cfg)
    n_actin, n_total = cell.n_actin, cell.n_total
    print(f"[dump] built cell: n_actin={n_actin:,} n_total={n_total:,}", flush=True)

    # ── PRE-STEP snapshot (build-time geometry): composed |F| + steric-only |F| ──────────────────────
    pos_pre = cell.pos_d.numpy().copy()                                  # (n_total, 3) combined
    rms0, rmin0, rmax0, com0 = _radii(pos_pre[:n_actin])
    res_start, f_total_pre_all = _residual_host(cell, cell.pos_d, cell.f_d)
    f_total_pre = f_total_pre_all[:n_actin]                              # per-actin-node composed |F| [pN]
    f_steric_pre = _steric_only_mag(cell)                               # per-actin-node steric-only |F| [pN]
    n_interp_pre = int(np.count_nonzero(f_steric_pre > 0.0))
    print(f"[dump] PRE-STEP  max|F|={res_start:.4g} pN  interpenetrating nodes(steric>0)={n_interp_pre:,}",
          flush=True)

    # ── WHERE THE RESIDUAL ACTUALLY IS, and why this block exists ───────────────────────────────────
    #
    # `res_start` is the max over ALL n_total nodes; `f_total_pre` above is sliced to `[:n_actin]`, and
    # everything downstream — the percentile table, the radial-shell medians, the "max/p50 = 1.50"
    # uniformity statistic that STATE.md's C-2 row cites to EXCLUDE construction — is computed from
    # that slice. Measured 2026-08-10 on a full native dump: the actin slice maxes at 0.0356 pN while
    # res_start is 47.41 pN, a factor of 1330. The residual is not on any actin node, so the statistic
    # that excluded construction was computed on the subset that does not carry it.
    #
    # This block adds no physics and changes no force. It partitions the SAME `f_total_pre_all` the
    # residual was taken from, so the answer cannot disagree with `res_start` by construction.
    n_myo = int(cell.ledger.get("n_myosin_particles", 0))
    n_nuc = int(cell.nucleus.n_verts) if cell.nucleus is not None else 0
    n_mem = int(cell.membrane.n_verts) if cell.membrane is not None else 0
    _blocks = [("actin", 0, n_actin), ("myosin", n_actin, n_actin + n_myo),
               ("nucleus", n_actin + n_myo, n_actin + n_myo + n_nuc),
               ("membrane", n_actin + n_myo + n_nuc, n_actin + n_myo + n_nuc + n_mem)]
    residual_by_block = {}
    for _name, _lo, _hi in _blocks:
        _hi = min(_hi, f_total_pre_all.size)
        if _hi <= _lo:
            continue
        _b = f_total_pre_all[_lo:_hi]
        residual_by_block[_name] = {
            "n": int(_b.size), "node_range": [int(_lo), int(_hi)],
            "max": float(_b.max()), "p50": float(np.percentile(_b, 50)),
            "p99": float(np.percentile(_b, 99)),
            "argmax_global_node": int(_lo + int(np.argmax(_b))),
            "fraction_of_sum_F2": float((_b ** 2).sum() / max((f_total_pre_all ** 2).sum(), 1e-300)),
        }
    _owner = max(residual_by_block, key=lambda k: residual_by_block[k]["max"])
    residual_by_block["_carries_the_maximum"] = _owner
    residual_by_block["_covered"] = int(sum(v["n"] for k, v in residual_by_block.items()
                                            if isinstance(v, dict)))
    residual_by_block["_n_total"] = int(f_total_pre_all.size)
    print("[dump] RESIDUAL BY BLOCK (max|F| pN): "
          + "  ".join(f"{k}={v['max']:.4g}" for k, v in residual_by_block.items() if isinstance(v, dict))
          + f"   -> carried by {_owner!r}"
          + (f"   ⚠ blocks cover {residual_by_block['_covered']:,} of "
             f"{residual_by_block['_n_total']:,} nodes"
             if residual_by_block["_covered"] != residual_by_block["_n_total"] else ""), flush=True)

    # ── ONE outer physical step (mirror driver.run_from_resting) ─────────────────────────────────────
    n_sub, iters, attempts = 0, 0, 0
    residual_candidate = res_start
    inner_converged = False
    inner_finite = True
    fluid_finite = True
    outer_accepted = False
    outer_rolled_back = True
    if cell.substrate is not None and cell.domain is not None:
        sched = PhysicalScheduler(
            substrate=cell.substrate, domain=cell.domain, membrane_bc=cell.membrane_bc,
            inner_solve=make_inner_solve(cell, n_inner, reshape_every, max_inner_retries),
            osmotic_difference=lambda t: PI_0_PA,
        )
        rep = sched.outer_step(dt_phys, biot_cfl_safety=biot_cfl_safety).readback()
        n_sub, iters, attempts = rep.n_biot_subcycles, rep.inner_iters, rep.inner_attempts
        residual_candidate = rep.inner_residual
        inner_converged = rep.inner_converged
        inner_finite = rep.inner_finite
        fluid_finite = rep.fluid_finite
        outer_accepted = rep.outer_accepted
        outer_rolled_back = rep.outer_rolled_back
    else:
        inner_solve = make_inner_solve(cell, n_inner, reshape_every, max_inner_retries)
        inner_d = inner_solve(dt_phys)
        inner_solve.commit_irreversible(inner_d.converged_d)  # type: ignore[attr-defined]
        wp.synchronize_device(cell.device)
        residual_candidate = float(inner_d.residual_d.numpy()[0])
        iters = int(inner_d.iters_d.numpy()[0])
        attempts = int(inner_d.attempts_d.numpy()[0])
        inner_converged = bool(inner_d.converged_d.numpy()[0])
        inner_finite = bool(inner_d.finite_d.numpy()[0])
        outer_accepted = inner_converged
        outer_rolled_back = not inner_converged

    # ── POST-STEP snapshot: the state of the cell that just ran --from-resting ───────────────────────
    pos_post = cell.pos_d.numpy().copy()
    rms1, rmin1, rmax1, com1 = _radii(pos_post[:n_actin])
    residual_committed, f_total_post_all = _residual_host(cell, cell.pos_d, cell.f_d)
    f_total_post = f_total_post_all[:n_actin]
    f_steric_post = _steric_only_mag(cell)
    n_interp_post = int(np.count_nonzero(f_steric_post > 0.0))
    com_drift = float(np.linalg.norm(com1 - com0))
    pos_finite = bool(np.all(np.isfinite(pos_post)))
    force_finite = bool(np.all(np.isfinite(f_total_post_all)))
    stable = bool(pos_finite and force_finite and inner_finite and fluid_finite
                  and inner_converged and outer_accepted and not outer_rolled_back
                  and np.isfinite(residual_committed)
                  and residual_committed <= max(2.0 * res_start, res_start + 1.0)
                  and com_drift < 0.5 and 0.5 * rms0 < rms1 < 2.0 * rms0 and rmax1 < 3.0 * rms0)
    print(f"[dump] POST-STEP candidate max|F|={residual_candidate:.4g} pN  "
          f"committed max|F|={residual_committed:.4g} pN  "
          f"accepted={outer_accepted} rolled_back={outer_rolled_back}  "
          f"interpenetrating nodes(steric>0)={n_interp_post:,}  "
          f"STABLE={stable}", flush=True)

    # ── myosin topology (global indices into the combined [actin|myosin] node array) ─────────────────
    if cell.myosin is not None:
        mf_bb = cell.myosin.backbone_bonds.numpy().astype(np.int32)      # (Kbb, 2)
        mf_hb = cell.myosin.head_bonds.numpy().astype(np.int32)          # (Khb, 2)
        mf_hn = cell.myosin.head_node.numpy().astype(np.int32)           # (2H,)
    else:
        mf_bb = np.zeros((0, 2), np.int32)
        mf_hb = np.zeros((0, 2), np.int32)
        mf_hn = np.zeros(0, np.int32)

    xl = cell.xl_d.numpy().astype(np.int32) if cell.n_xl else np.zeros((0, 2), np.int32)

    # ── deformable-mesh compartments (milestone-2): node ranges + LOCAL mesh faces + tether pairs ─────
    def _mesh_dump(comp, prefix):
        if comp is None:
            return {f"{prefix}_off": np.int64(-1), f"{prefix}_nverts": np.int64(0),
                    f"{prefix}_faces": np.zeros((0, 3), np.int32)}
        off = int(comp.node_off)
        return {f"{prefix}_off": np.int64(off), f"{prefix}_nverts": np.int64(comp.n_verts),
                f"{prefix}_faces": (comp.faces_d.numpy().astype(np.int32) - off)}   # LOCAL indices
    nuc_dump = _mesh_dump(cell.nucleus, "nucleus")
    mem_dump = _mesh_dump(cell.membrane, "membrane")
    # LINC tethers (nucleus node <-> cortex anchor, GLOBAL indices) for the coupling layer
    if cell.nucleus is not None and cell.nucleus.n_linc:
        linc = np.stack([cell.nucleus.linc_n_d.numpy(), cell.nucleus.linc_a_d.numpy()], axis=1).astype(np.int32)
    else:
        linc = np.zeros((0, 2), np.int32)

    # ── report + ledger (JSON-friendly) ─────────────────────────────────────────────────────────────
    def _pct(a):
        a = a[np.isfinite(a)]
        if not a.size:
            return {}
        return {"min": float(a.min()), "p50": float(np.percentile(a, 50)),
                "p99": float(np.percentile(a, 99)), "p999": float(np.percentile(a, 99.9)),
                "max": float(a.max())}

    report = {
        "source_note": "composed AssembledCell after one transactional --from-resting outer-step attempt",
        "n_actin": n_actin, "n_total": n_total,
        "n_myosin_particles": int(cell.ledger.get("n_myosin_particles", 0)),
        "n_fibers": int(cell.n_fibers), "n_crosslinks": int(cell.n_xl),
        "with_steric": cfg.with_steric, "with_myosin": cfg.with_myosin, "with_pressure": cfg.with_pressure,
        "dt_phys": dt_phys, "n_inner": n_inner, "max_inner_retries": max_inner_retries,
        "n_biot_subcycles": int(n_sub), "inner_iters": int(iters), "inner_attempts": int(attempts),
        "residual_start_pN": res_start,
        "residual_candidate_pN": float(residual_candidate),
        "residual_committed_pN": float(residual_committed),
        "residual_end_pN": float(residual_committed),
        "inner_converged": inner_converged, "inner_finite": inner_finite, "fluid_finite": fluid_finite,
        "outer_accepted": outer_accepted, "outer_rolled_back": outer_rolled_back,
        "com_drift_um": com_drift, "r_mean_start_um": rms0, "r_mean_end_um": rms1,
        "r_min_end_um": rmin1, "r_max_end_um": rmax1,
        "n_interpenetrating_pre": n_interp_pre, "n_interpenetrating_post": n_interp_post,
        # ⚠ f_total_*_pN are the ACTIN SLICE only. `residual_*_pN` is the max over ALL n_total nodes,
        # and the two differ by three orders of magnitude — see `residual_by_block`, which is the
        # partition that says which block actually carries it.
        "f_total_pre_pN": _pct(f_total_pre), "f_total_post_pN": _pct(f_total_post),
        "f_total_pre_scope": "actin slice [:n_actin] ONLY — not the array residual_start_pN is taken from",
        "residual_by_block": residual_by_block,
        "f_steric_pre_pN": _pct(f_steric_pre[f_steric_pre > 0]),
        "f_steric_post_pN": _pct(f_steric_post[f_steric_post > 0]),
        "pos_finite": pos_finite, "force_finite": force_finite, "stable": stable,
        "steric_sigma_um": float(cell.steric.sigma) if cell.steric else None,
        "steric_force_cap_pN": cfg.steric_force_cap,
        "wall_s": time.time() - t0,
    }

    # ── write .npz (f32 geometry/force; halves size, sub-nm precision irrelevant for rendering) ──────
    np.savez_compressed(
        out_npz,
        pos_pre=pos_pre.astype(np.float32), pos_post=pos_post.astype(np.float32),
        n_actin=np.int64(n_actin), n_total=np.int64(n_total),
        f_total_pre=f_total_pre.astype(np.float32), f_total_post=f_total_post.astype(np.float32),
        f_steric_pre=f_steric_pre.astype(np.float32), f_steric_post=f_steric_post.astype(np.float32),
        fiber_offsets=cell.foff_d.numpy().astype(np.int32),
        xl=xl, mf_backbone_bonds=mf_bb, mf_head_bonds=mf_hb, mf_head_node=mf_hn,
        linc=linc, **nuc_dump, **mem_dump,
        report_json=np.array(json.dumps(report), dtype=object),
        ledger_json=np.array(json.dumps(cell.ledger, default=float), dtype=object),
    )
    with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
        json.dump({"report": report, "ledger": cell.ledger}, fh, indent=2, default=float)
    print(f"[dump] wrote {out_npz}", flush=True)
    return report


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Dump composed Active Cell state after one --from-resting step.")
    p.add_argument("--from-resting", action="store_true", help="(default behaviour) run the resting outer step")
    p.add_argument("--n-filaments", type=int, default=70686, help="cortical F-actin count (full native=70686)")
    p.add_argument("--n-inner", type=int, default=60, help="inner mechanical iterations for the outer step")
    p.add_argument("--max-inner-retries", type=int, default=0,
                   help="additional same-time mechanics chunks after n-inner")
    p.add_argument("--dt-phys", type=float, default=0.05, help="outer physical step [s]")
    p.add_argument("--no-myosin", action="store_true", help="omit MyosinForce composition")
    p.add_argument("--no-steric", action="store_true", help="omit StericForce (DEFAULT ON — shows hotspots)")
    p.add_argument("--no-pressure", action="store_true", help="omit the Biot substrate + PressureCoupling")
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--out", type=str, required=True, help="output .npz path")
    return p


def main() -> None:
    args = _build_argparser().parse_args()
    cfg = CellConfig(
        n_filaments=args.n_filaments, with_myosin=not args.no_myosin, with_steric=not args.no_steric,
        with_pressure=not args.no_pressure, device=args.device)
    print(f"[dump] building composed cell: n_filaments={cfg.n_filaments} myosin={cfg.with_myosin} "
          f"steric={cfg.with_steric} pressure={cfg.with_pressure} device={cfg.device}", flush=True)
    dump_state(
        cfg, args.out, n_inner=args.n_inner, dt_phys=args.dt_phys,
        max_inner_retries=args.max_inner_retries,
    )


if __name__ == "__main__":
    main()
