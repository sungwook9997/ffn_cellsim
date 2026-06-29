#!/usr/bin/env python
"""KU-3.5-ACTIVE GATE-B: bipolar-COMPLETION root-cause diagnostic (2026-06-05).

The KU-3.5-active floor was localised (h3_ku35_stresslet) to GENERATION via
INCOMPLETE BIPOLAR RECRUITMENT: in the full-cell connected-mesh cortex only
~18-26 % of ENGAGED minifilaments form a COMPLETE bipolar pair (both a bound +
head and a bound − head, on DIFFERENT filaments); the remaining ~74-82 % are
single-sided drags whose reaction is absorbed by the backbone and generate no
contractile dipole. The complete pairs themselves are coherence ≈ −1 (perfectly
contractile) — so the floor is NOT aggregation-loss, it is the SCARCITY of
completed pairs.

This module asks WHY completion is low and quantifies the limiting factor. It
does NOT change any physics (measurement-only) and creates NO shared state.

Two layers:

(1) ROOT-CAUSE METRIC — ``antiparallel_partner_availability`` (pure geometry,
    HOOMD-free, unit-tested on synthetic lattices). For a minifilament whose +
    side grips actin filament A, a COMPLETE bipolar pair requires the − side to
    find an ANTIPARALLEL partner filament B — one whose minus-end-ward tangent
    m̂_B has the OPPOSITE sign of projection on the rod axis û to A's — within
    the − side's physical head reach (``head_off`` laterally across the backbone
    + ``capture_perp`` + the backbone axial half-span). The headline number is
    the FRACTION of engageable + heads that HAVE ≥1 such antiparallel partner
    reachable. If this fraction is ≈ the observed ``frac_complete_pairs``, the
    floor is OVERLAP-SCARCITY (antiparallel actin is simply not present within
    reach in the sparse ×40-coarse mesoscale cortex) — a STRUCTURAL limit, not a
    binding-kinetics or reach mis-set one.

(2) LIVE SWEEP — sweep the recruitment knobs (n_motors density, head-actin reach
    via ``reach_scale`` on capture_perp+max_bind_dist, filament density n_fil)
    over a few values each on the EXISTING connected-mesh cortex build and report
    how ``frac_complete_pairs`` (from h3_ku35_stresslet) and the
    antiparallel-availability move. Identifies WHICH knob limits completion.

Run:  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_completion_diag --self-test
      conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_completion_diag --sweep --quick
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_production_device,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- #
# (1) Pure-geometry antiparallel-partner availability (HOOMD-free, unit-tested)
# --------------------------------------------------------------------------- #
@dataclass(slots=True)
class CompletionGeometry:
    """Inputs for the antiparallel-partner-availability metric.

    All SI metres. ``positions`` is the tag-ordered cortex actin bead positions
    reshaped ``(F, N, 3)``; ``tangents`` is the per-filament PLUS-end-ward unit
    tangent ``(F, 3)`` (cortex convention: bead 0 = minus end, tangent points
    minus→plus, so the minus-end-ward tangent is ``m̂ = −tangents``).
    """

    positions: np.ndarray         # (F, N, 3) cortex actin bead positions
    tangents: np.ndarray          # (F, 3) plus-end-ward unit tangent per filament
    head_off: float               # head_rest_length [m] lateral head offset
    capture_perp: float           # head_actin_capture_perp [m] perp eligibility
    backbone_half: float          # backbone_length/2 [m] axial half-span of rod


def _minus_tangents(g: CompletionGeometry) -> np.ndarray:
    """Minus-end-ward unit tangent m̂ per filament (= −plus-tangent)."""
    m = -np.asarray(g.tangents, dtype=np.float64)
    n = np.linalg.norm(m, axis=1, keepdims=True)
    return m / np.maximum(n, 1.0e-30)


def antiparallel_partner_availability(
    g: CompletionGeometry,
    *,
    rng: np.random.Generator | None = None,
    n_probe: int | None = None,
    reach_scale: float = 1.0,
) -> dict:
    """Fraction of engageable + heads that have an ANTIPARALLEL partner in reach.

    Model of a candidate minifilament (one per probed actin bead):

    * It sits ON a "host" filament A at one of A's beads; the rod axis û is A's
      tangent (``generate_cortex_myosin_layout`` actin-aware placement: axis =
      local actin tangent). The + side grips A.
    * The − side projects laterally ``±head_off`` across the backbone, with
      perpendicular reach ``capture_perp`` and axial reach ``backbone_half`` — so
      the − head can bind any actin BEAD within an ellipsoid-ish reach of the rod
      centre. We use the conservative spherical reach
      ``R_reach = (head_off + capture_perp + backbone_half) · reach_scale``.
    * A bead on filament B is an ANTIPARALLEL partner iff
      ``sign(m̂_B · û) == −sign(m̂_A · û)`` (opposite polarity along the rod) and
      B ≠ A and B is within ``R_reach`` of the rod centre.

    A + head that grips A is "completable" iff ≥1 antiparallel partner bead is in
    reach. The returned ``availability`` is the mean over probed beads. This is
    the structural ceiling on ``frac_complete_pairs`` — kinetics/Bell-Evans can
    only LOWER it.

    Args:
        g: geometry.
        rng: probe sampler (default seed 7).
        n_probe: number of host beads to probe (default min(F, 400)).
        reach_scale: multiply the − side reach (sweep knob = the head-actin reach
            magnitude; ``head_actin_max_bind_dist``/``capture_perp`` proxy).

    Returns:
        dict with availability, mean_n_partners, frac_reach_empty (no actin at
        all in reach), n_probed, R_reach, plus the per-probe arrays for plots.
    """
    if rng is None:
        rng = np.random.default_rng(7)
    pos = np.asarray(g.positions, dtype=np.float64)
    F, N, _ = pos.shape
    flat = pos.reshape(F * N, 3)
    fil_of_bead = np.repeat(np.arange(F), N)
    mhat = _minus_tangents(g)                       # (F, 3)

    R_reach = float(
        (g.head_off + g.capture_perp + g.backbone_half) * reach_scale
    )

    from scipy.spatial import cKDTree
    tree = cKDTree(flat)

    n_probe = int(min(F, 400) if n_probe is None else n_probe)
    host_fils = rng.choice(F, size=min(n_probe, F), replace=False)

    avail = np.zeros(host_fils.size, dtype=bool)
    n_partners = np.zeros(host_fils.size, dtype=np.int64)
    n_actin_in_reach = np.zeros(host_fils.size, dtype=np.int64)
    for k, A in enumerate(host_fils):
        # Place the rod centre at a random bead of host A; axis û = A tangent.
        j = int(rng.integers(0, N))
        centre = pos[A, j]
        u_hat = np.asarray(g.tangents[A], dtype=np.float64)
        un = np.linalg.norm(u_hat)
        if un <= 0:
            continue
        u_hat = u_hat / un
        sA = np.sign(float(mhat[A] @ u_hat))
        if sA == 0.0:
            sA = 1.0
        idx = np.asarray(tree.query_ball_point(centre, r=R_reach), dtype=np.int64)
        if idx.size == 0:
            continue
        fb = fil_of_bead[idx]
        other = fb != A
        n_actin_in_reach[k] = int(other.sum())
        if not other.any():
            continue
        fb_other = np.unique(fb[other])
        # antiparallel: m̂_B·û has the OPPOSITE sign of A's.
        sB = np.sign(np.einsum("ij,j->i", mhat[fb_other], u_hat))
        is_anti = sB == -sA
        n_partners[k] = int(is_anti.sum())
        avail[k] = bool(is_anti.any())

    return dict(
        availability=float(avail.mean()) if avail.size else float("nan"),
        mean_n_partners=float(n_partners.mean()) if n_partners.size else float("nan"),
        frac_reach_empty=(
            float(np.mean(n_actin_in_reach == 0)) if n_actin_in_reach.size
            else float("nan")
        ),
        n_probed=int(host_fils.size),
        R_reach=R_reach,
        _avail=avail,
        _n_partners=n_partners,
        _n_actin_in_reach=n_actin_in_reach,
    )


# --------------------------------------------------------------------------- #
# (2) Live connected-mesh sweep adapter
# --------------------------------------------------------------------------- #
def _geometry_from_topology(topology, p_myo, *, reach_scale: float = 1.0):
    """Build a CompletionGeometry from a live CortexTopology + resolved myosin."""
    return CompletionGeometry(
        positions=np.asarray(topology.positions, dtype=np.float64),
        tangents=np.asarray(topology.tangents, dtype=np.float64),
        head_off=float(p_myo.head_rest_length),
        capture_perp=float(p_myo.head_actin_capture_perp) * reach_scale,
        backbone_half=0.5 * float(p_myo.backbone_length),
    )


def _build_cortex_for_sweep(*, n_fil, n_motors, reach_scale, seed, device="gpu",
                            n_ticks=40, allow_cpu_dev: bool = False):
    """Build the connected-mesh cortex + myosin and populate the bipolar
    bookkeeping by running a few myosin RECRUITMENT ticks (de-novo binding).

    Returns (sim, p_cortex, p_myo, topology, myosin_action). Uses the same
    build path as ``mcf7_fullcell_stage1`` (connected_mesh=True) so the sweep
    measures the REAL recruitment, not a sandbox.

    The integrator is NOT advanced: GATE-B measures RECRUITMENT (which heads can
    bind a partner) on the construction geometry, so we call the myosin
    updater's ``act()`` directly on the static frame. This isolates the bipolar
    completion question from any relaxation/M-SHAKE dynamics (and avoids the
    fresh-frame M-SHAKE cold-start that needs a warm-up), exactly the read we
    want: of the heads that engage, how many find an antiparallel partner.
    """
    from copy import deepcopy
    from dataclasses import replace as _replace

    import yaml
    import hoomd

    require_production_device(
        device, allow_cpu_dev=allow_cpu_dev, hoomd_module=hoomd
    )

    from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
    from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
    from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation

    cfg = deepcopy(yaml.safe_load(open(PKG / "configs" / "phase1_h3.yaml")))
    cfg["cortex"]["R_cell"] = 7.5e-6                 # MCF7
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = int(n_motors)
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    cfg["cortex"]["myosin"]["backbone_length"] = 700.0e-9
    # reach knob: scale the head-actin binding reach (capture_perp + bond ceiling).
    base_perp = float(cfg["cortex"]["myosin"]["head_actin_capture_perp"])
    base_max = float(cfg["cortex"]["myosin"]["head_actin_max_bind_dist"])
    cfg["cortex"]["myosin"]["head_actin_capture_perp"] = base_perp * reach_scale
    cfg["cortex"]["myosin"]["head_actin_max_bind_dist"] = base_max * reach_scale

    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        rng=np.random.default_rng(seed), connected_mesh=True,
        cm_z_struct=3.7, cm_bundle_mult=2,
    )
    sim = hw["sim"]
    myo_act = hw.get("myosin_action")
    sim.run(0)
    # A few myosin RECRUITMENT ticks (de-novo k_on binding) directly on the
    # static construction frame so frac_complete is populated, WITHOUT advancing
    # the integrator (GATE-B is a recruitment-geometry read). myo_act is a
    # hoomd.custom.Action with _sim_ref set at build (attach()); act() reads the
    # snapshot, binds, and rewrites bonds via set_snapshot.
    if myo_act is not None:
        if getattr(myo_act, "_sim_ref", None) is None:
            myo_act._sim_ref = sim
        for tick in range(int(n_ticks)):
            myo_act.act(tick)
    return sim, p, p_myo, hw["topology"], myo_act


def sweep(*, knob: str, values, n_fil, n_motors, reach_scale, seed=1,
          device="gpu", allow_cpu_dev: bool = False, verbose=True,
          n_ticks=40):
    """Sweep one recruitment knob; return per-value completion + availability.

    knob ∈ {"n_motors", "reach_scale", "n_fil"}. The other two are held at the
    passed defaults.
    """
    from ffn_sim.scripts.h3_ku35_stresslet import stresslet_ledger

    rows = []
    for v in values:
        nm, nf, rs = n_motors, n_fil, reach_scale
        if knob == "n_motors":
            nm = int(v)
        elif knob == "n_fil":
            nf = int(v)
        elif knob == "reach_scale":
            rs = float(v)
        else:
            raise ValueError(f"unknown knob {knob!r}")
        sim, p, p_myo, topology, myo_act = _build_cortex_for_sweep(
            n_fil=nf, n_motors=nm, reach_scale=rs, seed=seed, device=device,
            n_ticks=n_ticks, allow_cpu_dev=allow_cpu_dev)
        sled = stresslet_ledger(
            sim, p_myo=p_myo, myosin_action=myo_act,
            beads_per_filament=p.beads_per_filament)
        summ = sled["summary"]
        geo = _geometry_from_topology(topology, p_myo, reach_scale=rs)
        avail = antiparallel_partner_availability(
            geo, rng=np.random.default_rng(seed + 11), reach_scale=1.0)
        row = dict(
            knob=knob, value=v, n_motors=nm, n_fil=nf, reach_scale=rs,
            n_engaged=summ["n_engaged_motors"],
            frac_complete=summ["frac_complete_pairs"],
            frac_single=summ["frac_single_sided"],
            n_antiparallel=summ["n_antiparallel"],
            availability=avail["availability"],
            mean_n_partners=avail["mean_n_partners"],
            frac_reach_empty=avail["frac_reach_empty"],
            R_reach=avail["R_reach"],
        )
        rows.append(row)
        if verbose:
            fc = row["frac_complete"]
            fc_s = f"{fc*100:5.1f}%" if np.isfinite(fc) else "  n/a"
            print(
                f"  [{knob}={v}] engaged={row['n_engaged']:3d} "
                f"complete={fc_s} single={row['frac_single']*100:5.1f}% "
                f"| antiparallel-avail={row['availability']*100:5.1f}% "
                f"meanP={row['mean_n_partners']:.2f} R_reach={row['R_reach']*1e9:.0f}nm",
                flush=True,
            )
        del sim
    return rows


def make_figure(rows_by_knob: dict, out_png: Path, *, title: str = "") -> None:
    """frac_complete_pairs + antiparallel-availability vs the swept knob."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    knobs = [k for k in rows_by_knob if rows_by_knob[k]]
    n = len(knobs)
    if n == 0:
        return
    fig, axes = plt.subplots(1, n, figsize=(5.0 * n, 4.2), squeeze=False)
    for ax, knob in zip(axes[0], knobs):
        rows = rows_by_knob[knob]
        x = [r["value"] for r in rows]
        fc = [r["frac_complete"] * 100 for r in rows]
        av = [r["availability"] * 100 for r in rows]
        ax.plot(x, fc, "o-", color="C3", label="frac_complete_pairs (measured)")
        ax.plot(x, av, "s--", color="C0",
                label="antiparallel-availability (structural ceiling)")
        ax.set_xlabel(knob)
        ax.set_ylabel("percent of engaged + heads")
        ax.set_ylim(0, 100)
        ax.set_title(knob)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(title or "KU-3.5-active GATE-B: bipolar completion vs knobs")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Self-test (synthetic lattices → known availability), HOOMD-free
# --------------------------------------------------------------------------- #
def _lattice(F_side: int, spacing: float, length: float, polarity_mix: str,
             rng: np.random.Generator):
    """Build a square lattice of straight filaments in the z=0 plane.

    Each filament is ``N=4`` beads along its tangent. ``polarity_mix``:
      * "all_plus" — every filament tangent = +x (NO antiparallel partners).
      * "checker"  — alternate filaments flip tangent → antiparallel neighbours
        everywhere.
    """
    N = 4
    F = F_side * F_side
    pos = np.zeros((F, N, 3), dtype=np.float64)
    tang = np.zeros((F, 3), dtype=np.float64)
    f = 0
    for ix in range(F_side):
        for iy in range(F_side):
            x0 = ix * spacing
            y0 = iy * spacing
            if polarity_mix == "all_plus":
                sgn = 1.0
            elif polarity_mix == "checker":
                sgn = 1.0 if (ix + iy) % 2 == 0 else -1.0
            else:
                raise ValueError(polarity_mix)
            t = np.array([sgn, 0.0, 0.0])
            tang[f] = t
            for j in range(N):
                pos[f, j] = [x0 + (sgn * j) * (length / (N - 1)), y0, 0.0]
            f += 1
    return pos, tang


def self_test(verbose: bool = True) -> bool:
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    rng = np.random.default_rng(0)
    head_off = 200e-9
    capture_perp = 210e-9
    backbone_half = 350e-9       # 700 nm /2
    reach = head_off + capture_perp + backbone_half  # ~760 nm

    # (A) All-parallel lattice, neighbours WITHIN reach → availability MUST be 0
    #     (there is no antiparallel partner anywhere, only same-sign tangents).
    spacing = 0.4 * reach       # neighbours comfortably in reach
    pos, tang = _lattice(6, spacing, length=300e-9, polarity_mix="all_plus",
                         rng=rng)
    g = CompletionGeometry(pos, tang, head_off, capture_perp, backbone_half)
    rA = antiparallel_partner_availability(g, rng=np.random.default_rng(1))
    if verbose:
        print(f"\n[A] all-parallel, neighbours in reach: avail={rA['availability']:.3f} "
              f"reach_empty={rA['frac_reach_empty']:.3f}")
    check("all-parallel → availability == 0", rA["availability"] == 0.0)
    check("all-parallel → reach not empty (actin IS present)",
          rA["frac_reach_empty"] == 0.0)

    # (B) Checkerboard polarity, neighbours in reach → availability MUST be ~1
    #     (every filament has antiparallel neighbours within reach).
    pos2, tang2 = _lattice(6, spacing, length=300e-9, polarity_mix="checker",
                           rng=rng)
    g2 = CompletionGeometry(pos2, tang2, head_off, capture_perp, backbone_half)
    rB = antiparallel_partner_availability(g2, rng=np.random.default_rng(2))
    if verbose:
        print(f"[B] checker, neighbours in reach: avail={rB['availability']:.3f} "
              f"meanP={rB['mean_n_partners']:.2f}")
    check("checker → availability == 1.0", rB["availability"] == 1.0)
    check("checker → mean partners > 0", rB["mean_n_partners"] > 0.0)

    # (C) Checkerboard BUT spacing >> reach → availability collapses to 0
    #     (antiparallel actin EXISTS but is OUT OF REACH = overlap-scarcity).
    pos3, tang3 = _lattice(6, 5.0 * reach, length=300e-9, polarity_mix="checker",
                           rng=rng)
    g3 = CompletionGeometry(pos3, tang3, head_off, capture_perp, backbone_half)
    rC = antiparallel_partner_availability(g3, rng=np.random.default_rng(3))
    if verbose:
        print(f"[C] checker, spacing>>reach: avail={rC['availability']:.3f} "
              f"reach_empty={rC['frac_reach_empty']:.3f}")
    check("sparse checker → availability == 0 (out of reach)",
          rC["availability"] == 0.0)
    check("sparse → reach empty (no actin near rod centre)",
          rC["frac_reach_empty"] > 0.5)

    # (D) reach_scale MONOTONICITY: growing the reach on the sparse checker MUST
    #     recover availability (reach, not polarity, was the limiter there).
    rC_big = antiparallel_partner_availability(
        g3, rng=np.random.default_rng(3), reach_scale=6.0)
    if verbose:
        print(f"[D] sparse checker + reach×6: avail={rC_big['availability']:.3f}")
    check("reach×6 recovers availability on sparse checker",
          rC_big["availability"] > rC["availability"])

    if verbose:
        print(f"\n{'='*56}\nCOMPLETION-DIAG SELF-TEST "
              f"{'PASSED' if ok else 'FAILED'}\n{'='*56}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--sweep", action="store_true",
                    help="live connected-mesh recruitment-knob sweep")
    ap.add_argument("--quick", action="store_true",
                    help="fewer sweep points (2/knob) at production scale "
                         "(n_fil=1000, n_motors=100) for a fast pass")
    add_production_device_args(ap, default="gpu")
    ap.add_argument("--fig", action="store_true", help="write the sweep figure")
    args = ap.parse_args()

    if args.self_test:
        return 0 if self_test() else 1

    if args.sweep:
        validate_production_device_args(ap, args)
        # 1000 fil / 100 motors = production scale; 40 recruitment ticks gives a
        # populated frac_complete (~24 %, the diagnosed floor). quick = fewer
        # sweep points (same scale) so a fast pass still measures the real floor.
        n_fil0, n_motors0, n_ticks = 1000, 100, 40
        if args.quick:
            motor_vals = [50, 200]
            reach_vals = [1.0, 4.0]
            fil_vals = [1000, 2000]
        else:
            motor_vals = [50, 100, 200]
            reach_vals = [1.0, 2.0, 4.0]
            fil_vals = [500, 1000, 2000]
        print("=== KU-3.5-active GATE-B: bipolar-completion sweep "
              f"(connected-mesh, device={args.device}) ===", flush=True)
        rows_by_knob = {}
        print("\n-- knob: n_motors (motor density) --", flush=True)
        rows_by_knob["n_motors"] = sweep(
            knob="n_motors", values=motor_vals, n_fil=n_fil0, n_motors=n_motors0,
            reach_scale=1.0, device=args.device,
            allow_cpu_dev=args.allow_cpu_dev, n_ticks=n_ticks)
        print("\n-- knob: reach_scale (head-actin reach) --", flush=True)
        rows_by_knob["reach_scale"] = sweep(
            knob="reach_scale", values=reach_vals, n_fil=n_fil0,
            n_motors=n_motors0, reach_scale=1.0, device=args.device,
            allow_cpu_dev=args.allow_cpu_dev, n_ticks=n_ticks)
        print("\n-- knob: n_fil (filament density) --", flush=True)
        rows_by_knob["n_fil"] = sweep(
            knob="n_fil", values=fil_vals, n_fil=n_fil0, n_motors=n_motors0,
            reach_scale=1.0, device=args.device,
            allow_cpu_dev=args.allow_cpu_dev, n_ticks=n_ticks)
        if args.fig:
            out = PKG / "outputs" / "h3" / "figs" / "ku35_completion_diag_sweep.png"
            make_figure(rows_by_knob, out,
                        title="KU-3.5-active GATE-B: bipolar completion vs recruitment knobs")
            print(f"\n[FIG] wrote {out}", flush=True)
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
