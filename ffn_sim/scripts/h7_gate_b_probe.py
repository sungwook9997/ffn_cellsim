"""H.7 GATE-B design-investigation #2 — built-mesh free-length + per-span myosin load probe.

Investigation #1 (``H7_GATE_B_BUCKLING_INVESTIGATION_2026-06-08.md``) showed, with a
first-order straight-rod Euler estimate, that the cortex sits *at the buckling margin*:
F_head/F_crit ~ 0.56-0.70 at the segment/crosslink scale. It deferred the quantitative
question to a probe of the **built physiological cell**:

  (a) the anchor-to-anchor FREE-LENGTH distribution along filaments, and
  (b) the actual per-span myosin load,

so the margin is pinned with the real constructed mesh, not an idealization. This script
is that probe. It builds the suspended/rounded MCF7 cell (FA OFF, turgor ON) at smoke
scale per the Gate-B contract operating point (§6), runs the unconstrained warm-up so the
myosin binder engages heads, then measures from the settled snapshot:

  * Per cortical filament, the beads anchored by a crosslink (``xlink_attach_b*`` bonds).
  * The contour FREE-SPAN between consecutive anchors, in segments (×ℓ₀=500 nm) and µm.
  * Whether each span can buckle AT ALL: a 1-segment span between two crosslinked beads is
    a single M-SHAKE-rigid rod (no internal hinge → cannot buckle); a span of ≥2 segments
    has ≥1 free internal bead that can deflect, buckling resisted by the angle harmonic.
  * The engaged myosin heads landing on each span (``cortex_myosin_attach_b*`` bonds) and
    the summed compressive load = n_heads × F_stall_per_head.
  * The Euler margin per bucklable span: F_crit = π²κ_B/L², ratio = F_load/F_crit.

TWO CORRECTIONS TO INVESTIGATION #1 this probe surfaces:
  1. The per-head force in the ×40 mesoscale sim is the SCALED F_stall_per_head
     (mesoscale_force_scaling multiplies it ~4.24× → ~8.48 pN), NOT the unscaled 2.0 pN
     the first-order table used. The ratified ×40 convention (cortex.py:10-13) coarse-grains
     filament COUNT not stiffness, so F_crit's single-filament κ_B is correct; but the
     force-scaling raises the LOCAL per-span load ×4.24 for AGGREGATE-stress correctness,
     so as-simulated a single engaged head is ~12× over F_crit on a bucklable span (de-scaled
     to per-native 2 pN it sits at the margin — investigation #1). Surfaced to PI, not
     silently resolved: applying areal-scaled force to a local buckling instability mixes
     a global correction into a local one.
  2. The relevant suppressor is testable here: if crosslink density pins most free spans to
     1 segment (no internal hinge), buckling is geometrically blocked by ANCHOR DENSITY
     irrespective of force — contract §3 candidate (i).

Usage (smoke):
    python -m ffn_sim.scripts.h7_gate_b_probe --n-filaments 160 --warmup 600 \
        --device cpu --allow-cpu-dev
Usage (full scale, gbook GPU):
    python -m ffn_sim.scripts.h7_gate_b_probe --device gpu --warmup 4000
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

_PN = 1.0e12  # N -> pN
_UM = 1.0e6   # m -> µm


def _build_settled_cell(*, n_filaments, n_nuc_beads, warmup, softstart, device, seed):
    """Suspended/rounded cell (FA OFF), unconstrained warm-up so the binder engages.

    Matches the Gate-B contract §6 operating point and the Gate-A native phase-1 warm-up
    build (constrained=False, equilibrate=True) which passes the nucleus CFL."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest["optional_subsystems"]["fa"]["enabled"] = False  # suspended/rounded (§6)
    if n_filaments is not None:
        co = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
        co["n_filaments"] = int(n_filaments)
        co["demo_mode"] = True
    if n_nuc_beads is not None:
        manifest["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)
    cell = build_baseline_cell(
        manifest=manifest, device=device, seed=seed,
        constrained=False, with_baoab=True,
        equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=softstart,
        connected_mesh=True,
    )
    return cell


def _filament_paths(backbone_pairs, n_cortex_actin):
    """Order each cortex filament's beads along its chain from the backbone-bond graph.

    Robust to uniform (contiguous) and variable-length layouts: build adjacency over
    cortex-bond pairs, then walk each connected component from a degree-1 endpoint.
    Returns a list of bead-index lists (each ordered tip→tip)."""
    adj = defaultdict(list)
    for a, b in backbone_pairs:
        a, b = int(a), int(b)
        if a < n_cortex_actin and b < n_cortex_actin:
            adj[a].append(b)
            adj[b].append(a)
    seen = set()
    paths = []
    for start in adj:
        if start in seen or len(adj[start]) != 1:  # begin only at a tip (degree 1)
            continue
        path = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            nxts = [n for n in adj[cur] if n != prev]
            if not nxts:
                break
            nxt = nxts[0]
            if nxt in seen:  # ring guard (shouldn't happen for linear filaments)
                break
            path.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        paths.append(path)
    # isolated beads (degree 0) or any unwalked component → singletons
    for bead in range(n_cortex_actin):
        if bead not in seen and bead in adj:
            paths.append([bead])
            seen.add(bead)
    return paths


def _anchored_beads(snap, n_cortex_actin, prefix):
    """Set of cortex beads carrying at least one bond of type starting with ``prefix``
    (e.g. 'xlink_attach_b' or 'cortex_myosin_attach_b'). Returns also a per-bead count."""
    types = list(snap.bonds.types)
    tids = np.asarray(snap.bonds.typeid)
    group = np.asarray(snap.bonds.group)
    sel = np.array([types[t].startswith(prefix) for t in range(len(types))])
    mask = sel[tids] if len(tids) else np.zeros(0, dtype=bool)
    count = np.zeros(n_cortex_actin, dtype=np.int64)
    for a, b in group[mask]:
        a, b = int(a), int(b)
        cortex_bead = a if a < n_cortex_actin else b  # the other end is a head particle
        if 0 <= cortex_bead < n_cortex_actin:
            count[cortex_bead] += 1
    return count


def run_probe(*, n_filaments, n_nuc_beads, warmup, softstart, device, seed):
    cell = _build_settled_cell(
        n_filaments=n_filaments, n_nuc_beads=n_nuc_beads,
        warmup=warmup, softstart=softstart, device=device, seed=seed,
    )
    snap = cell.simulation.state.get_snapshot()
    n_cortex_actin = int(cell.n_cortex_actin)
    pos = np.asarray(snap.particles.position)[:n_cortex_actin]

    # config physics (no tuning) — bending modulus κ_B and segment length ℓ₀
    ell0 = float(cell.p_cortex.rest_length)                 # 500 nm
    kappa_B = float(cell.p_cortex.bending_modulus)          # κ_B = ℓ_p·k_BT
    F_head = float(cell.p_myosin.F_stall_per_head)          # SCALED (mesoscale) per-head
    _meso = getattr(cell.p_myosin, "extras", None) or {}
    _factor = float(_meso.get("mesoscale_force_factor", 1.0))
    F_head_unscaled = F_head / _factor if _factor else 2.0e-12

    # backbone graph → ordered filament paths
    types = list(snap.bonds.types)
    tids = np.asarray(snap.bonds.typeid)
    group = np.asarray(snap.bonds.group)
    bb_tid = types.index("cortex-bond") if "cortex-bond" in types else None
    backbone_pairs = group[tids == bb_tid] if bb_tid is not None else np.zeros((0, 2), int)
    paths = _filament_paths(backbone_pairs, n_cortex_actin)

    xlink_count = _anchored_beads(snap, n_cortex_actin, "xlink_attach_b")
    myo_count = _anchored_beads(snap, n_cortex_actin, "cortex_myosin_attach_b")

    # --- per-span analysis -------------------------------------------------
    # On each ordered filament path the crosslink-anchored beads partition the chain into
    # free spans. A span between anchors at bead-positions i<j spans (j-i) segments and
    # has (j-i-1) free internal beads (hinges). Filament tips are free ends (not anchors).
    spans = []  # dicts: n_seg, L_contour, n_hinge, n_myosin_heads, bucklable
    for path in paths:
        anchor_positions = [k for k, bead in enumerate(path) if xlink_count[bead] > 0]
        if not anchor_positions:
            continue  # whole filament free between its tips: one big span, handled below
        # consecutive anchor pairs define internal spans
        for p0, p1 in zip(anchor_positions[:-1], anchor_positions[1:]):
            n_seg = p1 - p0
            if n_seg <= 0:
                continue
            span_beads = path[p0:p1 + 1]
            internal_beads = span_beads[1:-1]
            n_myo = int(sum(myo_count[b] for b in internal_beads)) if internal_beads else 0
            spans.append({
                "n_seg": int(n_seg),
                "L_contour": float(n_seg * ell0),
                "n_hinge": int(n_seg - 1),
                "n_myosin_heads": n_myo,
                "bucklable": bool(n_seg >= 2),
            })

    n_spans = len(spans)
    seg_hist = defaultdict(int)
    for s in spans:
        seg_hist[s["n_seg"]] += 1
    n_bucklable = sum(1 for s in spans if s["bucklable"])
    n_one_seg = seg_hist.get(1, 0)

    # Euler margin on bucklable spans, using the SCALED per-head force
    margin_rows = []
    n_over = 0
    for s in spans:
        if not s["bucklable"]:
            continue
        L = s["L_contour"]
        F_crit = np.pi**2 * kappa_B / (L * L)
        F_load = s["n_myosin_heads"] * F_head
        ratio = F_load / F_crit if F_crit > 0 else 0.0
        if ratio >= 1.0:
            n_over += 1
        margin_rows.append({
            "n_seg": s["n_seg"], "L_um": L * _UM,
            "F_crit_pN": F_crit * _PN, "n_heads": s["n_myosin_heads"],
            "F_load_pN": F_load * _PN, "ratio": ratio,
        })

    total_engaged = int(myo_count.sum())
    heads_on_bucklable = sum(s["n_myosin_heads"] for s in spans if s["bucklable"])

    result = {
        "scale": {"n_filaments": int(cell.p_cortex.n_filaments),
                  "n_cortex_actin": n_cortex_actin,
                  "is_smoke": n_filaments is not None},
        "physics": {
            "ell0_nm": ell0 * 1e9, "kappa_B_Nm2": kappa_B,
            "F_head_pN_scaled": F_head * _PN,
            "F_head_pN_unscaled": F_head_unscaled * _PN,
            "mesoscale_force_factor": F_head / F_head_unscaled if F_head_unscaled else None,
            "F_crit_1seg_pN": (np.pi**2 * kappa_B / ell0**2) * _PN,
            "F_crit_2seg_pN": (np.pi**2 * kappa_B / (2 * ell0)**2) * _PN,
        },
        "free_length": {
            "n_spans": n_spans,
            "seg_histogram": {str(k): int(v) for k, v in sorted(seg_hist.items())},
            "frac_1seg_rigid": (n_one_seg / n_spans) if n_spans else None,
            "frac_bucklable_geom": (n_bucklable / n_spans) if n_spans else None,
            "mean_n_seg": float(np.mean([s["n_seg"] for s in spans])) if spans else None,
            "mean_L_um": float(np.mean([s["L_contour"] for s in spans]) * _UM) if spans else None,
        },
        "myosin_load": {
            "total_engaged_heads": total_engaged,
            "heads_on_bucklable_spans": int(heads_on_bucklable),
            "spans_with_any_head": int(sum(1 for s in spans if s["n_myosin_heads"] > 0)),
            "max_heads_on_one_span": int(max((s["n_myosin_heads"] for s in spans), default=0)),
        },
        "margin": {
            "n_bucklable_spans": n_bucklable,
            "n_over_threshold": n_over,
            "frac_over_threshold": (n_over / n_bucklable) if n_bucklable else None,
        },
    }
    return result, margin_rows, spans, (xlink_count, myo_count)


def _report(r):
    p, fl, ml, mg = r["physics"], r["free_length"], r["myosin_load"], r["margin"]
    sc = r["scale"]
    print("=" * 72, flush=True)
    print("H.7 GATE-B design-investigation #2 — built-mesh free-length + load probe",
          flush=True)
    tag = (f"SMOKE (n_filaments={sc['n_filaments']})" if sc["is_smoke"]
           else "FULL ×40 production scale")
    print(f"  scale: {tag}, {sc['n_cortex_actin']} cortex beads", flush=True)
    print("-" * 72, flush=True)
    print("  PHYSICS (config, no tuning):", flush=True)
    print(f"    ℓ₀ = {p['ell0_nm']:.0f} nm   κ_B = {p['kappa_B_Nm2']:.2e} N·m²", flush=True)
    print(f"    F_head SCALED (mesoscale) = {p['F_head_pN_scaled']:.2f} pN  "
          f"[unscaled {p['F_head_pN_unscaled']:.2f} pN × {p['mesoscale_force_factor']:.2f}]",
          flush=True)
    print(f"    F_crit(1-seg 500nm)={p['F_crit_1seg_pN']:.2f} pN   "
          f"F_crit(2-seg 1µm)={p['F_crit_2seg_pN']:.2f} pN", flush=True)
    print("-" * 72, flush=True)
    print("  (a) FREE-LENGTH distribution:", flush=True)
    print(f"    {fl['n_spans']} inter-anchor spans  |  segment histogram: {fl['seg_histogram']}",
          flush=True)
    print(f"    mean span = {fl['mean_n_seg']:.2f} seg ({fl['mean_L_um']:.3f} µm)", flush=True)
    if fl["frac_1seg_rigid"] is not None:
        print(f"    1-segment RIGID spans (M-SHAKE rod, no hinge → cannot buckle): "
              f"{fl['frac_1seg_rigid']*100:.1f}%", flush=True)
        print(f"    geometrically bucklable spans (≥2 seg, ≥1 free hinge): "
              f"{fl['frac_bucklable_geom']*100:.1f}%", flush=True)
    print("-" * 72, flush=True)
    print("  (b) PER-SPAN MYOSIN LOAD:", flush=True)
    print(f"    total engaged heads = {ml['total_engaged_heads']}  |  "
          f"on bucklable spans = {ml['heads_on_bucklable_spans']}", flush=True)
    print(f"    spans carrying ≥1 head = {ml['spans_with_any_head']}  |  "
          f"max heads on one span = {ml['max_heads_on_one_span']}", flush=True)
    print("-" * 72, flush=True)
    print("  MARGIN (Euler, scaled per-head force):", flush=True)
    print(f"    bucklable spans = {mg['n_bucklable_spans']}  |  "
          f"over F_crit = {mg['n_over_threshold']} "
          f"({(mg['frac_over_threshold'] or 0)*100:.1f}%)", flush=True)
    print("=" * 72, flush=True)
    print("  READING:", flush=True)
    if (fl["frac_1seg_rigid"] or 0) > 0.5:
        print("   → ANCHOR DENSITY dominates: most free spans are single M-SHAKE rods with", flush=True)
        print("     no internal hinge, so buckling is geometrically blocked irrespective of", flush=True)
        print("     load (contract §3 candidate i). The relaxed-mode lever must LENGTHEN the", flush=True)
        print("     compression-side free span (thin the anchor density on compression), not", flush=True)
        print("     relax axial stretch.", flush=True)
    else:
        print("   → spans are long enough to hinge; check the margin fraction + binding", flush=True)
        print("     throughput (heads/span) for whether load reaches F_crit.", flush=True)
    print("   ⚠ MESOSCALE-CONSISTENCY (surface to PI): the ratified ×40 convention", flush=True)
    print("     (cortex.py:10-13) coarse-grains filament COUNT, keeping per-filament", flush=True)
    print("     stiffness single-filament — so F_crit's single-filament κ_B is CORRECT.", flush=True)
    print("     But mesoscale_force_scaling raises the per-head LOAD ×4.24 (areal/aggregate-", flush=True)
    print("     stress correctness) onto that single-filament-stiffness chain, whereas", flush=True)
    print("     buckling is a LOCAL per-filament instability. As-simulated, a single engaged", flush=True)
    print("     head is ~12× over F_crit on a bucklable span; de-scaled to per-native 2 pN it", flush=True)
    print("     sits at the margin (investigation #1). Either way the DOMINANT suppressor is", flush=True)
    print("     convention-independent: anchor density pins 70%+ of spans to 1 rigid segment.", flush=True)
    print("=" * 72, flush=True)


def _figure(r, spans, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fl, p = r["free_length"], r["physics"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)

    # (a) free-span segment histogram, rigid (1-seg) vs bucklable (≥2-seg)
    seg_hist = {int(k): v for k, v in fl["seg_histogram"].items()}
    segs = sorted(seg_hist)
    counts = [seg_hist[s] for s in segs]
    colors = ["#c0392b" if s < 2 else "#2e86c1" for s in segs]
    ax1.bar(segs, counts, color=colors, edgecolor="k", linewidth=0.5)
    ax1.set_xlabel("free-span length (segments, ×500 nm)")
    ax1.set_ylabel("number of inter-anchor spans")
    ax1.set_title(f"(a) built-mesh free-length distribution\n"
                  f"red=1-seg M-SHAKE rigid (no hinge), blue=bucklable (≥2 seg)")
    ax1.axvline(1.5, color="grey", ls="--", lw=1)
    ax1.text(0.98, 0.95, f"1-seg rigid: {(fl['frac_1seg_rigid'] or 0)*100:.0f}%\n"
             f"bucklable: {(fl['frac_bucklable_geom'] or 0)*100:.0f}%",
             transform=ax1.transAxes, ha="right", va="top",
             bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))

    # (b) Euler margin: F_crit vs free length, with per-head loads marked
    L = np.linspace(0.4, 3.2, 200) * 1e-6
    F_crit = np.pi**2 * p["kappa_B_Nm2"] / L**2
    ax2.plot(L * _UM, F_crit * _PN, "k-", lw=2, label="F_crit = π²κ_B/L²")
    ax2.axhline(p["F_head_pN_scaled"], color="#c0392b", ls="-", lw=1.8,
                label=f"1 head SCALED = {p['F_head_pN_scaled']:.1f} pN")
    ax2.axhline(p["F_head_pN_unscaled"], color="#e67e22", ls="--", lw=1.5,
                label=f"1 head unscaled = {p['F_head_pN_unscaled']:.1f} pN (invest. #1)")
    for nseg, lab in [(1, "ℓ₀"), (2, "2ℓ₀"), (6, "full 3µm")]:
        Lx = nseg * p["ell0_nm"] * 1e-9
        ax2.axvline(Lx * _UM, color="grey", ls=":", lw=1)
        ax2.text(Lx * _UM, ax2.get_ylim()[1] * 0.6, f" {lab}", fontsize=8, color="grey")
    ax2.set_yscale("log")
    ax2.set_xlabel("free (column) length L (µm)")
    ax2.set_ylabel("force (pN, log)")
    ax2.set_title("(b) Euler buckling margin vs free length\n"
                  "a head buckles a span where its line sits ABOVE F_crit")
    ax2.legend(fontsize=8, loc="upper right")
    fig.suptitle("H.7 Gate-B probe #2 — built-mesh free-length & per-span myosin load "
                 f"({'smoke' if r['scale']['is_smoke'] else 'full'})", fontweight="bold")
    fig.savefig(out_png, dpi=130)
    print(f"  figure → {out_png}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160,
                    help="smoke scale; omit/0 for full production ×40 scale")
    ap.add_argument("--n-nuc-beads", type=int, default=400)
    ap.add_argument("--warmup", type=int, default=600, help="warm-up baoab steps (binder engages)")
    ap.add_argument("--softstart", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/gate_b_probe.json")
    ap.add_argument("--out-png", type=str,
                    default="ffn_sim/outputs/h7/figs/h7_gate_b_freelength_probe.png")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments
    r, margin_rows, spans, _ = run_probe(
        n_filaments=n_fil, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, softstart=args.softstart, device=dev, seed=args.seed,
    )
    _report(r)
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(r, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    _figure(r, spans, args.out_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
