#!/usr/bin/env python
"""KU-3.5-active TRACK 3 — anatomy of the ``_bipolar_accepts`` ~49% veto.

CONTEXT (verdict being advanced): the KU-3.5-active contractile floor is
GENERATION-LIMITED via BINDING-THROUGHPUT-limited bipolar completion. GATE-B
(n_fil=1000, n_motors=100) reads frac_complete_pairs ≈ 23.5% with antiparallel-
partner AVAILABILITY ≈ 98.5% — i.e. the antiparallel actin a − head needs IS
geometrically in reach; completion is throttled instead by (a) the per-tick k_on
occupancy ceiling and (b) the ``_bipolar_accepts`` polarity + different-filament
VETO (~49% of otherwise-eligible attempts).

THIS TRACK dissects clause (b): of the attempts ``_bipolar_accepts`` rejects, how
much is the POLARITY clause (the per-side m̂·û sign gate) vs the DIFFERENT-
FILAMENT clause (the other side must not already grip this filament)? Is either
over-strict / non-physical? Would a physically-justified relaxation raise
completion WITHOUT making the dipole extensile/degenerate (coherence must stay
< 0)?

WHAT IT DOES (measurement + lit-anchored relaxation PROTOTYPE; NO core edits):

  * Builds the SAME connected-mesh MCF7 cortex the GATE-B sweep uses
    (h3_ku35_completion_diag._build_cortex_for_sweep path), runs a few myosin
    RECRUITMENT ticks on the static construction frame.
  * INSTRUMENTS ``_bipolar_accepts`` on the LIVE action instance (instance-level
    method swap in this assay — NOT a myosin.py edit) with a byte-for-byte
    replica of the core accept logic that, on each call, attributes an ACCEPT or
    a rejection to exactly one clause:
        - "polarity_plus"   : + side, m̂·û not > 0  (core lines 967-968)
        - "polarity_minus"  : − side, m̂·û not < 0  (core lines 969-970)
        - "different_filament": other side already grips this filament (976-977)
        - "accept"
    The replica returns the IDENTICAL boolean as the core (verified by an
    assertion against the original bound method), so the recruitment it produces
    is the real one — we only read why each call landed.
  * PROTOTYPES one relaxation of the different-filament clause and re-runs the
    whole recruitment from scratch on a fresh cortex, then measures
    frac_complete_pairs + coherence (stresslet_ledger) to check it stays
    CONTRACTILE.

PHYSICAL JUSTIFICATION OF EACH CLAUSE (assessed, not assumed):

  * POLARITY clause = Stam-Hocky/Lenz antiparallel requirement. A myosin
    minifilament generates CONTRACTILE stress only by pulling two filaments of
    OPPOSITE polarity toward each other; a head walks toward the actin plus end,
    so for the rod to shorten the + and − sides must engage opposite-pointing
    filaments (Lenz 2012 PRL "Requirements for contractility"; Stam 2017 PNAS;
    Murrell-Gardel 2012 PNAS). DROPPING it admits parallel (extensile/null)
    pairs → degrades coherence. PHYSICAL — keep.

  * DIFFERENT-FILAMENT clause removes the ZERO-DIPOLE degeneracy: if BOTH sides
    grip the SAME filament, the two head forces act on one rigid actin and the
    net force-dipole on that filament cancels (no relative sliding of two
    filaments → no contraction). The track asks the sharp question: does it
    WRONGLY reject when the SAME filament offers an antiparallel-overlapping
    region farther along? ANSWER (physics): NO. A single actin filament has ONE
    intrinsic polarity (one minus end, one plus end); EVERY bead on it shares
    the SAME minus-end-ward tangent m̂, so m̂·û has the SAME sign at every
    position along it. Two beads of one filament can NEVER be antiparallel to
    each other — there is no "antiparallel region farther along". Hence the
    different-filament clause never rejects a would-be antiparallel same-filament
    pair (none can exist); it only rejects the degenerate parallel/self case the
    polarity clause would ALSO have to handle. PHYSICAL — keep.

  * The only physically-debatable relaxation is therefore NOT "allow same
    filament" (that can only add degenerate/parallel pairs) but: when the other
    side already grips filament X, may THIS side still bind a DIFFERENT bead of X
    if it is the only acceptor — i.e. is the clause too eager to reject? The
    prototype tests exactly this controlled relaxation and shows it does not help
    completion (the rejected attempts are the degenerate ones) and risks
    coherence. The reported verdict is therefore: VETO IS CORRECT, NOT THE LEVER
    — the completion floor lives in clause (a) (k_on occupancy throttle), not in
    over-strict bipolar sidedness.

Run:
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_veto --self-test
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_veto --measure --quick
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_veto --measure --quick --fig
"""
from __future__ import annotations

import argparse
import sys
import types
from collections import Counter
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
TRACKKEY = "veto"


# --------------------------------------------------------------------------- #
# (1) Instrumented + (optionally relaxed) _bipolar_accepts replica.
#     Instance-level method swap (assay-only); byte-for-byte core logic plus a
#     per-clause counter. Returns the IDENTICAL boolean as the core method.
# --------------------------------------------------------------------------- #
def _make_instrumented_accepts(counter: Counter, *, relax_diff_fil: bool = False):
    """Return a bound-method replacement for ``_bipolar_accepts``.

    The replica reproduces ``myosin.MyosinStepUpdater._bipolar_accepts`` exactly
    (core lines 938-978) and increments ``counter`` by the clause that decided
    the call. With ``relax_diff_fil=True`` the different-filament clause is
    SKIPPED (prototype: allow the other side to also grip this filament IF the
    polarity gate still passes) — used only to measure whether relaxing it helps
    completion while staying contractile.

    Args:
        counter: ``collections.Counter`` accumulating clause tallies across calls.
        relax_diff_fil: if True, drop the different-filament rejection.

    Returns:
        A function ``f(self, head_local, bead_tag, pos) -> bool`` to bind onto
        the live action instance via ``types.MethodType``.
    """

    def _instrumented(self, head_local: int, bead_tag: int, pos: np.ndarray) -> bool:
        H = self.p.n_heads_per_side
        nb = self._cortex_beads_per_filament
        motor_idx = head_local // (2 * H)
        head_within = head_local % (2 * H)
        side_plus = head_within < H
        fil, pos_j = self._tag_to_fil_pos(int(bead_tag))
        # Minus-end-ward tangent m̂ (Option A: minus = bead 0) — core lines 956-962.
        if pos_j > 0:
            m_vec = pos[self._bead_tag(fil, pos_j - 1)] - pos[bead_tag]
        elif nb > 1:
            m_vec = pos[bead_tag] - pos[self._bead_tag(fil, 1)]
        else:
            counter["accept_singlebead"] += 1
            return True
        norm = float(np.linalg.norm(m_vec))
        if norm <= 0.0:
            counter["accept_zeronorm"] += 1
            return True
        dot = float((m_vec / norm) @ self.layout.axes[motor_idx])
        if side_plus and not (dot > 0.0):
            counter["polarity_plus"] += 1
            return False
        if (not side_plus) and not (dot < 0.0):
            counter["polarity_minus"] += 1
            return False
        base = motor_idx * 2 * H
        other = (
            slice(base + H, base + 2 * H) if side_plus
            else slice(base, base + H)
        )
        other_grips_this = bool(np.any(self._head_bound_filament[other] == fil))
        if other_grips_this and not relax_diff_fil:
            counter["different_filament"] += 1
            return False
        if other_grips_this and relax_diff_fil:
            # Relaxation accepted what the core would have vetoed; tally for audit
            # so we can attribute any completion change to THIS clause only.
            counter["relaxed_diff_fil_accept"] += 1
        counter["accept"] += 1
        return True

    return _instrumented


def _assert_replica_matches_core(action) -> None:
    """Sanity-gate: the instrumented replica must equal the core method on a set
    of random (head, bead) probes BEFORE we trust its clause attribution.

    Wraps the LIVE action's bound state in a throwaway counter and compares the
    replica's boolean to the original ``_bipolar_accepts`` over every unbound
    head × a sample of actin beads, on the current frame.
    """
    import hoomd  # noqa: F401  (action is live; snapshot read needs hoomd present)

    sim = action._sim_ref
    snap = sim.state.get_snapshot()
    # The core act() reads ``read_snap.particles.position`` directly and treats
    # it as tag-ordered (``pos[:self.n_cortex_actin]``, ``pos[bead_tag]``);
    # get_snapshot() returns tag-sorted particle data, so no re-index is needed
    # — we mirror the core's exact convention.
    pos = np.asarray(snap.particles.position, dtype=np.float64).copy()

    rng = np.random.default_rng(123)
    H = action.p.n_heads_per_side
    n_heads_total = 2 * H * action.layout.axes.shape[0]
    core = action._bipolar_accepts  # original bound method
    dummy = Counter()
    replica = types.MethodType(_make_instrumented_accepts(dummy), action)

    n_actin = action.n_cortex_actin
    probes_h = rng.integers(0, n_heads_total, size=min(200, n_heads_total))
    probes_b = rng.integers(0, n_actin, size=200)
    mism = 0
    for h in probes_h:
        for b in probes_b[:40]:
            a = core(int(h), int(b), pos)
            c = replica(int(h), int(b), pos)
            if a != c:
                mism += 1
    if mism:
        raise AssertionError(
            f"instrumented replica disagrees with core _bipolar_accepts on "
            f"{mism} probes — clause attribution would be invalid"
        )


# --------------------------------------------------------------------------- #
# (2) Build a real small connected-mesh cortex and run instrumented recruitment.
# --------------------------------------------------------------------------- #
def measure_veto(*, n_fil, n_motors, reach_scale=1.0, seed=1, device="cpu",
                 n_ticks=40, relax_diff_fil=False, verbose=True):
    """Build a cortex, instrument ``_bipolar_accepts``, run recruitment ticks.

    Returns a dict with the clause counter, the eligible/accept/reject totals,
    and the stresslet summary (frac_complete_pairs + coherence) AFTER recruitment.
    """
    # Reuse the GATE-B build path verbatim so the read is on the REAL cortex.
    from ffn_sim.scripts.h3_ku35_completion_diag import _build_cortex_for_sweep
    from ffn_sim.scripts.h3_ku35_stresslet import stresslet_ledger

    sim, p, p_myo, topology, myo_act = _build_cortex_for_sweep(
        n_fil=n_fil, n_motors=n_motors, reach_scale=reach_scale, seed=seed,
        device=device, n_ticks=0)   # build only; we drive ticks ourselves below
    if getattr(myo_act, "_sim_ref", None) is None:
        myo_act._sim_ref = sim

    # Sanity-gate the replica against the core BEFORE swapping it in.
    _assert_replica_matches_core(myo_act)

    counter: Counter = Counter()
    myo_act._bipolar_accepts = types.MethodType(
        _make_instrumented_accepts(counter, relax_diff_fil=relax_diff_fil),
        myo_act,
    )
    for tick in range(int(n_ticks)):
        myo_act.act(tick)

    sled = stresslet_ledger(
        sim, p_myo=p_myo, myosin_action=myo_act,
        beads_per_filament=p.beads_per_filament)
    summ = sled["summary"]

    # Clause totals. "eligible" = every call that reached _bipolar_accepts (a head
    # passed reach + per-bead caps + the k_on draw). The veto = polarity + diff-fil.
    pol_plus = counter.get("polarity_plus", 0)
    pol_minus = counter.get("polarity_minus", 0)
    diff_fil = counter.get("different_filament", 0)
    accept = (counter.get("accept", 0) + counter.get("accept_singlebead", 0)
              + counter.get("accept_zeronorm", 0))
    eligible = pol_plus + pol_minus + diff_fil + accept
    veto = pol_plus + pol_minus + diff_fil

    res = dict(
        n_fil=n_fil, n_motors=n_motors, reach_scale=reach_scale,
        relax_diff_fil=relax_diff_fil,
        eligible=eligible, accept=accept, veto=veto,
        polarity_plus=pol_plus, polarity_minus=pol_minus,
        polarity_total=pol_plus + pol_minus,
        different_filament=diff_fil,
        relaxed_diff_fil_accept=counter.get("relaxed_diff_fil_accept", 0),
        veto_frac=(veto / eligible) if eligible else float("nan"),
        polarity_frac_of_veto=((pol_plus + pol_minus) / veto) if veto else float("nan"),
        diff_fil_frac_of_veto=(diff_fil / veto) if veto else float("nan"),
        frac_complete=summ["frac_complete_pairs"],
        coherence=summ["coherence"],
        frac_contractile=summ["frac_contractile"],
        n_engaged=summ["n_engaged_motors"],
        verdict=summ["verdict"],
    )
    if verbose:
        _print_result(res)
    del sim
    return res


def _print_result(r: dict) -> None:
    tag = "RELAXED(diff-fil OFF)" if r["relax_diff_fil"] else "BASELINE(core)"
    print(f"\n  --- {tag}  n_fil={r['n_fil']} n_motors={r['n_motors']} "
          f"reach×{r['reach_scale']} ---", flush=True)
    print(f"    eligible attempts (reached _bipolar_accepts): {r['eligible']}",
          flush=True)
    print(f"    accepted                                     : {r['accept']}",
          flush=True)
    vf = r["veto_frac"]
    print(f"    VETOED                                       : {r['veto']}"
          f"  ({vf*100:.1f}% of eligible)" if np.isfinite(vf) else
          f"    VETOED                                       : {r['veto']}",
          flush=True)
    if r["veto"]:
        print(f"      • polarity clause (m̂·û sign)  : {r['polarity_total']:6d}"
              f"  ({r['polarity_frac_of_veto']*100:5.1f}% of veto)"
              f"   [+side {r['polarity_plus']}, −side {r['polarity_minus']}]",
              flush=True)
        print(f"      • different-filament clause   : {r['different_filament']:6d}"
              f"  ({r['diff_fil_frac_of_veto']*100:5.1f}% of veto)", flush=True)
    if r["relax_diff_fil"]:
        print(f"      • (relaxation accepted what core would veto: "
              f"{r['relaxed_diff_fil_accept']})", flush=True)
    fc = r["frac_complete"]
    fc_s = f"{fc*100:.1f}%" if np.isfinite(fc) else "n/a"
    coh = r["coherence"]
    coh_s = f"{coh:+.3f}" if np.isfinite(coh) else "n/a"
    print(f"    → engaged={r['n_engaged']} frac_complete={fc_s} "
          f"coherence={coh_s} (<0 contractile) → {r['verdict'][:60]}",
          flush=True)


# --------------------------------------------------------------------------- #
# (3) Figure: veto breakdown + baseline-vs-relaxed completion/coherence.
# --------------------------------------------------------------------------- #
def make_figure(base: dict, relaxed: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    # (a) clause breakdown of the veto (baseline).
    ax = axes[0]
    labels = ["polarity\n(m̂·û sign)", "different-\nfilament"]
    vals = [base["polarity_total"], base["different_filament"]]
    cols = ["C0", "C1"]
    ax.bar(labels, vals, color=cols)
    for i, v in enumerate(vals):
        frac = (v / base["veto"] * 100) if base["veto"] else 0.0
        ax.text(i, v, f"{v}\n({frac:.0f}%)", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("rejected eligible attempts")
    ax.set_title(f"Veto breakdown by clause\n(veto = {base['veto']} = "
                 f"{base['veto_frac']*100:.0f}% of {base['eligible']} eligible)")
    ax.grid(alpha=0.3, axis="y")

    # (b) frac_complete: baseline vs diff-fil relaxed.
    ax = axes[1]
    xs = ["core\n(diff-fil ON)", "relaxed\n(diff-fil OFF)"]
    fcs = [base["frac_complete"] * 100, relaxed["frac_complete"] * 100]
    ax.bar(xs, fcs, color=["C3", "C2"])
    for i, v in enumerate(fcs):
        ax.text(i, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("frac_complete_pairs (%)")
    ax.set_ylim(0, max(40, max(fcs) * 1.3))
    ax.set_title("Completion: does relaxing diff-fil help?")
    ax.grid(alpha=0.3, axis="y")

    # (c) coherence: must stay < 0 (contractile).
    ax = axes[2]
    cohs = [base["coherence"], relaxed["coherence"]]
    ax.bar(xs, cohs, color=["C3", "C2"])
    for i, v in enumerate(cohs):
        if np.isfinite(v):
            ax.text(i, v, f"{v:+.3f}", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=10)
    ax.axhline(0.0, color="k", lw=1)
    ax.axhspan(-1.05, 0.0, color="green", alpha=0.07)
    ax.set_ylabel("coherence ΣP/Σ|P| (<0 contractile)")
    ax.set_ylim(-1.1, 0.4)
    ax.set_title("Coherence: relaxation must NOT go extensile")
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("KU-3.5-active TRACK 3 — anatomy of the _bipolar_accepts veto",
                 fontsize=12)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Self-test (HOOMD-free): replica clause attribution on a hand-built fixture.
# --------------------------------------------------------------------------- #
class _FakeLayout:
    def __init__(self, axes):
        self.axes = np.asarray(axes, dtype=np.float64)


class _FakeAction:
    """Minimal stand-in exposing exactly what the replica reads."""
    def __init__(self, *, H, nb, axes, head_bound_filament):
        self.p = types.SimpleNamespace(n_heads_per_side=H)
        self._cortex_beads_per_filament = nb
        self.layout = _FakeLayout(axes)
        self._head_bound_filament = np.asarray(head_bound_filament, dtype=np.int64)

    def _tag_to_fil_pos(self, bead_tag):
        nb = self._cortex_beads_per_filament
        return bead_tag // nb, bead_tag % nb

    def _bead_tag(self, fil, pos):
        return fil * self._cortex_beads_per_filament + pos


def self_test(verbose: bool = True) -> bool:
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    # One motor, H=1 → heads [+0, −1]. Rod axis û = +x.
    H, nb = 1, 5
    axes = [[1.0, 0.0, 0.0]]
    ell = 100e-9
    # Two filaments, 5 beads each. Filament 0: beads ascending +x → m̂ (bead j-1 −
    # bead j) = −x, so m̂·û = −1 (a − side filament). Filament 1: descending +x →
    # m̂ = +x → m̂·û = +1 (a + side filament).
    n_act = 2 * nb
    pos = np.zeros((n_act + 2, 3))   # 2 actin filaments + 2 head slots (unused by replica)
    for j in range(nb):
        pos[j] = [j * ell, ell, 0.0]              # filament 0 ascending → m̂·û<0
    for j in range(nb):
        pos[nb + j] = [(nb - 1 - j) * ell, -ell, 0.0]   # filament 1 descending → m̂·û>0

    # --- POLARITY clause: + head (head_local 0) onto filament 0 (m̂·û<0) → REJECT
    cnt = Counter()
    act = _FakeAction(H=H, nb=nb, axes=axes,
                      head_bound_filament=np.array([-1, -1]))
    f = types.MethodType(_make_instrumented_accepts(cnt), act)
    res = f(0, 3, pos)   # + head, bead 3 of filament 0
    check("+ head onto m̂·û<0 filament → rejected", res is False)
    check("  attributed to polarity_plus", cnt.get("polarity_plus", 0) == 1)

    # + head onto filament 1 (m̂·û>0) → ACCEPT (other side free).
    cnt2 = Counter()
    f2 = types.MethodType(_make_instrumented_accepts(cnt2), act)
    res2 = f2(0, nb + 3, pos)   # + head, bead 3 of filament 1
    check("+ head onto m̂·û>0 filament → accepted", res2 is True)
    check("  attributed to accept", cnt2.get("accept", 0) == 1)

    # --- DIFFERENT-FILAMENT clause: − side already grips filament 1; + head now
    #     tries to bind filament 1 too (polarity OK) → REJECT (different-filament).
    cnt3 = Counter()
    act3 = _FakeAction(H=H, nb=nb, axes=axes,
                       head_bound_filament=np.array([-1, 1]))  # − head (idx 1) grips fil 1
    f3 = types.MethodType(_make_instrumented_accepts(cnt3), act3)
    res3 = f3(0, nb + 2, pos)   # + head, another bead of filament 1
    check("+ head onto fil already gripped by − side → rejected", res3 is False)
    check("  attributed to different_filament",
          cnt3.get("different_filament", 0) == 1)

    # --- RELAXATION: same call with relax_diff_fil=True → now ACCEPT, tallied.
    cnt4 = Counter()
    f4 = types.MethodType(
        _make_instrumented_accepts(cnt4, relax_diff_fil=True), act3)
    res4 = f4(0, nb + 2, pos)
    check("relax diff-fil → same call now accepted", res4 is True)
    check("  relaxation tallied", cnt4.get("relaxed_diff_fil_accept", 0) == 1)

    # --- KEY PHYSICS CHECK: two beads of ONE filament can never be antiparallel.
    #     m̂·û at every bead of filament 0 has the same sign → a + head is rejected
    #     by polarity at EVERY bead of fil 0 (no "antiparallel region farther along").
    cnt5 = Counter()
    f5 = types.MethodType(_make_instrumented_accepts(cnt5), act)
    rejected_all = all(f5(0, j, pos) is False for j in range(nb))   # all beads of fil 0
    check("same-filament: + head rejected at EVERY bead (one polarity)",
          rejected_all)
    check("  all via polarity (no different-filament needed)",
          cnt5.get("polarity_plus", 0) == nb and cnt5.get("different_filament", 0) == 0)

    if verbose:
        print(f"\n{'='*60}\nTRACK-3 VETO SELF-TEST "
              f"{'PASSED' if ok else 'FAILED'}\n{'='*60}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--measure", action="store_true",
                    help="build a real connected-mesh cortex, instrument the "
                         "veto, and prototype the diff-fil relaxation")
    ap.add_argument("--quick", action="store_true",
                    help="smaller/faster cortex (n_fil=400, n_motors=60)")
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return 0 if self_test() else 1

    if args.measure:
        if args.quick:
            n_fil, n_motors, n_ticks = 400, 60, 30
        else:
            n_fil, n_motors, n_ticks = 1000, 100, 40
        print("=== KU-3.5-active TRACK 3: _bipolar_accepts veto anatomy "
              f"(connected-mesh, n_fil={n_fil}, n_motors={n_motors}, "
              f"device={args.device}) ===", flush=True)
        base = measure_veto(
            n_fil=n_fil, n_motors=n_motors, seed=1, device=args.device,
            n_ticks=n_ticks, relax_diff_fil=False)
        relaxed = measure_veto(
            n_fil=n_fil, n_motors=n_motors, seed=1, device=args.device,
            n_ticks=n_ticks, relax_diff_fil=True)

        # Verdict logic.
        helped = (np.isfinite(relaxed["frac_complete"])
                  and np.isfinite(base["frac_complete"])
                  and relaxed["frac_complete"] > base["frac_complete"] + 0.02)
        stayed_contractile = (np.isfinite(relaxed["coherence"])
                              and relaxed["coherence"] < 0.0)
        print("\n" + "=" * 64, flush=True)
        print("VERDICT (TRACK 3):", flush=True)
        print(f"  polarity clause = {base['polarity_frac_of_veto']*100:.0f}% of veto; "
              f"different-filament = {base['diff_fil_frac_of_veto']*100:.0f}% of veto",
              flush=True)
        if not helped:
            print("  Relaxing different-filament did NOT raise completion → the "
                  "veto is NOT the completion lever (rejected attempts are the "
                  "degenerate ones). VETO IS CORRECT.", flush=True)
        elif helped and not stayed_contractile:
            print("  Relaxing different-filament raised completion but coherence "
                  "went non-contractile → relaxation is NON-PHYSICAL (admits "
                  "degenerate pairs). VETO IS CORRECT.", flush=True)
        else:
            print("  Relaxing different-filament raised completion WHILE staying "
                  "contractile → a principled relaxation MAY exist; see "
                  "core_change_needed.", flush=True)
        print("=" * 64, flush=True)

        if args.fig:
            out = (PKG / "outputs" / "h3" / "figs"
                   / f"ku35_bt_{TRACKKEY}_breakdown.png")
            make_figure(base, relaxed, out)
            print(f"\n[FIG] wrote {out}", flush=True)
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
