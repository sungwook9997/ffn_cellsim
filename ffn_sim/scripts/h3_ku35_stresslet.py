#!/usr/bin/env python
"""KU-3.5-ACTIVE per-minifilament BIPOLAR STRESSLET diagnostic (2026-06-05).

The aggregation ledger (`h3_ku35_aggregation.py`) localised the open KU-3.5-active
floor to GENERATION (the coherent ceiling is itself ~840× below band, while the
*aggregation* of the bonds that exist is ~94% efficient). The literature
(`docs/ACTIN_ARCHITECTURE_NOTES.md` §16 + synthesis: Lenz2012, Murrell2012,
Stam2017, Ennomani2016, Miyazaki2015) says disordered/cortex actomyosin only
generates contractile stress when each myosin minifilament acts as a *bipolar
contractile stresslet* — a force-dipole that pulls two ANTIPARALLEL actin
filaments toward each other (symmetry-broken; not a single head dragging one
filament with its reaction absorbed by the backbone).

This module treats each minifilament as a POTENTIAL bipolar stresslet and
measures, per minifilament:
  * which actin filaments its + side and − side heads bind,
  * the polarity relation of those filaments to the rod axis û
    (ANTIPARALLEL → contractile-intended / PARALLEL → extensile / SINGLE-sided
    → no dipole / INCOHERENT),
  * the signed scalar force-DIPOLE P along û (P < 0 = contractile, P > 0 =
    extensile), P = Σ_bound (f_on_bead·û)·((r_bead − rod_center)·û),
  * per-head F/F_stall (bond-frame; Hill validity),
  * engagement depth (grip s).

The headline question — *does local head-walking become a coherent cortex-scale
bipolar contractile stress?* — is answered by two aggregate numbers:
  * **frac_complete_pairs** — fraction of engaged minifilaments with BOTH a bound
    + head and a bound − head on DIFFERENT filaments (a real dipole, not a
    single-sided drag),
  * **coherence = ΣP / Σ|P|** over complete pairs — −1 = perfectly coherent
    contractile, 0 = washed out (random signs), +1 = extensile.

If frac_complete is low and/or coherence ≈ 0, the per-head forces never organise
into a net contractile stresslet → GENERATION-limited *at the stresslet level*.

MEASUREMENT-ONLY — changes no physics. Fixture-validated (`--self-test` +
`tests/test_ku35_dipole_gate.py`): known geometry → known polarity/dipole sign.

Run:  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_stresslet --self-test
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]

# polarity-relation labels
POL_ANTIPARALLEL = "antiparallel"   # + and − sides bind opposite-polarity filaments → contractile dipole intended
POL_PARALLEL = "parallel"           # same-sign polarity → extensile/degenerate
POL_SINGLE = "single"               # only one side bound → no dipole (drag, reaction into backbone)
POL_INCOHERENT = "incoherent"       # bound both sides but polarity signs don't form a clean dipole
POL_FREE = "free"                   # no heads bound


def classify_polarity(dot_plus: float | None, dot_minus: float | None) -> str:
    """Classify a minifilament's bipolar polarity from the per-side m̂·û means.

    Args:
        dot_plus: mean of (minus-end-ward tangent m̂ · rod axis û) over BOUND
            + side heads, or None if no + head bound.
        dot_minus: same over BOUND − side heads, or None.

    Returns:
        One of POL_ANTIPARALLEL / POL_PARALLEL / POL_SINGLE / POL_FREE. A proper
        bipolar contractile dipole (Stam-Hocky) has the + side on m̂·û > 0 and
        the − side on m̂·û < 0 (opposite signs = ANTIPARALLEL).
    """
    if dot_plus is None and dot_minus is None:
        return POL_FREE
    if dot_plus is None or dot_minus is None:
        return POL_SINGLE
    if dot_plus > 0.0 > dot_minus:
        return POL_ANTIPARALLEL
    if dot_plus < 0.0 < dot_minus:
        return POL_ANTIPARALLEL          # mirror (û sign convention) still antiparallel
    return POL_PARALLEL                  # same sign on both sides → not contractile-bipolar


# --------------------------------------------------------------------------- #
# Core per-minifilament stresslet (operates on extracted arrays / action)
# --------------------------------------------------------------------------- #
def minifilament_stresslets(
    *,
    pos_byTag: np.ndarray,
    action,
    p_myo,
    beads_per_filament: int,
    F_stall: float | None = None,
) -> list[dict]:
    """Per-minifilament bipolar-stresslet records from a live MyosinStepUpdater.

    Args:
        pos_byTag: ``(N, 3)`` tag-ordered positions.
        action: ``MyosinStepUpdater`` with grip-walk bound state
            (``_head_bound_to_actin``, ``_head_bound_filament``,
            ``_head_bound_bead_pos``, ``_head_grip_s``, ``layout``).
        p_myo: ``ResolvedCortexMyosin`` (k_head_actin, F_stall_per_head, sizes).
        beads_per_filament: cortex actin beads per filament (for m̂ tangent).
        F_stall: stall force [N]; defaults to ``p_myo.F_stall_per_head`` (scaled).

    Returns:
        One dict per ENGAGED minifilament (≥1 bound head), with keys: motor_idx,
        n_plus, n_minus, complete (bool), plus_fils, minus_fils, polarity,
        dipole_P (signed N·m, <0 contractile), mean_F_over_Fstall, sum_F_N.
    """
    bound = np.asarray(action._head_bound_to_actin)
    fil_of = np.asarray(action._head_bound_filament)
    grip = np.asarray(action._head_grip_s)
    H = int(p_myo.n_heads_per_side)
    N = int(p_myo.n_backbone)
    per_motor = int(p_myo.n_particles_per_motor)
    mt0 = int(action.layout.motor_tag_start)
    k_ha = float(p_myo.k_head_actin)
    k_hs = float(p_myo.k_head_spring)
    k_series = 1.0 / (1.0 / k_ha + 1.0 / k_hs)   # Hill-stall series stiffness
    r0_eps = 0.0  # grip_walk attach r0 ≈ 0; F = k·r
    F_stall = float(F_stall if F_stall is not None else p_myo.F_stall_per_head)
    nb = int(beads_per_filament)

    def bead_tag(fil, pos):
        return fil * nb + pos

    n_motors = int(action.layout.axes.shape[0])
    out = []
    for m in range(n_motors):
        base = m * 2 * H
        heads = np.arange(base, base + 2 * H)
        bnd = heads[bound[heads] >= 0]
        if bnd.size == 0:
            continue
        # current rod center + axis from the backbone end beads (matches act()).
        bb0 = mt0 + m * per_motor + 0
        bbN = mt0 + m * per_motor + (N - 1)
        bb_tags = np.arange(mt0 + m * per_motor, mt0 + m * per_motor + N)
        center = pos_byTag[bb_tags].mean(axis=0)
        u_vec = pos_byTag[bbN] - pos_byTag[bb0]
        un = float(np.linalg.norm(u_vec))
        u_hat = u_vec / un if un > 0 else np.array([1.0, 0.0, 0.0])

        plus_dots, minus_dots = [], []
        plus_fils, minus_fils = set(), set()
        P = 0.0
        sumF = 0.0
        ratios = []          # bond-frame  F_bond/F_stall = k_ha·r / F_stall
        ratios_series = []   # Hill-stall  F_series/F_stall = k_series·min(s,r) / F_stall
        for h in bnd:
            side_plus = (int(h) - base) < H
            bt = int(bound[h])
            fil = int(fil_of[h]) if fil_of[h] >= 0 else bt // nb
            pos_j = bt - fil * nb
            head_tag = action._head_global_tag(int(h))
            r_bead = pos_byTag[bt]
            r_head = pos_byTag[head_tag]
            d = r_head - r_bead
            r = float(np.linalg.norm(d))
            F = k_ha * max(r - r0_eps, 0.0)
            s_grip = float(grip[int(h)])
            F_ser = k_series * max(min(s_grip, r), 0.0)   # Hill-stall load
            sumF += F
            ratios.append(F / F_stall if F_stall > 0 else np.nan)
            ratios_series.append(F_ser / F_stall if F_stall > 0 else np.nan)
            # force ON the bead = toward the head (attach bond pulls them together)
            f_on_bead = (F * (d / r)) if r > 0 else np.zeros(3)
            x = float((r_bead - center) @ u_hat)         # bead offset along rod
            P += float(f_on_bead @ u_hat) * x            # dipole moment along û
            # minus-end-ward tangent m̂ of the bound filament (Option A: minus=0)
            if pos_j > 0:
                m_vec = pos_byTag[bead_tag(fil, pos_j - 1)] - r_bead
            elif nb > 1:
                m_vec = r_bead - pos_byTag[bead_tag(fil, 1)]
            else:
                m_vec = None
            dot = None
            if m_vec is not None:
                mn = float(np.linalg.norm(m_vec))
                if mn > 0:
                    dot = float((m_vec / mn) @ u_hat)
            if side_plus:
                if dot is not None:
                    plus_dots.append(dot)
                plus_fils.add(fil)
            else:
                if dot is not None:
                    minus_dots.append(dot)
                minus_fils.add(fil)

        dp = float(np.mean(plus_dots)) if plus_dots else None
        dm = float(np.mean(minus_dots)) if minus_dots else None
        complete = bool(plus_fils) and bool(minus_fils) and len(plus_fils & minus_fils) == 0
        polarity = classify_polarity(dp, dm)
        if not complete and polarity == POL_ANTIPARALLEL:
            polarity = POL_SINGLE if not (plus_fils and minus_fils) else POL_INCOHERENT
        out.append(dict(
            motor_idx=m,
            n_plus=int(sum(1 for h in bnd if (int(h) - base) < H)),
            n_minus=int(sum(1 for h in bnd if (int(h) - base) >= H)),
            complete=complete,
            plus_fils=sorted(plus_fils),
            minus_fils=sorted(minus_fils),
            polarity=polarity,
            dipole_P=float(P),                 # < 0 contractile, > 0 extensile
            mean_F_over_Fstall=float(np.nanmean(ratios)) if ratios else float("nan"),
            mean_Fseries_over_Fstall=float(np.nanmean(ratios_series)) if ratios_series else float("nan"),
            max_Fbond_over_Fstall=float(np.nanmax(ratios)) if ratios else float("nan"),
            frac_superstall_bond=float(np.mean(np.asarray(ratios) > 1.0)) if ratios else float("nan"),
            sum_F_N=float(sumF),
        ))
    return out


def stresslet_summary(stresslets: list[dict]) -> dict:
    """Aggregate the per-minifilament records into the headline coherence metrics."""
    n_eng = len(stresslets)
    if n_eng == 0:
        return dict(n_engaged_motors=0, frac_complete_pairs=float("nan"),
                    frac_single_sided=float("nan"), n_antiparallel=0,
                    n_parallel=0, sumP=0.0, sum_absP=0.0, coherence=float("nan"),
                    frac_contractile=float("nan"), verdict="NO-ENGAGED-MOTORS")
    complete = [s for s in stresslets if s["complete"]]
    single = [s for s in stresslets if s["polarity"] == POL_SINGLE]
    P = np.array([s["dipole_P"] for s in complete]) if complete else np.array([])
    sumP = float(P.sum()) if P.size else 0.0
    sum_absP = float(np.abs(P).sum()) if P.size else 0.0
    coherence = (sumP / sum_absP) if sum_absP > 0 else float("nan")  # <0 contractile-coherent
    frac_contractile = float(np.mean(P < 0)) if P.size else float("nan")
    frac_complete = len(complete) / n_eng
    frac_single = len(single) / n_eng
    n_anti = sum(1 for s in stresslets if s["polarity"] == POL_ANTIPARALLEL)
    n_par = sum(1 for s in stresslets if s["polarity"] == POL_PARALLEL)

    if frac_complete < 0.15:
        verdict = (f"STRESSLET-INCOMPLETE: only {frac_complete*100:.0f}% of engaged "
                   "minifilaments form a bipolar pair — most are single-sided drags "
                   "(reaction absorbed by backbone, no contractile dipole)")
    elif not np.isfinite(coherence) or coherence > -0.2:
        verdict = (f"STRESSLET-INCOHERENT: complete pairs exist but coherence="
                   f"{coherence:.2f} (≈0 = contractile/extensile cancel) → no NET "
                   "contractile stress")
    else:
        verdict = (f"STRESSLET-CONTRACTILE: {frac_complete*100:.0f}% complete pairs, "
                   f"coherence={coherence:.2f} (<0 = net contractile) — dipoles DO "
                   "organise into contractile stress")
    return dict(
        n_engaged_motors=n_eng,
        frac_complete_pairs=frac_complete,
        frac_single_sided=frac_single,
        n_antiparallel=n_anti, n_parallel=n_par,
        sumP=sumP, sum_absP=sum_absP, coherence=coherence,
        frac_contractile=frac_contractile,
        verdict=verdict,
    )


# --------------------------------------------------------------------------- #
# Live-sim adapter
# --------------------------------------------------------------------------- #
def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position, dtype=np.float64)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def stresslet_ledger(sim, *, p_myo, myosin_action, beads_per_filament: int) -> dict:
    """Compute the bipolar-stresslet ledger for a live HOOMD snapshot.

    Returns ``{"summary": ..., "stresslets": [...]}``; an empty/no-motor sim
    yields the NO-ENGAGED-MOTORS summary.
    """
    if myosin_action is None or getattr(myosin_action, "_head_bound_to_actin", None) is None:
        return dict(summary=stresslet_summary([]), stresslets=[])
    pos = _tagpos(sim)
    sl = minifilament_stresslets(
        pos_byTag=pos, action=myosin_action, p_myo=p_myo,
        beads_per_filament=beads_per_filament)
    return dict(summary=stresslet_summary(sl), stresslets=sl)


def format_stresslet(ledger: dict, *, indent: str = "    ") -> str:
    s = ledger["summary"]
    if s["n_engaged_motors"] == 0:
        return f"{indent}[STRESSLET] no engaged motors"
    return (
        f"{indent}[STRESSLET] engaged_motors={s['n_engaged_motors']} "
        f"complete_pairs={s['frac_complete_pairs']*100:.0f}% "
        f"single_sided={s['frac_single_sided']*100:.0f}% "
        f"(anti={s['n_antiparallel']} par={s['n_parallel']})\n"
        f"{indent}      ΣP={s['sumP']:.3e} Σ|P|={s['sum_absP']:.3e} "
        f"coherence={s['coherence']:.3f} (<0=contractile) "
        f"frac_contractile={s['frac_contractile']*100:.0f}%\n"
        f"{indent}      → {s['verdict']}"
    )


# --------------------------------------------------------------------------- #
# Self-test (known geometry → known polarity / dipole sign), HOOMD-free
# --------------------------------------------------------------------------- #
class _FakeLayout:
    def __init__(self, axes, motor_tag_start):
        self.axes = np.asarray(axes, dtype=np.float64)
        self.motor_tag_start = int(motor_tag_start)


class _FakeAction:
    """Minimal stand-in exposing the fields minifilament_stresslets reads."""
    def __init__(self, *, bound, fil, grip, axes, motor_tag_start, H, N, per_motor):
        self._head_bound_to_actin = np.asarray(bound, dtype=np.int64)
        self._head_bound_filament = np.asarray(fil, dtype=np.int64)
        self._head_grip_s = np.asarray(grip, dtype=np.float64)
        self.layout = _FakeLayout(axes, motor_tag_start)
        self._H, self._N, self._per = H, N, per_motor

    def _head_global_tag(self, h):
        # mirror MyosinStepUpdater._head_global_tag
        H, N, per = self._H, self._N, self._per
        motor_idx = h // (2 * H)
        within = h % (2 * H)
        return self.layout.motor_tag_start + motor_idx * per + N + within


class _P:
    pass


def self_test(verbose: bool = True) -> bool:
    """One minifilament, two ANTIPARALLEL filaments, both heads bound → the
    diagnostic must report complete + antiparallel + CONTRACTILE (P<0)."""
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    # geometry: rod along +x centered at origin; H=1, N=3, per_motor=5.
    H, N = 1, 3
    per_motor = N + 2 * H
    nb = 5            # beads per filament
    ell0 = 500e-9
    k_ha = 1e-6
    F_stall = 5e-13
    p = _P(); p.n_heads_per_side = H; p.n_backbone = N
    p.n_particles_per_motor = per_motor; p.k_head_actin = k_ha
    p.k_head_spring = 1e-6
    p.F_stall_per_head = F_stall

    # actin: filament 0 along +x (minus=bead0 at x=0 → plus end +x), at y=+ell0;
    #        filament 1 ANTIPARALLEL: minus end at +x (so m̂ points +x... we build
    #        bead positions so its minus-ward tangent is opposite filament 0).
    # tags: actin beads 0..2*nb-1, then motor particles.
    motor_tag_start = 2 * nb
    pos = np.zeros((2 * nb + per_motor, 3))
    # filament 0: beads along +x, y=+ell0 ; minus-ward tangent (bead j-1 - bead j) = -x
    for j in range(nb):
        pos[j] = [j * ell0, +ell0, 0.0]
    # filament 1: beads along +x BUT polarity reversed — place so minus end is at
    # high x. We just reverse the x-order: bead0 at +4ell0 ... bead4 at 0, y=-ell0.
    for j in range(nb):
        pos[nb + j] = [(nb - 1 - j) * ell0, -ell0, 0.0]
    # minifilament backbone along +x centered at x=2ell0, y=0
    bb = motor_tag_start
    for i in range(N):
        pos[bb + i] = [2 * ell0 + (i - (N - 1) / 2) * ell0, 0.0, 0.0]
    # + head bound to filament 0 bead 3 (x=3ell0, +x of center) ; sits at bead, pulled
    plus_head = bb + N
    minus_head = bb + N + 1
    # heads sit INBOARD of their bound beads (rod center is x=2ell0), so the
    # attach force pulls each bead toward the centre along û = a CONTRACTILE
    # dipole (force has a real axial component, not perpendicular).
    pos[plus_head] = [2 * ell0, +ell0 * 0.5, 0.0]   # inboard of filament-0 bead3 (x=3ell0)
    pos[minus_head] = [2 * ell0, -ell0 * 0.5, 0.0]  # inboard of filament-1 bead (x=1ell0)

    # bound state: + head -> filament 0 (tag = 0*nb + 3 = 3), bead_pos 3
    #              − head -> filament 1 (tag = nb + (nb-1 - 1)=nb+3 -> x=1ell0)
    bound = np.full(2 * H, -1, dtype=np.int64)
    fil = np.full(2 * H, -1, dtype=np.int64)
    grip = np.zeros(2 * H)
    # filament 0 bead at x=3ell0 is index 3 -> global tag 3
    bound[0] = 3; fil[0] = 0
    # filament 1: x=1ell0 is at j where (nb-1-j)=1 -> j=3 -> global tag nb+3=8; its
    # filament index = 8 // nb = 1, pos within = 8 - 1*nb = 3
    bound[1] = nb + 3; fil[1] = 1

    act = _FakeAction(bound=bound, fil=fil, grip=grip, axes=[[1.0, 0, 0]],
                      motor_tag_start=motor_tag_start, H=H, N=N, per_motor=per_motor)
    sl = minifilament_stresslets(pos_byTag=pos, action=act, p_myo=p,
                                 beads_per_filament=nb, F_stall=F_stall)
    if verbose:
        print("\n[1] one minifilament, two antiparallel filaments, both heads bound:")
        for s in sl:
            print("   ", {k: s[k] for k in ("complete", "polarity", "dipole_P", "n_plus", "n_minus")})
    check("one engaged minifilament", len(sl) == 1)
    s0 = sl[0]
    check("complete bipolar pair", s0["complete"] is True)
    check("polarity = antiparallel", s0["polarity"] == POL_ANTIPARALLEL)
    check("dipole P < 0 (contractile)", s0["dipole_P"] < 0)

    summ = stresslet_summary(sl)
    check("summary frac_complete = 1.0", summ["frac_complete_pairs"] == 1.0)
    check("summary coherence < 0 (contractile)", summ["coherence"] < 0)
    check("verdict CONTRACTILE", "CONTRACTILE" in summ["verdict"])

    # ---- single-sided case: only + head bound → SINGLE, not a dipole ----
    bound2 = np.array([3, -1], dtype=np.int64); fil2 = np.array([0, -1], dtype=np.int64)
    act2 = _FakeAction(bound=bound2, fil=fil2, grip=grip, axes=[[1.0, 0, 0]],
                       motor_tag_start=motor_tag_start, H=H, N=N, per_motor=per_motor)
    sl2 = minifilament_stresslets(pos_byTag=pos, action=act2, p_myo=p,
                                  beads_per_filament=nb, F_stall=F_stall)
    if verbose:
        print("\n[2] single-sided (only + head bound):")
        print("   ", {k: sl2[0][k] for k in ("complete", "polarity")})
    check("single-sided → not complete", sl2[0]["complete"] is False)
    check("single-sided → polarity SINGLE", sl2[0]["polarity"] == POL_SINGLE)
    summ2 = stresslet_summary(sl2)
    check("single-sided summary verdict INCOMPLETE", "INCOMPLETE" in summ2["verdict"])

    if verbose:
        print(f"\n{'='*56}\nSTRESSLET SELF-TEST {'PASSED' if ok else 'FAILED'}\n{'='*56}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return 0 if self_test() else 1
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
