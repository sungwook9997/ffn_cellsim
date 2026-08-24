"""KU-3.5-active binding-throughput SYNTHESIS — one principled lever decision.

Track-unique synthesis assay (read-only; NO core edits). Consolidates the three
binding-throughput investigation tracks into one auditable arithmetic that
decides whether the KU-3.5-active g_soft floor is reachable by any
literature-anchored binding-kinetics lever.

The three tracks established (all verified against the repo + KB before this
synthesis was written):

  TRACK 1 (duty / off-rate):   head_actin_k_off0 = 10/s is ~29x FASTER than the
    NMIIA actin-detachment reference off-rate 0.35/s (Tam 2021 bioRxiv
    2021.02.23.432588 p37, quoting Stam-Hocky 2015 citing Wang 2003 JBC +
    Kovacs 2003 JBC 278:38132: "the reference off-rate koff(0) for non-muscle
    myosin is 0.35 s^-1 (IIA) and 1.71 s^-1 (IIB)... we adopt IIA"). Verified
    verbatim in paper_chunks [Tam2021_Jour p37 line 844-846]. Re-anchoring
    10 -> 0.35 raises equilibrium duty 0.833 -> 0.993 AND relaxes the D2-batch
    CFL clamp (batch_dt * k_off0 <= 1e-3, myosin.py:330-341), lifting per-tick
    p_bind ~28x. Prototyped effect: bound-frac 11.4% -> 38.4%, completion
    41.6% -> 77.8%, coherence -0.82 -> -0.96. g_soft rose only 3.17x and stays
    ~20,800x UNDER band. NEGATIVE at the g_soft level.

  TRACK 2 (run window):  the GATE-B 23.5% completion is a 40-tick (0.17 tau_eq)
    sub-equilibrium artifact; the equilibrated steady-state completion at
    physiological kinetics is ~70% (reached in ~4 tau_eq ~ 70 ms). A
    MEASUREMENT-window fix, not a kinetics change. Removes a confound but does
    not change the generation budget.

  TRACK 3 (bipolar veto):  the ~49% veto is ~100% the mandatory POLARITY clause
    (Lenz-Gardel-Dinner 2012 NJP 14:033037; Murrell 2015 NRMCB antiparallel
    requirement) and ~0% the different-filament clause (dormant; dropping it
    accepts 0 extra binds). NO principled lever here; relaxing polarity would
    wash out coherence (extensile) = a regression.

THE DECIDING ARITHMETIC (this file): g_soft does NOT scale with completion or
duty because the per-head FORCE-GENERATION budget is structurally capped far
below the band. Each bound head is Hill-bounded at F <= F_stall = 0.5 pN. The
ABSOLUTE coherent ceiling = F_stall * n_motors * n_heads_per_side. The band
requires a cut force F_band = gamma_lo * 2*pi*R. Their ratio is the verdict.

This is the per-head x N "GENERATION-LIMITED" decomposition that
scripts/h3_ku35_aggregation.py already classifies (gamma_ceiling < band ->
GENERATION-LIMITED); this file states it in closed form so the synthesis is
self-contained and the magic-number/band-tuning posture is auditable.

Run:
    conda run -n ffn_sim python -m aleph.scripts.h3_ku35_bt_synthesis
    conda run -n ffn_sim python -m aleph.scripts.h3_ku35_bt_synthesis --self-test

NO band-tuning: every number below is either a config literal or a literature
anchor (cited inline). Nothing here is chosen to make the band pass; the point
of the file is to show the band is UNREACHABLE at physiological kinetics and to
quantify by how much, so the floor is reported HONESTLY as real.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, asdict


# ---------------------------------------------------------------------------
# Repo-pinned constants (config literals from aleph/configs/phase1_h3.yaml
# cortex.myosin + KU-3.5 band from scripts/h3_ku35_aggregation.py BAND).
# NONE of these is introduced or tuned here; they are READ from the model.
# ---------------------------------------------------------------------------
BAND_LO_NPM = 0.35e-3          # N/m  KU-3.5 active cortical-tension band lower bound
BAND_HI_NPM = 0.65e-3          # N/m  KU-3.5 active cortical-tension band upper bound
F_STALL_PER_HEAD_N = 0.5e-12   # N    config: F_stall_per_head = 0.5 pN (Kovacs 2003)
N_MOTORS = 100                 # config: n_motors_per_cell (Salbreux 2012 3/um^2)
N_HEADS_PER_SIDE = 10          # config: n_heads_per_side (brief literal)
R_CELL_MCF7_M = 7.5e-6         # m    MCF7 radius (Wagner 2011; GATE-B scale)
R_CELL_GENERIC_M = 10.0e-6     # m    generic R_cell (config: R_cell, KU-3.17)

# --- binding-kinetics anchors (config current vs literature) ---
K_OFF0_MODEL = 10.0            # 1/s  config: head_actin_k_off0 (mislabelled "NMII duty / Veigel 2002")
K_OFF0_LIT_NMIIA = 0.35        # 1/s  NMIIA reference off-rate (Tam 2021 p37; Stam-Hocky 2015; Wang 2003; Kovacs 2003)
K_ON_MODEL = 50.0             # 1/s  config: head_actin_k_on
CFL_EVENTS_CEILING = 1.0e-3    # -    D2-batch CFL: batch_dt * k_off0 <= 1e-3 (myosin.py:330)


@dataclass
class GenerationVerdict:
    """Closed-form GENERATION-LIMITED decomposition for one (R, scenario)."""

    scenario: str
    R_cell_m: float
    # band-required cut force
    F_band_required_N: float
    # absolute coherent ceiling (every head bound, fully stalled, perfectly coherent)
    F_ceiling_absolute_N: float
    ceiling_over_band: float
    # realistic operating point (bound_frac x coherence x typical F/F_stall)
    bound_frac: float
    coherence: float
    F_over_Fstall_typ: float
    F_realistic_N: float
    realistic_over_band: float
    g_soft_realistic_Npm: float
    under_band_factor: float          # band_lo / g_soft_realistic
    is_generation_limited: bool       # ceiling < band  ->  True


def band_required_cut_force(R_cell_m: float, band_lo: float = BAND_LO_NPM) -> float:
    """F_band = gamma_lo * 2*pi*R  (the cut force a great-circle plane must carry).

    Mirrors scripts/h3_ku35_aggregation.py:250 exactly.
    """
    return band_lo * 2.0 * math.pi * R_cell_m


def absolute_coherent_ceiling(
    n_motors: int = N_MOTORS,
    n_heads: int = N_HEADS_PER_SIDE,
    f_stall: float = F_STALL_PER_HEAD_N,
) -> float:
    """Σ|F| if EVERY head is bound, fully stalled (F=F_stall), perfectly coherent.

    This is the physically-unachievable best case — the Hill force-velocity law
    caps each head at F_stall, so no binding-kinetics lever can exceed it.
    """
    return f_stall * n_motors * n_heads


def realistic_generation(
    R_cell_m: float,
    bound_frac: float,
    coherence: float,
    f_over_fstall: float,
    n_motors: int = N_MOTORS,
    n_heads: int = N_HEADS_PER_SIDE,
    f_stall: float = F_STALL_PER_HEAD_N,
) -> float:
    """Σ|F| at a realistic operating point (the measured GATE-B / lit-anchored read)."""
    F_typ = f_over_fstall * f_stall
    return F_typ * n_motors * n_heads * bound_frac * coherence


def verdict_for(
    scenario: str,
    R_cell_m: float,
    bound_frac: float,
    coherence: float,
    f_over_fstall: float,
) -> GenerationVerdict:
    F_band = band_required_cut_force(R_cell_m)
    F_ceiling = absolute_coherent_ceiling()
    F_real = realistic_generation(R_cell_m, bound_frac, coherence, f_over_fstall)
    g_real = F_real / (2.0 * math.pi * R_cell_m)
    return GenerationVerdict(
        scenario=scenario,
        R_cell_m=R_cell_m,
        F_band_required_N=F_band,
        F_ceiling_absolute_N=F_ceiling,
        ceiling_over_band=F_ceiling / F_band,
        bound_frac=bound_frac,
        coherence=coherence,
        F_over_Fstall_typ=f_over_fstall,
        F_realistic_N=F_real,
        realistic_over_band=F_real / F_band,
        g_soft_realistic_Npm=g_real,
        under_band_factor=BAND_LO_NPM / g_real if g_real > 0 else float("inf"),
        is_generation_limited=(F_ceiling < F_band),
    )


def cfl_batch_dt_for(k_off0: float, dt_s: float = 6.98e-7) -> float:
    """Max batch_dt the D2-batch CFL allows for a given k_off0 (myosin.py:333)."""
    return CFL_EVENTS_CEILING / k_off0


def equilibrium_duty(k_on: float, k_off0: float) -> float:
    """Single-head zero-load duty p_bound = k_on/(k_on+k_off0)."""
    return k_on / (k_on + k_off0)


def synthesize() -> dict:
    """Run the full closed-form synthesis and return an audit dict."""
    # The three operating points we care about (all from the verified tracks):
    #   - model:    measured GATE-B read  (bound 11.4%, coh 0.82, F/Fst 0.48)
    #   - lit_max:  TRACK-1 best principled lever (lit off-rate + CFL-max batch_dt;
    #               bound 38.4%, coh 0.96, F/Fst 0.34), TRACK-2 equilibrated window
    #   - ideal:    physically-impossible 100% bound / fully stalled / coherent=1
    verdicts = [
        verdict_for("model_GATE-B (10/s, 40-tick)", R_CELL_MCF7_M,
                    bound_frac=0.114, coherence=0.82, f_over_fstall=0.48),
        verdict_for("lit_off-rate + CFL-max + equilibrated", R_CELL_MCF7_M,
                    bound_frac=0.384, coherence=0.96, f_over_fstall=0.34),
        verdict_for("IDEAL (100% bound, fully stalled, coh=1)", R_CELL_MCF7_M,
                    bound_frac=1.0, coherence=1.0, f_over_fstall=1.0),
        verdict_for("IDEAL @ generic R=10um", R_CELL_GENERIC_M,
                    bound_frac=1.0, coherence=1.0, f_over_fstall=1.0),
    ]

    duty_model = equilibrium_duty(K_ON_MODEL, K_OFF0_MODEL)
    duty_lit = equilibrium_duty(K_ON_MODEL, K_OFF0_LIT_NMIIA)

    return dict(
        anchors=dict(
            band_Npm=(BAND_LO_NPM, BAND_HI_NPM),
            k_off0_model=K_OFF0_MODEL,
            k_off0_lit_NMIIA=K_OFF0_LIT_NMIIA,
            off_rate_overset_factor=K_OFF0_MODEL / K_OFF0_LIT_NMIIA,
            equilibrium_duty_model=duty_model,
            equilibrium_duty_lit=duty_lit,
            cfl_batch_dt_model_s=cfl_batch_dt_for(K_OFF0_MODEL),
            cfl_batch_dt_lit_s=cfl_batch_dt_for(K_OFF0_LIT_NMIIA),
            cfl_relax_factor=cfl_batch_dt_for(K_OFF0_LIT_NMIIA)
            / cfl_batch_dt_for(K_OFF0_MODEL),
        ),
        verdicts=[asdict(v) for v in verdicts],
        decision=dict(
            generation_limited=all(
                v.is_generation_limited for v in verdicts if "IDEAL" in v.scenario
            ),
            ideal_ceiling_over_band=verdicts[2].ceiling_over_band,
            best_principled_lever="lit off-rate (10->0.35/s, Tam 2021) + "
            "CFL-relaxed batch_dt + equilibrated measurement window",
            best_principled_under_band_factor=verdicts[1].under_band_factor,
            floor_real_at_physiological_kinetics=True,
        ),
    )


def _print_report(audit: dict) -> None:
    a = audit["anchors"]
    print("=" * 74)
    print("KU-3.5-active BINDING-THROUGHPUT SYNTHESIS — generation-limited verdict")
    print("=" * 74)
    print(f"  band (KU-3.5 active)        : {a['band_Npm'][0]*1e3:.2f}-"
          f"{a['band_Npm'][1]*1e3:.2f} mN/m")
    print(f"  k_off0 model vs lit (NMIIA) : {a['k_off0_model']} vs "
          f"{a['k_off0_lit_NMIIA']} /s  ({a['off_rate_overset_factor']:.0f}x over-set)")
    print(f"  equilibrium duty model/lit  : {a['equilibrium_duty_model']:.3f} -> "
          f"{a['equilibrium_duty_lit']:.3f}")
    print(f"  CFL batch_dt relax (lever)  : {a['cfl_batch_dt_model_s']:.2e} -> "
          f"{a['cfl_batch_dt_lit_s']:.2e} s  ({a['cfl_relax_factor']:.0f}x headroom)")
    print("-" * 74)
    for v in audit["verdicts"]:
        print(f"  [{v['scenario']}]  R={v['R_cell_m']*1e6:.1f}um")
        print(f"     F_band required        = {v['F_band_required_N']:.3e} N")
        print(f"     absolute ceiling Sum|F|= {v['F_ceiling_absolute_N']:.3e} N "
              f"({v['ceiling_over_band']:.4f}x band)")
        print(f"     realistic Sum|F|       = {v['F_realistic_N']:.3e} N "
              f"({v['realistic_over_band']:.4f}x band)")
        print(f"     g_soft realistic       = {v['g_soft_realistic_Npm']*1e3:.3e} mN/m "
              f"-> {v['under_band_factor']:.0f}x UNDER band_lo")
        print(f"     GENERATION-LIMITED     = {v['is_generation_limited']}")
        print()
    d = audit["decision"]
    print("-" * 74)
    print(f"  DECISION: floor is GENERATION-LIMITED (per-head F x N), not "
          f"binding-throughput.")
    print(f"     IDEAL coherent ceiling is only {d['ideal_ceiling_over_band']:.4f}x "
          f"the band even at 100% bound / fully stalled / coherence=1.")
    print(f"     best principled lever : {d['best_principled_lever']}")
    print(f"     ...still {d['best_principled_under_band_factor']:.0f}x UNDER band.")
    print(f"     floor real at physiological kinetics : "
          f"{d['floor_real_at_physiological_kinetics']}")
    print("=" * 74)


def _self_test() -> None:
    """Sanity gate — invariants that must hold for the synthesis to be trustworthy."""
    audit = synthesize()
    a = audit["anchors"]
    v = {x["scenario"]: x for x in audit["verdicts"]}

    # 1. off-rate over-set factor is ~29x (10 / 0.35).
    assert abs(a["off_rate_overset_factor"] - 10.0 / 0.35) < 1e-9, "off-rate factor"

    # 2. Re-anchoring off-rate RAISES duty (TRACK-1 falsification of "duty too low").
    assert a["equilibrium_duty_lit"] > a["equilibrium_duty_model"], "duty should rise"
    assert a["equilibrium_duty_lit"] > 0.99, "lit duty toward saturation"

    # 3. CFL relax factor == off-rate factor (the CFL clamp is linear in 1/k_off0).
    assert abs(a["cfl_relax_factor"] - a["off_rate_overset_factor"]) < 1e-6, "CFL relax"

    # 4. THE DECIDING INVARIANT: the IDEAL coherent ceiling is BELOW the band.
    ideal = v["IDEAL (100% bound, fully stalled, coh=1)"]
    assert ideal["is_generation_limited"], "ideal must be generation-limited"
    assert ideal["ceiling_over_band"] < 1.0, "ideal ceiling must be < band"
    # quantitatively ~0.03x (5e-10 / 1.65e-8).
    assert 0.01 < ideal["ceiling_over_band"] < 0.1, "ideal ceiling ~0.03x band"

    # 5. The best principled lever stays UNDER band (honest floor).
    lit = v["lit_off-rate + CFL-max + equilibrated"]
    assert lit["realistic_over_band"] < 1.0, "best lever still under band"
    assert lit["under_band_factor"] > 10.0, "best lever far under band"

    # 6. Even relaxing to the larger generic radius keeps it generation-limited
    #    (band-required force grows with R, so a bigger cell is no easier).
    assert v["IDEAL @ generic R=10um"]["is_generation_limited"], "generic R still GL"

    # 7. g_soft does NOT scale 1:1 with bound_frac: model->lit duty rises ~3.4x
    #    but g_soft ratio is bounded by (bound x coh x F/Fst) ratio, NOT the band.
    model = v["model_GATE-B (10/s, 40-tick)"]
    g_ratio = lit["g_soft_realistic_Npm"] / model["g_soft_realistic_Npm"]
    assert g_ratio < 10.0, "g_soft moves <1 decade while completion ~doubles"

    print("SELF-TEST PASS — all 7 synthesis invariants hold.")
    print(f"  IDEAL ceiling/band      = {ideal['ceiling_over_band']:.4f}x "
          f"(generation-limited even at 100% / stalled / coherent)")
    print(f"  best-lever under-band   = {lit['under_band_factor']:.0f}x")
    print(f"  g_soft(lit)/g_soft(model)= {g_ratio:.2f}x  (completion ~doubled)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true",
                    help="run the synthesis sanity gate and exit")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    _print_report(synthesize())


if __name__ == "__main__":
    main()
