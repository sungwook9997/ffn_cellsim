---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# MCF7 resting turgor Π₀ — provenance, units-explicit derivation, recommendation

**Date:** 2026-07-24 · **Scope:** read-only analysis + recommendation for the single most load-bearing
baseline in the model, `assemble.py:128 PI_0_PA = 40.0` (also `fluid/params_i0b1.yaml` Q3
`delta_p_hydrostatic_rest`). **Trigger:** `PARAM_PROVENANCE_AUDIT_2026-07-24.md` row 1 flagged 40 Pa as the
worst CONVENIENCE/proxy default — a Fischer-Friedrich **HeLa** value used for an **MCF7** baseline. **No sim
code was changed and nothing was committed;** this memo is the sourcing analysis a PI evidence-card
ratification needs.

---

## Headline

**Recommended MCF7 resting turgor: Π₀ ≈ 72 Pa (range 48–107 Pa).**
Central 72 Pa; low 48; high 107. This is **~1.8× the current 40 Pa**, which sits *below* the entire MCF7 range.

- **Provenance:** derived by Young–Laplace `ΔP = 2·γ_eff/R` from the **direct MCF7** effective cortical
  tension γ_eff = **0.27 mN/m** (IQR 0.18–0.40, n=27), Hosseini 2020 (Adv Sci 7:2001276,
  doi:10.1002/advs.202001276), at the MCF7 suspended radius **R = 7.5 µm** (Wagner 2011, PMC3147247).
  KB anchor: **KB-6.1.7 / gate VG-gamma-MCF7-suspended** (P2 Card 3).
- **Confidence: MODERATE.** It is anchored to a *real, direct MCF7 measurement*, but γ_eff is figure-digitized
  (KB confidence Low, wide IQR) and the γ→ΔP inversion assumes Laplace equilibrium. There is **no direct MCF7
  or carcinoma turgor-pressure measurement** anywhere in the KB — turgor is hard to measure directly (it needs
  a rounding-pressure / micropipette assay; FF 2014 did it only for mitotic HeLa).
- **This is NOT a hand-swap config change.** Per CLAUDE.md (no gate-loosening / no changing a sourced value off
  its source without PI) and the existing yaml policy, the MCF7 Π₀ must be landed by the **PI-authored
  pressure-consistency gate** over {R, γ_cortex, γ_mem, ΔP} = **P2 Card 3, currently PI-pending**. Recommend PI
  ratify Card 3 (KB-6.1.7) and then land 72 Pa in ONE documented commit that *also* raises the resting-myosin
  setpoint and re-runs GATE-A/GATE-B (see §4). Until then, 40 Pa stays a labeled HeLa diagnostic proxy.

---

## 1. How Π₀ enters the mechanics (confirmed operative relation)

`assemble.py`:
- `R_CELL_UM = 7.5` — MCF7 suspended outer (membrane) radius (Wagner 2011, audit row 27, SOURCED).
- `PI_0_PA = 40.0` — resting turgor Π₀, engine comment *"1 Pa == 1 pN/µm² so 40 in engine pressure units"*.
- Seeded as the Biot substrate mean pressure `p_bar = p = Π₀` (assemble.py:32); a uniform field has zero bulk
  `grad(p)` but exerts the boundary traction `p·n` on the live membrane via `MembranePressureTraction`.

**Operative force balance (confirmed).** For a pressurized thin shell the interior excess pressure ΔP is
balanced by the surface tension γ via **Young–Laplace**:

```
    ΔP = 2 γ / R          ⇔          γ = ΔP · R / 2            (thin sphere, pressure inside)
```

This is the relation the codebase uses throughout: `params_i0b1.yaml` (`gamma = dP*R/2`, lines 172/199),
`bleb_experiment.py` (`a_crit = 2 γ_mem/ΔP`), and the `resting_bound_myosin` docstring in `assemble.py:262`
(*"γ_cortex = ΔP·R/2 − γ_mem the ERM transmits"*). So **Π₀ sets the ENTIRE resting cortical-tension budget**:
the total effective tension the shell must carry is `γ_total = ΔP·R/2`, split into the membrane bilayer tension
γ_mem (sourced 10 pN/µm = 0.01 mN/m, KB-3.B1.1) + the actomyosin cortex tension γ_cortex.

---

## 2. Sourced γ and direct-turgor values (with DOIs / KB-ids)

| Quantity | Value | Cell | Source / KB-id | Status |
|---|---|---|---|---|
| **Direct turgor (only one on record)** | ΔP = **40 ± 30 Pa** (γ_surf 0.17 ± 0.13 mN/m) | **HeLa** (mitotic-rounded) | Fischer-Friedrich 2014, Sci Rep 4:6213 (PMC4148660); **KB-6.1.6** (PI-ratified as proxy) / KB-DRAFT-3.B-30 | HeLa proxy; *"NO MCF7 datum"*. This is the current 40 Pa. |
| **MCF7 effective cortical tension (primary)** | γ_eff = **0.27 mN/m** (IQR 0.18–0.40, n=27) | **MCF7** (suspended interphase, AFM parallel-plate) | Hosseini 2020, Adv Sci 7:2001276, doi:10.1002/advs.202001276; **KB-6.1.7 / VG-gamma-MCF7-suspended** | Direct MCF7; figure-digitized, KB confidence Low. **γ_eff is TOTAL/effective (γ_mem already folded in).** |
| MCF7 active cortical tension (upper cross-check) | γ ≈ 0.41 mN/m | MCF7 suspended | Hosseini 2021, Biophys J (PMC8391033) | Upper cross-check, protocol reconciliation pending. |
| Generic cortex band (NOT MCF7) | γ ≈ 0.5–1 mN/m | generic | Salbreux 2012 (via KB-3.5) | Too generic to anchor Π₀ (see §3 sensitivity). |
| Modeling-contract target (NOT a measurement) | γ = 0.35–0.65 mN/m | — | KU-3.5 / CORTICAL_TENSION_RECORD_2026-06-30 | A gate target, not literature; Hosseini 0.27 sits below it. |

TAG confirmed there is **no per-line MCF7 cortical-tension value stored as a verdicted constant** — the KB flags
it "figure-locked/scarce"; the 0.27 mN/m is the digitized project record behind KB-6.1.7. And there is **no
direct MCF7/carcinoma turgor** — the only turgor anchor is the HeLa FF 2014 value.

---

## 3. Units-explicit reconciliation (this is where a slip hides)

**Unit facts (verified):**
- `1 mN/m = 1e-3 N/m`. `1 pN/µm = 1e-12 N / 1e-6 m = 1e-6 N/m = 1e-3 mN/m` ⇒ **1 mN/m = 1000 pN/µm**.
  (So γ_eff = 0.27 mN/m = **270 pN/µm**; γ_mem = 10 pN/µm = 0.01 mN/m.)
- `1 Pa = 1 N/m² = 1e-12 N / 1e-12 m² = 1 pN/µm²` ⇒ the engine comment is correct.
- R = 7.5 µm = 7.5e-6 m.

**Route (b): DERIVE Π₀ = 2γ_eff/R from the MCF7 γ — all in SI, no shortcuts:**

```
central:  ΔP = 2 · 0.27e-3 N/m / 7.5e-6 m = 0.54e-3 / 7.5e-6 = 72.0 Pa
low:      ΔP = 2 · 0.18e-3      / 7.5e-6   = 0.36e-3 / 7.5e-6 = 48.0 Pa
high:     ΔP = 2 · 0.40e-3      / 7.5e-6   = 0.80e-3 / 7.5e-6 = 106.7 Pa  ≈ 107 Pa
```

→ **Π₀ = 72 Pa (48–107).** Cross-check in engine units: γ_eff = 270 pN/µm, R = 7.5 µm ⇒
ΔP = 2·270/7.5 = **72 pN/µm² = 72 Pa.** ✓ (identical — the two unit systems agree, no slip.)

**Double-count check (resolved).** Hosseini's γ_eff is the *effective/total* surface tension of the living,
myosin-active MCF7 cell — γ_mem is already inside it. So the correct **single reading is `ΔP = 2·γ_eff/R`;
do NOT re-add γ_mem** (`ΔP = 2(γ_cortex+γ_mem)/R` would over-count the membrane). Numerically it barely matters
here — γ_mem = 0.01 mN/m is only ~4% of 0.27, i.e. 72 vs ~77 Pa — but the *reading* matters for §4. This
matches the already-PI-ratified note in `params_i0b1.yaml:185` (`hosseini2020_gamma_link`).

**Is 40 Pa consistent with any sourced MCF7 γ?** No.
`40 Pa ⇒ γ_total = ΔP·R/2 = 40 · 7.5e-6 / 2 = 1.5e-4 N/m = 0.15 mN/m` — **below the bottom of the Hosseini IQR
(0.18)**, i.e. below the entire MCF7 range, ~45% under the central. 40 Pa is self-consistent only with the
**HeLa** γ_surf = 0.17 mN/m at HeLa geometry. It is "sourced" for the wrong cell — exactly the
unphysical-baseline → meaningless-comparison failure the CLAUDE.md HARD rule names.

**Sensitivity to WHICH γ (the honest uncertainty of route b):**

| γ used | mN/m | Π₀ = 2γ/R at R=7.5 µm | Comment |
|---|---|---|---|
| HeLa (FF 2014) | 0.17 | **45 Pa** | ≈ current 40 Pa; wrong cell. |
| **MCF7 Hosseini 2020 (recommended)** | **0.27** | **72 Pa** | Direct MCF7, effective tension. |
| MCF7 Hosseini 2021 (cross-check) | 0.41 | 109 Pa | Upper bound. |
| Generic Salbreux band | 0.5–1.0 | 133–267 Pa | Too generic; overshoots — do NOT use. |
| KU-3.5 modeling target | 0.35–0.65 | 93–173 Pa | A gate target, not a measurement. |

The spread shows the derived route is only as good as the γ it uses; the MCF7-specific Hosseini value is the
right anchor, and it is what gives 72 Pa. (Note: the naive example in the tasking, γ=0.5 mN/m → 133 Pa, uses
the *generic* band, not the MCF7 γ — that is exactly the trap: use the cell-specific γ.)

---

## 4. Coupling to the myosin (flag — Π₀ and the resting-myosin setpoint are NOT independent)

The resting cortical tension is **not** `γ_total = γ_turgor + γ_myosin` (parallel additive contributors).
Turgor is the outward **load**; the tension (membrane + actomyosin) is the **response that balances it** — they
sit on opposite sides of the Laplace force balance and are equal-and-opposite at equilibrium:

```
    ΔP · R / 2   =   γ_mem  +  γ_cortex(actomyosin)          (equilibrium)
    γ_cortex_required  =  ΔP·R/2 − γ_mem
```

| Π₀ | γ_total = ΔP·R/2 | γ_mem | **γ_cortex (myosin) required** |
|---|---|---|---|
| 40 Pa (current) | 150 pN/µm = 0.15 mN/m | 10 pN/µm | 140 pN/µm = **0.14 mN/m** |
| **72 Pa (recommended)** | 270 pN/µm = 0.27 mN/m | 10 pN/µm | 260 pN/µm = **0.26 mN/m** |

**So raising Π₀ 40→72 raises the actomyosin target ~1.86× — the myosin contribution does NOT drop, it must
RISE.** The intuition "raise turgor so myosin can relax" is the wrong framing: a higher-turgor cell needs a
*higher* cortical tension to contain it. Consequences to sequence in ONE commit, not piecemeal:

1. The `resting_bound_myosin_fraction × resting_bound_myosin_force_pn` setpoint (currently a PI-GAP, both
   `None`) that delivers γ_cortex must be scaled to ~0.26 mN/m, not 0.14. Setting Π₀=72 while leaving the
   myosin at the 40-Pa level leaves the cortex out of equilibrium (or dumps the imbalance onto the membrane —
   the `realization_gap` already noted in `params_i0b1.yaml:195`).
2. Consistency, not double-counting: because Hosseini's γ_eff = γ_mem + γ_cortex = ΔP·R/2, setting Π₀=72 AND a
   myosin setpoint that delivers ~0.26 mN/m is *self-consistent* — you are matching the same measured γ_eff
   from both sides. The genuine double-count risk is only if someone derived Π₀ from the full γ_eff and *then*
   also added γ_mem on top, or added myosin *on top of* a turgor already meant to represent the total tension.
3. **GATE-B interaction:** the emergent dynamic γ was validated against the 40-Pa (0.14 mN/m) target. At 72 Pa
   the whole resting tension budget rises ~1.8×, so the GATE-B emergent-γ target rises correspondingly and the
   milestone should be re-stated on the 72-Pa baseline (already caveated in the param audit §"GATE-A/GATE-B").

---

## 5. Recommendation

- **Value:** `PI_0_PA` should be **72 Pa** (uncertainty band **48–107 Pa**), from Hosseini 2020 MCF7 γ_eff =
  0.27 mN/m (IQR 0.18–0.40) via `ΔP = 2·γ_eff/R` at R = 7.5 µm, **γ_mem NOT re-added**.
- **Confidence:** MODERATE. Anchored to a direct MCF7 measurement, but γ_eff is figure-digitized (KB Low, wide
  IQR) and the inversion assumes Laplace equilibrium; there is no direct MCF7 turgor measurement to
  cross-check. Report Π₀ with its 48–107 band, never as a bare point value.
- **Config change or PI ratification first? → PI ratification first.** This is the single most load-bearing
  baseline; it also shifts the myosin target (§4) and the already-banked GATE-A/GATE-B γ magnitudes. CLAUDE.md
  forbids swapping a sourced value off-source without PI, and the yaml already records the standing PI policy:
  *"the MCF7 dP (48–107 Pa) is landed ONLY by the PI-authored pressure-consistency gate over
  {R, γ_cortex, γ_mem, dP_hyd}, NEVER a hand-swap."* The mechanism already exists as **P2 Card 3 /
  KB-6.1.7 / VG-gamma-MCF7-suspended (PI-pending).**
- **Concrete path:** (1) PI ratifies Card 3 (register the Hosseini 2020 KnowledgeClaim KB-6.1.7 + the
  pressure-consistency gate). (2) In ONE documented commit: set `PI_0_PA = 72.0` with the Hosseini/Laplace
  provenance, set the resting-myosin setpoint to deliver γ_cortex ≈ 0.26 mN/m, engage radial ERM pairing, and
  re-run GATE-A/GATE-B on the 72-Pa baseline. (3) Until (1), keep 40 Pa as the explicitly labeled HeLa
  diagnostic proxy — do not silently edit it.

---

## Appendix — files / anchors

- `aleph/components/incumbent/assemble.py:117,128,262` — R_CELL_UM, PI_0_PA, resting-myosin γ_cortex=ΔP·R/2−γ_mem.
- `aleph/components/fluid/params_i0b1.yaml:152–204` — Q3 `delta_p_hydrostatic_rest`: FF-HeLa source, the
  `hosseini2020_gamma_link` (72 Pa reading, γ_mem-not-re-added), `mcf7_consistency_gate`, `realization_gap`.
- `aleph/components/incumbent/compartments.py:121` — γ_mem = 10 pN/µm (KB-3.B1.1).
- `aleph/docs/v2_audit/AC_DECISION_CARDS_2026-07-22.md` Card 3 — Hosseini γ registration + double-count check.
- `aleph/docs/v2_audit/PARAM_PROVENANCE_AUDIT_2026-07-24.md` row 1 — the trigger.
- KB: KB-6.1.6 (FF HeLa proxy, PI-ratified) · KB-DRAFT-3.B-30 (turgor 40±30 Pa) · KB-6.1.7 /
  VG-gamma-MCF7-suspended (Hosseini 2020 MCF7 γ_eff, PI-pending) · KB-3.B1.1 (γ_mem) · Wagner 2011 (R=7.5 µm).
