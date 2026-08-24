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

# PI Decision Queue — 2026-05-30 (autonomous-build staging)

> Items the autonomous /loop session staged but did **not** apply, because they touch governed surfaces (gate-contract / frozen-integrator / `ffn/foundation` push / magic-number trigger). CLAUDE.md hard rules forbid silent override. Each has a recommendation; PI ratifies on wake.
> Source: mechanism audit [`MECHANISM_AUDIT_2026-05-30.md`](MECHANISM_AUDIT_2026-05-30.md). Autonomous-applied (additive, non-governed) work is tracked in the session todo + commits.

## A. Gate-contract changes (no-gate-loosening rule → PI only)

| # | Item | Current | Proposed | Rationale | Audit ref |
|---|---|---|---|---|---|
| A1 | **Filamin slip→catch** + revise the slip-sign Sanity-Gate test | `test_crosslinkers.py` asserts `dk_off/dF>0` (slip) for filamin | two-pathway catch-slip; test → catch contract | Filamin (70% of cortex pool) is a documented catch bond (Ehrlicher 2011; Rognoni 2012; Gieseke 2013). A pure-slip cortex sheds crosslinks under load where the real cortex stiffens → plausible KU-3.5 contributor | §2, §3 |
| A2 | **Myosin binned-r0 Hill ratchet → physical stepping** + add n_bins grid-invariance gate | Hill ratchets a binned bond rest-length; emergent γ depends on (n_bins, bin_width) | translate bound anchor by v·batch_dt; γ invariant under n_bins sweep | Grid-dependent runtime construct directly feeding the failing KU-3.5 gate | §2, §4#9 |
| A3 | **KU-5.1 density band** ratification | placeholder [50,150] /µm² (vis decoration only) | literature-anchored band ~100/µm² (Abraham 1999; Bieling 2016) | Needs a real band before KU-5.1 can PASS/FAIL | §6 |
| A4 | **KU-3.5 interphase vs metaphase** | band [0.35,0.65] mN/m | clarify regime | enclosed-volume pressure (co-req) implies metaphase γ~1.6 mN/m > current ceiling → KU-3.5 may be interphase-specific | §5#2 |
| A5 | **Arp2/3 debranching** breaks the "monotonic bond count" invariant | bond count grows monotonically | allow force-dependent debranching (Pandit 2020) | Bounds n_total + adds real physics; current invariant would fail | §5#5 |
| A6 | **D4 ECM rebanding** re-derivation | G_0 [1,50] Pa (PI-authorized 2026-05-21) | optionally tighten via re-derived percolation prefactor | transparency; already ratified, low priority | §2 |

| A7 | **S5 tag-space unification (FA + myosin coexistence)** | FA path raises `NotImplementedError` for FA + {myosin, xlink, lamellipodium} | generalize `IntegrinBondUpdater` tag→row lookup so integrins need not be tags `[0,n_int)`; FA block in natural append order | **EMPIRICALLY CONFIRMED 2026-05-30** by the S0–S2 wiring (commit `0f655f6`): reused `IntegrinBondUpdater` assumes integrins at global tags `[0,n_int)`, colliding with cortex-actin tag bookkeeping. This is the gating dependency for KU-3.5 v4 (cortex+myosin+FA). A subagent is implementing it (additive, no physics change). | §6 S5 |

> **Note on A7:** this is additive code (no governed-surface edit — no gate band, no frozen integrator, no Pereverzev/Kong param change), so it is being IMPLEMENTED autonomously, not held for PI. It is listed here only because it is the empirical confirmation of the §6 S5 design step and the gate that unblocks KU-3.5 v4. The genuinely PI-gated co-dependencies for v4 remain **B1** (global dt) + **B2** (equilibration prelude) + **A3/A4** (KU-5.1 band / KU-3.5 regime).

## B. Frozen-integrator changes (`integrator/` freeze → PI only)

| # | Item | Proposed | Rationale | Audit ref |
|---|---|---|---|---|
| B1 | **Global `dt = min(τ)` reconciliation** | single CFL across ECM+cortex+FA; assert SHAKE convergence vs dt | ECM stiffer bonds under cortex dt may violate CFL once wired | §1#15 |
| B2 | **`equilibrate_no_shear` prelude in Cell.build** | mandatory soft-start before BAOAB production loop | prevents step-0 force explosion when ECM/FA add construction overlap | §1#11 |
| B3 | **Production dt co-tune** (path B, PI 2026-05-30) | sweep max-stable dt on the FA-contracting v4 system | the biggest wall-time lever; must be tuned on the real target, not floating-shell v3 | opt-decision |
| B4 | **lazy-Fixman + banded-inv** (FIXMAN_LAZY_EVAL_DESIGN A2) | ~1.2–1.4×, bit-for-bit within Sanity-Gate tolerance | CPU-side speedup; needs parity re-validation | FIXMAN doc |
| B5 | **SHAKE λ warm-start** (new) | seed Newton solve with previous step's multipliers | targets the #1 measured hotspot (SHAKE 0.46 ms/step); biggest payoff in the loaded regime FA creates | profile 2026-05-30 |

## C. Parameter changes (magic-number / literature corrections)

Additive-documented in code comments where possible; the ones that change runtime physics are staged here.

| # | Param | Current | Corrected | Source | Note |
|---|---|---|---|---|---|
| C1 | Integrin catch F* | 6.99 pN (Pereverzev analytic) | ~30 pN (Kong 2009 two-state) | Kong 2009 | self-flagged in-config; do at FA wiring |
| C2 | Capping δ_cap | 0.3 nm (mislabeled "0.3 pN") | 2.7 nm (= δ_elong) | Li/Bieling eLife 2022 | single most useful KU-5.1 capping fix |
| C3 | Myosin v0 | 1.0 µm/s | ~0.2 µm/s NMIIA | Kovács 2003 | or document as deliberate mesoscale effective value |
| C4 | α-actinin k_off0 | 1.0 /s | ~0.066 /s (single-molecule) | Ferrer 2008 | or justify bulk value |
| C5 | Arp2/3 branch t0 | 72° | 70° (68±9°) | Fäßler 2020 | + derive k_angle from the distribution (replaces "TBD" magic number) |
| C6 | WAVE k_wave_pin | 1e-4 N/m (CFL magic number) | replace with explicit membrane OR document as CFL-bounded confinement | — | resolved by S7 membrane |

**Citation fixes (doc-only).** Pending a precise per-occurrence pass (NOT a blanket replace — "Funk 2022" is a phantom paper that maps to TWO real papers by context): capping (`k_cap_0`, `delta_cap`, "D1 capping") → **Li/Bieling eLife 2022**; abortive-branching threshold (`abortive_pressure_Pa`) → **Funk 2021 Nat Commun** (the real branching/CP-NPF paper). Locations: `configs/phase1_h5.yaml` (3×: lines ~59/61/63) + `cell/lamellipodium.py` (7×). Also: Furuike 2001 re-cite for filamin (wrong observable — measures Ig unfolding, not bond off-rate); k_ERM 1e-4 provenance note (it is the *correct* per-linker value per Alert 2015, not a softening error). Deferred from the 2026-05-30 autonomous pass to keep citation attribution exact.

## D. Governance

| # | Item | Status |
|---|---|---|
| D1 | `ffn/foundation` push (77 commits behind: H.2/H.3/H.4/H.5) | PI sign-off required |
| D2 | Stale branch cleanup (`phase1/h4-fa-clutch`, `phase1/h5-lamellipodium` are stale pointers) | PI confirm before deletion |
| D3 | `.sync-conflict-*` output shadows (10 cosmetic, output-only) | leave / PI confirm cleanup |

## Recommended order on PI wake
1. Ratify **A3 + A4** (KU-5.1 band + KU-3.5 regime) — unblocks the production targets.
2. Approve **B1 + B2** (dt reconciliation + equilibration prelude) — required for the FA wiring to run at all.
3. Approve **C1–C5** as a batch (literature corrections) — applied at FA wiring.
4. **A1 + A2** (filamin catch + myosin stepping) — correctness, can follow first v4 smoke.
5. Defer **B3–B5** (speed) to after a correct v4; **D1–D3** (governance) at session closeout.
