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

# DCM real-timescale compaction — cadherin load-sharing cluster + maturation-gated rearrangement (DESIGN)

**Date:** 2026-07-09 · **Owner:** Lead session · **Branch:** `dcm/main` · **Status:** DESIGN (PI reopened decision (c); design before implementing per CLAUDE.md)
**Prereq / supersedes-scope:** `DCM_CADHERIN_MATURATION_2026-07-09.md` (maturation built but does not engage at lit τ),
`project-dcm-compaction-turgor-blocked` (aggregate-σ compaction is drag-limited, cadherin-independent),
`project-dcm-ipc-large-dt` (aggregate-σ end-state valid; RATE wrong for the wrong reason).

---

## 0. Why we are here (PI decision, and the finding that forces this design)

PI stopped the DCM line at decision (c): the aggregate-σ (Foty-Steinberg liquid-drop) compaction is valid as a
**mechanical END-STATE only** — dense sphere, virial stress, porosity end-value, condition comparisons — but its
**RATE is wrong for the wrong reason**: ~100 s = drag-limited radial pull, not the real 24–48 h set by cadherin
junction remodelling. Two decision tests nailed why the current model can't produce the real rate:

1. **Freezing the junction leaves aggregate-σ compaction unchanged** (porosity 0.236→0.234). The liquid-drop pull
   densifies a loose aggregate by gap-closing + turgor-limited cell deformation — it never needs a cadherin bond to
   break. So compaction is cadherin-**independent** → cadherin can't be its rate-limiter.
2. **Maturation never engages at the lit timescale.** A junction that persists should mature (E-cadherin
   clustering + α-catenin/vinculin reinforcement, KB-4.3) and its lifetime should rise to 5–30 min (KB-4.11). But
   the current bond **ruptures in ~0.036–0.74 s**, long before `τ_mature=600 s` → `m=1−exp(−age/τ_mature)≈0` → the
   maturation blend is a no-op.

**Root cause (read from the code, 2026-07-09):** `dcm/dcm_cadherin_gpu.py::cad_break_kernel` evaluates the catch-slip
`k_off` at the correct **per-molecule** load `F₁ = k_trans·(L−r0)` (parallel springs at a common extension all bear
`F₁`, independent of the count) — that part is right. The flaw is the **survival semantics**: a single Bernoulli draw
`p_break = 1−exp(−k_off(F₁)·dt)` removes the **entire bundle of `bundle_n` molecules at once**. A real cadherin cluster
is `N_b` parallel trans-dimers that **share the load**; one molecule unbinding leaves `N_b−1` still bonded (which
re-share the load) and the empty slot **rebinds** at `k_on` — the junction dies **only when all `N_b` are
simultaneously unbound (`m→0`)**. Collective cluster lifetime `T(N_b) ≫ 1/k_off`. This is the memory's "bundle FORCE
is ×N but BREAK is single-molecule": the current model is a **lumped proxy** (whole-bundle-at-single-rate) for a
fine-grained cluster — exactly the abstraction the CLAUDE.md hard rule forbids.

## 1. Target physics — parallel-bond cluster with load sharing (Erdmann–Schwarz; Bell/Rakshit per molecule)

A node-pair junction is a **cluster of `N_b` parallel E-cadherin trans-dimers** gripping a shared extension `L`:

- **Per-molecule load** `F₁ = k_trans·max(0, L−r0)` — the same for all engaged molecules (parallel springs, common
  extension). The transmitted junction force is `F_junc = m·k_trans·(L−r0)` with `m` = the **currently engaged count**
  (NOT a fixed `bundle_n`).
- **Per-molecule unbinding** at the faithful Rakshit-2012 catch-slip rate `ε(F₁) = effective_k_off(F₁)`
  (`validation/cadherin_sliding_rebinding.py`, exact X-dimer generator, unchanged — this stays the molecular oracle).
- **Per-empty-slot rebinding** at `k_on` (rest-symmetric `k_on = k_off(0) ≈ 28/s`; a mature junction reinforces this).
- **Cluster state** `m ∈ {0,…,N_b}` is a birth–death chain: death `d_m = m·ε(F₁)`, birth `b_m = (N_b−m)·k_on`.
  The junction is removed only at `m=0`. Mean first-passage time `T(N_b)` = the **emergent junction lifetime**.

This is the standard adhesion-cluster dynamics (Erdmann & Schwarz 2004 *PRL* 92:108102; Bell 1978). It is the
fine-grained realization the hard rule mandates (individual molecules bind/unbind/share-load; the junction lifetime is
**emergent**, not a lumped single rate).

## 2. ⚠️ The load-sharing FINDING that reshapes the design (surface to PI)

Computing `T(N_b)` analytically (BD-MFPT, `effective_k_off` per molecule, `k_on=27.96/s`) at physiological
per-molecule loads gives:

| per-molecule load | ε(F₁) [1/s] | ρ=k_on/ε | T(N_b=5) | T(N_b=10) | T(N_b=20) | T(N_b=50) | T(N_b=100) |
|---|---|---|---|---|---|---|---|
| 0 pN   | 27.96 | 1.00 | 0.31 s | 4.2 s  | **33 min** | 2.6e4 yr | 1.4e19 yr |
| 10 pN  | 51.6  | 0.54 | 0.09 s | 0.37 s | 11.7 s | 22 d | 7.6e7 yr |
| 20 pN  | 65.0  | 0.43 | 0.06 s | 0.19 s | 2.7 s  | 12 h | 4e4 yr |
| 29 pN (catch peak) | 18.4 | 1.52 | 0.89 s | 40 s | **54 h** | 2.5e9 yr | — |

**Two things this shows.**
1. **The mechanism WORKS in principle**: load-sharing lifts the junction from the single-molecule ~0.036 s into the
   KB-4.11 5–30 min band — which a single-rate break can never do. The band is crossed near `N_b ≈ 15–25`.
2. **But `T(N_b)` is super-exponentially sensitive to `N_b` and load near the `ρ~1` critical transition.** The KB
   density anchor `N_b=100` (ρ_cad·A_junction, KB-4.1/4.11/4.17, the fidelity-remediation `BUNDLE=100`) gives
   `T ~ 10²⁶ s` = **permanent bonds → tissue frozen → NO rearrangement → NO compaction**. Picking `N_b` to land on
   30 min would be **tuning to a target = a magic-number-rule violation**. So the raw cluster lifetime is **not a
   robust, derivable rate-limiter** — it is poised at a critical point.

**Resolution (the design pivot).** The robust, **lit-anchored** rate-limiter is the **maturation timescale
`τ_mature` (KB-4.11, 5–30 min)**, not the fragile cluster lifetime:

- A **nascent** junction has a **small** engaged capacity `N_b^nascent` (few molecules) → short cluster lifetime →
  it **turns over fast → allows neighbour exchange (rearrangement)**.
- If the junction **persists under load**, it **matures** over `τ_mature`: recruited capacity grows
  `N_b(age): N_b^nascent → N_b^mature` (E-cadherin clustering) → cluster lifetime rises **steeply** → the junction
  **locks**.
- So the rate at which junctions lock = `1/τ_mature` ≈ 5–30 min = **robust + lit-anchored** (independent of the
  fragile `T(N_b)` value). Compaction proceeds only while junctions are young/unlocked; the cumulative
  maturation-gated rearrangement sets the **min–hr → 24–48 h** tissue timescale.
- **Crucially, (a) enables (b):** load-sharing is what lets a nascent junction survive the ~0.036 s single-molecule
  death long enough to *age toward `τ_mature`*. Without (a), maturation is dead on arrival (the current bug).

`N_b^nascent`/`N_b^mature` are anchored to nascent-vs-mature cadherin cluster sizes (KB-4.1/4.3); the **rate** the
model reports is `τ_mature` (measured/lit), and the emergent `T(N_b)` is a **cross-check** that nascent turns over
(< s) and mature locks (≫ min), NOT a tuned number. This keeps the no-magic-number contract: we anchor `τ_mature` and
the two capacities to literature and **measure** the compaction rate — we never pick a value to hit a rate.

## 3. Mechanism — concrete

Per node-pair junction bond `b`, carry two per-bond fields alongside the existing `age`:
- `m_b` — currently engaged molecule count (float/int in `[0, N_b(age)]`).
- `N_b(age_b) = N_nascent + (N_mature − N_nascent)·(1 − exp(−age_b/τ_mature))` — maturing capacity (default off →
  `N_b ≡ bundle_n`, back-compat).

Each KMC sub-step `δt` (the existing `_n_subcycle` cadence, unchanged):
1. `F₁ = k_trans·max(0, L−r0)` (per molecule) → `ε = effective_k_off(F₁)` (existing table).
2. `n_off ~ Binomial(m_b, 1−exp(−ε·δt))` (each engaged molecule unbinds independently — exact Bernoulli sum in-kernel).
3. `n_on ~ Binomial(N_b−m_b, 1−exp(−k_on·δt))` (each empty slot rebinds).
4. `m_b ← m_b − n_off + n_on`; `age_b ← age_b + δt`.
5. **Junction removed iff `m_b == 0`** (absorbing; the pair returns to the free pool and may re-form a fresh nascent
   junction via the existing mutual-nearest `cad_form`/`cad_partner`).
6. **Force** applied that step: `F_junc = m_b·k_trans·(L−r0)` (engaged-count-scaled), replacing the fixed
   `bundle_n·k_trans`. `f_contract` likewise scales with `m_b`.

**Reuse unchanged:** the Rakshit `effective_k_off` table + `build_koff_device`; `cad_partner_kernel` /
`cad_form_kernel` (mutual-nearest formation — a new junction starts at `m = N_nascent`, `age=0`); `_n_subcycle`
(the sub-cycling cadence is if anything *less* critical now that the junction lifetime is long, but keep it for the
fast nascent population); the aggregate-σ / turgor / IPC stack. **New:** the cluster birth–death break kernel
(`cad_break_cluster_kernel`) + `m`/`N_b` per-bond arrays. **Modify:** the force launch (engaged-count `m_b`, not
`bundle_n`) at `dcm_warp_decohesion.py:1217`.

## 4. Validation gates (write BEFORE implementing; do not loosen; PI-authored contracts)

- **G1 — cluster-lifetime oracle (analytic):** the in-kernel stochastic cluster mean lifetime `T(N_b, F₁)` matches
  the analytic BD-MFPT (§2) within Monte-Carlo error, at `N_b ∈ {1,5,10,20}` and `F₁ ∈ {0,10,20,29} pN`.
  `N_b=1` reproduces the single-molecule `1/effective_k_off` (byte-continuity with the old model).
- **G2 — nascent turns over, mature locks:** at the lit `N_nascent`/`N_mature`, emergent nascent lifetime `< 1 s`
  (rearrangement-permissive) and mature lifetime `≫ τ_mature` (locking) — the qualitative separation, not a tuned value.
- **G3 — maturation now ENGAGES:** with load-sharing on, the mean junction age distribution reaches `~τ_mature`
  (junctions survive to mature), i.e. `⟨m·1[age>τ_mature]⟩ > 0` — directly fixing the `DCM_CADHERIN_MATURATION`
  no-engage finding.
- **G4 — compaction RATE is maturation-gated (the payoff):** in a compacting run, the compaction half-time scales
  with `τ_mature` (sweep `τ_mature` → half-time tracks it), NOT with `accel_dt`/drag. Freezing maturation
  (`N_b≡const`) recovers a different (faster, drag-limited) rate — i.e. cadherin now **does** gate the rate (inverts
  the decision-test result). *This is the gate that certifies the redesign; it is where (b) is proven.*
- **G5 — catch peak preserved:** `catch_peak_force` still ≈ 29 pN per molecule (the molecular oracle is untouched).
- **G6 — back-compat:** `cluster=False` (or `N_b=1`) → byte-identical to the current break (same RNG draw order).
- **G7 — end-state unchanged:** the mechanical end-state (porosity/σ_vm/shape) at full compaction is consistent with
  the aggregate-σ result (we change the RATE and its cause, not the validated end-state) — unless (b) reveals the
  end-state itself was rearrangement-dependent, in which case surface to PI.
- **G8 — NATIVE (A5000):** `N=400`/`N=2000` on cuda:0 — stable (no NaN), maturation-gated rate, pen-free.

## 5. Staging (isolate-first, like the FF treadmill)

- **S1 — cluster break, HOST path + G1/G6 oracle.** Add `m`/`N_b` to `CadherinBondHost._tick` (numpy Binomial); gate
  the emergent lifetime vs the analytic BD-MFPT; `N_b=1` byte-identical. **[this session]**
- **S2 — maturation-capacity `N_b(age)` + G2/G3.** Grow `N_b` over `τ_mature`; show nascent turns over / mature locks
  / maturation engages.
- **S3 — GPU kernel (`cad_break_cluster_kernel`) + engaged-count force + G5/G6 parity** (CPU==GPU, catch peak, back-compat).
- **S4 — compaction coupling + G4 (the payoff).** Run loose→compact with cluster+maturation; show the rate is
  `τ_mature`-gated. **Decision point (b):** if the aggregate-σ pull still bypasses junctions (G4 fails because
  compaction doesn't require bond breaking), escalate to explicit **cadherin-tension→neighbour-exchange (T1)** so the
  compaction path *requires* rearrangement — surface to PI before that deeper change.
- **S5 — NATIVE (A5000) + G8.**

## 6. Risks / open questions (surface to PI)

1. **The `T(N_b)` critical sensitivity (§2) is the headline risk.** The design deliberately makes `τ_mature`
   (robust, lit) the rate-limiter and treats `T(N_b)` as a cross-check — but if the emergent behaviour still hinges on
   the exact `N_nascent`/`N_mature`, that is a magic-number exposure → surface, do not tune.
2. **Does compaction actually require rearrangement?** Decision-test 1 says the loose→dense aggregate-σ path is pure
   gap-closing (no bond breaking) → long-lived junctions may NOT slow it. If so, G4 fails and (b) needs the explicit
   T1/neighbour-exchange path (S4 decision point). This is the crux of whether (a) alone suffices or (b) is a separate build.
3. **Binomial in-kernel** for large `N_b` — Bernoulli-sum is O(N_b) per bond; fine at N_b≤~100, but if mature N_b is
   huge, use a Gaussian/Poisson approx above a threshold (guard the approximation error).
4. **`k_on` for an empty slot** may be < the rest `k_on` if the membranes locally separate on unbinding; keep the
   rest-symmetric `k_on` as the anchor and flag if the cluster is unrealistically stable.
5. **`f_contract`/`bundle_n` semantics change** (force now scales with engaged `m_b`, not fixed `bundle_n`) — audit
   every call-site that read `bundle_n` as a force multiplier (`dcm_warp_decohesion.py:1219,1254,1491,1780`).

## 6b. S2 RESULT (2026-07-09) — maturation must be a CONTACT property, not a bond property

Building S2 surfaced a second finding: **bond-level age does NOT make maturation engage, even with the S1
load-sharing cluster.** A nascent cluster's lifetime (~0.24 s at rest, even load-sharing) is « τ_mature=600 s,
so an individual trans-dimer bond ruptures long before its `age` reaches τ_mature — the capacity `N_b(age)` barely
grows (`1−exp(−0.24/600)≈4e-4`) before the bond dies. This is the *same* `DCM_CADHERIN_MATURATION` no-engage wall,
now shown to persist under load sharing.

**Fix (implemented):** maturation is a property of the sustained **CONTACT (apposition)**, not a single bond's
uninterrupted lifetime. Track per-**node** `contact_age` (incremented while the node is apposed to / bonded with
another cell, reset when it leaves contact); a bond's capacity `N_b` derives from its endpoints' `contact_age`. A
persisting contact then accumulates maturation **through bond turnover** (break→reform, which the host already does
via mutual-nearest re-formation) → `N_b` climbs nascent→mature over τ_mature → the junction **locks**, at the robust
rate 1/τ_mature. **Gates PASS** (`test_cadherin_cluster.py`, 29 total): G2 nascent(<1 s)/mature(≫τ) separation,
**G3 maturation ENGAGES** (a sustained contact matures to n_mature and its turnover collapses; single-molecule
P(survive to τ)~e⁻⁶¹≈0), G3b contact-loss resets maturation. Figure `dcm_cadherin_maturation_engages.png`
(N_b 4→25, cluster 2 breaks vs single-molecule 141 over 10 s). This closes the "maturation never engages" problem;
`n_nascent` remains a controlled variable flagged for a firm KB-4.3 anchor.

## 6c. S4 RESULT (2026-07-10) — G4: aggregate-σ compaction is NOT maturation-rate-limited → (b) needed

Native N=400 loose-fcc, σ=5 mN/m, cluster ON, 48 s (STEPS=6000) on gbook A5000, 3 conditions sharing σ and
differing only in maturation: **no maturation** / **τ_mature=30 s (engages within the run)** / **τ_mature=600 s
(lit; barely engages)**. Figure `dcm_cadherin_s4_compaction_rate.png` (porosity(t) + Rg(t)).

| condition | porosity 0→48 s | saturation |
|---|---|---|
| cluster, no maturation | 0.631 → **0.392** | saturates by ~23 s |
| cluster + mat τ=30 s | 0.631 → **0.341** | still declining at 48 s |
| cluster + mat τ=600 s | 0.631 → **0.351** | still declining at 48 s |

**Decisive: the compaction RATE is NOT gated by τ_mature.** The τ=30 s and τ=600 s trajectories are **nearly
identical** (0.341 vs 0.351, curves overlap) despite a **20× difference in τ_mature** — if maturation set the
timescale, τ=600 s would compact ~20× slower; it does not. Two phases: (1) the first ~20 s is IDENTICAL across all
three = fast **drag-limited gap-closing** (the σ liquid-drop pulling loose cells together — needs no bond breaking,
confirming decision-test ①); (2) the maturation conditions then continue past where no-maturation saturates, but
that extra compaction is driven by junction STRENGTH/nucleation (nascent m=4 « full m=20 → weaker, more-dynamic
adhesion lets σ densify further), on the σ/drag timescale — NOT on τ_mature (else τ=30 s ≠ τ=600 s).

**G4 VERDICT: FAIL — aggregate-σ compaction bypasses the junction rate-limiter (as decision-test ① predicted, now
confirmed with the full cluster+maturation redesign at native scale).** The S1–S3 mechanism is correct and does
change the junction physics (lifetime emergent, maturation engages, force ∝ engaged m), but the COMPACTION DRIVER
(Foty liquid-drop σ) densifies by radial gap-closing on the drag timescale, so no cadherin property can make its
RATE ~τ_mature. **→ the (b) escalation is required** (design §5 S4 decision point): to get a genuinely
maturation-rate-limited (min–hr) compaction, the σ liquid-drop drive must be SUBORDINATED to a mechanism where
densification REQUIRES junction rearrangement (explicit cadherin-tension → T1 neighbour-exchange), so long-lived
matured junctions actually throttle the rate. This is a PI-scoped deeper change (surfaced with this data).

## 7. One-line summary
Replace the lumped "whole cadherin bundle breaks at the single-molecule rate" with a **fine-grained load-sharing
parallel-bond cluster** (`m` engaged molecules, birth–death, junction dies only at `m→0`) so junctions survive long
enough to **mature**; make the **maturation timescale `τ_mature` (lit-anchored, robust) the compaction rate-limiter**
(nascent turns over → mature locks), with the super-exponentially-sensitive raw cluster lifetime used only as a
cross-check — turning the drag-limited ~100 s end-state into a maturation-gated min–hr rearrangement rate.
