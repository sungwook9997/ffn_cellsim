# SF tension generation-limit — synthesis & γ-floor unification (2026-06-09)

**Session (ii) "ACCEPT THE FLOOR"** · branch `h7/sf-generation-floor` (from
`h7/compartment-platform`). This is the honest fallback half of the PI session-split
(PLATFORM_PI_QUEUE §"SESSION SPLIT"): document the ventral-stress-fiber (SF) tension as
**generation-limited**, unify it with the cortical active-γ floor as **one ½·n·f·ℓ
generation-limit story with one missing MCF7 density datum**, and surface the gate-reframe
options to PI **without loosening any gate**. The sibling session (i) "CLOSE THE FLOOR"
chases the SF NMII density datum that could close Kumar via Route B; if it lands a
defensible datum this conclusion updates toward it (see §6 Reconciliation).

Cross-links: `H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09.md` (cortical γ-floor, the sibling
limit), `H7_MYOSIN_OVERLAP_MECHANISM_DESIGN_2026-06-09.md` §9–§11 (per-head fix + SF-2c
sarcomeric result), `BASAL_MESH_DESIGN_2026-06-09.md` (the 2-layer apparatus this measures),
`H7_MYOSIN_DENSITY_LITERATURE_2026-06-09.md` (the density deep-research). Tool:
`scripts/h7_basal_sf_force_budget.py` (now emits the labeled decomposition below).

> **✅ RECONCILED with session (i) 2026-06-09 (commit `f6b956c`, branch
> `h7/sf-nmii-forcescale`, doc `H7_SF_NMII_FORCESCALE_RESULT_2026-06-09.md`).** Session (i)
> ran the density deep-research (104-agent, 3-vote verify) and **REFUTE'd** Kumar-closability:
> the molecular content is solid (~30 heads/side, ~300 nm; Billington 2013 et al.) but the
> **parallel minifilament count per SF cross-section — the one multiplier Route B needs — has
> no usable datum and is GEOMETRICALLY IMPOSSIBLE** at the required magnitude (band needs
> ~590–1760; a 50–250 nm-radius cross-section holds O(5–15)). So the generation-bound
> conclusion below is **no longer provisional — both sessions independently confirm it and
> jointly HALT→PI.** Session (i) also lands a **decisive gate-reframe** (Kassianidou/Kumar 2017
> PNAS): the Kumar 10–30 nN is mostly **network/prestress**; single-fiber **active** myosin is
> ≈5–6 nN, and in our platform that network/prestress is supplied by the **basal mesh + FA
> anchoring + turgor**, not `sf_myosin_` generation — so the active component is ~7× under even
> the 6 nN active target, irreducible without the impossible count. Folded into §5–§7 below.

---

## 1. One-paragraph answer

The single-SF tension is **motor-generation-limited**, exactly like the cortical active-γ
floor — same ½·n·f·ℓ aggregate-motor-force budget, same missing MCF7 NMII density datum.
The quasi-static force budget (`h7_basal_sf_force_budget.py`, REFUTE — confirmed decisive by
session (i), no longer provisional) gives a raw
mesoscale single-SF tension of **≈5.0 pN** (φ≈0.99 engaged × 10 heads × 0.5 pN brief-literal
stall), **~2014× under the Kumar 2006 10–30 nN band**. This is **N_filaments-INDEPENDENT**:
N_filaments / k_actin set the bundle's *strain* (ε = T/μ_SF) at a given tension, **not the
tension ceiling**, which is set by motor generation. The 2014× gap factors cleanly into named
physical causes (§3): a small per-minifilament parameter fidelity factor — **~3.4×** against
session (i)'s directly-anchored per-minifilament stall (17 pN, Stachowiak 2009; the defensible
value) up to **~11×** against the heads×per-head dipole upper bound (56 pN) — and a ~100×–~1760×
**cross-sectional NMII count**: how many native minifilaments must act in parallel across ONE SF
cross-section, i.e. the Route-B native:effective force-scale factor = the **missing MCF7 SF-NMII
density datum**. Session (i) showed this count datum **does not exist and is geometrically
impossible** at the required magnitude (a 50–250 nm-radius cross-section holds O(5–15)
minifilaments; the band needs ~590–1760) — the gap is therefore **irreducible**, not a tuning
artifact. And the budget is only an **upper bound
on a perfectly-organised bundle**: the §11 SF-2c result (random mixed-polarity bundle →
**−61 pN, slackening, not contraction**) shows that even with the correct per-motor force and
count, tension is not realised without **sarcomeric polarity organisation** — the SF analogue
of the cortical "tension is overlap/architecture-set, not count-set" finding.

---

## 2. What is established (build) vs what is generation-limited

The 2-layer basal contractile apparatus **BUILD is complete and gated** (B1–B4 + S1 +
bending; `BASAL_MESH_DESIGN`): connected basal mesh (giant 0.963, z 3.078, L/lc 6.22),
flat ventral surface positioning FA, F-layer cables anchored FA→FA force-free, `sf_myosin_`
NMII placed with the bipolar/Hill gate satisfiable, cortical-γ no-contamination held, and
the CFL wall solved (bending wired → bend-before-stretch ~1.5e5 → AFINES soft-stretch →
explicit BAOAB feasible). **The only remaining blocker to the live Kumar gate is the
generation-limit.** That is the subject of this doc.

---

## 3. The generation gap, decomposed (labeled — mandate ①)

`decompose_generation_gap` (`h7_basal_sf_force_budget.py`) factors the raw 2014× gap so each
multiplier is attributable, not opaque. **All factors are evaluated at the UPPER engagement
φ≈0.99 and assume coherent axial summation (cos = 1) — both optimistic (§5).**

| factor | value | physical meaning / source | lever |
|---|---|---|---|
| raw single-SF tension | ≈ 5.0 pN | φ 0.99 × 10 heads × 0.5 pN (brief literal D5); NO mesoscale scaling, NO continuous_stroke | — |
| **per-minifilament fidelity** | **~11×** = ~2.8× × ~4× | brief-literal minifilament (10 heads × 0.5 pN ≈ 5 pN) → **literature dipole 56 pN** (28 heads × 2 pN; Billington 2013 × Chugh 2017). The 11× = (28/10 heads) × (2/0.5 pN per-head) — both **parameter** choices, NOT delivery (the budget uses F_stall directly, already assuming the §9 fix) | brief-vs-literature minifilament parameterisation (same category/magnitude as the §9 per-head recovery, distinct mechanism) |
| **cross-sectional NMII count** | **~180× (floor) / ~360× (centre)** | after the per-minifilament fix, how many **literature minifilaments must act coherently across ONE SF cross-section** to reach Kumar = Route-B native:effective factor | the **missing MCF7 SF-NMII density datum** (session (i)) |
| engagement realism | **~9× WORSE** | budget uses upper φ≈0.99; the fine-grained cortex realises only ~11% (geometric/availability limited, `h7_engagement_saturation.json`) | structural (availability), not a credit |

Consistency: per-minifilament fidelity × cross-sectional count == the raw factor-to-band
(11.2 × 360 ≈ 4028 to centre; 11.2 × 180 ≈ 2014 to floor); locked by
`test_sf_force_budget_decomposition.py`.

**Reading.** Of the 2014×, only ~11× is a brief-vs-literature minifilament *parameter* choice
(2.8× heads-per-minifilament + 4× per-head stall — the budget already assumes ideal §9-style
F_stall delivery, so this is NOT the delivery/stiffness gap); the dominant ~180× is a
*count/density* quantity with **no MCF7 datum**.
A Kumar SF needs **~180–360 literature minifilaments coherent per cross-section** — an
order-of-magnitude that session (i) must verify against native ventral-SF NMII content; it is
plausibly high (a NMII stack is ~10–30 minifilaments; ~180–360 implies ~6–12 stacks bundled
per cross-section), which is precisely why the honest fallback (this doc) may stand.

---

## 4. Unification with the cortical active-γ floor

The SF tension floor and the cortical active-γ floor are **the same physics**: aggregate
motor force = (geometric motor density) × (per-minifilament dipole force), both bottoming out
on per-head stall × engaged heads, both missing the same MCF7 density anchor.

| | cortical active-γ floor | SF tension floor |
|---|---|---|
| budget form | γ_active = ½ · n₂D · f · ℓ (surface tension; force/length) | T = N_cross · f (line tension; force) |
| per-minifilament dipole f | **56 pN** (28 heads × 2 pN) | **56 pN** (literature ref); brief-literal model uses 5 pN (→ §3 11× fidelity) |
| geometric density | areal **n₂D** [µm⁻²] | cross-sectional **N_cross** [per SF] |
| target | MCF7 active γ ≈ 0.39–0.41 mN/m (Hosseini/Fischer-Friedrich 2021) | Kumar 2006 single-SF 10–30 nN |
| gap at lit/proxy density | **~40–80× under** | **~180× under** (after per-minifilament fix) |
| density needed | n₂D ≈ 16–47 µm⁻² (~27–80× the HeLa proxy 0.6 µm⁻²) | N_cross ≈ 180–360 lit-minifilaments/SF |
| **MCF7 datum** | **NONE** (Nie 2015 HeLa proxy only) | **NONE** (session (i) deep-research target) |
| architecture caveat | tension is **overlap/length-set, independent of myosin count** (Chugh 2017; Truong Quang 2021) | tension needs **sarcomeric polarity organisation**; random polarity → no rectification (§11; Hotulainen-Lappalainen 2006) |

**The one story.** Both floors are **generation-bound at the literature parameters**, both
are gated by a **missing MCF7 native-NMII density datum** (areal for the cortex, cross-sectional
for the SF), and both carry the same caveat that raising raw motor *count* is the **wrong
lever** — tension is set by **organisation** (cortical actin overlap; SF sarcomeric polarity),
which *rectifies* generated motor force into sustained tension. This is consistent with the
project's authoritative record (`CORTICAL_TENSION_RECORD_2026-06-04`: active channel "open
force-generation problem") and is a **model-fidelity / scope conclusion, not a parameter to
tune to pass a gate**.

---

## 5. Honest bounds (do not read the budget optimistically)

The §3 budget is a **ceiling**, not a prediction of the current construction:

1. **Engagement.** Uses the zero-load upper φ≈0.99. The realised engagement on the
   fine-grained model is ~11% → the realistic budget is **~9× worse** than tabulated.
2. **Coherence / rectification (decisive).** The budget assumes cos = 1 (all motor force adds
   along the cable). The §11 SF-2c experiment (same-seed paired, coherent signed-traction,
   §9-corrected motor, n_fil=24, 20 minifilaments, 68 engaged heads) measured the
   force-ON − force-OFF differential at **−61.45 ± 27.74 pN — significant and NEGATIVE
   (slackening)**: a random mixed-polarity end-anchored bundle does **not** rectify bipolar
   sliding into end-ward contraction. So the budget is realised **only** with sarcomeric
   polarity organisation; the current random-polarity construction realises ~0 (or slightly
   negative) net traction even with the per-head fix.
3. **Per-head value itself is unsettled across lines.** The SF budget uses 0.5 pN/head (brief
   literal D5); the cortical ½·n·f·ℓ uses 2 pN/head (Billington/Chugh dipole); the §9
   stiffness analysis uses F_stall ≈ 8.48 pN/head (single-myosin stall, Finer 1994 / Kaya-
   Higuchi 2010). These differ by ~4–17× and should be reconciled into one anchored per-head
   value before any deliverable cite — flagged as a **cross-line consistency item for PI**
   (it does not change the qualitative generation-bound conclusion, only the exact factor).

---

## 6. Reconciliation with session (i) — LANDED 2026-06-09 (REFUTE; conclusion DECISIVE)

Session (i) "CLOSE THE FLOOR" (`f6b956c`, `h7/sf-nmii-forcescale`,
`H7_SF_NMII_FORCESCALE_RESULT_2026-06-09.md`) ran the density deep-research (104-agent,
3-vote verify) and **REFUTE'd** — Kumar is **not** closable from a measured native-NMII
datum via Route B. The two conclusions agree; the generation-bound finding is now **decisive,
jointly HALT→PI**. (i)'s findings sharpen this doc on three axes:

1. **The molecular correction is reconciled and SMALL.** (i)'s directly-anchored
   per-minifilament stall is **fs ≈ 17 pN** (Stachowiak 2009 muscle extrapolation, ORDER,
   self-flagged not a direct NMII measurement) → a **~3.4×** correction over the model's 5 pN
   minifilament. This is the *defensible* value; my §3 "~11×" used the heads×per-head dipole
   (56 pN) as an **upper bound** (it assumes all ~28 heads at full stall simultaneously). Both
   are small vs the count factor — the §5.3 per-head spread (0.5 / 1.7 / 2 / 8.48 pN across
   brief / Stachowiak-extrap / dipole / single-myosin-stall) resolves toward (i)'s 17 pN as the
   per-minifilament anchor, with the others being per-*head* values on different bases.

2. **The cross-sectional count is GEOMETRICALLY IMPOSSIBLE — the gap is irreducible.** The only
   literature parallel-count estimate (~50/cross-section) was **REFUTED 0–3**; back-solving the
   band (10–30 nN ÷ 17 pN) needs **~590–1760** minifilaments, but a 50–250 nm-radius SF
   cross-section with ~10–30 actin filaments holds **O(5–15)** — the band exceeds geometry by
   ~40–350×. So the count datum the §3 requirement names **does not exist and cannot fit**. This
   strengthens this doc's earlier "plausibly high" to **physically impossible**: Route B cannot
   close Kumar, and **the cortical γ-floor cannot get its density datum from the SF side either**
   (the SF count datum is itself missing).

3. **The gate target itself is mostly NOT single-fiber active generation** (the decisive
   reframe, HIGH confidence; Kassianidou, Brand, Schwarz & **Kumar** 2017, PNAS 114:2622, U2OS
   active-cable-network analysis). Single-fiber aggregate motor stall ≈ **6 nN**; a connecting SF
   adds only ~5 nN of **active** myosin force; the ~25 nN center / 10–30 nN band is largely
   **network connectivity / prestress**. In our platform that prestress is supplied by the
   **connected basal mesh + FA anchoring + whole-cell turgor**, *not* `sf_myosin_` generation
   alone. Against the ~6 nN single-fiber **active** target, even the refuted-generous 50 × 17 pN
   ≈ 0.85 nN is **~7× under** — so the active component is consistent with the model's order, and
   the band is not a pure-active target.

**Net.** Both sessions converge: the SF active-generation channel is generation-limited at the
literature parameters, the missing/impossible cross-section count makes the gap irreducible, and
the Kumar band bundles network/prestress the active motors do not (and need not) supply. The
honest verdict is **REFUTE on Kumar-as-pure-active-generation + HALT→PI on the gate scope** (§7).

> **Integration note (for the Lead).** Session (i) also extended
> `scripts/h7_basal_sf_force_budget.py` on its own branch (molecular correction + missing-count
> accounting + active/network split). This branch (`h7/sf-generation-floor`) extended the *same*
> file (`decompose_generation_gap`). The two are complementary but will need a manual merge at
> integration — they touch the same module on different branches (1-session-1-worktree, as
> designed).

---

## 7. Gate-reframe options — SURFACED to PI, NOT applied (mandate ③; no gate-loosening)

The Kumar 10–30 nN gate and the cortical active-γ band are validation contracts; this session
**does not edit them**. These are options for PI to consider (parallel to the cortical PI items
in `H7_ACTIVE_GAMMA_SYNTHESIS` §"PI decision items"):

1. **Active-fraction / network-vs-single-fiber scope (now literature-anchored by (i)).** Kumar
   10–30 nN is the *total* single-SF tension, and session (i)'s Kassianidou/Kumar 2017 reframe
   shows it is **mostly network/prestress**: single-fiber **active** myosin force is ≈**5–6 nN**,
   the remainder is network connectivity / prestress (and the SF-2c run independently shows a
   large passive pre-tension, force-OFF coherent traction ≈ +3476 pN from EV/bundling). In our
   platform the network/prestress component is supplied by the **connected basal mesh + FA
   anchoring + turgor**, not `sf_myosin_` generation — so the gate for the *active* channel
   should be the **~6 nN single-fiber active target**, not the 10–30 nN composite. (Even against
   ~6 nN the active component is ~7× under, so this re-scope does not "rescue" the gate — it
   makes the comparison *correct*, which is the point.) Mirrors cortical item #2. **Recommended
   primary reframe.**

2. **Tool-role scope.** Is the absolute Kumar band the right gate for the *fine-grained*
   single-cell tool, or is the fine-grained tool's role to **supply** the active-stress
   coefficient (per-minifilament dipole f, density-scaling ∝ ρ_NMII) to the coarser layers,
   with the absolute band validated at the layer that owns the density datum? (Mirrors cortical
   item #3.) This is the natural home for the "generation-bound, density-datum-gated" result.

3. **Density-conditioned band.** State the gate as conditional on the (currently unknown) MCF7
   SF-NMII cross-sectional count: until that datum exists, the SF tension is a **provisional
   prediction**, not a pass/fail — the same status the cortical density-datum carries.

None of these is a gate edit; any contract change is a PI sign-off item.

---

## 8. What this is NOT

- **Not** an N_filaments / k_actin problem — those set strain, not the tension ceiling (the
  whole point of the quasi-static budget; locked N-independent).
- **Not** a code bug in placement (B4 gate PASS), construction (force-free, gated), or the
  budget arithmetic (decomposition consistency-tested).
- **Not** a license to raise density, per-head stall, or engagement to pass the gate.
- **Not** the final SF word: the **sarcomeric-polarity construction** (graded +ends→FA /
  −ends→centre + periodic α-actinin; §11 NEXT) is the PI-gateable build that would test whether
  traction emerges from {corrected motor} × {sarcomeric structure} — the missing organisation
  factor, distinct from the count datum (i) chases.

---

## 9. Deliverables (this session)

- `scripts/h7_basal_sf_force_budget.py` — added `decompose_generation_gap` (labeled,
  N-independent factors + cortical-γ unification string + upper-bound caveat); CLI prints the
  decomposition; JSON gains a `decomposition` block.
- `tests/test_sf_force_budget_decomposition.py` — 5 consistency gates (product recovers raw
  gap; lit minifilament == cortical 56 pN; fidelity ≈ 11×; cross-section count order-hundreds;
  engagement is a penalty not a credit).
- `docs/v2_audit/SF_GENERATION_LIMIT_SYNTHESIS_2026-06-09.md` — this doc.
- `scripts/h7_sf_generation_limit_vis.py` → `outputs/h7/figs/h7_sf_generation_limit.png`.

## Figures

- **`figs/h7_sf_generation_limit.png`** — (left) labeled waterfall of the single-SF force
  budget: raw ≈5 pN brief-literal minifilament → ×~11 per-minifilament fidelity → ×~180
  (floor) / ~360 (centre) cross-sectional NMII count, reaching the overlaid Kumar 2006
  10–30 nN band (log y-axis, noted); annotated with the coherent-bundle upper-bound /
  §11-slackening / engagement-realism caveats. (right) the cortex ∥ SF unification — same
  ½·n·f·ℓ budget, same f≈56 pN dipole, same missing-MCF7-density-datum gap, same
  organisation-not-count lever.
