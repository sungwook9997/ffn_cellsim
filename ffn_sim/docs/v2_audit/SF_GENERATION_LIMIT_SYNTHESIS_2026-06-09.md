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

---

## 1. One-paragraph answer

The single-SF tension is **motor-generation-limited**, exactly like the cortical active-γ
floor — same ½·n·f·ℓ aggregate-motor-force budget, same missing MCF7 NMII density datum.
The quasi-static force budget (`h7_basal_sf_force_budget.py`, REFUTE provisional) gives a raw
mesoscale single-SF tension of **≈5.0 pN** (φ≈0.99 engaged × 10 heads × 0.5 pN brief-literal
stall), **~2014× under the Kumar 2006 10–30 nN band**. This is **N_filaments-INDEPENDENT**:
N_filaments / k_actin set the bundle's *strain* (ε = T/μ_SF) at a given tension, **not the
tension ceiling**, which is set by motor generation. The 2014× gap factors cleanly into named
physical causes (§3): ~11× per-minifilament parameter fidelity (the brief-literal 5 pN
minifilament vs the literature 56 pN dipole — the *same* per-head/stiffness correction the
cortical line already executed in §9), and ~180× (to band floor) / ~360× (to band centre)
**cross-sectional NMII count** — how many literature minifilaments must act coherently across
ONE SF cross-section, i.e. the Route-B native:effective force-scale factor = the **missing
MCF7 SF-NMII density datum** (session (i)'s target). And the budget is only an **upper bound
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

## 6. Reconciliation with session (i)

Session (i) "CLOSE THE FLOOR" is deriving the native ventral-SF NMII content (minifilaments
per cross-section / per µm + per-minifilament stall) and the Route-B mesoscale force-scale
factor. As of this writing the sibling worktree (`ffn_cellsim-nmii`, `h7/sf-nmii-forcescale`)
has not yet committed a result.

- **If (i) lands a defensible MCF7/native SF-NMII datum that gives N_cross ≈ 180–360** (the §3
  requirement), Route B closes Kumar, the SF line is unblocked, **and the cortical γ-floor
  gets its density datum from the SF side** — this doc updates to "generation-limit closed by
  density datum (i)" and the gate-reframe (§7) becomes moot for the count axis (the
  organisation/rectification caveat §5.2 would still need the sarcomeric construction).
- **If (i) finds no datum giving ~180–360**, this generation-bound conclusion stands: the
  fine-grained model, at the best available density, predicts sub-band single-SF tension —
  honestly reported, gate-reframe surfaced to PI.

Either way, whichever session lands the result updates PLATFORM_PI_QUEUE + the Dev-Logs board.

---

## 7. Gate-reframe options — SURFACED to PI, NOT applied (mandate ③; no gate-loosening)

The Kumar 10–30 nN gate and the cortical active-γ band are validation contracts; this session
**does not edit them**. These are options for PI to consider (parallel to the cortical PI items
in `H7_ACTIVE_GAMMA_SYNTHESIS` §"PI decision items"):

1. **Active-fraction scope.** Kumar 10–30 nN is the *total* single-SF tension. The SF-2c run
   shows a large **passive** pre-tension (force-OFF coherent traction ≈ +3476 pN from EV /
   bundling). The blebbistatin-sensitive (active) fraction is the true myosin-generation
   target; if a substantial part of Kumar is passive actin/EA pre-tension, the *active* target
   is below the nominal band and the generation gap shrinks. (Mirrors cortical item #2.)

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
