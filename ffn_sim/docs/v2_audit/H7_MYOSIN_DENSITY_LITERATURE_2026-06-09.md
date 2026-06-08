# H.7 active-γ — myosin-density literature reconciliation (deep-research, 2026-06-09)

**Branch** `h7/full-cell-integration`. PI authorized acting on the recommendation + continuing.
Resolves the 3 PI decision items from `H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09.md` against the
primary literature via a verified deep-research sweep (101 agents, 19 primary sources fetched,
25 claims adversarially 3-vote-verified, 22 confirmed / 3 killed). This is a DATUM/literature
finding, not a model change.

## The question
The active-γ floor is generation/density-bound: the active-gel envelope at the model's params
is ~18–36× under band. Item #1 asked — is the model's myosin density (0.6/µm², HeLa proxy) too
low, i.e. is there a defensible MCF7/cortical density ~16–21/µm² that would close the gap?

## Verified findings

1. **There is NO higher cortical minifilament density in the literature.** The *only* directly-
   quantified cortical NMII minifilament areal density is **Nie et al. 2015** (Cytoskeleton,
   PMC4361371): 40±20 MRLC-GFP foci / 64 µm² = **~0.63/µm², interphase adherent HeLa medial
   cortex** — *exactly the model's value*. The widely-cited super-resolution NMII studies
   (Fenix 2016 PMC4850034, "Beach 2014") quantify the **leading edge / cleavage furrow**, not
   the medial cortex, and report filament-class % + stack lengths, **not** a calibrated
   minifilaments/µm². ⇒ **The model's 0.6/µm² is the only available cortical datum; there is no
   measured value to "correct up" to.** (Confidence: high, 3-0.)

2. **Nie 2015 is a diffraction-limited confocal intensity ESTIMATE** (one MRLC-GFP focus is not
   established to equal one bipolar minifilament; super-res shows minifilaments stack), so it
   **likely undercounts** — but no paper quantifies by how much. This is a soft lower-bound
   argument, not a hard datum for a higher number. (3-0.)

3. **MCF7-specific cortical ACTIVE tension IS measured** — a genuine new MCF7 datum, better than
   the suspended-total Hosseini 2020 0.27: **Hosseini, Frenzel & Fischer-Friedrich 2021** (Biophys
   J 120(16):3516, DOI 10.1016/j.bpj.2021.05.006, PMID 34022239, PMC8391033), AFM parallel-plate
   confinement (γ = F/[A_con(1/R1+1/R2)] — the same observable the model computes):
   **control MCF-7 γ_act = 0.39–0.41 mN/m interphase, 0.72–0.78 mN/m mitosis** (Tables 1–2; mitosis
   ≈1.85–1.9× interphase). γ_act is the *active/steady-state* component (distinct from the area
   modulus K_h and viscosity in the same rows). (3-0; the 0.39 datum voted 2-1 but both numbers
   verified verbatim.) ⇒ **the active target is MCF7-specific ~0.4 mN/m**, not 0.27.

4. **The active component is ~half of total tension** (resolves item #2): blebbistatin drops
   cortical tension >50% (Tinevez 2009 PNAS, PMC2765453); L929 total 0.41 mN/m → active ~0.2;
   interphase HeLa total 0.17±0.13 (FF2014). So the active tension is the SAME ORDER as total —
   the model's 0.005 mN/m is **~40–80× too low** against the 0.4 MCF7 active target. (3-0.)

5. **The active-gel relation has a factor-2 convention spread in the primary literature**
   (resolves my loop-6 vs re-correction flip-flop — BOTH are published):
   - **Tinevez 2009** (the canonical Salbreux/Joanny/Paluch bleb paper): *"cortical tension is
     equal to ζΔμh/2"* — i.e. γ = ½·σ_active·h **(½ included).**
   - **Fischer-Friedrich 2014** (PMC4148660): σ = γ/h (≈8 kPa metaphase at h≈200 nm) **(no ½).**
   ⇒ the envelope is honestly a **range 0.005–0.010 mN/m** (½ → 1 convention), NOT a single value.
   The framework/formula is correct (3-0); the open issue is purely the inputs + this convention.

6. **f_minifil (~56 pN) and ℓ (~301 nm) are well-anchored and NOT the gap source.** Billington
   2013 (PMID 24072716): NM IIA = 29 molecules/filament (58 heads), contour 301±24 nm; Melli
   2018 (eLife, PMC5829915): ~30 molecules, ~300 nm. Per-minifilament force ~50–56 pN. (3-0.)
   The config's `n_heads_per_side=28` (≈28–29 heads/side, full-stall one-sided force ≈56 pN) is
   consistent with Billington. (The deep-research "2× overcount" caveat came from imprecise
   wording in my query — "molecules/side" — not the config; the config is fine.)

7. **The resolution is ARCHITECTURE/OVERLAP, not density** (the decisive reframe). Chugh 2017
   (Nat Cell Biol, PMC5536221) + Truong Quang 2021 (Nat Commun, PMC8586027), both verified 3-0:
   - tension is *maximal at intermediate actin filament length* (non-monotonic in thickness);
   - *"cortical tension can be modulated independently of myosin motor density and activity, by
     fine-tuning actin filament structural properties"*; depleting filament-length regulators
     drops tension *as much as* myosin inhibition;
   - interphase→mitosis tension rises **>3× despite LOWER cortical myosin** (Truong Quang);
   - tension correlates with **actin-myosin OVERLAP**; incompletely-overlapping myosin generates
     stress inefficiently.

## Reconciliation — what this means for the model

- **Item #1 (density datum): RESOLVED — no fix available.** There is no measured cortical density
  higher than the model's 0.6/µm²; it is the only cortical datum and is itself a likely-undercount
  confocal estimate. **Do not change the production density** — there is nothing defensible to
  change it to. (Adopting ~16–47/µm² to close γ would be inventing a number = magic-number
  violation, and per Chugh/Truong Quang would *misattribute* architecture-dependent tension to raw
  count.)
- **Item #2 (active target): RESOLVED — update the anchor.** The MCF7-specific active target is
  **γ_act ≈ 0.39–0.41 mN/m interphase (Hosseini/FF 2021)**, not the 0.27 suspended-total. The
  active fraction is ~½ of total, and γ_act IS the active component directly. This *widens* the
  model gap to **~40–80×** (0.4 ÷ 0.005–0.010), and to reach it at the anchored f/ℓ would need
  n_2D ≈ 24–47/µm² (~40–80× the proxy) — confirming the gap is genuine and large.
- **Item #3 (mechanism/scope): RESOLVED — the gap is architecture/overlap/transmission, not count.**
  The literature is explicit that cortical tension is set substantially by actin architecture and
  actin-myosin overlap *independently of myosin number*. This **vindicates the transmission/lever
  framing** (Gate-A) and the buckling/condensation lit-study (Miyazaki/Lenz,
  [[project-actin-architecture-litstudy]]), and explains the model's floor: the mean-field
  "tension = number × force × length" is incomplete, and the model's rigid M-SHAKE backbone +
  low actin-myosin overlap don't realize the architecture-dependent stress that carries most of
  the real cortical tension. **Raising density is the wrong lever; the right lever is the
  architecture/overlap/transmission physics** — which Gate-B's M-SHAKE-relax buckling did NOT
  capture (per-filament load < Euler threshold), so the architecture lever is overlap/length/
  connectivity, a distinct and unbuilt mechanism.

## Net verdict (closes the active-γ saga)
The active-myosin cortical-γ floor is **real and large (~40–80× under the MCF7 active target
0.4 mN/m)**, and it is **NOT a density-datum error** (the model uses the only measured cortical
density). The gap is **architecture/overlap/transmission-rooted** per primary literature
(Chugh 2017, Truong Quang 2021) — i.e. the mean-field motor-count generation the model and the
½·n·f·ℓ envelope capture is genuinely insufficient because real cortical tension depends on
actin-filament-length / actin-myosin-overlap architecture that the current rigid-backbone build
does not realize. This is a model-FIDELITY / scope conclusion to surface to PI, not a parameter
to tune.

## SE-registration candidate (PI-gated, do not auto-register)
- **Hosseini2021_BiophysJ** — Hosseini, Frenzel & Fischer-Friedrich 2021, Biophys J 120(16):
  3516-3526, DOI 10.1016/j.bpj.2021.05.006, PMID 34022239, PMC8391033. MCF-7 cortical γ_act
  0.39–0.41 mN/m (interphase) / 0.72–0.78 (mitosis), AFM parallel-plate confinement. The
  MCF7-specific ACTIVE cortical-tension anchor for KU-3.5-ACTIVE. (Companion to the existing
  Hosseini2020_AdvSci total-tension datum.)
- Supporting: Chugh2017_NatCellBiol (PMC5536221), TruongQuang2021_NatCommun (PMC8586027) —
  architecture/overlap regulation of cortical tension; Tinevez2009_PNAS (PMC2765453) — γ=ζΔμh/2 +
  blebbistatin >50%; Nie2015_Cytoskeleton (PMC4361371) — the 0.63/µm² cortical density (already used).

## Provenance
Deep-research run wf_7847bc41-6d9 (transcript under subagents/workflows/). Full verified result:
`outputs/.../tasks/wny8b2o1m.output` (mirrored reasoning here). 22/25 claims confirmed 3-0/2-1;
3 killed (a duplicate-weaker-source density restatement, a thickness-tension anti-correlation
over-claim, a leading-edge motor-group claim).
