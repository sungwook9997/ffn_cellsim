# CFD-Transport Research Program — planning dossier (2026-07-16)

Planning documents produced ahead of the **CFD-based engine rebuild** (fluid-first Biot poroelastic
substrate + transported G-actin monomer field + emergent actomyosin network). Simulation is on hold
(`ff/mech-hierarchy`); this package captures the literature, hypotheses, and experimental plan that
motivate and will validate the new engine.

**I0-A PI contract (2026-07-16):** the new engine executes only Warp-CUDA GPU kernels—HOOMD is never run.
Production NMII is explicit backbone+individual heads with Hill FV. First baseline:
MCF7×collagen×α2β1, 70,686 cortical F-actin plus separately density×geometry-derived non-cortical populations.

## ★ Manuscript

**[manuscript/MANUSCRIPT.pdf](manuscript/MANUSCRIPT.pdf)** (25 pp, Korean, bilingual title) — the full
research-program manuscript: *"세포골격의 필드 작동과 세포내 수송 / Field Actuation and Intracellular
Transport of the Cytoskeleton"*. 11 sections + Abstract + 13 original figures (schematics + physics
plots, `manuscript/make_figures.py`) + novelty + honest limits + 66 DOI-verified references + appendix.
Source: [manuscript/MANUSCRIPT.md](manuscript/MANUSCRIPT.md). Thesis: external E/B fields cannot drive
the cell interior (Re≈10⁻¹⁰ turbulence forbidden · membrane shielding 10³–10⁴× · reaction-limited SF
growth); the modelable survivor is motor-driven reaction-advection-diffusion of a mass-conserving
G-actin field. Regenerate figures: `python manuscript/make_figures.py`.

## Contents

| File | What |
|---|---|
| [**ENGINE_ARCHITECTURE_PLAN.md**](ENGINE_ARCHITECTURE_PLAN.md) | 🆕 The new-engine build plan: three-layer **living cortex** — fine-grained actomyosin (discrete ③) ⊗ Biot poroelastic cytosol (continuum fluid ①) ⊗ mass-conserving G-actin monomer field (continuum solute ②), two-way coupled in one coordinate frame using a conservative field grid plus a separate neighbour HashGrid. Frozen-cortex motivation, the v1-guardrail (continuum only for genuine continua, filaments stay discrete), governing equations per layer, coupling map, P0–P4 phasing with gates, risks, KB gaps. Ties together the `v2_audit` fluid/dynamic-fiber roadmaps and promotes the solute-transport layer to first-class. |
| [LITERATURE_DOSSIER.md](LITERATURE_DOSSIER.md) | 66 papers across 15 topics — actin retrograde flow + molecular clutch, stress-fiber assembly, formin/polymerization kinetics, myosin-II, cytoplasm rheology + transport-limited actin, low-Reynolds physics, cytoplasmic streaming/advection, galvanotaxis, E/B field effects on filaments, membrane electrical shielding, magnetic actuation, amoeboid/bleb flow. **All DOIs machine-verified** (crossref/PubMed). 🆕 = newly KB-registered. |
| [HYPOTHESES.md](HYPOTHESES.md) | H1–H12, each with rationale, quantitative prediction, null, and falsification route. |
| [WETLAB_EXPERIMENTAL_PLAN.md](WETLAB_EXPERIMENTAL_PLAN.md) · [**PDF**](WETLAB_EXPERIMENTAL_PLAN.pdf) | 7 bench aims (galvanotaxis rig, magnetic actuation, monomer-transport FRAP, SF-at-FA assembly, retrograde-flow↔clutch stiffness sweep, in-vitro electro-orientation, **TFM cross-cutting traction readout**) with readouts, controls, a go/no-go decision matrix, timeline, and risk register. TFM (Aim 7) supplies the signaling-vs-direct-force discriminator via traction on/off kinetics. |

## KB registration (Notion SoT — done this session)

Registered to the Contract-Graph via the Notion API (`scratchpad/register_kb.py`, 2-pass with Evidence relations):
- **52 SourceEvidence** rows (the papers NOT already in the KB) — Source Type mapped (Direct measurement 41 / Review consensus 10 / Model-derived 1); all DOI-verified.
- **13 KnowledgeClaim** rows (`KB-DRAFT-7-01 … 7-13`, **status=draft**) capturing the quantitative session findings, each linked (Evidence relation) to its backing SourceEvidence rows.

⚠️ KnowledgeClaim `KB-DRAFT-*` ids are **DRAFT** — the PI assigns the final `KB-x.y` ids and promotes `draft → verified/PI-ratified`. No ValidationGate / ModelContract rows were auto-created (PI-authored per CLAUDE.md). Refresh the read layers with `outputs/tag_kb/refresh.sh` + `outputs/obsidian_rag_full/refresh.sh`.

## The one-line physics thread

Across the session's field-actuation ideas the same wall recurs — **"drive the cell interior directly with an
E/B field" fails on low-Re + membrane-shielding + reaction-limited kinetics** — and the physically-correct,
modelable survivor is **motor-driven flow and conservative transport of a mass-conserving G-actin field**
in the new CFD engine. See `~/.claude/.../memory/project-field-actuation-and-actin-gf-gap.md` and the
`H1–H12` map.
