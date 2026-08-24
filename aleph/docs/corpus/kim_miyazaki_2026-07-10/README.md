# Kim Taeyoon + Miyazaki lab corpus — absorbed knowledge layer (2026-07-10)

Durable, page-by-page absorption of the 77-paper Kim/Miyazaki corpus (raw PDFs live under
`aleph/references/2026_07_10/`, which is gitignored; they are also in the TAG content layer
`paper_chunks`/BM25). This directory is the **derived knowledge layer** — the actual mechanistic content
mined for engine improvement — persisted in the repo so it survives session end.

## Layout
- `dossiers/pNN.md` — one complete faithful dossier per paper: bibliography + lab, model & methods, every
  governing equation transcribed, force laws, full parameter table, quantitative results (validation
  targets), per-figure notes, and a per-page coverage log proving full read.
- `lenses/pNN.{dcm,ff,tag}.md` — three independent engine lenses per paper:
  - `.dcm` — what it offers the Deformable Cell Model (surface-mesh cell) engine.
  - `.ff` — what it offers the Filament-FEM (fine-grained Cytosim-physics) engine.
  - `.tag` — validation-gap analysis vs the existing TAG gates/params (what new validation data it adds).
- `digests/digest_{dcm,ff,tag}.md` — per-engine structured summaries of all 77 papers, sorted by
  relevance score (top_mechanisms + headline per paper). The synthesis input.
- `worklist.json` — pNN ↔ file / DOI / n_pages / title map.
- `kb_validation_state_snapshot.md` — the KB gates/params/contracts as of the absorption (the TAG-lens
  reference).

## The deliverables built from this layer
- `../ENGINE_ROADMAP_KIM_MIYAZAKI_2026-07-10.md` — prioritized DCM + FF improvement roadmaps,
  new-engine judgment (verdict: no new engine; add an FF microtubule/kinesin compartment), and how the
  corpus resolves three standing open problems (γ-floor, aggregate-compaction rate/T1, crawl grid-drag).
- `../../outputs/tag_kb/SE_REGISTRATION_CANDIDATES_2026-07-10_kim-miyazaki.md` — Notion-ready
  SourceEvidence / KnowledgeClaim / Parameter / ValidationGate candidates (gates/contracts staged for
  PI ratification), conflicts to surface, and an SI-fetch backlog.

## Provenance / integrity
- 308 subagent passes (77 dossier + 231 lens); the 15 tail lenses that hit the account session-limit
  were reconstructed by the Lead from their full dossiers (marked in-file). Verified 77/77 complete.
- Experimental values are overlay/cross-check, never fit targets. Simulation-origin constants are
  flagged, not registered. Force constants behind SI tables (Mulla T1/T2, Fan T5, Nam S1, Matsuda S1)
  are on an SI-fetch backlog before KB registration.
- Corpus: 51 Kim-lab, 24 Miyazaki-lab (labels via first/corresponding-author affiliation).
