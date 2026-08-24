# SourceEvidence registration candidates — 2026-06-03

After the TAG re-index (kb.duckdb rebuilt 2026-06-03: 144 citation_keys / 6292 chunks /
150 PDFs; **94 linked to a Notion SourceEvidence by DOI, 56 NOT linked**), these papers are
**searchable in TAG (BM25) but absent from the Obsidian graph** because Obsidian
SourceEvidence nodes come only from the Notion SourceEvidence DB (SoT). To appear in the
graph + resolve citations cleanly they need a Notion SourceEvidence row (which also replaces
the auto-slugified `citation_key` with a proper Author-Year key).

> ⚠️ Per the citation-integrity hard rule, **verify each before creating its SE row** (the
> prior audit found 3 hallucinated sources). SourceEvidence is edited in Notion, not here.

## A. Newly ingested this session (19) — the actionable new set

### Directly project-relevant — register first
| TAG citation_key (auto) | Real paper | Relevance |
|---|---|---|
| `molecular-determinants-of-cadherin-ideal-dc7506` | **Manibog et al. 2016**, *Molecular determinants of cadherin ideal bond formation* | **L2.5** cadherin catch-bond |
| `a-structure-based-sliding-rebinding-mech-a31a18` | **Lou & Zhu 2007** (Biophys J), *A structure-based sliding-rebinding mechanism for catch bonds* | **L2.5** catch-bond mechanism |
| `supporting-information-c2d14d` | **Rakshit et al. 2012 SI** (*Ideal, catch and slip bonds in cadherin adhesion*) | SI of already-indexed `Rakshit2012_PNAS` — link to that SE row, don't double-register |
| `Stephens2017_MBoC` *(clean key, already SE-linked ✓)* | Stephens et al. 2017, *Chromatin and lamin A determine two mechanical regimes of the nucleus* | H.9 nucleus — **already linked, no action** |
| `Murrell2015_NRMCB` *(clean key, already SE-linked ✓)* | Murrell et al. 2015 NRMCB, actomyosin contractility review | **already linked, no action** |

### Actomyosin-network-mechanics cluster (H.3 cortex / semiflexible networks) — register as a batch
| TAG citation_key (auto) | Real paper |
|---|---|
| `actin-cortex-architecture-regulates-cell-3ddcda` | Chugh et al. 2017 Nat Cell Biol, *Actin cortex architecture regulates cell surface tension* |
| `f-actin-buckling-coordinates-contractili-a399a1` | Murrell & Gardel 2012 PNAS, *F-actin buckling coordinates contractility and severing* |
| `architecture-and-connectivity-govern-act-af52c5` | Ennomani et al. 2016 Curr Biol, *Architecture and connectivity govern actin network contractility* |
| `contractile-units-in-disordered-actomyos-2c287c` | Lenz et al., *Contractile units in disordered actomyosin bundles arise from F-actin buckling* |
| `protein-friction-and-filament-bending-fa-71b702` | *Protein friction and filament bending facilitate contraction of disordered actomyosin networks* |
| `filament-rigidity-and-connectivity-tune--d0fab4` | *Filament rigidity and connectivity tune the deformation modes of active biopolymer networks* |
| `a-theory-that-predicts-behaviors-of-diso-baeb74` | Belmonte et al., *A theory that predicts behaviors of disordered cytoskeletal networks* |
| `a-versatile-framework-for-simulating-the-b443c7` | *A versatile framework for simulating the dynamic mechanical structure of cytoskeletal networks* (Cytosim-class) |
| `actomyosin-pulsation-and-ows-in-an-activ-be51eb` | *Actomyosin pulsation and flows in an active elastomer with turnover and network remodeling* |
| `received-8-feb-2016-accepted-13-jul-2016-468ea9` | Nat Commun 2016, *Disordered actomyosin…* (key-extraction failed — fix on registration) |
| `research-article-89ab03` | *Filament turnover tunes both force generation and dissipation to control…* |
| `research-article-b69991` | **MEDYAN** — *Mechanochemical simulations of contraction and polarity alignment…* |
| `royalsocietypublishing-org-journal-rsif-b67b03` | J. R. Soc. Interface (Juel Pørtner, Mularski et al.) |

### Off-topic — should NOT be in references/ (remove rather than register)
| TAG citation_key | Paper |
|---|---|
| `the-in-uence-of-interdomain-interactions-d32903` | *Influence of interdomain interactions on intradomain motions in yeast phosphoglycerate kinase* — unrelated protein-dynamics paper, noise |

## B. Older unlinked backlog (~37 more)

The other ~37 of the 56 SE-unlinked papers pre-date this session (membrane/lipid, AUC,
protein-dynamics, and some genuinely-relevant: `moving-cell-boundaries-drive-nuclear-sha`,
`feeling-for-filaments-quanti-cation-of-t`, `motility-driven-glass-and-jamming-transi`,
`modeling-semi-exible-polymer-networks-590f71` (151 chunks), `nonlinear-elasticity-in-biological-gels`).
A chunk of these look like off-topic ingests (`spida-surveys`, `synaptobrevin-transmembrane`,
`prion-protein-antibody`, `thermodynamics-of-micelle-formation`) that may not belong in
`references/` at all. Triage separately.

> Full live list: `python outputs/tag_kb/... ` — `select citation_key,count(*) from paper_chunks where se_uid is null or se_uid='' group by 1`.
