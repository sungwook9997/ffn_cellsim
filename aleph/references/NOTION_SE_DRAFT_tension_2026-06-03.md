# Notion SourceEvidence DRAFT — cortical/membrane/tissue-tension batch (2026-06-03)

> **For the Lead to create as Notion SourceEvidence rows at closeout.** One block per paper.
> Source: 11 PI-supplied tension PDFs (dropped 2026-06-03), triaged in
> `aleph/docs/CORTICAL_TENSION_TRIAGE_2026-06-03.md` (§6 bibliographic, §1–§4 classification).
> All 11 appended to `aleph/references/tag_corpus.json` (`source:"new"`,
> `linked_to_source_evidence:false`, `se_uid:null`, `n_chunks:0` — Lead reruns refresh to embed).
>
> ⚠️ **INTEGRITY FLAG (hard).** De Belly, Keren, and Lüchtefeld measure **MEMBRANE tension** — a
> DIFFERENT physical quantity from actomyosin **CORTICAL** tension. **None gives an absolute pN/µm
> cortical value.** No KB entry may imply they anchor the KU-3.5 cortical band. Only **Chugh** (and
> the bridge papers, which USE cortical tension as a term) concern the KU-3.5 cortical quantity.

---

## 1. Chugh2017_NatureCellBiology  ⭐
- **citation_key:** Chugh2017_NatureCellBiology
- **title:** Actin cortex architecture regulates cell surface tension
- **authors:** Chugh P., Clark A.G., Smith M.B., Cassani D.A.D., Dierkes K., Ragab A., Roux P.P., Charras G., Salbreux G., Paluch E.K.
- **year:** 2017
- **venue:** Nature Cell Biology 19(6):689–697
- **DOI:** 10.1038/ncb3525
- **path:** emss-72183.pdf
- **relevance:** KU-3.5 magnitude oracle — mechanistic cortex peaks at ~0.37 mN/m = band floor; T₀=230 pN/µm; filament-length optimum ~400–500 nm; tension = stall-force × connectivity, not myosin count.
- **classification:** acceptance oracle (magnitude) + architecture constant + mechanism corroboration
- **⚠️ membrane-vs-cortical:** **CORTICAL** (actomyosin) — this IS the KU-3.5 quantity.

## 2. DeBelly2023_Cell
- **citation_key:** DeBelly2023_Cell
- **title:** Cell protrusions and contractions generate long-range membrane tension propagation
- **authors:** De Belly H. et al.
- **year:** 2023
- **venue:** Cell 186(14):3049–3061
- **DOI:** 10.1016/j.cell.2023.05.014
- **path:** 1-s2.0-S0092867423005330-main.pdf
- **relevance:** Long-range tension transmission near-unattenuated (inter-point delay 1.2±1.2 s) iff force engages a continuous cortex. STAGE-2 already ruled transmission OUT as our floor → confirmatory. Lasting value: membrane↔cortex coupling oracle (3-tier model) for future H.9 membrane layer.
- **classification:** confirmatory (transmission not the floor) + H.9 membrane-coupling oracle (deferred)
- **⚠️ membrane-vs-cortical:** **MEMBRANE** tension — different quantity; **no absolute pN/µm cortical value**; MUST NOT anchor KU-3.5.

## 3. Keren2023_Cell
- **citation_key:** Keren2023_Cell
- **title:** Membrane tension goes the distance (preview/commentary on De Belly 2023)
- **authors:** Keren K.
- **year:** 2023
- **venue:** Cell 186(14):2956–2958
- **DOI:** 10.1016/j.cell.2023.05.033
- **path:** 1-s2.0-S0092867423005858-main.pdf
- **relevance:** 3-page Cell preview of De Belly 2023; same membrane-tension-transmission story, no new data.
- **classification:** preview/commentary (confirmatory context)
- **⚠️ membrane-vs-cortical:** **MEMBRANE** tension — different quantity; **no cortical value**; MUST NOT anchor KU-3.5.

## 4. Luchtefeld2024_NatureMethods
- **citation_key:** Luchtefeld2024_NatureMethods
- **title:** Dynamic monitoring of membrane tension changes (membrane, not cortical)
- **authors:** Lüchtefeld I. et al.
- **year:** 2024
- **venue:** Nature Methods 21:1063–1073
- **DOI:** 10.1038/s41592-024-02277-8
- **path:** s41592-024-02277-8.pdf
- **relevance:** O(0.1–0.5 mN/m) membrane-tension CHANGES on HFF; cytoskeleton CONFINES membrane tension (confinement length ~0.7–1.7 µm). Membrane-tension perturbation-scale anchor only.
- **classification:** membrane-tension perturbation-scale anchor (EXTEND/H.9 context)
- **⚠️ membrane-vs-cortical:** **MEMBRANE** tension — different quantity; the O(0.1–0.5 mN/m) is a membrane CHANGE, NOT a cortical band; MUST NOT conflate with / anchor KU-3.5.

## 5. NishitaniMiura2025_arXiv
- **citation_key:** NishitaniMiura2025_arXiv
- **title:** Quantitative analysis of cell membrane tension via the curvature-velocity law V=σκ
- **authors:** Nishitani W., Miura T.
- **year:** 2025
- **venue:** arXiv:2504.14887
- **DOI:** 10.48550/arXiv.2504.14887
- **path:** Quantitative_Analysis_of_Cell_Membrane_Tension_in_.pdf
- **relevance:** Effective surface tension from a shape time-series via V=σκ (+ volume term), regressed on retraction regions only (MDCK, 2D). Mirrorable now from our GSD trajectories. ⚠️ extracted σ has units m²/s (mobility×tension), NOT mN/m → RELATIVE/comparative observable, not an absolute band anchor; take absolute magnitude from Chugh.
- **classification:** observable protocol (relative σ) — NOT a magnitude anchor
- **⚠️ membrane-vs-cortical:** titled "membrane tension" but is a shape-derived effective surface tension (m²/s, relative); neither a cortical nor a membrane absolute anchor — observable only.

## 6. Fastabend2026_PhysicalReviewE  ⭐
- **citation_key:** Fastabend2026_PhysicalReviewE
- **title:** Cortical tension links curvature to tissue growth in the cellular Potts model
- **authors:** Fastabend et al.
- **year:** 2026
- **venue:** Physical Review E 113:024403
- **DOI:** 10.1103/z152-x4l1
- **path:** z152-x4l1.pdf
- **relevance:** The published mechanism for L2 `A/A₀ = a + b/R + c/R²`: cortical tension = CPM perimeter term `λ_P(P−P₀)²`; 2D Young-Laplace `R=λ/σ`; growth ∝ local interfacial curvature (∝1/R), ΔP as proliferation signal.
- **classification:** acceptance oracle + bridge-relation (single-cell γ → curvature → growth → A/A₀(R))
- **⚠️ membrane-vs-cortical:** **CORTICAL** — uses actomyosin cortical tension as the CPM perimeter term (not membrane tension).

## 7. Thiticharoentam2026_bioRxiv
- **citation_key:** Thiticharoentam2026_bioRxiv
- **title:** 3D vertex model of interfacial-tension modulation (Γ = cortical contractility − adhesion)
- **authors:** Thiticharoentam, Fukamachi, Horiguchi, Okuda
- **year:** 2026
- **venue:** bioRxiv 2026.03.17.712503
- **DOI:** 10.1101/2026.03.17.712503
- **path:** 2026.03.17.712503v1.full.pdf
- **relevance:** `Γ_interface = cortical contractility − adhesion` (verbatim) → our γ minus E-cadherin catch-bond cohesion; monolayer→3D-cap transition once free-surface (cortical) tension > ~0.2 κ₀ (underwrites L2.6); virial tissue-stress coarse-grain (observable); differential-Γ sorting (Steinberg DAH).
- **classification:** acceptance oracle + constant-anchor (Γ=cortical−adhesion; 0.2 κ₀ cap threshold) + observable (virial stress)
- **⚠️ membrane-vs-cortical:** **CORTICAL** — Γ built from actomyosin cortical contractility minus adhesion (not membrane tension).

## 8. Roffay2021_Development  ⭐
- **citation_key:** Roffay2021_Development
- **title:** Inferring cell junction tension and pressure from cell geometry
- **authors:** Roffay C., Chan C.J., Guirao B., Hiiragi T., Graner F.
- **year:** 2021
- **venue:** Development 148:dev192773
- **DOI:** 10.1242/dev.192773
- **path:** DownloadCombinedArticleAndSupplmentPdf.pdf (article + supplement combined)
- **relevance:** Vertex force balance `Σtᵢ=0`; Young-Laplace 3D `ΔP = t·K`, `K=1/R+1/R'` (pressure↔tension↔curvature for a real 3D cap); outer (cell-medium) tension ≈ 1.6–2× interior cell-cell tension (pipette-validated) ⇒ aggregate surface is cortex-derived. ⚠️ mouse embryo — ratio transfers, absolute a.u. do not.
- **classification:** acceptance oracle (Young-Laplace 3D) + observable (vertex-tension inference) + constant-anchor (outer/inner ratio)
- **⚠️ membrane-vs-cortical:** **CORTICAL/junctional** (actomyosin cortex sets junction tension), not membrane tension.

## 9. Runser2023_bioRxiv  (peripheral)
- **citation_key:** Runser2023_bioRxiv
- **title:** 3D simulation of tissue mechanics with cell polarization
- **authors:** Runser S., Vetter R., Iber D.
- **year:** 2023
- **venue:** bioRxiv 2023.03.28.534574 (ETH Library 10.3929/ethz-b-000653764)
- **DOI:** 10.3929/ethz-b-000653764
- **path:** 2023.03.28.534574v1.full.pdf
- **relevance:** 3D tissue-mechanics simulation with cell polarization; tissue-sim context, peripheral to γ/cohesion. Scan at incorporation. NOTE: a duplicate copy `2023.03.28.534574v1.full (1).pdf` exists in references/ — register ONE.
- **classification:** peripheral (tissue-scale sim context)
- **⚠️ membrane-vs-cortical:** N/A — tissue-mechanics sim, not a single-cell tension measurement.

## 10. Tao2020_BiophysicalJournal  (peripheral)
- **citation_key:** Tao2020_BiophysicalJournal
- **title:** Tuning cell motility via cell tension with a mechanochemical cell migration model
- **authors:** Tao K., Wang J., Kuang X., Wang W., Liu F., Zhang L.
- **year:** 2020
- **venue:** Biophysical Journal (published version of bioRxiv 847046)
- **DOI:** 10.1016/j.bpj.2020.07.030
- **path:** main (1).pdf  (= published version)
- **relevance:** Mechanochemical migration model; cell tension tunes motility. Peripheral to γ/cohesion (migration-focused). ⚠️ DUPLICATE: `847046v1.full.pdf` is the bioRxiv PREPRINT of the SAME paper (registered the published `main (1).pdf` as canonical; do not create a 2nd SE row for the preprint).
- **classification:** peripheral (cell-migration model)
- **⚠️ membrane-vs-cortical:** "cell tension" is a lumped model parameter (not resolved membrane vs cortical); not a band anchor.

## 11. DeVecchis2021_BiophysicalJournal_mmc4  (supplement)
- **citation_key:** DeVecchis2021_BiophysicalJournal_mmc4
- **title:** Molecular dynamics simulations of Piezo1 channel opening by increases in membrane tension (supplement mmc4)
- **authors:** De Vecchis D., Beech D.J., Kalli A.C.
- **year:** 2021
- **venue:** Biophysical Journal (supplemental file mmc4 of the article)
- **DOI:** 10.1016/j.bpj.2021.02.006
- **path:** mmc4.pdf
- **relevance:** Piezo1 membrane-tension MD; supplement-level only (this PDF is the SUPPLEMENTARY material `mmc4`, not the main article). Membrane MD, peripheral.
- **classification:** peripheral / supplement-level (membrane MD)
- **⚠️ membrane-vs-cortical:** **MEMBRANE** tension (MD of a mechanosensitive channel); not a cortical anchor; supplement, not a primary tension source.

---

## Dedup / supplement notes (for the Lead)
- **`main (1).pdf` (#10, published Tao 2020) duplicates `847046v1.full.pdf` (bioRxiv preprint).** Same paper, two versions. Registered the published version in tag_corpus; create only ONE SE row.
- **`2023.03.28.534574v1.full (1).pdf` duplicates `2023.03.28.534574v1.full.pdf` (#9, Runser).** Register ONE.
- **`mmc4.pdf` (#11) is a SUPPLEMENT** (Piezo1 MD supplementary material), not a standalone primary paper.

## INTEGRITY summary (carry into every KB entry)
De Belly (#2, S0092867423005330), Keren (#3, S0092867423005858), and Lüchtefeld (#4,
s41592-024-02277-8) measure **MEMBRANE tension** — a DIFFERENT quantity from actomyosin CORTICAL
tension. None gives an absolute pN/µm cortical value. **No KB entry may imply they anchor the KU-3.5
cortical band.** The only KU-3.5 cortical magnitude anchor in this batch is **Chugh 2017** (#1).
