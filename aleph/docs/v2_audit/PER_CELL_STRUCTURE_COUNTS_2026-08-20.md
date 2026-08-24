# Decision 1 — the five per-cell structure counts. Two are derivable, three are not, and one KB row does not hold.

**Status:** SEARCH RESULT. PI-requested 2026-08-20. Selects no value. Everything below is either read
out of this repository's own KB or retrieved from **PubMed** during this audit, with DOIs.

**Why it blocks:** session A2 built five arena populations — microtubule, intermediate filament,
filopodium, lamellipodium, stress fibre — and **none of the five has a sourced per-cell structure
count**, so all five are `explicit` counts carrying `source_class=PI_GAP`. Half the node budget rests on
them, and the rendered PHASE 1 cell looks sparse for exactly this reason: the counts are placeholders.

---

## 1. What was searched, and what it returned

| surface | result |
|---|---|
| KB `parameter` table | **0** rows for a per-cell count of any of the five |
| KB `knowledge_claim` | densities, stiffnesses and assembly routes — **0** per-cell counts |
| KB `paper_chunks`, **12,810 chunks**, regex over the whole corpus | MT count **4** hits, all unrelated (plant cortical arrays; bundle length). Filopodia **1** hit, qualitative. Stress fibre, IF, lamellipodial density: **0** |
| PubMed | see §3 |

**Session A2's report is confirmed by an exhaustive search of the corpus**, not by agreement.

---

## 2. Two ARE derivable — the KB already carries a verified density

### Lamellipodium — `KB-3.7`, status `verified`

> *"Width 5-20 µm; thickness ~200 nm (150-200 nm sheet); **fil density ~100 per µm leading edge**"*

A count follows from the leading-edge width, which is a geometry choice rather than a missing datum:

| leading edge | filaments |
|---:|---:|
| 5 µm | 500 |
| 10 µm | 1,000 |
| 20 µm | 2,000 |
| *placeholder in the build* | *200* |

So the standing cell is **2.5–10x under** on this population, and closing it needs a width, not a count.

### Filopodium — `KB-3.8`, status `verified` … but see §4

> *"width 100-300 nm; length 1-10 um; **10-30 fil/bundle**; **r_filo~5 um**"*

If `r_filo` is the SPACING between filopodia, the count follows from covering the cell surface:

| spacing | filopodia on 4πR² = 707 µm² |
|---:|---:|
| 3 µm | 25 |
| 5 µm | 9 |
| 8 µm | 4 |
| *placeholder in the build* | *50* |

⚠ **That reading is not established — see §4.** The bundle count 10–30 is separately sourced and the
build's 20 sits inside it.

---

## 3. Three are genuinely absent

**Microtubule.** No per-cell count anywhere in the corpus. The only figure in this repository is a
docstring in `components/solid/microtubule.py` — *"~250–600 MTs/epithelial cell (order-of-magnitude)"* —
whose source is itself. The shipped `n_mt=40` is 6–15x under even that.

⚠ **And a proxy may be the wrong instrument here.** According to PubMed, Lingle, Lutz, Ingle, Maihle &
Salisbury (1998), *"Centrosome hypertrophy in human breast tumors"*, PNAS 95(6):2950-5 —
[10.1073/pnas.95.6.2950](https://doi.org/10.1073/pnas.95.6.2950) — report that breast **adenocarcinoma**
centrosomes show *"increased microtubule nucleating capacity in comparison to centrosomes of normal
breast epithelial and stromal tissues"*. MCF7 is a breast adenocarcinoma line, so an epithelial proxy is
likely a floor rather than a centre. Using it is a decision, not a default.

**Intermediate filament.** The KB has persistence lengths and the keratin/vimentin identity split
(`KB-3.26`, `KB-DRAFT-3.B-13`: MCF7 keratin `l_p` 0.3–0.5 µm) but **no count and no density**.

**Stress fibre.** Both KB claims (`KB-DRAFT-7-07`, `-08`) are about mechanics and assembly routes.
**No count.**

---

## 4. ⚠ A KB row does not hold up, and it is the one §2 leans on

`KB-3.8` is `verified` with `confidence: High` and cites **Mogilner2005_BiophysJ; Pronk2008_PNAS;
MattilaLappalainen2008_NRMCB**. Of those three, **only Mattila & Lappalainen has any text in the
corpus** — the other two have zero chunks.

Searching all 41 Mattila chunks: **no `r_filo`, no filopodia spacing, no bundle filament count, no
width.** What it does support is the length scale — *"fibroblast filopodia and nerve growth-cone
filopodia rarely exceed 10 µm in length"* — and, for count, only *"several filopodia per cell explore
the environment"*, in macrophages, qualitatively.

**So `r_filo ~5 µm` has no retrievable support in this repository.** It may be correct and sourced to
Mogilner 2005, which is simply not ingested. Until that PDF is in the corpus, §2's filopodium
derivation rests on an unverifiable row, and `KB-3.8`'s `verified` status is broader than its evidence.

⚠ Recorded as a KB defect, not a physics finding: the row's OTHER entries may be fine. What is wrong is
that a `verified/High` row cannot be checked against its own citations.

---

## 5. What the PI is being asked

| # | population | the decision | why it is that shape |
|---|---|---|---|
| 1 | lamellipodium | **a leading-edge width** (5–20 µm) | count then follows from `KB-3.7`'s verified 100/µm |
| 2 | filopodium | **ingest Mogilner 2005**, or declare a spacing as a test axis | `r_filo` is unverifiable today (§4) |
| 3 | microtubule | **declare 250–600 as a test axis**, and decide whether a normal-epithelial proxy may stand for an adenocarcinoma line | only figure in the repo, and PubMed says the direction of the error is known |
| 4 | intermediate filament | **declare a test axis** — no band exists | — |
| 5 | stress fibre | **declare a test axis** — no band exists | — |

Under standing ruling #1 every one of these is an axis with a band and a scope, not a value. Declaring
them lets the build proceed with the artifacts saying `PI_GAP` in their own records, which is what
A2's builders already do.
