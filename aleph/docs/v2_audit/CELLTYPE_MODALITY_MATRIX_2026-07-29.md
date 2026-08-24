# Cell-type x measurement-method matrix — extracted from `kb.duckdb`

274 numeric `knowledge_claim` rows scanned · **22 cell types** matched · 13 methods probed

A cell means the KB says something NUMERIC about that (cell, method) pair. It does not mean the value is usable: check the claim's own status and its citation's `source_audit` verdict.


## Coverage per cell type

| cell type | lineage | claims | methods named |
|---|---|---:|---|
| **MCF7** | breast · non-invasive carcinoma | 22 | AFM (colloidal/bead), AFM (parallel-plate), AFM (sharp-tip), AFM (unspecified), MRS (nanowire), TFM (traction) |
| **MDA-MB-231** | breast · metastatic | 12 | AFM (colloidal/bead), AFM (sharp-tip), AFM (unspecified), micropipette, tether pulling |
| **MCF10A** | breast · normal epithelial | 6 | AFM (colloidal/bead), AFM (sharp-tip), AFM (unspecified), MRS (nanowire), TFM (traction), micropipette |
| **SKOV3** | ovary · carcinoma | 6 | AFM (colloidal/bead), AFM (unspecified) |
| **neutrophil** | blood · primary — NON-ADHERENT | 5 | AFM (unspecified), microwell |
| **HeLa** | cervix · carcinoma | 4 | AFM (unspecified), tether pulling |
| **IOSE364** | ovary · normal surface epithelial | 2 | AFM (unspecified) |
| **PC-3** | prostate · metastatic | 2 | TFM (traction) |
| **blastomere** | embryo · ARCHETYPE | 2 | — |
| **hMSC** | stem · mesenchymal | 2 | AFM (unspecified), micropipette |
| **A375M** | melanoma · metastatic+ (stage 4/4) | 1 | SICM |
| **A375P** | melanoma · metastatic (stage 3/4) | 1 | SICM |
| **C. elegans** | model organism | 1 | — |
| **HL60** | blood · AML — NON-ADHERENT, cortex-only | 1 | AFM (unspecified), microwell |
| **Jurkat** | blood · ALL — NON-ADHERENT, cortex-only | 1 | AFM (unspecified), microwell |
| **K562** | blood · CML — NON-ADHERENT, cortex-only | 1 | AFM (unspecified), microwell |
| **MDCK** | kidney · epithelial (canine) | 1 | — |
| **ME10538** | melanoma · VGP (stage 2/4) | 1 | SICM |
| **PANC-1** | pancreas · carcinoma | 1 | — |
| **RBC** | blood · anucleate ARCHETYPE | 1 | — |
| **WM35** | melanoma · RGP (stage 1/4) | 1 | SICM |
| **fibroblast** | connective · primary | 1 | tether pulling |

## The matrix

| cell type | AFM (sharp-tip) | AFM (colloidal/bead) | AFM (parallel-plate) | AFM (unspecified) | micropipette | MRS (nanowire) | TFM (traction) | SICM | microwell | tether pulling | (method not stated) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **MCF7** | `KB-6.1.1` | `KB-6.1.1` | `KB-6.1.7` | `KB-PIV-8` | · | `KB-3.B3.1` | `KB-6.1.3` | · | · | · | `KB-3.33` `KB-3.26` `KB-DRAFT-7-02` `KB-DRAFT-3-21` `KB-DRAFT-3.B-30` `KB-DRAFT-3.B-29` `KB-DRAFT-3.B-26` `KB-DRAFT-3.B-24` `KB-DRAFT-3.B-19` `KB-DRAFT-3.B-15` `KB-DRAFT-3.B-13` `KB-DRAFT-3.B-05` `KB-5.18` `KB-5.20` `KB-PIV-5` `KB-PIV-7` `KB-6.1.2` |
| **MDA-MB-231** | `KB-6.1.1` | `KB-6.1.1` | · | `KB-6.1.6` `KB-6.1.5` | `KB-6.1.5` | · | · | · | · | `KB-6.1.6` | `KB-3.26` `KB-DRAFT-3.B-13` `KB-DRAFT-3.B-05` `KB-6.1.4` `KB-6.2.2` `KB-6.3.1` `KB-6.2.1` `KB-PIV-7` `KB-3.B1.5` |
| **MCF10A** | `KB-6.1.1` | `KB-6.1.1` | · | `KB-PIV-8` | `KB-3.32` | `KB-3.B3.1` | `KB-6.1.3` | · | · | · | `KB-PIV-7` |
| **SKOV3** | · | `KB-6.4.3` | · | `KB-6.4.4` `KB-6.4.2` `KB-6.4.1` `KB-6.4.5` | · | · | · | · | · | · | `KB-3.B1.5` |
| **neutrophil** | · | · | · | `KB-6.2.4` | · | · | · | · | `KB-6.2.4` | · | `KB-3.B2.3` `KB-1.V.3.3` `KB-3.B1.2` `KB-1.8` |
| **HeLa** | · | · | · | `KB-6.1.6` | · | · | · | · | · | `KB-6.1.6` | `KB-3.30` `KB-DRAFT-3.B-30` `KB-3.B3.2` |
| **IOSE364** | · | · | · | `KB-6.4.4` | · | · | · | · | · | · | `KB-3.B1.5` |
| **PC-3** | · | · | · | · | · | · | `KB-6.2.1` | · | · | · | `KB-3.B1.5` |
| **blastomere** | · | · | · | · | · | · | · | · | · | · | `KB-3.1` `KB-3.5` |
| **hMSC** | · | · | · | `KB-6.3.1` | `KB-6.3.1` | · | · | · | · | · | `KB-6.3.4` |
| **A375M** | · | · | · | · | · | · | · | `KB-6.2.2` | · | · | · |
| **A375P** | · | · | · | · | · | · | · | `KB-6.2.2` | · | · | · |
| **C. elegans** | · | · | · | · | · | · | · | · | · | · | `KB-DRAFT-3-21` |
| **HL60** | · | · | · | `KB-6.2.4` | · | · | · | · | `KB-6.2.4` | · | · |
| **Jurkat** | · | · | · | `KB-6.2.4` | · | · | · | · | `KB-6.2.4` | · | · |
| **K562** | · | · | · | `KB-6.2.4` | · | · | · | · | `KB-6.2.4` | · | · |
| **MDCK** | · | · | · | · | · | · | · | · | · | · | `KB-3.B3.2` |
| **ME10538** | · | · | · | · | · | · | · | `KB-6.2.2` | · | · | · |
| **PANC-1** | · | · | · | · | · | · | · | · | · | · | `KB-3.B1.5` |
| **RBC** | · | · | · | · | · | · | · | · | · | · | `KB-3.B1.2` |
| **WM35** | · | · | · | · | · | · | · | `KB-6.2.2` | · | · | · |
| **fibroblast** | · | · | · | · | · | · | · | · | · | `KB-3.B1.4` | · |

## Methods that appear with NO cell type attached

These are generic/model claims — a method the KB discusses without binding it to a line.

- **AFM (unspecified)** — KB-1.14, KB-1.V.1.2, KB-1.V.4.1, KB-1.V.4.3, KB-6.2.3, KB-6.3.2
- **SICM** — KB-6.2.3
- **TFM (traction)** — KB-2.19
- **cone-plate rheology** — KB-1.12, KB-1.14, KB-1.32, KB-DRAFT-7-13
- **tether pulling** — KB-3.B1.1, KB-DRAFT-3.B-02, KB-DRAFT-3.B-04, KB-DRAFT-3.B-09, KB-DRAFT-3.B-20

## Citation audit of the sources behind these claims

`source_audit` over 555 citation keys: **OK** 318, **CHECK** 166, **NO_DOI_FOUND** 49, **DOI_MISMATCH** 15, **DOI_DEAD** 7

`CLAUDE.md`: confirm a source's verdict is `OK` before citing it in a deliverable. `STATE.md` (c) 11 records that `CHECK` is ~95% a year-parsing regex artifact rather than a citation-quality statement, but that is a reason to look, not a reason to skip the look.
