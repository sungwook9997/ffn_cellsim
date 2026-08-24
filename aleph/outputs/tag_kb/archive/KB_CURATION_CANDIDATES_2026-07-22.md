# KB Curation Candidates — 2026-07-22

Generated during the TAG/Obsidian tidy pass. These are PI-decision items:
relations that CANNOT be auto-derived from a SoT field (auto-wiring would
fabricate provenance, violating the citation-integrity rule). Fix in Notion
(SoT), then run refresh.sh. Never hand-edit the vault or duckdb.

- SourceEvidence duplicate citation_keys: 7
- SourceEvidence duplicate DOIs: 10
- Isolated SE, empty claims (unattached refs): 82
- Isolated KB-DRAFT claims (promote or retire): 41
- Isolated ValidationGates (no links): 3
- Isolated RunResult/R_FF (ops-harvest territory): 42
- Isolated CodeMapping (ops-harvest territory): 10

## 1. SourceEvidence duplicates — DEDUP (keep claim-linked row, migrate, delete other)

| citation_key | rows (uid) |
|---|---|
| Belmonte2017_MolSystBiol | SE546+SE287 |
| Bi2016_NatPhys | SE113+SE5 |
| Bieling2016_Cell | SE97+SE7 |
| ChughPaluch2018_JCS | SE83+SE2 |
| Murrell2015_NRMCB | SE84+SE4 |
| Park2015_NatMater | SE114+SE6 |
| Salbreux2012_TCB | SE82+SE1 |

### Same DOI, different/duplicate citation_key
| DOI | rows |
|---|---|
| https://doi.org/10.1016/j.cell.2015.11.057 | Bieling2016_Cell(SE97) | Bieling2016_Cell(SE7) |
| https://doi.org/10.1016/j.cub.2014.05.069 | Bovellan2014_CurrBiol(SE545) | Bovellan2014(SE430) |
| https://doi.org/10.1016/j.tcb.2012.07.001 | Salbreux2012_TCB(SE82) | Salbreux2012_TCB(SE1) |
| https://doi.org/10.1038/nmat4357 | Park2015_NatMater(SE114) | Park2015_NatMater(SE6) |
| https://doi.org/10.1038/nphys3471 | Bi2016_NatPhys(SE113) | Bi2016_NatPhys(SE5) |
| https://doi.org/10.1038/nrm4012 | Murrell2015_NRMCB(SE84) | Murrell2015_NRMCB(SE4) |
| https://doi.org/10.1038/s41598-025-04841-4 | SciRep2025_GliomaAFM(SE189) | Masud2025_SciRep(SE190) |
| https://doi.org/10.1103/physrevlett.91.108102 | HeadLevineMacKintosh2003_PRL(SE44) | HeadLevineMacKintosh2003_PRE(SE45) |
| https://doi.org/10.1242/jcs.186254 | ChughPaluch2018_JCS(SE83) | ChughPaluch2018_JCS(SE2) |
| https://doi.org/10.15252/msb.20177796 | Belmonte2017_MolSystBiol(SE546) | Belmonte2017_MolSystBiol(SE287) |

## 2. Isolated SourceEvidence — empty claims (attach to a KnowledgeClaim, or mark reference-only)

- Alert2016 (SE415) — doi=https://doi.org/10.1103/physrevlett.116.068101 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Ballestrem2001 (SE478) — doi=https://doi.org/10.1083/jcb.200107107 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Balog2007_BiophysicalJournal (SE327) — doi=https://doi.org/10.1529/biophysj.106.093195 anchor=ORPHAN + OFF-DOMAIN flag 2026-06-11 (AI review): topic = yeast phosphoglycerate kinase enzyme intradomain motions — unrelated to cytoskeleton/cell-mechanics KB. Likely mis-ingested. No claim linked. PI review for relevance/removal.
- Billington2013 (SE426) — doi=https://doi.org/10.1074/jbc.m113.499848 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Biot1941 (SE457) — doi=https://doi.org/10.1063/1.1712886 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Biot1941-2 (SE466) — doi=https://doi.org/10.1016/j.semcdb.2008.01.008 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Blanchoin2014 (SE433) — doi=https://doi.org/10.1152/physrev.00018.2013 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Block2017 (SE448) — doi=https://doi.org/10.1103/physrevlett.118.048101 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Block2018 (SE450) — doi=https://doi.org/10.1126/sciadv.aat1161 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Bois2011 (SE434) — doi=https://doi.org/10.1103/physrevlett.106.028103 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Bovellan2014 (SE430) — doi=https://doi.org/10.1016/j.cub.2014.05.069 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Brangwynne2006 (SE438) — doi=https://doi.org/10.1083/jcb.200601060 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Brinkman1949 (SE459) — doi=https://doi.org/10.1007/bf02120313 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Chandran2004 (SE481) — doi=https://doi.org/10.1115/1.1688774 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Claessens2006 (SE437) — doi=https://doi.org/10.1038/nmat1718 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Dejardin2020 (SE451) — doi=https://doi.org/10.1083/jcb.201908036 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Dogterom1997 (SE439) — doi=https://doi.org/10.1126/science.278.5339.856 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Enrico2009 (SE482) — doi=https://doi.org/10.1103/physrevlett.102.088102 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Erdmann2004 (SE476) — doi=https://doi.org/10.1103/physrevlett.92.108102 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- ErdmannAlbertSchwarz2013 (SE422) — doi=https://doi.org/10.1063/1.4827497 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Fassler2020_EMBOJ (SE335) — doi=https://doi.org/10.15252/embj.2019104254 anchor=registered 2026-06-30 (FF Stage 6L inc2)
- Figard2014 (SE418) — doi=https://doi.org/10.4161/bioa.29069 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Fischer2020 (SE462) — doi=https://doi.org/10.3389/fcell.2020.00393 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Fiucci2002 (SE419) — doi=https://doi.org/10.1038/sj.onc.1205300 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Furuike2001 (SE454) — doi=https://doi.org/10.1016/s0014-5793(01)02497-8 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Gillespie2001 (SE484) — doi=https://doi.org/10.1063/1.1378322 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Gruening2021_AppliedSciences (SE329) — doi=https://doi.org/10.3390/app11125689 anchor=ORPHAN + peripheral flag 2026-06-11 (AI review): topic = automatic actin-filament quantification + cell-shape modeling of OSTEOBLASTS (imaging/method, not MCF7 cortex). No clean single-cell KB claim match. Left unlinked. PI review.
- Guo2013 (SE445) — doi=https://doi.org/10.1016/j.bpj.2013.08.037 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Guo2017 (SE465) — doi=https://doi.org/10.1073/pnas.1705179114 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Happel1983 (SE458) — doi=https://doi.org/10.1007/978-94-009-8352-6 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Helfrich1973-2 (SE432) — doi=https://doi.org/10.1515/znc-1973-11-1209 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Hill1938 (SE428) — doi=https://doi.org/10.1098/rspb.1938.0050 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Hoffmann2009 (SE467) — doi=https://doi.org/10.1152/physrev.00037.2007 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- HonerkampSmith2013 (SE414) — doi=https://doi.org/10.1103/physrevlett.111.038103 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Hotulainen2006 (SE436) — doi=https://doi.org/10.1083/jcb.200511093 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Iturri2020 (SE480) — doi=https://doi.org/10.3390/cells9040935 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Jiang2013 (SE488) — doi=https://doi.org/10.1016/j.bpj.2013.06.021 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kassianidou2017 (SE435) — doi=https://doi.org/10.1073/pnas.1606649114 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kay2017 (SE472) — doi=https://doi.org/10.3389/fcell.2017.00041 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kedem1958 (SE468) — doi=https://doi.org/10.1016/0006-3002(58)90330-5 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kojima1994_PNAS (SE334) — doi=https://doi.org/10.1073/pnas.91.26.12962 anchor=registered 2026-06-30 (FF Stage 6c)
- Kollman2011 (SE443) — doi=https://doi.org/10.1038/nrm3209 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kong2013 (SE474) — doi=https://doi.org/10.1016/j.molcel.2013.01.015 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kovcs2003 (SE425) — doi=https://doi.org/10.1074/jbc.m305453200 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Kovcs2007 (SE424) — doi=https://doi.org/10.1073/pnas.0701181104 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Laan2012 (SE444) — doi=https://doi.org/10.1016/j.cell.2012.01.007 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Lammerding2011 (SE460) — doi=https://doi.org/10.1002/cphy.c100038 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Lavie1998 (SE420) — doi=https://doi.org/10.1074/jbc.273.49.32380 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Lipowsky1991 (SE421) — doi=https://doi.org/10.1038/349475a0 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Lombardi2011 (SE464) — doi=https://doi.org/10.1074/jbc.m111.233700 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Lymn1971 (SE453) — doi=https://doi.org/10.1021/bi00801a004 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Maître2012 (SE477) — doi=https://doi.org/10.1126/science.1225399 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Mitchison1984 (SE442) — doi=https://doi.org/10.1038/312237a0 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Mücke2004 (SE447) — doi=https://doi.org/10.1016/j.jmb.2003.11.038 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Nie2015_Cytoskeleton (SE337) — doi=https://doi.org/10.1002/cm.21207 anchor=registered 2026-06-30 (gamma-floor; PI-approved Decision 2)
- Nédélec2007 (SE440) — doi=https://doi.org/10.1088/1367-2630/9/11/427 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Nöding2014 (SE452) — doi=https://doi.org/10.1016/j.bpj.2014.09.050 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Peskin2002 (SE486) — doi=https://doi.org/10.1017/s0962492902000077 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Pollard1986 (SE483) — doi=https://doi.org/10.1083/jcb.103.6.2747 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Raucher1999 (SE417) — doi=https://doi.org/10.1016/s0006-3495(99)77040-2 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Roll-Mecak2010 (SE446) — doi=https://doi.org/10.1016/j.ceb.2009.11.001 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Rotne1969 (SE456) — doi=https://doi.org/10.1063/1.1670977 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Saffman1975 (SE413) — doi=https://doi.org/10.1073/pnas.72.8.3111 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Schafer1996 (SE487) — doi=https://doi.org/10.1083/jcb.135.1.169 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Scriven1960 (SE412) — doi=https://doi.org/10.1016/0009-2509(60)87003-0 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Serwas2022 (SE431) — doi=https://doi.org/10.1016/j.devcel.2022.04.012 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Smith2007 (SE479) — doi=https://doi.org/10.1371/journal.pbio.0050268 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- StamAlbertsGardelMunro2015 (SE423) — doi=https://doi.org/10.1016/j.bpj.2015.03.030 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Stephens2018 (SE461) — doi=https://doi.org/10.1091/mbc.e17-06-0410 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Strychalski2016 (SE490) — doi=https://doi.org/10.1016/j.bpj.2016.01.012 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Suarez2011 (SE485) — doi=https://doi.org/10.1016/j.cub.2011.03.064 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Suzuki2017 (SE489) — doi=https://doi.org/10.1073/pnas.1616001114 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Svitkina1996 (SE449) — doi=https://doi.org/10.1083/jcb.135.4.991 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Syeda2016 (SE469) — doi=https://doi.org/10.1016/j.cell.2015.12.031 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Truong2021 (SE429) — doi=https://doi.org/10.1038/s41467-021-26611-2 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Tseng2004 (SE463) — doi=https://doi.org/10.1242/jcs.01073 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Veigel2003 (SE427) — doi=https://doi.org/10.1038/ncb1060 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Vicente-Manzanares2009 (SE455) — doi=https://doi.org/10.1038/nrm2786 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Vinzenz2012_JCS (SE336) — doi=none anchor=registered 2026-06-30 (FF Stage 6L inc2) - DOI NEEDS-VERIFY
- Wagner2011_FreeRadicBiolMed (SE544) — doi=https://doi.org/10.1016/j.freeradbiomed.2011.05.024 anchor=OK (DOI verified live via PubMed 2026-07-21; PMC open-access)
- Walker1988 (SE441) — doi=https://doi.org/10.1083/jcb.107.4.1437 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify
- Yao2016 (SE473) — doi=https://doi.org/10.1038/ncomms11966 anchor=auto-ingest 2026-07-15 (AUDIT-2026-07-15 missing-inventory audit); Source Type + Claims unclassified — PI to classify

## 3. Isolated KB-DRAFT claims — promote to KB-x.y or retire

- KB-DRAFT-1-01 — status=draft unit=Unit 1 ECM — Excluded volume / steric fiber-fiber repulsion
- KB-DRAFT-1-02 — status=draft unit=Unit 1 ECM — Compressed-fiber buckling (tension-compression asymmetry)
- KB-DRAFT-1-04 — status=draft unit=Unit 1 ECM — Fiber-length polydispersity in ECM build
- KB-DRAFT-3-01 — status=draft unit=Unit 3 Cell — Stiff polar beam bending κ_MT (EI = k_BT·L_p)
- KB-DRAFT-3-02 — status=draft unit=Unit 3 Cell — Euler buckling F_crit = π²·EI/L² as EMERGENT compression strut
- KB-DRAFT-3-03 — status=draft unit=Unit 3 Cell — MT tip ↔ cortex soft excluded-volume contact
- KB-DRAFT-3-04 — status=draft unit=Unit 3 Cell — MT length monodisperse (all arms L=6µm)
- KB-DRAFT-3-05 — status=draft unit=Unit 3 Cell — Anisotropic slender-body (cylinder) drag for MT nodes
- KB-DRAFT-3-06 — status=draft unit=Unit 3 Cell — MT polarity (plus-end / minus-end)
- KB-DRAFT-3-07 — status=draft unit=Unit 3 Cell — Dynamic instability (growth / catastrophe / rescue / shrink)
- KB-DRAFT-3-08 — status=draft unit=Unit 3 Cell — γ-tubulin nucleation from MTOC (aster generation / MT number regulatio
- KB-DRAFT-3-09 — status=draft unit=Unit 3 Cell — Kinesin / dynein directed transport ALONG microtubules
- KB-DRAFT-3-10 — status=draft unit=Unit 3 Cell — Cortical-dynein pulling / MTOC-aster centering + nucleus positioning
- KB-DRAFT-3-11 — status=draft unit=Unit 3 Cell — MT↔actin / MT↔cortex distributed crosslinks (MAPs, +TIPs, spectraplaki
- KB-DRAFT-3-12 — status=draft unit=Unit 3 Cell — MT severing (katanin / spastin)
- KB-DRAFT-3-15 — status=draft unit=Unit 3 Cell — Filamin catch-slip crosslinker (Pereverzev) not in the runtime — produ
- KB-DRAFT-3-20 — status=draft unit=Unit 3 Cell — Stress-fiber prestress magnitude unreachable (5-30 nN target; twin of 
- KB-DRAFT-3-22 — status=draft unit=Unit 3 Cell — Filopodium tip complex + membrane sheath + mechanosensing readout abse
- KB-DRAFT-3-23 — status=draft unit=Unit 3 Cell — Fascin / espin / fimbrin / villin bundlers UNSOURCED — filopodium & mi
- KB-DRAFT-3-24 — status=draft unit=Unit 3 Cell — SF graded antiparallel polarity + sarcomeric Z-body/band placement (ge
- KB-DRAFT-3-25 — status=draft unit=Unit 3 Cell — KMC event loop / tau-leaping scheduler (topology-count driver)
- KB-DRAFT-3-26 — status=draft unit=Unit 3 Cell — Polymerize/depolymerize as TOPOLOGY change (bead insertion/removal)
- KB-DRAFT-3-27 — status=draft unit=Unit 3 Cell — Cofilin/ADF severing (creates new barbed + pointed ends)
- KB-DRAFT-3-28 — status=draft unit=Unit 3 Cell — Barbed-end capping-protein + formin processive elongation (state machi
- KB-DRAFT-3.B-01 — status=draft unit=Unit 3.B Compartment — Nuclear envelope as a coherent 2D shell (Helfrich bending kappa_NE + l
- KB-DRAFT-3.B-02 — status=draft unit=Unit 3.B Compartment — Chromatin as an internal crosslinked polymer network (short-strain ela
- KB-DRAFT-3.B-03 — status=draft unit=Unit 3.B Compartment — Nucleoplasm as a distinct viscous (poroelastic) medium
- KB-DRAFT-3.B-04 — status=draft unit=Unit 3.B Compartment — LINC complex — nucleus<->cytoskeleton force transfer (nesprin/SUN bond
- KB-DRAFT-3.B-05 — status=draft unit=Unit 3.B Compartment — Nuclear stiffness setpoint E_nuc (MCF7 in-situ 399 Pa) — magnitude anc
- KB-DRAFT-3.B-17 — status=draft unit=Unit 3.B Compartment — Engaged (myosin–actin overlap) fraction
- KB-DRAFT-3.B-19 — status=draft unit=Unit 3.B Compartment — RPY / Oseen hydrodynamic mobility (long-range solvent-mediated node co
- KB-DRAFT-3.B-20 — status=draft unit=Unit 3.B Compartment — Brinkman interface term at bleb / sub-membrane (Darcy porous ↔ free St
- KB-DRAFT-3.B-22 — status=draft unit=Unit 3.B Compartment — In-plane 2-D membrane FLUID — Scriven-Boussinesq surface viscosity eta
- KB-DRAFT-3.B-23 — status=draft unit=Unit 3.B Compartment — Bleb nucleation + inflation (emergent membrane detachment under hydros
- KB-DRAFT-3.B-24 — status=draft unit=Unit 3.B Compartment — Area reservoir / wrinkle-fold unfolding (Raucher-Sheetz tension platea
- KB-DRAFT-3.B-25 — status=draft unit=Unit 3.B Compartment — Spontaneous curvature C0 (leaflet-asymmetry / protein-induced curvatur
- KB-DRAFT-3.B-26 — status=draft unit=Unit 3.B Compartment — Kedem-Katchalsky hydraulic drainage (water flux J_v, membrane clock ta
- KB-DRAFT-3.B-28 — status=draft unit=Unit 3.B Compartment — Regulatory Volume change (RVI/RVD) via dynamic solute state N(t)
- KB-DRAFT-3.B-29 — status=draft unit=Unit 3.B Compartment — Aquaporin-emergent hydraulic conductivity Lp
- KB-DRAFT-3.B-32 — status=draft unit=Unit 3.B Compartment — Scriven-Boussinesq in-plane membrane lipid flow (2-D surface viscosity
- KB-DRAFT-3.B-33 — status=draft unit=Unit 3.B Compartment — HAZARD-B6: explicit ECM fiber net + continuum (Winkler) substrate in s

## 4. Isolated ValidationGates — wire to contract/run/code (PI-authored)

- VG-231-G0 — type=Verification status=not-started — G0 — HeLa rounded-cell absolute-mechanics oracle (ΔP & γ jointly)
- VG-231-G5 — type=Validation status=not-started — G5 — lineage generalization: NIH/3T3 migration-traction
- VG-231-G6 — type=Validation status=not-started — G6 — final 3D invasion challenge: HT1080 MT1-MMP tunnel

## 5. Isolated RunResult / CodeMapping — closed by harvest_ops.py --apply

RunResult/R_FF (42): RUN-h1_baoab_freeze-closeout, RUN-h3-gpu-native38k-feas, RUN-h3-gpu-pilot-principled, RUN-h3-gpu-probe-v50, RUN-h3-ku35-degcap-verify, RUN-h3-ku35-smoke-seed1-fanin, RUN-h7-closeout, RUN-h7-closeout (Run-3951), RUN-h7-force-stack-profile, RUN-h7-force-stack-profile-sweep, RUN-h7-h7-basal-apparatus-gate, RUN-h7-h7-basal-mesh-gate, RUN-h7-h7-basal-nmii-gate, RUN-h7-h7-basal-sf-force-budget, RUN-h7-h7-basal-surface-gate, RUN-h7-h7-cadherin-junction-activation-gate, RUN-h7-h7-compartment-gpu-port, RUN-h7-h7-dt-push-test, RUN-h7-h7-intermediate-filaments-activation-gate, RUN-h7-h7-junctional-actin-activation-gate, RUN-h7-h7-linc-activation-gate, RUN-h7-h7-membrane-reservoir-activation-gate, RUN-h7-h7-microtubules-activation-gate, RUN-h7-h7-native-fullcell-go, RUN-h7-h7-no-nucleus-dt, RUN-h7-h7-nucleus-contribution, RUN-h7-h7-osmotic-activation-gate, RUN-h7-h7-stress-fibers-activation-gate, RUN-h7-hotloop-scaling-smoke, RUN-h_0_2-closeout, RUN-h_dcm_active-closeout, RUN-h_dcm_gpu-closeout, RUN-h_dcm_gpu_lod-closeout, RUN-h_dcm_native-closeout, RUN-h_dcm_two_stage-closeout, R_FF_STAGE6J_GPU_NATIVE_2026-06-30, R_FF_STAGE6K_TURNOVER_2026-06-30, R_FF_STAGE6L_UNIFIED_ARCHITECTURE_2026-06-30, R_FF_STAGE6M_CROSSLINK_STIFFNESS_2026-06-30, R_FF_STAGE6N_CROSSLINK_TURNOVER_2026-06-30, R_FF_STAGE6O_NATIVE_GAMMA_2026-07-01, R_FF_STAGE6P_CORTICAL_TENSION_DEFINITION_2026-07-01

CodeMapping (10): CM-common-compartments, CM-common-filament_math, CM-common-gsd_traj, CM-common-integrity, CM-common-production_policy, CM-common-sim_realtime, CM-common-surface_manifold, CM-integrator-baoab, CM-integrator-baoab_device, CM-integrator-constrained_baoab

