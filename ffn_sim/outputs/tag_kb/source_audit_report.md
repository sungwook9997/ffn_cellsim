# SourceEvidence hallucination audit (CrossRef)

Audited **243** SourceEvidence rows.

## Summary

| verdict | n | meaning |
|---|---|---|
| DOI_DEAD | 1 | DOI does not resolve — likely fabricated/wrong |
| NO_DOI_NOMATCH | 8 | no DOI and no confident match — unverifiable |
| CHECK | 23 | partial match (author XOR year) — review |
| NO_DOI_FOUND | 102 | real paper found — just missing DOI (suggested) |
| OK | 109 | DOI resolves, author+year match — verified |

**9 high-suspicion rows** (top of list).

## All rows (suspicion-ranked)

| verdict | citation_key | DOI | suggested | CrossRef note |
|---|---|---|---|---|
| DOI_DEAD | Funk2021_eLife | https://doi.org/10.7554/eLife.72860 |  | DOI did not resolve on CrossRef |
| NO_DOI_NOMATCH | Andreu2021_NatCommun |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Bi2015_PRX |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | K562_BBRC2019 |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Mui2016_PNAS |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Nam2016_PNAS |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Ray2017_NatCommun |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Yang_NRMCB |  |  | no CrossRef author match in top hits |
| NO_DOI_NOMATCH | Yao2011_NatCommun |  |  | no CrossRef author match in top hits |
| CHECK | AllenTildesley |  | 10.1093/oso/9780198803195.001.0001 | author-only match (year differs); candidate: allen 2017 — Computer Simulation of Liquids |
| CHECK | CalzadoMartin2016_ACSNano | https://doi.org/10.1021/acsnano.5b07162 |  | CrossRef: calzado-martín 2016 — Effect of Actin Organization on the Stiffness of Living Breast Cancer  |
| CHECK | Cox_NatCommun |  | 10.1016/j.bpj.2014.11.3079 | author-only match (year differs); candidate: cox 2015 — Probing the Mechanosensitivity of Piezo1 Channels |
| CHECK | Friedl2011_Cell | https://doi.org/10.1083/jcb.200909003 |  | CrossRef: friedl 2009 — Plasticity of cell migration: a multiscale tuning model |
| CHECK | Geiger_NRMCB |  | 10.1038/nrm2593 | author-only match (year differs); candidate: geiger 2009 — Environmental sensing through focal adhesions |
| CHECK | Hosseini2021_BiophysJ | https://doi.org/10.1016/j.bpj.2021.05.00 |  | CrossRef: maciuba 2021 — Facile tethering of stable and unstable proteins for optical tweezers  |
| CHECK | HuertaLopez2024_SciAdv |  | 10.1016/s0091-6749(00)90871-6 | author-only match (year differs); candidate: lopez-gonzalez 2000 — 442 Ulcerative colitis and food allergy |
| CHECK | KhalilFriedl_TCB |  | 10.1039/c0ib00052c | author-only match (year differs); candidate: khalil 2010 — Determinants of leader cells in collective cell migration |
| CHECK | Liew2024_CMBE | https://doi.org/10.1007/s12195-024-00802 |  | CrossRef: rubenstein 2024 — Perspectives on Incorporating a Large Language Model into a Cellular a |
| CHECK | Lindstrom2010_Biomaterials | https://doi.org/10.1103/PhysRevE.82.0519 |  | CrossRef: lindström 2010 — Biopolymer network geometries: Characterization, regeneration, and ela |
| CHECK | MegeIshiyama_JBCReview |  | 10.1101/2022.12.08.519566 | author-only match (year differs); candidate: osman 2022 — Modulation of the E-cadherin in human cells infected
        |
| CHECK | NanoConvergence2021_Glioma |  | 10.1021/acsnano.5c20682 | author-only match (year differs); candidate: li 2026 — Active Mechanics Governs Three-Dimensional Cell-in-Cell Form |
| CHECK | Panzetta2019_PNAS |  | 10.1007/978-1-4613-2045-6_7 | author-only match (year differs); candidate: maschio 1987 — Nutrition in dialysis patients |
| CHECK | Parsons_FAReview |  | 10.1038/sj.onc.1203877 | author-only match (year differs); candidate: parsons 2000 — Focal Adhesion Kinase: a regulator of focal adhesion dynamic |
| CHECK | PolacheckChen_ARBE |  | 10.1146/annurev-bioeng-110220-031722 | author-only match (year differs); candidate: chen 2022 — Current Developments and Challenges of mRNA Vaccines |
| CHECK | RocaCusachs_IntegrinReviews |  | 10.1016/j.isci.2020.100907 | author-only match (year differs); candidate: lerche 2020 — Integrin Binding Dynamics Modulate Ligand-Specific Mechanose |
| CHECK | SciRep2025_GliomaAFM | https://doi.org/10.1038/s41598-025-04841 |  | CrossRef: masud 2025 — Exploring the heterogeneity in glioblastoma cellular mechanics using i |
| CHECK | SensPlastino_NRMCB |  | 10.1088/0953-8984/27/27/273103 | author-only match (year differs); candidate: sens 2015 — Membrane tension and cytoskeleton organization in cell motil |
| CHECK | Smelser2015_BMMB | https://doi.org/10.1007/s10237-015-0680- |  | CrossRef: burd 2015 — Finite element implementation of a multiscale model of the human lens  |
| CHECK | TrepatSahai_NRC |  | 10.1038/ncb2548 | author-only match (year differs); candidate: friedl 2012 — Classifying collective cancer cell invasion |
| CHECK | Vallotton2003_BiophysJ |  | 10.70675/f72da841zbffcz4d68z9607zd6e9ada27f8d | author-only match (year differs); candidate: vallotton  — La décision publique et la crise |
| CHECK | YapKovacs_JCS |  | 10.3410/f.738614181.793585959 | author-only match (year differs); candidate: yap 2021 — Faculty Opinions recommendation of Adherens junction regulat |
| CHECK | deGennes_LiquidCrystals |  | 10.1093/oso/9780198520245.001.0001 | author-only match (year differs); candidate: gennes 1993 — The Physics of Liquid Crystals |
| NO_DOI_FOUND | Adar2020_JCS |  | 10.1073/pnas.1918203117 | candidate: adar 2020 — Active volume regulation in adhered cells |
| NO_DOI_FOUND | AlonsoMatilla2023_BiophysJ |  | 10.1101/2023.11.13.566882 | candidate: alonso-matilla 2023 — Cell intrinsic mechanical regulation of plasma membrane accu |
| NO_DOI_FOUND | Banerjee2015_PRL |  | 10.1103/physrevlett.114.228101 | candidate: banerjee 2015 — Propagating Stress Waves During Epithelial Expansion |
| NO_DOI_FOUND | Bangasser2017_NCB |  | 10.3791/56219 | candidate: bangasser 2017 — Touchscreen Sustained Attention Task (SAT) for Rats |
| NO_DOI_FOUND | Beaune2014_PNAS |  | 10.4000/books.pur.49676 | candidate: beaune 2014 — Conclusions |
| NO_DOI_FOUND | Bertet2004_Nature |  | 10.1038/nature02590 | candidate: bertet 2004 — Myosin-dependent junction remodelling controls planar cell i |
| NO_DOI_FOUND | Bobrowska2021_Cells |  | 10.31338/uw.9788323548799.pp.90-102 | candidate: bobrowska 2021 — Retoryka współczesnej estetyki i sztuki: figura i narracja w |
| NO_DOI_FOUND | Borghi2012_PNAS |  | 10.1073/pnas.1204390109 | candidate: borghi 2012 — E-cadherin is under constitutive actomyosin-generated tensio |
| NO_DOI_FOUND | Brodland2002 |  | 10.1115/1.1449491 | candidate: brodland 2002 — The Differential Interfacial Tension Hypothesis (DITH): A Co |
| NO_DOI_FOUND | BroederszMacKintosh2014_RMP |  | 10.1103/revmodphys.86.995 | candidate: broedersz 2014 — Modeling semiflexible polymer networks |
| NO_DOI_FOUND | Cadart2019_NatPhysics |  | 10.1038/s41567-019-0629-y | candidate: cadart 2019 — The physics of cell-size regulation across timescales |
| NO_DOI_FOUND | Cai2014_Cell |  | 10.1201/b17354 | candidate: cai 2014 — Ambient Diagnostics |
| NO_DOI_FOUND | Cai2022_FCDB |  | 10.3389/fcell.2021.817104 | candidate: cai 2022 — Interplay Between Iron Overload and Osteoarthritis: Clinical |
| NO_DOI_FOUND | Carlsson2018_JPhysCondMat |  | 10.2210/pdb5v7f/pdb | candidate: carlsson 2018 — T4 lysozyme Y18Ymi |
| NO_DOI_FOUND | CavalcantiAdam2007_BiophysJ |  | 10.1039/b614008d | candidate: girard 2007 — Cellular chemomechanics at interfaces: sensing, integration  |
| NO_DOI_FOUND | Cavey2008_Cell |  | 10.1038/nature06953 | candidate: cavey 2008 — A two-tiered mechanism for stabilization and immobilization  |
| NO_DOI_FOUND | Chaudhuri2015_NatCommun |  | 10.1038/ncomms7365 | candidate: chaudhuri 2015 — Substrate stress relaxation regulates cell spreading |
| NO_DOI_FOUND | Conklin2011 |  | 10.1515/9781575066288 | candidate: conklin 2011 — Oath Formulas in Biblical Hebrew |
| NO_DOI_FOUND | Doss2020_PNAS |  | 10.1201/9780429344206 | candidate: ramamoorthy 2020 — Biology, Chemistry, and Applications of Apocarotenoids |
| NO_DOI_FOUND | Douezan2011_PNAS |  | 10.1073/pnas.1018057108 | candidate: douezan 2011 — Spreading dynamics and wetting transition of cellular aggreg |
| NO_DOI_FOUND | EloseguiArtola2016_NatMater |  | 10.1016/j.bpj.2017.05.020 | candidate: elosegui-artola 2017 — Amoebae as Mechanosensitive Tanks |
| NO_DOI_FOUND | Evans2007_RBP |  | 10.1353/scu.2007.0015 | candidate: ferris 2007 — Walker Evans, 1974 |
| NO_DOI_FOUND | EvansRitchie1997_BiophysJ |  | 10.1016/s0006-3495(97)78802-7 | candidate: evans 1997 — Dynamic strength of molecular adhesion bonds |
| NO_DOI_FOUND | FarooquiFenteany2005_JCS |  | 10.1242/jcs.01577 | candidate: farooqui 2005 — Multiple rows of cells behind an epithelial wound edge exten |
| NO_DOI_FOUND | Foty1996_PNAS |  | 10.1242/dev.122.5.1611 | candidate: foty 1996 — Surface tensions of embryonic tissues predict their mutual e |
| NO_DOI_FOUND | Franz2023_NatCommun |  | 10.5194/egusphere-2024-217 | candidate: franz 2024 — Reply to Comment on Franz et al. (2023): A reinterpretation  |
| NO_DOI_FOUND | Fredberg2015_FASEB |  | 10.1096/fasebj.29.1_supplement.85.5 | candidate: krishnan 2015 — Force Chains And Gap Formation in Thrombin‐induced Endotheli |
| NO_DOI_FOUND | Fritzsche2013_MBC |  | 10.1007/978-3-642-37495-1_2 | candidate: fritzsche 2013 — Lebesgue-Theorie |
| NO_DOI_FOUND | Fu2024_PNAS |  | 10.20944/preprints202406.1065.v2 | candidate: fu 2024 — An Introduction to Cosmos Thermodynamics |
| NO_DOI_FOUND | Garcia2015_PNAS |  | 10.1073/pnas.1510973112 | candidate: garcia 2015 — Physics of active jamming during collective cellular motion  |
| NO_DOI_FOUND | Goldmann2002 |  | 10.1006/cbir.2002.0900 | candidate: goldmann 2002 — p56<sup>lck</sup> CONTROLS PHOSPHORYLATION OF FILAMIN (ABP‐2 |
| NO_DOI_FOUND | GonzalezRodriguez2012_Science |  | 10.1126/science.1226418 | candidate: gonzalez-rodriguez 2012 — Soft Matter Models of Developing Tissues and Tumors |
| NO_DOI_FOUND | Guo2017_PNAS |  | 10.1145/3145690.3145698 | candidate: guo 2017 — Importance sampling measured BRDFs based on second order sph |
| NO_DOI_FOUND | Han2017_PNAS |  | 10.2172/1477879 | candidate: han 2017 — Closeout Report for CTEQ Summer School 2017 |
| NO_DOI_FOUND | Han2021_eLife |  | 10.4211/hs.7c5e032bdc7648a4a8c863a2c175e5b6 | candidate: han 2021 — DataShare for Han et al., 2021 Ecological Indicators |
| NO_DOI_FOUND | Hannezo2014_PNAS |  | 10.1073/pnas.1312076111 | candidate: hannezo 2013 — Theory of epithelial sheet morphology in three dimensions |
| NO_DOI_FOUND | HeadLevineMacKintosh2003_PRE |  | 10.1103/physrevlett.91.108102 | candidate: head 2003 — Deformation of Cross-Linked Semiflexible Polymer Networks |
| NO_DOI_FOUND | HeadLevineMacKintosh2003_PRL |  | 10.1103/physrevlett.91.108102 | candidate: head 2003 — Deformation of Cross-Linked Semiflexible Polymer Networks |
| NO_DOI_FOUND | Helfrich1973 |  | 10.1515/znc-1973-11-1209 | candidate: helfrich 1973 — Elastic Properties of Lipid Bilayers: Theory and Possible Ex |
| NO_DOI_FOUND | Jansen2018_BiophysJ |  | 10.1016/j.jas.2018.02.016 | candidate: jansen 2018 — On the use of Cu isotope signatures in archaeometallurgy: A  |
| NO_DOI_FOUND | Jetta2023_FCDB |  | 10.3389/fcell.2023.1198109 | candidate: jetta 2023 — Epithelial cells sense local stiffness via Piezo1 mediated c |
| NO_DOI_FOUND | Kage2017_NatCommun |  | 10.1038/ncomms14832 | candidate: kage 2017 — FMNL formins boost lamellipodial force generation |
| NO_DOI_FOUND | Kim2021_PNAS |  | 10.1145/3508259.3508271 | candidate: kim 2021 — Color Separated Restoration for Lightweight Single Image Sup |
| NO_DOI_FOUND | Kothari2018_JApplMech |  | 10.2514/6.2018-0530 | candidate: rustagi 2018 — Gyroscopic Stabilization of Flying Wing Aircraft |
| NO_DOI_FOUND | Krause2014_NRMCB |  | 10.1038/nrm3861 | candidate: krause 2014 — Steering cell migration: lamellipodium dynamics and the regu |
| NO_DOI_FOUND | Kunda2008_CurrBiol |  | 10.1016/j.cub.2007.12.051 | candidate: kunda 2008 — Moesin Controls Cortical Rigidity, Cell Rounding, and Spindl |
| NO_DOI_FOUND | Labernadie2017_NCB |  | 10.1038/ncb3478 | candidate: labernadie 2017 — A mechanically active heterotypic E-cadherin/N-cadherin adhe |
| NO_DOI_FOUND | LeDuc2010_JCB |  | 10.1083/jcb.201001149 | candidate: le duc 2010 — Vinculin potentiates E-cadherin mechanosensing and is recrui |
| NO_DOI_FOUND | LeimkuhlerMatthews2013 |  | 10.1063/1.4802990 | candidate: leimkuhler 2013 — Robust and efficient configurational molecular sampling via  |
| NO_DOI_FOUND | LewisGrandl2015_eLife |  | 10.7554/elife.12088 | candidate: lewis 2015 — Mechanical sensitivity of Piezo1 ion channels can be tuned b |
| NO_DOI_FOUND | Lindstrom2013_SoftMatter |  | 10.1093/obo/9780199766567-0108 | candidate: lindstrom 2013 — Cargo Cults |
| NO_DOI_FOUND | Liu2019_ProstateMech |  | 10.1007/978-981-13-6962-9 | candidate: liu 2019 — Deuteride Materials |
| NO_DOI_FOUND | Liu2025_NatPhysics |  | 10.22541/au.176275808.84444425/v1 | candidate: liu 2025 — Comment on Bandyopadhyay et al. |
| NO_DOI_FOUND | Maitre2012_Nature |  | 10.3917/lav.daler.2012.01.0096 | candidate: bubrovszky 2012 — Pathologies schizophréniques |
| NO_DOI_FOUND | Maitre2015_NCB |  | 10.1051/shsconf/20152000001 | candidate: colón de carvajal 2015 — Préface |
| NO_DOI_FOUND | Malinova2021_NatCommun |  | 10.1017/nps.2020.87 | candidate: malinova 2021 — Politics of Memory and Nationalism |
| NO_DOI_FOUND | Malinverno2016_NatMater |  | 10.1002/2016gl070096 | candidate: nole 2016 — Short‐range, overpressure‐driven methane migration in coarse |
| NO_DOI_FOUND | MalyBorisy2001_PNAS |  | 10.1073/pnas.181338798 | candidate: maly 2001 — Self-organization of a propulsive actin network as an evolut |
| NO_DOI_FOUND | Manibog2014_NatCommun |  | 10.1038/ncomms4941 | candidate: manibog 2014 — Resolving the molecular mechanism of cadherin catch bond for |
| NO_DOI_FOUND | MarantanMahadevan2018_AmJPhys |  | 10.1119/1.5003376 | candidate: marantan 2018 — Mechanics and statistics of the worm-like chain |
| NO_DOI_FOUND | Mogilner2002_2003 |  | 10.1002/3527601503.ch14 | candidate: scholey 2002 — Mitotic Spindle Motors |
| NO_DOI_FOUND | MogilnerOster2003_BJ |  | 10.1016/s0006-3495(03)74969-8 | candidate: mogilner 2003 — Force Generation by Actin Polymerization II: The Elastic Rat |
| NO_DOI_FOUND | Munster2013_PNAS |  | 10.7551/mitpress/8982.001.0001 | candidate: munster 2013 — An Aesthesia of Networks |
| NO_DOI_FOUND | Nakamura2007 |  | 10.1083/jcb.200707073 | candidate: nakamura 2007 — Structural basis of filamin A functions |
| NO_DOI_FOUND | Notbohm2016_BiophysJ |  | 10.1016/j.bpj.2016.05.019 | candidate: notbohm 2016 — Cellular Contraction and Polarization Drive Collective Cellu |
| NO_DOI_FOUND | Odijk1995 |  | 10.1037/e493132004-001 | candidate: odijk 1995 — A syntactic condition on definite descriptions |
| NO_DOI_FOUND | Ozawa2020_JCB |  | 10.1007/978-981-15-4190-2_7 | candidate: ozawa 2020 — Comprehensive Registry in Japan |
| NO_DOI_FOUND | Pereverzev2005_BiophysJ |  | 10.1007/s11018-005-0098-9 | candidate: pereverzev 2005 — Measurement of the spectral density of broadband radio inter |
| NO_DOI_FOUND | PerezGonzalez2019_NatPhys |  | 10.2172/1597056 | candidate: perez-gonzalez 2019 — Dirac and Majorana Neutrino Signatures of Primordial Black H |
| NO_DOI_FOUND | Petitjean2010_BiophysJ |  | 10.1016/j.bpj.2010.01.030 | candidate: petitjean 2010 — Velocity Fields in a Collectively Migrating Epithelium |
| NO_DOI_FOUND | PollardBorisy2003_Cell |  | 10.1016/s0092-8674(03)00120-x | candidate: pollard 2003 — Cellular Motility Driven by Assembly and Disassembly of Acti |
| NO_DOI_FOUND | Proestaki2019_ExpMech |  | 10.1007/s11340-018-00453-4 | candidate: proestaki 2019 — Modulus of Fibrous Collagen at the Length Scale of a Cell |
| NO_DOI_FOUND | Pronk2008_PNAS |  | 10.1057/dev.2008.28 | candidate: pronk 2008 — Climate, Scarcities and Development |
| NO_DOI_FOUND | Prost2015_NatPhys |  | 10.1038/nphys3224 | candidate: prost 2015 — Active gel physics |
| NO_DOI_FOUND | Provenzano2006_BMCMed |  | 10.1186/1741-7015-4-38 | candidate: provenzano 2006 — Collagen reorganization at the tumor-stromal interface facil |
| NO_DOI_FOUND | Rakshit2012_PNAS |  | 10.1073/pnas.1208349109 | candidate: rakshit 2012 — Ideal, catch, and slip bonds in cadherin adhesion |
| NO_DOI_FOUND | Ray2021_COCellBiol |  | 10.1016/j.ceb.2021.05.004 | candidate: ray 2021 — Aligned forces: Origins and mechanisms of cancer disseminati |
| NO_DOI_FOUND | Reffay2014_NCB |  | 10.1038/ncb2917 | candidate: reffay 2014 — Interplay of RhoA and mechanical forces in collective cell m |
| NO_DOI_FOUND | Riahi2015_NatCommun |  | 10.1038/ncomms7556 | candidate: riahi 2015 — Notch1–Dll4 signalling and mechanical force regulate leader  |
| NO_DOI_FOUND | Ruscone2024_Bioinformatics |  | 10.46471/gigabyte.136 | candidate: noël 2024 — PhysiMeSS - a new physiCell addon for extracellular matrix m |
| NO_DOI_FOUND | Saraswathibhatla2023_NRMCB |  | 10.1038/s41580-023-00583-1 | candidate: saraswathibhatla 2023 — Cell–extracellular matrix mechanotransduction in 3D |
| NO_DOI_FOUND | SchwarzSafran2013_RMP |  | 10.1103/revmodphys.85.1327 | candidate: schwarz 2013 — Physics of adherent cells |
| NO_DOI_FOUND | SerraPicamal2012_NatPhys |  | 10.1038/nphys2355 | candidate: serra-picamal 2012 — Mechanical waves during tissue expansion |
| NO_DOI_FOUND | Steinberg1996_CurrTopDevBiol |  | 10.1387/ijdb.8877443 | candidate: drawbridge 1996 — Morphogenesis of the axolotl pronephric duct: a model system |
| NO_DOI_FOUND | Stylianopoulos2012_PNAS |  | 10.1073/pnas.1213353109 | candidate: stylianopoulos 2012 — Causes, consequences, and remedies for growth-induced solid  |
| NO_DOI_FOUND | TapiaRojo2019_SciAdv |  | 10.1063/1.5126071 | candidate: tapia-rojo 2019 — Thermal versus mechanical unfolding in a model protein |
| NO_DOI_FOUND | TrepatFredberg2011_TCB |  | 10.1016/j.tcb.2011.06.006 | candidate: trepat 2011 — Plithotaxis and emergent dynamics in collective cellular mig |
| NO_DOI_FOUND | Trichet2012_PNAS |  | 10.3917/bupsy.520.0365 | candidate: trichet 2012 — La notion de kakon . Histoire et enjeux psychopathologiques |
| NO_DOI_FOUND | Varma2016_BiophysJ |  | 10.4135/9781526429629 | candidate: varma 2016 — Swiss Roll |
| NO_DOI_FOUND | Vassalli2023_Cancers |  | 10.3917/mem.086.0017 | candidate: vassalli 2023 — L’hospitalité citoyenne : bien plus qu’un hébergement |
| NO_DOI_FOUND | Venturini2020_Science |  | 10.32614/cran.package.dmbc | candidate: venturini 2020 — dmbc: Model Based Clustering of Binary Dissimilarity Measure |
| NO_DOI_FOUND | Vishwakarma2018_NatCommun |  | 10.1038/s41467-018-05927-6 | candidate: vishwakarma 2018 — Mechanical interactions among followers determine the emerge |
| NO_DOI_FOUND | Wahlsten2023_ActaBiomater |  | 10.1016/j.actbio.2023.08.030 | candidate: wahlsten 2023 — Multiscale mechanical analysis of the elastic modulus of ski |
| NO_DOI_FOUND | Wang1997_BiophysJ |  | 10.1115/97-aa-050 | candidate: wang 1997 — Vibration of Skew Sandwich Plates With Laminated Facings |
| NO_DOI_FOUND | Wei2024_ActaBiomater |  | 10.26434/chemrxiv-2023-8c9vh-v2 | candidate: wei 2024 — Persistent Topological Laplacians  -- a Survey |
| NO_DOI_FOUND | Wenger2007 |  | 10.1097/mlr.0b013e31815b97bf | candidate: wenger 2007 — Comorbidity and Quality of Care: Clarification From Author |
| NO_DOI_FOUND | WirtzFriedl2011 |  | 10.1038/nrc3080 | candidate: wirtz 2011 — The physics of cancer: the role of physical interactions and |
| NO_DOI_FOUND | Yamada2019_NRMCB |  | 10.1038/s41580-019-0172-9 | candidate: yamada 2019 — Mechanisms of 3D cell migration |
| NO_DOI_FOUND | YangKaufman2009_BiophysJ |  | 10.1016/j.bpj.2008.10.063 | candidate: yang 2009 — Rheology and Confocal Reflectance Microscopy as Probes of Me |
| NO_DOI_FOUND | Yao2014_NatCommun |  | 10.1002/chin.201448205 | candidate: yao 2014 — ChemInform Abstract: New Polyesters from Talaromyces flavus. |
| NO_DOI_FOUND | Yao2014_SciRep |  | 10.1038/srep04960 | candidate: yao 2014 — Photoacoustic computed microscopy |
| NO_DOI_FOUND | Yen2020_BBRC |  | 10.1016/j.bbrc.2020.03.146 | candidate: yen 2020 — Alteration of Young’s modulus in mesenchymal stromal cells d |
| OK | BakerChen2012_JCS | https://doi.org/10.1242/jcs.079509 |  | CrossRef: baker 2012 — Deconstructing the third dimension – how 3D culture microenvironments  |
| OK | Bangasser2013_BiophysJ | https://doi.org/10.1016/j.bpj.2013.06.02 |  | CrossRef: bangasser 2013 — Determinants of Maximal Force Transmission in a Motor-Clutch Model of  |
| OK | Bastatas2012_BBA | https://doi.org/10.1016/j.bbagen.2012.02 |  | CrossRef: bastatas 2012 — AFM nano-mechanics and calcium dynamics of prostate cancer cells with  |
| OK | Bell1978_Science | https://doi.org/10.1126/science.347575 |  | CrossRef: bell 1978 — Models for the Specific Adhesion of Cells to Cells |
| OK | Bi2016_NatPhys | https://doi.org/10.1038/nphys3471 |  | CrossRef: bi 2015 — A density-independent rigidity transition in biological tissues |
| OK | Bi2016_NatPhys | https://doi.org/10.1038/nphys3471 |  | CrossRef: bi 2015 — A density-independent rigidity transition in biological tissues |
| OK | Bieling2016_Cell | https://doi.org/10.1016/j.cell.2015.11.0 |  | CrossRef: bieling 2016 — Force Feedback Controls Motor Activity and Mechanical Properties of Se |
| OK | Bieling2016_Cell | https://doi.org/10.1016/j.cell.2015.11.0 |  | CrossRef: bieling 2016 — Force Feedback Controls Motor Activity and Mechanical Properties of Se |
| OK | Bredfeldt2014_JBiomedOpt | https://doi.org/10.1117/1.JBO.19.1.01600 |  | CrossRef: bredfeldt 2014 — Computational segmentation of collagen fibers from second-harmonic gen |
| OK | BrochardWyart2006_PNAS | https://doi.org/10.1073/pnas.0602012103 |  | CrossRef: brochard-wyart 2006 — Hydrodynamic narrowing of tubes extruded from cells |
| OK | Buckley2014_Science | https://doi.org/10.1126/science.1254211 |  | CrossRef: buckley 2014 — The minimal cadherin-catenin complex binds to actin filaments under fo |
| OK | Bustamante1994_Science | https://doi.org/10.1126/science.8079175 |  | CrossRef: bustamante 1994 — Entropic Elasticity of λ-Phage DNA |
| OK | Caille2002_JBiomech | https://doi.org/10.1016/S0021-9290(01)00 |  | CrossRef: caille 2002 — Contribution of the nucleus to the mechanical properties of endothelia |
| OK | Cambria2024_NatRevCancer | https://doi.org/10.1038/s41568-023-00656 |  | CrossRef: cambria 2024 — Linking cell mechanical memory and cancer metastasis |
| OK | ChanOdde2008_Science | https://doi.org/10.1126/science.1163595 |  | CrossRef: chan 2008 — Traction Dynamics of Filopodia on Compliant Substrates |
| OK | Charras2008_BJ | https://doi.org/10.1529/biophysj.107.113 |  | CrossRef: charras 2008 — Life and Times of a Cellular Bleb |
| OK | Chaudhuri2015_NatMater | https://doi.org/10.1038/nmat4489 |  | CrossRef: chaudhuri 2015 — Hydrogels with tunable stress relaxation regulate stem cell fate and a |
| OK | Chen2020_BeilsteinJNanotechnol | https://doi.org/10.3762/bjnano.11.45 |  | CrossRef: chen 2020 — Examination of the relationship between viscoelastic properties and th |
| OK | ChughPaluch2018_JCS | https://doi.org/10.1242/jcs.186254 |  | CrossRef: chugh 2018 — The actin cortex at a glance |
| OK | ChughPaluch2018_JCS | https://doi.org/10.1242/jcs.186254 |  | CrossRef: chugh 2018 — The actin cortex at a glance |
| OK | Coceano2016_Nanotechnology | https://doi.org/10.1088/0957-4484/27/6/0 |  | CrossRef: coceano 2015 — Investigation into local cell mechanics by atomic force microscopy map |
| OK | Conklin2011_AmJPathol | https://doi.org/10.1016/j.ajpath.2010.11 |  | CrossRef: conklin 2011 — Aligned Collagen Is a Prognostic Signature for Survival in Human Breas |
| OK | Darling2008_JBiomech | https://doi.org/10.1016/j.jbiomech.2007. |  | CrossRef: darling 2008 — Viscoelastic properties of human mesenchymally-derived stem cells and  |
| OK | DelRio2009_Science | https://doi.org/10.1126/science.1162912 |  | CrossRef: del rio 2009 — Stretching Single Talin Rod Molecules Activates Vinculin Binding |
| OK | DenaisRaab2016_Science | https://doi.org/10.1126/science.aad7297 |  | CrossRef: denais 2016 — Nuclear envelope rupture and repair during cancer cell migration |
| OK | Derenyi2002_PRL | https://doi.org/10.1103/PhysRevLett.88.2 |  | CrossRef: derényi 2002 — Formation and Interaction of Membrane Tubes |
| OK | Dessard2024_NanoscaleAdv | https://doi.org/10.1039/D4NA00003J |  | CrossRef: dessard 2024 — Cytoplasmic viscosity is a potential biomarker for metastatic breast c |
| OK | Dimova2014_ACIS | https://doi.org/10.1016/j.cis.2014.03.00 |  | CrossRef: dimova 2014 — Recent developments in the field of bending rigidity measurements on m |
| OK | Discher2005_Science | https://doi.org/10.1126/science.1116995 |  | CrossRef: discher 2005 — Tissue Cells Feel and Respond to the Stiffness of Their Substrate |
| OK | DizMunoz2013_TCB | https://doi.org/10.1016/j.tcb.2012.09.00 |  | CrossRef: diz-muñoz 2013 — Use the force: membrane tension as an organizer of cell shape and moti |
| OK | DoyleYamada2009_JCB | https://doi.org/10.1083/jcb.200810041 |  | CrossRef: doyle 2009 — One-dimensional topography underlies three-dimensional fibrillar cell  |
| OK | Dupont2011_Nature | https://doi.org/10.1038/nature10137 |  | CrossRef: dupont 2011 — Role of YAP/TAZ in mechanotransduction |
| OK | EloseguiArtola2017_Cell | https://doi.org/10.1016/j.cell.2017.10.0 |  | CrossRef: elosegui-artola 2017 — Force Triggers YAP Nuclear Entry by Regulating Transport across Nuclea |
| OK | Engler2006_Cell | https://doi.org/10.1016/j.cell.2006.06.0 |  | CrossRef: engler 2006 — Matrix Elasticity Directs Stem Cell Lineage Specification |
| OK | FengSun2021_FCDB | https://doi.org/10.3389/fcell.2021.71883 |  | CrossRef: fan 2021 — Substrate Stiffness Modulates the Growth, Phenotype, and Chemoresistan |
| OK | GalWeihs2012_CBB | https://doi.org/10.1007/s12013-012-9356- |  | CrossRef: gal 2012 — Intracellular Mechanics and Activity of Breast Cancer Cells Correlate  |
| OK | Gittes1993_JCB | https://doi.org/10.1083/jcb.120.4.923 |  | CrossRef: gittes 1993 — Flexural rigidity of microtubules and actin filaments measured from th |
| OK | Guilak2000_BBRC | https://doi.org/10.1006/bbrc.2000.2360 |  | CrossRef: guilak 2000 — Viscoelastic Properties of the Cell Nucleus |
| OK | Guo2020_CellRegen | https://doi.org/10.1186/s13619-020-00054 |  | CrossRef: guo 2020 — Consistent apparent Young’s modulus of human embryonic stem cells and  |
| OK | Hammerick2010_TissueEngA | https://doi.org/10.1089/ten.tea.2010.021 |  | CrossRef: hammerick 2011 — Elastic Properties of Induced Pluripotent Stem Cells |
| OK | Harada2014_JCB | https://doi.org/10.1083/jcb.201308029 |  | CrossRef: harada 2014 — Nuclear lamin stiffness is a barrier to 3D migration, but softness can |
| OK | Hochmuth1996_BiophysJ | https://doi.org/10.1016/S0006-3495(96)79 |  | CrossRef: hochmuth 1996 — Deformation and flow of membrane into tethers extracted from neuronal  |
| OK | Hosseini2020_AdvSci | https://doi.org/10.1002/advs.202001276 |  | CrossRef: hosseini 2020 — EMT‐Induced Cell‐Mechanical Changes Enhance Mitotic Rounding Strength |
| OK | Hynes2002_Cell | https://doi.org/10.1016/s0092-8674(02)00 |  | CrossRef: hynes 2002 — Integrins |
| OK | Isambert1995_JBC | https://doi.org/10.1074/jbc.270.19.11437 |  | CrossRef: isambert 1995 — Flexibility of Actin Filaments Derived from Thermal Fluctuations |
| OK | Janmey2019_PhysiolRev | https://doi.org/10.1152/physrev.00013.20 |  | CrossRef: janmey 2020 — Stiffness Sensing by Cells |
| OK | Kanchanawong2010_Nature | https://doi.org/10.1038/nature09621 |  | CrossRef: kanchanawong 2010 — Nanoscale architecture of integrin-based cell adhesions |
| OK | Kong2009_Nature | https://doi.org/10.1083/jcb.200810002 |  | CrossRef: kong 2009 — Demonstration of catch bonds between an integrin and its ligand |
| OK | KraningRush2012_PLoSOne | https://doi.org/10.1371/journal.pone.003 |  | CrossRef: kraning-rush 2012 — Cellular Traction Stresses Increase with Increasing Metastatic Potenti |
| OK | Lammerding2004_JCI | https://doi.org/10.1172/JCI19670 |  | CrossRef: lammerding 2004 — Lamin A/C deficiency causes defective nuclear mechanics and mechanotra |
| OK | Levental2009_Cell | https://doi.org/10.1016/j.cell.2009.10.0 |  | CrossRef: levental 2009 — Matrix Crosslinking Forces Tumor Progression by Enhancing Integrin Sig |
| OK | Li2008_BBRC | https://doi.org/10.1016/j.bbrc.2008.07.0 |  | CrossRef: li 2008 — AFM indentation study of breast cancer cells |
| OK | Licup2015_PNAS | https://doi.org/10.1073/pnas.1504258112 |  | CrossRef: licup 2015 — Stress controls the mechanics of collagen networks |
| OK | Liu2015_Cell | https://doi.org/10.1016/j.cell.2015.01.0 |  | CrossRef: liu 2015 — Confinement and Low Adhesion Induce Fast Amoeboid Migration of Slow Me |
| OK | Lo2000_BiophysJ | https://doi.org/10.1016/S0006-3495(00)76 |  | CrossRef: lo 2000 — Cell Movement Is Guided by the Rigidity of the Substrate |
| OK | Lomakin2020_Science | https://doi.org/10.1126/science.aba2894 |  | CrossRef: lomakin 2020 — The nucleus acts as a ruler tailoring cell responses to spatial constr |
| OK | Marchetti2013_RMP | https://doi.org/10.1103/revmodphys.85.11 |  | CrossRef: marchetti 2013 — Hydrodynamics of soft active matter |
| OK | MarkoSiggia1995_Macromolecules | https://doi.org/10.1021/ma00130a008 |  | CrossRef: marko 1995 — Stretching DNA |
| OK | MattilaLappalainen2008_NRMCB | https://doi.org/10.1038/nrm2406 |  | CrossRef: mattila 2008 — Filopodia: molecular architecture and cellular functions |
| OK | McKenzie2018_SciRep | https://doi.org/10.1038/s41598-018-25589 |  | CrossRef: mckenzie 2018 — The mechanical microenvironment regulates ovarian cancer cell morpholo |
| OK | Moeendarbary2013_NatMater | https://doi.org/10.1038/nmat3517 |  | CrossRef: moeendarbary 2013 — The cytoplasm of living cells behaves as a poroelastic material |
| OK | Mogilner2005_BiophysJ | https://doi.org/10.1529/biophysj.104.056 |  | CrossRef: mogilner 2005 — The Physics of Filopodial Protrusion |
| OK | MogilnerOster1996_BiophysJ | https://doi.org/10.1016/S0006-3495(96)79 |  | CrossRef: mogilner 1996 — Cell motility driven by actin polymerization |
| OK | Molter2022_FCDB | https://doi.org/10.3389/fcell.2022.93251 |  | CrossRef: molter 2022 — Prostate cancer cells of increasing metastatic potential exhibit diver |
| OK | MotteKaufman2013_Biopolymers | https://doi.org/10.1002/bip.22133 |  | CrossRef: motte 2012 — Strain stiffening in collagen I networks |
| OK | Mueller2017_Cell | https://doi.org/10.1016/j.cell.2017.07.0 |  | CrossRef: mueller 2017 — Load Adaptation of Lamellipodial Actin Networks |
| OK | Murrell2015_NRMCB | https://doi.org/10.1038/nrm4012 |  | CrossRef: murrell 2015 — Forcing cells into shape: the mechanics of actomyosin contractility |
| OK | Murrell2015_NRMCB | https://doi.org/10.1038/nrm4012 |  | CrossRef: murrell 2015 — Forcing cells into shape: the mechanics of actomyosin contractility |
| OK | Ofek2009_JBiomech | https://doi.org/10.1016/j.jbiomech.2009. |  | CrossRef: ofek 2009 — In situ mechanical properties of the chondrocyte cytoplasm and nucleus |
| OK | Omidvar2014_JBiomech | https://doi.org/10.1016/j.jbiomech.2014. |  | CrossRef: omidvar 2014 — Atomic force microscope-based single cell force spectroscopy of breast |
| OK | Oria2017_Nature | https://doi.org/10.1038/nature24662 |  | CrossRef: oria 2017 — Force loading explains spatial sensing of ligands by cells |
| OK | Otto2015_NatMethods | https://doi.org/10.1038/nmeth.3281 |  | CrossRef: otto 2015 — Real-time deformability cytometry: on-the-fly cell mechanical phenotyp |
| OK | Pajerowski2007_PNAS | https://doi.org/10.1073/pnas.0702576104 |  | CrossRef: pajerowski 2007 — Physical plasticity of the nucleus in stem cell differentiation |
| OK | Palchesko2012_PLoSOne | https://doi.org/10.1371/journal.pone.005 |  | CrossRef: palchesko 2012 — Development of Polydimethylsiloxane Substrates with Tunable Elastic Mo |
| OK | Park2015_NatMater | https://doi.org/10.1038/nmat4357 |  | CrossRef: park 2015 — Unjamming and cell shape in the asthmatic airway epithelium |
| OK | Park2015_NatMater | https://doi.org/10.1038/nmat4357 |  | CrossRef: park 2015 — Unjamming and cell shape in the asthmatic airway epithelium |
| OK | PathakKumar2012_PNAS | https://doi.org/10.1073/pnas.1118073109 |  | CrossRef: pathak 2012 — Independent regulation of tumor cell migration by matrix stiffness and |
| OK | PaulKonstantopoulos2017_NatRevCancer | https://doi.org/10.1038/nrc.2016.123 |  | CrossRef: paul 2016 — Cancer cell motility: lessons from migration in confined spaces |
| OK | Pillarisetti2011_CellReprogram | https://doi.org/10.1089/cell.2011.0028 |  | CrossRef: pillarisetti 2011 — Mechanical Phenotyping of Mouse Embryonic Stem Cells: Increase in Stif |
| OK | Plodinec2012_NatNanotechnol | https://doi.org/10.1038/nnano.2012.167 |  | CrossRef: plodinec 2012 — The nanomechanical signature of breast cancer |
| OK | Plotnikov2012_Cell | https://doi.org/10.1016/j.cell.2012.11.0 |  | CrossRef: plotnikov 2012 — Force Fluctuations within Focal Adhesions Mediate ECM-Rigidity Sensing |
| OK | PontesGauthier2017_SCDB | https://doi.org/10.1016/j.semcdb.2017.08 |  | CrossRef: pontes 2017 — Membrane tension: A challenging but universal physical parameter in ce |
| OK | Rawicz2000_BiophysJ | https://doi.org/10.1016/S0006-3495(00)76 |  | CrossRef: rawicz 2000 — Effect of Chain Length and Unsaturation on Elasticity of Lipid Bilayer |
| OK | Rosenbluth2006_BiophysJ | https://doi.org/10.1529/biophysj.105.067 |  | CrossRef: rosenbluth 2006 — Force Microscopy of Nonadherent Cells: A Comparison of Leukemia Cell D |
| OK | Ruprecht2015_Cell | https://doi.org/10.1016/j.cell.2015.01.0 |  | CrossRef: ruprecht 2015 — Cortical Contractility Triggers a Stochastic Switch to Fast Amoeboid C |
| OK | Saarinen2015_PLoSOne | https://doi.org/10.1371/journal.pone.014 |  | CrossRef: saarinen 2015 — Differential Predictive Roles of A- and B-Type Nuclear Lamins in Prost |
| OK | Salbreux2012_TCB | https://doi.org/10.1016/j.tcb.2012.07.00 |  | CrossRef: salbreux 2012 — Actin cortex mechanics and cellular morphogenesis |
| OK | Salbreux2012_TCB | https://doi.org/10.1016/j.tcb.2012.07.00 |  | CrossRef: salbreux 2012 — Actin cortex mechanics and cellular morphogenesis |
| OK | Sharma2014_IntegrBiol | https://doi.org/10.1039/c3ib40246k |  | CrossRef: sharma 2014 — The role of Rho GTPase in cell stiffness and cisplatin resistance in o |
| OK | Stephens2017_MBoC | https://doi.org/10.1091/mbc.e16-09-0653 |  | CrossRef: stephens 2017 — Chromatin and lamin A determine two different mechanical response regi |
| OK | Stewart2011_Nature | https://doi.org/10.1038/nature09642 |  | CrossRef: stewart 2011 — Hydrostatic pressure and the actomyosin cortex drive mitotic cell roun |
| OK | Storm2005_Nature | https://doi.org/10.1038/nature03521 |  | CrossRef: storm 2005 — Nonlinear elasticity in biological gels |
| OK | Stroka2014_Cell | https://doi.org/10.1016/j.cell.2014.02.0 |  | CrossRef: stroka 2014 — Water Permeation Drives Tumor Cell Migration in Confined Microenvironm |
| OK | Swaminathan2011_CancerRes | https://doi.org/10.1158/0008-5472.CAN-11 |  | CrossRef: swaminathan 2011 — Mechanical Stiffness Grades Metastatic Potential in Patient Tumor Cell |
| OK | Swift2013_Science | https://doi.org/10.1126/science.1240104 |  | CrossRef: swift 2013 — Nuclear Lamin-A Scales with Tissue Stiffness and Enhances Matrix-Direc |
| OK | Tambe2011_NatMater | https://doi.org/10.1038/nmat3025 |  | CrossRef: tambe 2011 — Collective cell guidance by cooperative intercellular forces |
| OK | Thiery2009_Cell | https://doi.org/10.1016/j.cell.2009.11.0 |  | CrossRef: thiery 2009 — Epithelial-Mesenchymal Transitions in Development and Disease |
| OK | Tinevez2009_PNAS | https://doi.org/10.1073/pnas.0903353106 |  | CrossRef: tinevez 2009 — Role of cortical tension in bleb growth |
| OK | TitushkinCho2007_BiophysJ | https://doi.org/10.1529/biophysj.107.107 |  | CrossRef: titushkin 2007 — Modulation of Cellular Mechanics during Osteogenic Differentiation of  |
| OK | Trepat2009_NatPhys | https://doi.org/10.1038/nphys1269 |  | CrossRef: trepat 2009 — Physical forces during collective cell migration |
| OK | Tsujita2021_NatCommun | https://doi.org/10.1038/s41467-021-26156 |  | CrossRef: tsujita 2021 — Homeostatic membrane tension constrains cancer cell dissemination by c |
| OK | VincentEngler2013_BiotechnolJ | https://doi.org/10.1002/biot.201200205 |  | CrossRef: vincent 2013 — Mesenchymal stem cell durotaxis depends on substrate stiffness gradien |
| OK | WalcottSun2010_PNAS | https://doi.org/10.1073/pnas.0912739107 |  | CrossRef: walcott 2010 — A mechanical model of actin stress fiber formation and substrate elast |
| OK | Wolf2013_JCB | https://doi.org/10.1083/jcb.201210152 |  | CrossRef: wolf 2013 — Physical limits of cell migration: Control by ECM space and nuclear de |
| OK | Woodcock2023_Cells | https://doi.org/10.3390/cells12192401 |  | CrossRef: woodcock 2023 — Measuring Melanoma Nanomechanical Properties in Relation to Metastatic |
| OK | XuSulchek2012_PLoSOne | https://doi.org/10.1371/journal.pone.004 |  | CrossRef: xu 2012 — Cell Stiffness Is a Biomarker of the Metastatic Potential of Ovarian C |
| OK | YangAnseth2014_NatMater | https://doi.org/10.1038/nmat3889 |  | CrossRef: yang 2014 — Mechanical memory and dosing influence stem cell fate |
| OK | Yonemura2010_NCB | https://doi.org/10.1038/ncb2055 |  | CrossRef: yonemura 2010 — α-Catenin as a tension transducer that induces adherens junction devel |
| OK | Zbiral2023_IJMS | https://doi.org/10.3390/ijms241512208 |  | CrossRef: zbiral 2023 — Characterization of Breast Cancer Aggressiveness by Cell Mechanics |