# SourceEvidence hallucination audit (CrossRef)

Audited **332** SourceEvidence rows.

## Summary

| verdict | n | meaning |
|---|---|---|
| DOI_DEAD | 6 | DOI does not resolve — likely fabricated/wrong |
| DOI_MISMATCH | 7 | DOI resolves to a DIFFERENT paper (author+year both off) |
| CHECK | 96 | partial match (author XOR year) — review |
| NO_DOI_FOUND | 46 | real paper found — just missing DOI (suggested) |
| OK | 177 | DOI resolves, author+year match — verified |

**13 high-suspicion rows** (top of list).

## All rows (suspicion-ranked)

| verdict | citation_key | DOI | suggested | CrossRef note |
|---|---|---|---|---|
| DOI_DEAD | 2023-published-by-the-company-of-biologi-8c227c | https://doi.org/10.1242/dev.201099: |  | DOI did not resolve on CrossRef |
| DOI_DEAD | Funk2021_eLife | https://doi.org/10.7554/eLife.72860 |  | DOI did not resolve on CrossRef |
| DOI_DEAD | NishitaniMiura2025_arXiv | https://doi.org/10.48550/arxiv.2504.1488 |  | DOI did not resolve on CrossRef |
| DOI_DEAD | Runser2023_bioRxiv | https://doi.org/10.3929/ethz-b-000653764 |  | DOI did not resolve on CrossRef |
| DOI_DEAD | Thiticharoentam2026_bioRxiv | https://doi.org/10.1101/2026.03.17.71250 |  | DOI did not resolve on CrossRef |
| DOI_DEAD | j-indian-inst-sci-vol-101-1-97-112-janua-639fda | https://doi.org/10.1038/s4156 |  | DOI did not resolve on CrossRef |
| DOI_MISMATCH | Bcher2021_FrontPhys | https://doi.org/10.3389/fphy.2021.753230 |  | CrossRef: bächer 2021 — A Three-Dimensional Numerical Model of an Active Cell Cortex in the Vi |
| DOI_MISMATCH | Brningk2019_SciRep | https://doi.org/10.1038/s41598-019-54117 |  | CrossRef: brüningk 2019 — A cellular automaton model for spheroid response to radiation and hype |
| DOI_MISMATCH | DeBelly2023_Cell | https://doi.org/10.1016/j.cell.2023.05.0 |  | CrossRef: de belly 2023 — Cell protrusions and contractions generate long-range membrane tension |
| DOI_MISMATCH | DeVecchis2021_BiophysicalJournal | https://doi.org/10.1016/j.bpj.2021.02.00 |  | CrossRef: de vecchis 2021 — Molecular dynamics simulations of Piezo1 channel opening by increases  |
| DOI_MISMATCH | JuelPrtner2025_JRSocInterface | https://doi.org/10.1098/rsif.2024.0885 |  | CrossRef: juel pørtner 2025 — Viscoelastic differences between isolated and live MCF7 cancer cell nu |
| DOI_MISMATCH | Lchtefeld2024_NatMethods | https://doi.org/10.1038/s41592-024-02277 |  | CrossRef: lüchtefeld 2024 — Dissecting cell membrane tension dynamics and its effect on Piezo1-med |
| DOI_MISMATCH | MerinoCasallo2022_CellAdhesionAmpMigration | https://doi.org/10.1080/19336918.2022.20 |  | CrossRef: merino-casallo 2022 — Unravelling cell migration: defining movement from the cell surface |
| CHECK | AllenTildesley |  | 10.1093/oso/9780198803195.001.0001 | author-only match (year differs); candidate: allen 2017 — Computer Simulation of Liquids |
| CHECK | Arruda2026_NpjSystBiolAppl | https://doi.org/10.1038/s41540-026-00648 |  | CrossRef: arruda 2026 — Simulation-based inference of cell migration dynamics in complex spati |
| CHECK | Arslan2024_CurrentBiology | https://doi.org/10.1016/j.cub.2023.11.06 |  | CrossRef: arslan 2024 — Adhesion-induced cortical flows pattern E-cadherin-mediated cell conta |
| CHECK | Audoin2022_SciRep | https://doi.org/10.1038/s41598-022-18950 |  | CrossRef: audoin 2022 — Tumor spheroids accelerate persistently invading cancer cells |
| CHECK | Balog2007_BiophysicalJournal | https://doi.org/10.1529/biophysj.106.093 |  | CrossRef: balog 2007 — The Influence of Interdomain Interactions on the Intradomain Motions i |
| CHECK | Banerjee2017_NatCommun | https://doi.org/10.1038/s41467-017-01130 |  | CrossRef: banerjee 2017 — Actomyosin pulsation and flows in an active elastomer with turnover an |
| CHECK | Belmonte2017_MolSystBiol | https://doi.org/10.15252/msb.20177796 |  | CrossRef: belmonte 2017 — A theory that predicts behaviors of disordered cytoskeletal networks |
| CHECK | Berens2015_JoVE | https://doi.org/10.3791/53409 |  | CrossRef: berens 2015 — A Cancer Cell Spheroid Assay to Assess Invasion in a 3D Setting |
| CHECK | Betorz2023_JournalOfTheMechanicsAnd | https://doi.org/10.1016/j.jmps.2023.1053 |  | CrossRef: betorz 2023 — A computational model for early cell spreading, migration, and competi |
| CHECK | Bohec2025_Jour | https://doi.org/10.64898/2025.12.15.6944 |  | CrossRef: bohec 2025 — Control of cellular cortical tension and shape by RhoGTPase signalling |
| CHECK | Boot2021_AdvancesInPhysicsX | https://doi.org/10.1080/23746149.2021.19 |  | CrossRef: boot 2021 — Spheroid mechanics and implications for cell invasion |
| CHECK | Botticelli2025_FrontBioengBiotechnol | https://doi.org/10.3389/fbioe.2025.15159 |  | CrossRef: botticelli 2025 — A hybrid computational model of cancer spheroid growth with ribose-ind |
| CHECK | Bustamante2021_Biofabrication | https://doi.org/10.1088/1758-5090/abe025 |  | CrossRef: bustamante 2021 — Biofabrication of spheroids fusion-based tumor models: computational s |
| CHECK | CalzadoMartin2016_ACSNano | https://doi.org/10.1021/acsnano.5b07162 |  | CrossRef: calzado-martín 2016 — Effect of Actin Organization on the Stiffness of Living Breast Cancer  |
| CHECK | Chandrasekaran2024_NatCommun | https://doi.org/10.1038/s41467-024-46726 |  | CrossRef: chandrasekaran 2024 — Kinetic trapping organizes actin filaments within liquid-like protein  |
| CHECK | Chen2017_MBE | https://doi.org/10.3934/mbe.2018016 |  | CrossRef: chen 2017 — A multiscale model for heterogeneous tumor spheroid in vitro |
| CHECK | Chen2024_NatPhys | https://doi.org/10.1038/s41567-024-02626 |  | CrossRef: chen 2024 — Energy partitioning in the cell cortex |
| CHECK | Chugh2017_NatCellBiol | https://doi.org/10.1038/ncb3525 |  | CrossRef: chugh 2017 — Actin cortex architecture regulates cell surface tension |
| CHECK | Cornell2025_NatCellBiol | https://doi.org/10.1038/s41556-025-01807 |  | CrossRef: cornell 2025 — Target cell cortical tension regulates macrophage trogocytosis |
| CHECK | Cox_NatCommun | https://doi.org/10.1038/ncomms10366 |  | CrossRef: cox 2016 — Removal of the mechanoprotective influence of the cytoskeleton reveals |
| CHECK | De2018_CommunBiol | https://doi.org/10.1038/s42003-018-0084- |  | CrossRef: de 2018 — A general model of focal adhesion orientation dynamics in response to  |
| CHECK | Dhandapani2023_AdvHealthcareMaterials | https://doi.org/10.1002/adhm.202300164 |  | CrossRef: dhandapani 2023 — In Vitro 3D Spheroid Model Preserves Tumor Microenvironment of Hot and |
| CHECK | Dmitrieff2017_ProcNatlAcadSciUSA | https://doi.org/10.1073/pnas.1618041114 |  | CrossRef: dmitrieff 2017 — Balance of microtubule stiffness and cortical tension determines the s |
| CHECK | Ennomani2016_CurrentBiology | https://doi.org/10.1016/j.cub.2015.12.06 |  | CrossRef: ennomani 2016 — Architecture and Connectivity Govern Actin Network Contractility |
| CHECK | Fang2016_PhysRevE | https://doi.org/10.1103/physreve.93.0424 |  | CrossRef: fang 2016 — Modeling the mechanics of cells in the cell-spreading process driven b |
| CHECK | Fastabend2026_PhysRevE | https://doi.org/10.1103/z152-x4l1 |  | CrossRef: fastabend 2026 — Cortical tension links curvature to tissue growth in the cellular Pott |
| CHECK | Feng2025_BiochemicalEngineeringJo | https://doi.org/10.1016/j.bej.2024.10956 |  | CrossRef: feng 2025 — The impact of 3D tumor spheroid maturity on cell migration and invasio |
| CHECK | Flormann2024_ProcNatlAcadSciUSA | https://doi.org/10.1073/pnas.2320372121 |  | CrossRef: flormann 2024 — The structure and mechanics of the cell cortex depend on the location  |
| CHECK | Freedman2017_BiophysicalJournal | https://doi.org/10.1016/j.bpj.2017.06.00 |  | CrossRef: freedman 2017 — A Versatile Framework for Simulating the Dynamic Mechanical Structure  |
| CHECK | Friedl2011_Cell | https://doi.org/10.1083/jcb.200909003 |  | CrossRef: friedl 2009 — Plasticity of cell migration: a multiscale tuning model |
| CHECK | Fritzsche2016_SciAdv | https://doi.org/10.1126/sciadv.1501337 |  | CrossRef: fritzsche 2016 — Actin kinetics shapes cortical network structure and mechanics |
| CHECK | Fritzsche2017_NatCommun | https://doi.org/10.1038/ncomms14347 |  | CrossRef: fritzsche 2017 — Self-organizing actin patterns shape membrane architecture but not cel |
| CHECK | Fu2019_ColloidsAndSurfacesBBioi | https://doi.org/10.1016/j.colsurfb.2018. |  | CrossRef: fu 2019 — Spontaneous formation of tumor spheroid on a hydrophilic filter paper  |
| CHECK | Garlick2022_SciRep | https://doi.org/10.1038/s41598-022-06702 |  | CrossRef: garlick 2022 — Simple methods for quantifying super-resolved cortical actin |
| CHECK | Geiger2022_PLoSONE | https://doi.org/10.1371/journal.pone.026 |  | CrossRef: geiger 2022 — Directed invasion of cancer cell spheroids inside 3D collagen matrices |
| CHECK | Geiger_NRMCB |  | 10.1038/nrm2593 | author-only match (year differs); candidate: geiger 2009 — Environmental sensing through focal adhesions |
| CHECK | Gowrishankar2012_Cell | https://doi.org/10.1016/j.cell.2012.05.0 |  | CrossRef: gowrishankar 2012 — Active Remodeling of Cortical Actin Regulates Spatiotemporal Organizat |
| CHECK | Gruening2021_AppliedSciences | https://doi.org/10.3390/app11125689 |  | CrossRef: gruening 2021 — Automatic Actin Filament Quantification and Cell Shape Modeling of Ost |
| CHECK | Herold2023_PLoSComputBiol | https://doi.org/10.1371/journal.pcbi.101 |  | CrossRef: herold 2023 — Development of a scoring function for comparing simulated and experime |
| CHECK | Hou2018_SciRep | https://doi.org/10.1038/s41598-018-25337 |  | CrossRef: hou 2018 — TASI: A software tool for spatial-temporal quantification of tumor sph |
| CHECK | HuertaLopez2024_SciAdv | https://doi.org/10.1126/sciadv.adf9758 |  | CrossRef: huerta-lópez 2024 — Cell response to extracellular matrix viscous energy dissipation outwe |
| CHECK | Jin2020_Jour | https://doi.org/10.1101/2020.12.06.41385 |  | CrossRef: jin 2020 — Mathematical model of tumour spheroid experiments with real-time cell  |
| CHECK | Jo2025_ActaBiomaterialia | https://doi.org/10.1016/j.actbio.2025.05 |  | CrossRef: jo 2025 — Reciprocal folding dynamics in cellular networks at the stroma-basemen |
| CHECK | K562_BBRC2019 | https://doi.org/10.1016/j.bbrc.2019.06.0 |  | CrossRef: takahashi 2019 — Inhibition of EP2/EP4 prostanoid receptor-mediated signaling suppresse |
| CHECK | Kadzik2026_Jour | https://doi.org/10.64898/2026.05.24.7275 |  | CrossRef: kadzik 2026 — Rapid actin filament turnover maintains cortical connectivity while al |
| CHECK | Keren2023_Cell | https://doi.org/10.1016/j.cell.2023.05.0 |  | CrossRef: keren 2023 — Effective membrane tension: A long-range integrator of cellular dynami |
| CHECK | KhalilFriedl_TCB |  | 10.1039/c0ib00052c | author-only match (year differs); candidate: khalil 2010 — Determinants of leader cells in collective cell migration |
| CHECK | Kim2012_IntegrBiol | https://doi.org/10.1039/c2ib20159c |  | CrossRef: kim 2012 — Integrating focal adhesion dynamics, cytoskeleton remodeling, and acti |
| CHECK | Kim2023_Jour | https://doi.org/10.21203/rs.3.rs-3736468 |  | CrossRef: kim 2023 — 3D spheroid model of adipose-derived stem cell and breast cancer cell  |
| CHECK | Lenz2012_PhysRevLett | https://doi.org/10.1103/physrevlett.108. |  | CrossRef: lenz 2012 — Contractile Units in Disordered Actomyosin Bundles Arise from F-Actin  |
| CHECK | Li2022_BiophysicalJournal | https://doi.org/10.1016/j.bpj.2022.09.03 |  | CrossRef: li 2022 — Network dynamics of the nonlinear power-law relaxation of cell cortex |
| CHECK | Lindstrom2010_Biomaterials | https://doi.org/10.1103/PhysRevE.82.0519 |  | CrossRef: lindström 2010 — Biopolymer network geometries: Characterization, regeneration, and ela |
| CHECK | Linsmeier2016_NatCommun | https://doi.org/10.1038/ncomms12615 |  | CrossRef: linsmeier 2016 — Disordered actomyosin networks are sufficient to produce cooperative a |
| CHECK | Lou2007_BiophysicalJournal | https://doi.org/10.1529/biophysj.106.097 |  | CrossRef: lou 2007 — A Structure-Based Sliding-Rebinding Mechanism for Catch Bonds |
| CHECK | Mangani2025_Cancers | https://doi.org/10.3390/cancers17213512 |  | CrossRef: mangani 2025 — Spheroid-Based 3D Models to Decode Cell Function and Matrix Effectors  |
| CHECK | Manibog2016_ProcNatlAcadSciUSA | https://doi.org/10.1073/pnas.1604012113 |  | CrossRef: manibog 2016 — Molecular determinants of cadherin ideal bond formation: Conformation- |
| CHECK | McFadden2017_PLoSComputBiol | https://doi.org/10.1371/journal.pcbi.100 |  | CrossRef: mcfadden 2017 — Filament turnover tunes both force generation and dissipation to contr |
| CHECK | MegeIshiyama_JBCReview | https://doi.org/10.1101/cshperspect.a028 |  | CrossRef: mège 2017 — Integration of Cadherin Adhesion and Cytoskeleton at
                  |
| CHECK | Miroshnikova2018_NatCellBiol | https://doi.org/10.1038/s41556-017-0005- |  | CrossRef: miroshnikova 2017 — Adhesion forces and cortical tension couple cell proliferation and dif |
| CHECK | Miyazaki2015_NatCellBiol | https://doi.org/10.1038/ncb3142 |  | CrossRef: miyazaki 2015 — Cell-sized spherical confinement induces the spontaneous formation of  |
| CHECK | Moazzeni2021_PhysRevE | https://doi.org/10.1103/physreve.103.032 |  | CrossRef: moazzeni 2021 — Single-cell mechanical analysis and tension quantification via electro |
| CHECK | Murrell2012_ProcNatlAcadSciUSA | https://doi.org/10.1073/pnas.1214753109 |  | CrossRef: murrell 2012 — F-actin buckling coordinates contractility and severing in a biomimeti |
| CHECK | Mykuliak2020_BiophysicalJournal | https://doi.org/10.1016/j.bpj.2020.07.03 |  | CrossRef: mykuliak 2020 — Mechanical Unfolding of Proteins—A Comparative Nonequilibrium Molecula |
| CHECK | Nayak2023_Cancers | https://doi.org/10.3390/cancers15194846 |  | CrossRef: nayak 2023 — Three-Dimensional In Vitro Tumor Spheroid Models for Evaluation of Ant |
| CHECK | Nieto2026_EuropeanJournalOfPharmac | https://doi.org/10.1016/j.ejps.2025.1073 |  | CrossRef: nieto 2026 — Integrating simulation and experimental validation of nutrient-limited |
| CHECK | Odenthal2013_PLoSComputBiol | https://doi.org/10.1371/journal.pcbi.100 |  | CrossRef: odenthal 2013 — Analysis of Initial Cell Spreading Using Mechanistic Contact Formulati |
| CHECK | Parsons_FAReview |  | 10.1038/sj.onc.1203877 | author-only match (year differs); candidate: parsons 2000 — Focal Adhesion Kinase: a regulator of focal adhesion dynamic |
| CHECK | PolacheckChen_ARBE |  | 10.1146/annurev-bioeng-110220-031722 | author-only match (year differs); candidate: chen 2022 — Current Developments and Challenges of mRNA Vaccines |
| CHECK | Popov2016_PLoSComputBiol | https://doi.org/10.1371/journal.pcbi.100 |  | CrossRef: popov 2016 — MEDYAN: Mechanochemical Simulations of Contraction and Polarity Alignm |
| CHECK | RocaCusachs_IntegrinReviews |  | 10.1016/j.isci.2020.100907 | author-only match (year differs); candidate: lerche 2020 — Integrin Binding Dynamics Modulate Ligand-Specific Mechanose |
| CHECK | Rodenhizer2018_AdvHealthcareMaterials | https://doi.org/10.1002/adhm.201701174 |  | CrossRef: rodenhizer 2018 — The Current Landscape of 3D In Vitro Tumor Models: What Cancer Hallmar |
| CHECK | Roffay2021_Development | https://doi.org/10.1242/dev.192773 |  | CrossRef: roffay 2021 — Inferring cell junction tension and pressure from cell geometry |
| CHECK | Rosenbauer2023_JPhysChemB | https://doi.org/10.1021/acs.jpcb.2c08114 |  | CrossRef: rosenbauer 2023 — Multiscale Modeling of Spheroid Tumors: Effect of Nutrient Availabilit |
| CHECK | Sakamoto2024_CellReportsPhysicalScien | https://doi.org/10.1016/j.xcrp.2024.1023 |  | CrossRef: sakamoto 2024 — Substrate geometry and topography induce F-actin reorganization and ch |
| CHECK | Samarage2015_DevelopmentalCell | https://doi.org/10.1016/j.devcel.2015.07 |  | CrossRef: samarage 2015 — Cortical Tension Allocates the First Inner Cells of the Mammalian Embr |
| CHECK | SciRep2025_GliomaAFM | https://doi.org/10.1038/s41598-025-04841 |  | CrossRef: masud 2025 — Exploring the heterogeneity in glioblastoma cellular mechanics using i |
| CHECK | SensPlastino_NRMCB |  | 10.1088/0953-8984/27/27/273103 | author-only match (year differs); candidate: sens 2015 — Membrane tension and cytoskeleton organization in cell motil |
| CHECK | Serwas2021_Jour | https://doi.org/10.1101/2021.06.28.45026 |  | CrossRef: serwas 2021 — Actin force generation in vesicle formation: mechanistic insights from |
| CHECK | Shah2025_Cells | https://doi.org/10.3390/cells14100732 |  | CrossRef: shah 2025 — Modeling Tumor Microenvironment Complexity In Vitro: Spheroids as Phys |
| CHECK | Slater2021_SoftMatter | https://doi.org/10.1039/d0sm01911a |  | CrossRef: slater 2021 — Transient mechanical interactions between cells and viscoelastic extra |
| CHECK | Stam2017_ProcNatlAcadSciUSA | https://doi.org/10.1073/pnas.1708625114 |  | CrossRef: stam 2017 — Filament rigidity and connectivity tune the deformation modes of activ |
| CHECK | Tam2021_Jour | https://doi.org/10.1101/2021.02.23.43258 |  | CrossRef: tam 2021 — Protein Friction and Filament Bending Facilitate Contraction of Disord |
| CHECK | Tao2019_Jour | https://doi.org/10.1101/847046 |  | CrossRef: tao 2019 — Tuning cell motility via cell tension with a mechanochemical cell migr |
| CHECK | Tao2020_BiophysicalJournal | https://doi.org/10.1016/j.bpj.2020.04.03 |  | CrossRef: tao 2020 — Tuning Cell Motility via Cell Tension with a Mechanochemical Cell Migr |
| CHECK | Thiticharoentam2026_Jour | https://doi.org/10.64898/2026.03.17.7125 |  | CrossRef: thiticharoentam 2026 — Computational Design for Engineering Layered Tissue Architectures via  |
| CHECK | TrepatSahai_NRC |  | 10.1038/ncb2548 | author-only match (year differs); candidate: friedl 2012 — Classifying collective cancer cell invasion |
| CHECK | Vallotton2003_BiophysJ |  | 10.70675/f72da841zbffcz4d68z9607zd6e9ada27f8d | author-only match (year differs); candidate: vallotton  — La décision publique et la crise |
| CHECK | Wang2016_PLoSONE | https://doi.org/10.1371/journal.pone.016 |  | CrossRef: wang 2016 — Spheroid Formation of Hepatocarcinoma Cells in Microwells: Experiments |
| CHECK | Wang2021_IEEESensorsJ | https://doi.org/10.1109/jsen.2020.304859 |  | CrossRef: wang 2021 — Quantification of Single-Cell Cortical Tension Using Multiple Constric |
| CHECK | Warmt2021_NewJPhys | https://doi.org/10.1088/1367-2630/ac254e |  | CrossRef: warmt 2021 — Differences in cortical contractile properties between healthy epithel |
| CHECK | Winklbauer2015_JournalOfCellScience | https://doi.org/10.1242/jcs.174623 |  | CrossRef: winklbauer 2015 — Cell adhesion strength from cortical tension – an integration of conce |
| CHECK | Xiong2007_NatPrec | https://doi.org/10.1038/npre.2007.62.2 |  | CrossRef: xiong 2007 — A three-dimensional stochastic spatio-temporal model of cell spreading |
| CHECK | Yang_NRMCB | https://doi.org/10.1038/s41580-020-0237- |  | CrossRef: yang 2020 — Guidelines and definitions for research on epithelial–mesenchymal tran |
| CHECK | Zhu2022_Organoids | https://doi.org/10.3390/organoids1020012 |  | CrossRef: zhu 2022 — 3D Tumor Spheroid and Organoid to Model Tumor Microenvironment for Can |
| CHECK | Zhuo2018_LightSciAppl | https://doi.org/10.1038/s41377-018-0001- |  | CrossRef: zhuo 2018 — Quantitative analysis of focal adhesion dynamics using photonic resona |
| CHECK | deGennes_LiquidCrystals |  | 10.1093/oso/9780198520245.001.0001 | author-only match (year differs); candidate: gennes 1993 — The Physics of Liquid Crystals |
| NO_DOI_FOUND | AlonsoMatilla2023_BiophysJ |  | 10.1101/2023.11.13.566882 | candidate: alonso-matilla 2023 — Cell intrinsic mechanical regulation of plasma membrane accu |
| NO_DOI_FOUND | Bangasser2017_NCB |  | 10.3791/56219 | candidate: bangasser 2017 — Touchscreen Sustained Attention Task (SAT) for Rats |
| NO_DOI_FOUND | Beaune2014_PNAS |  | 10.4000/books.pur.49676 | candidate: beaune 2014 — Conclusions |
| NO_DOI_FOUND | Bobrowska2021_Cells |  | 10.31338/uw.9788323548799.pp.90-102 | candidate: bobrowska 2021 — Retoryka współczesnej estetyki i sztuki: figura i narracja w |
| NO_DOI_FOUND | Cai2014_Cell |  | 10.1201/b17354 | candidate: cai 2014 — Ambient Diagnostics |
| NO_DOI_FOUND | Cai2022_FCDB |  | 10.3389/fcell.2021.817104 | candidate: cai 2022 — Interplay Between Iron Overload and Osteoarthritis: Clinical |
| NO_DOI_FOUND | Carlsson2018_JPhysCondMat |  | 10.2210/pdb5v7f/pdb | candidate: carlsson 2018 — T4 lysozyme Y18Ymi |
| NO_DOI_FOUND | CavalcantiAdam2007_BiophysJ |  | 10.1039/b614008d | candidate: girard 2007 — Cellular chemomechanics at interfaces: sensing, integration  |
| NO_DOI_FOUND | Conklin2011 |  | 10.1515/9781575066288 | candidate: conklin 2011 — Oath Formulas in Biblical Hebrew |
| NO_DOI_FOUND | Doss2020_PNAS |  | 10.1201/9780429344206 | candidate: ramamoorthy 2020 — Biology, Chemistry, and Applications of Apocarotenoids |
| NO_DOI_FOUND | EloseguiArtola2016_NatMater |  | 10.1016/j.bpj.2017.05.020 | candidate: elosegui-artola 2017 — Amoebae as Mechanosensitive Tanks |
| NO_DOI_FOUND | Fritzsche2013_MBC |  | 10.1007/978-3-642-37495-1_2 | candidate: fritzsche 2013 — Lebesgue-Theorie |
| NO_DOI_FOUND | Fu2024_PNAS |  | 10.20944/preprints202406.1065.v2 | candidate: fu 2024 — An Introduction to Cosmos Thermodynamics |
| NO_DOI_FOUND | Guo2017_PNAS |  | 10.1145/3145690.3145698 | candidate: guo 2017 — Importance sampling measured BRDFs based on second order sph |
| NO_DOI_FOUND | Han2017_PNAS |  | 10.2172/1477879 | candidate: han 2017 — Closeout Report for CTEQ Summer School 2017 |
| NO_DOI_FOUND | Han2021_eLife |  | 10.4211/hs.7c5e032bdc7648a4a8c863a2c175e5b6 | candidate: han 2021 — DataShare for Han et al., 2021 Ecological Indicators |
| NO_DOI_FOUND | Helfrich1973 |  | 10.1515/znc-1973-11-1209 | candidate: helfrich 1973 — Elastic Properties of Lipid Bilayers: Theory and Possible Ex |
| NO_DOI_FOUND | Jansen2018_BiophysJ |  | 10.1016/j.jas.2018.02.016 | candidate: jansen 2018 — On the use of Cu isotope signatures in archaeometallurgy: A  |
| NO_DOI_FOUND | Kim2021_PNAS |  | 10.1145/3450507.3457434 | candidate: lee 2021 — Isle of reflections |
| NO_DOI_FOUND | Kothari2018_JApplMech |  | 10.2514/6.2018-0530 | candidate: rustagi 2018 — Gyroscopic Stabilization of Flying Wing Aircraft |
| NO_DOI_FOUND | Lindstrom2013_SoftMatter |  | 10.1093/obo/9780199766567-0108 | candidate: lindstrom 2013 — Cargo Cults |
| NO_DOI_FOUND | Liu2019_ProstateMech |  | 10.1007/978-981-13-6962-9 | candidate: liu 2019 — Deuteride Materials |
| NO_DOI_FOUND | Liu2025_NatPhysics |  | 10.22541/au.176275808.84444425/v1 | candidate: liu 2025 — Comment on Bandyopadhyay et al. |
| NO_DOI_FOUND | Maitre2012_Nature |  | 10.3917/lav.daler.2012.01.0096 | candidate: bubrovszky 2012 — Pathologies schizophréniques |
| NO_DOI_FOUND | Maitre2015_NCB |  | 10.1051/shsconf/20152000001 | candidate: colón de carvajal 2015 — Préface |
| NO_DOI_FOUND | Malinova2021_NatCommun |  | 10.34076/20713797_2021_2_29 | candidate: malinova 2021 — «Interest» as a legal concept: problems of doctrinal definit |
| NO_DOI_FOUND | MarantanMahadevan2018_AmJPhys |  | 10.1119/1.5003376 | candidate: marantan 2018 — Mechanics and statistics of the worm-like chain |
| NO_DOI_FOUND | Mogilner2002_2003 |  | 10.1002/3527601503.ch14 | candidate: scholey 2002 — Mitotic Spindle Motors |
| NO_DOI_FOUND | Munster2013_PNAS |  | 10.7551/mitpress/8982.001.0001 | candidate: munster 2013 — An Aesthesia of Networks |
| NO_DOI_FOUND | Odijk1995 |  | 10.1037/e493132004-001 | candidate: odijk 1995 — A syntactic condition on definite descriptions |
| NO_DOI_FOUND | Ozawa2020_JCB |  | 10.1007/978-981-15-4190-2_7 | candidate: ozawa 2020 — Comprehensive Registry in Japan |
| NO_DOI_FOUND | Pereverzev2005_BiophysJ |  | 10.1007/s11018-005-0098-9 | candidate: pereverzev 2005 — Measurement of the spectral density of broadband radio inter |
| NO_DOI_FOUND | PerezGonzalez2019_NatPhys |  | 10.2172/1597056 | candidate: perez-gonzalez 2019 — Dirac and Majorana Neutrino Signatures of Primordial Black H |
| NO_DOI_FOUND | Pronk2008_PNAS |  | 10.1057/dev.2008.28 | candidate: pronk 2008 — Climate, Scarcities and Development |
| NO_DOI_FOUND | Ray2021_COCellBiol |  | 10.1016/j.ceb.2021.05.004 | candidate: ray 2021 — Aligned forces: Origins and mechanisms of cancer disseminati |
| NO_DOI_FOUND | Steinberg1996_CurrTopDevBiol |  | 10.1387/ijdb.8877443 | candidate: drawbridge 1996 — Morphogenesis of the axolotl pronephric duct: a model system |
| NO_DOI_FOUND | TapiaRojo2019_SciAdv |  | 10.1101/816801 | candidate: tapia-rojo 2019 — Thermal versus Mechanical Unfolding in a Model Protein |
| NO_DOI_FOUND | Trichet2012_PNAS |  | 10.3917/bupsy.520.0365 | candidate: trichet 2012 — La notion de kakon . Histoire et enjeux psychopathologiques |
| NO_DOI_FOUND | Varma2016_BiophysJ |  | 10.4135/9781526429629 | candidate: varma 2016 — Swiss Roll |
| NO_DOI_FOUND | Vassalli2023_Cancers |  | 10.3917/mem.086.0017 | candidate: vassalli 2023 — L’hospitalité citoyenne : bien plus qu’un hébergement |
| NO_DOI_FOUND | Venturini2020_Science |  | 10.32614/cran.package.dmbc | candidate: venturini 2020 — dmbc: Model Based Clustering of Binary Dissimilarity Measure |
| NO_DOI_FOUND | Wang1997_BiophysJ |  | 10.1115/97-aa-050 | candidate: wang 1997 — Vibration of Skew Sandwich Plates With Laminated Facings |
| NO_DOI_FOUND | Wei2024_ActaBiomater |  | 10.26434/chemrxiv-2023-8c9vh-v2 | candidate: wei 2024 — Persistent Topological Laplacians  -- a Survey |
| NO_DOI_FOUND | Wenger2007 |  | 10.1097/mlr.0b013e31815b97bf | candidate: wenger 2007 — Comorbidity and Quality of Care: Clarification From Author |
| NO_DOI_FOUND | WirtzFriedl2011 |  | 10.1038/nrc3080 | candidate: wirtz 2011 — The physics of cancer: the role of physical interactions and |
| NO_DOI_FOUND | Yao2014_NatCommun |  | 10.1038/ncomms6114 | candidate: yao 2014 — Functional annotation of colon cancer risk SNPs |
| OK | Adar2020_JCS | https://doi.org/10.1073/pnas.1918203117 |  | CrossRef: adar 2020 — Active volume regulation in adhered cells |
| OK | Andreu2021_NatCommun | https://doi.org/10.1038/s41467-021-24383 |  | CrossRef: andreu 2021 — The force loading rate drives cell mechanosensing through both reinfor |
| OK | BakerChen2012_JCS | https://doi.org/10.1242/jcs.079509 |  | CrossRef: baker 2012 — Deconstructing the third dimension – how 3D culture microenvironments  |
| OK | Banerjee2015_PRL | https://doi.org/10.1103/physrevlett.114. |  | CrossRef: banerjee 2015 — Propagating Stress Waves During Epithelial Expansion |
| OK | Bangasser2013_BiophysJ | https://doi.org/10.1016/j.bpj.2013.06.02 |  | CrossRef: bangasser 2013 — Determinants of Maximal Force Transmission in a Motor-Clutch Model of  |
| OK | Bastatas2012_BBA | https://doi.org/10.1016/j.bbagen.2012.02 |  | CrossRef: bastatas 2012 — AFM nano-mechanics and calcium dynamics of prostate cancer cells with  |
| OK | Bell1978_Science | https://doi.org/10.1126/science.347575 |  | CrossRef: bell 1978 — Models for the Specific Adhesion of Cells to Cells |
| OK | Bertet2004_Nature | https://doi.org/10.1038/nature02590 |  | CrossRef: bertet 2004 — Myosin-dependent junction remodelling controls planar cell intercalati |
| OK | Bi2015_PRX | https://doi.org/10.1103/physrevx.6.02101 |  | CrossRef: bi 2016 — Motility-Driven Glass and Jamming Transitions in Biological Tissues |
| OK | Bi2016_NatPhys | https://doi.org/10.1038/nphys3471 |  | CrossRef: bi 2015 — A density-independent rigidity transition in biological tissues |
| OK | Bi2016_NatPhys | https://doi.org/10.1038/nphys3471 |  | CrossRef: bi 2015 — A density-independent rigidity transition in biological tissues |
| OK | Bieling2016_Cell | https://doi.org/10.1016/j.cell.2015.11.0 |  | CrossRef: bieling 2016 — Force Feedback Controls Motor Activity and Mechanical Properties of Se |
| OK | Bieling2016_Cell | https://doi.org/10.1016/j.cell.2015.11.0 |  | CrossRef: bieling 2016 — Force Feedback Controls Motor Activity and Mechanical Properties of Se |
| OK | Borghi2012_PNAS | https://doi.org/10.1073/pnas.1204390109 |  | CrossRef: borghi 2012 — E-cadherin is under constitutive actomyosin-generated tension that is  |
| OK | Bredfeldt2014_JBiomedOpt | https://doi.org/10.1117/1.JBO.19.1.01600 |  | CrossRef: bredfeldt 2014 — Computational segmentation of collagen fibers from second-harmonic gen |
| OK | BrochardWyart2006_PNAS | https://doi.org/10.1073/pnas.0602012103 |  | CrossRef: brochard-wyart 2006 — Hydrodynamic narrowing of tubes extruded from cells |
| OK | Brodland2002 | https://doi.org/10.1115/1.1449491 |  | CrossRef: brodland 2002 — The Differential Interfacial Tension Hypothesis (DITH): A Comprehensiv |
| OK | BroederszMacKintosh2014_RMP | https://doi.org/10.1103/revmodphys.86.99 |  | CrossRef: broedersz 2014 — Modeling semiflexible polymer networks |
| OK | Buckley2014_Science | https://doi.org/10.1126/science.1254211 |  | CrossRef: buckley 2014 — The minimal cadherin-catenin complex binds to actin filaments under fo |
| OK | Bustamante1994_Science | https://doi.org/10.1126/science.8079175 |  | CrossRef: bustamante 1994 — Entropic Elasticity of λ-Phage DNA |
| OK | Cadart2019_NatPhysics | https://doi.org/10.1038/s41567-019-0629- |  | CrossRef: cadart 2019 — The physics of cell-size regulation across timescales |
| OK | Caille2002_JBiomech | https://doi.org/10.1016/S0021-9290(01)00 |  | CrossRef: caille 2002 — Contribution of the nucleus to the mechanical properties of endothelia |
| OK | Cambria2024_NatRevCancer | https://doi.org/10.1038/s41568-023-00656 |  | CrossRef: cambria 2024 — Linking cell mechanical memory and cancer metastasis |
| OK | Cavey2008_Cell | https://doi.org/10.1038/nature06953 |  | CrossRef: cavey 2008 — A two-tiered mechanism for stabilization and immobilization of E-cadhe |
| OK | ChanOdde2008_Science | https://doi.org/10.1126/science.1163595 |  | CrossRef: chan 2008 — Traction Dynamics of Filopodia on Compliant Substrates |
| OK | Charras2008_BJ | https://doi.org/10.1529/biophysj.107.113 |  | CrossRef: charras 2008 — Life and Times of a Cellular Bleb |
| OK | Chaudhuri2015_NatCommun | https://doi.org/10.1038/ncomms7365 |  | CrossRef: chaudhuri 2015 — Substrate stress relaxation regulates cell spreading |
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
| OK | Douezan2011_PNAS | https://doi.org/10.1073/pnas.1018057108 |  | CrossRef: douezan 2011 — Spreading dynamics and wetting transition of cellular aggregates |
| OK | DoyleYamada2009_JCB | https://doi.org/10.1083/jcb.200810041 |  | CrossRef: doyle 2009 — One-dimensional topography underlies three-dimensional fibrillar cell  |
| OK | Dupont2011_Nature | https://doi.org/10.1038/nature10137 |  | CrossRef: dupont 2011 — Role of YAP/TAZ in mechanotransduction |
| OK | EloseguiArtola2017_Cell | https://doi.org/10.1016/j.cell.2017.10.0 |  | CrossRef: elosegui-artola 2017 — Force Triggers YAP Nuclear Entry by Regulating Transport across Nuclea |
| OK | Engler2006_Cell | https://doi.org/10.1016/j.cell.2006.06.0 |  | CrossRef: engler 2006 — Matrix Elasticity Directs Stem Cell Lineage Specification |
| OK | Evans2007_RBP | https://doi.org/10.1353/scu.2007.0015 |  | CrossRef: ferris 2007 — Walker Evans, 1974 |
| OK | EvansRitchie1997_BiophysJ | https://doi.org/10.1016/s0006-3495(97)78 |  | CrossRef: evans 1997 — Dynamic strength of molecular adhesion bonds |
| OK | FarooquiFenteany2005_JCS | https://doi.org/10.1242/jcs.01577 |  | CrossRef: farooqui 2005 — Multiple rows of cells behind an epithelial wound edge extend cryptic  |
| OK | FengSun2021_FCDB | https://doi.org/10.3389/fcell.2021.71883 |  | CrossRef: fan 2021 — Substrate Stiffness Modulates the Growth, Phenotype, and Chemoresistan |
| OK | Ferrer2008_PNAS | https://doi.org/10.1073/pnas.0706124105 |  | CrossRef: ferrer 2008 — Measuring molecular rupture forces between single actin filaments and  |
| OK | Foty1996_PNAS | https://doi.org/10.1242/dev.122.5.1611 |  | CrossRef: foty 1996 — Surface tensions of embryonic tissues predict their mutual envelopment |
| OK | Franz2023_NatCommun | https://doi.org/10.5194/egusphere-2024-2 |  | CrossRef: franz 2024 — Reply to Comment on Franz et al. (2023): A reinterpretation of the 1.5 |
| OK | Fredberg2015_FASEB | https://doi.org/10.1096/fasebj.29.1_supp |  | CrossRef: krishnan 2015 — Force Chains And Gap Formation in Thrombin‐induced Endothelial Permeab |
| OK | GalWeihs2012_CBB | https://doi.org/10.1007/s12013-012-9356- |  | CrossRef: gal 2012 — Intracellular Mechanics and Activity of Breast Cancer Cells Correlate  |
| OK | Garcia2015_PNAS | https://doi.org/10.1073/pnas.1510973112 |  | CrossRef: garcia 2015 — Physics of active jamming during collective cellular motion in a monol |
| OK | Gittes1993_JCB | https://doi.org/10.1083/jcb.120.4.923 |  | CrossRef: gittes 1993 — Flexural rigidity of microtubules and actin filaments measured from th |
| OK | Goldmann2002 | https://doi.org/10.1006/cbir.2002.0900 |  | CrossRef: goldmann 2002 — p56<sup>lck</sup> CONTROLS PHOSPHORYLATION OF FILAMIN (ABP‐280) AND RE |
| OK | GonzalezRodriguez2012_Science | https://doi.org/10.1126/science.1226418 |  | CrossRef: gonzalez-rodriguez 2012 — Soft Matter Models of Developing Tissues and Tumors |
| OK | Guilak2000_BBRC | https://doi.org/10.1006/bbrc.2000.2360 |  | CrossRef: guilak 2000 — Viscoelastic Properties of the Cell Nucleus |
| OK | Guo2020_CellRegen | https://doi.org/10.1186/s13619-020-00054 |  | CrossRef: guo 2020 — Consistent apparent Young’s modulus of human embryonic stem cells and  |
| OK | Hammerick2010_TissueEngA | https://doi.org/10.1089/ten.tea.2010.021 |  | CrossRef: hammerick 2011 — Elastic Properties of Induced Pluripotent Stem Cells |
| OK | Hannezo2014_PNAS | https://doi.org/10.1073/pnas.1312076111 |  | CrossRef: hannezo 2013 — Theory of epithelial sheet morphology in three dimensions |
| OK | Harada2014_JCB | https://doi.org/10.1083/jcb.201308029 |  | CrossRef: harada 2014 — Nuclear lamin stiffness is a barrier to 3D migration, but softness can |
| OK | HeadLevineMacKintosh2003_PRE | https://doi.org/10.1103/physrevlett.91.1 |  | CrossRef: head 2003 — Deformation of Cross-Linked Semiflexible Polymer Networks |
| OK | HeadLevineMacKintosh2003_PRL | https://doi.org/10.1103/physrevlett.91.1 |  | CrossRef: head 2003 — Deformation of Cross-Linked Semiflexible Polymer Networks |
| OK | Hochmuth1996_BiophysJ | https://doi.org/10.1016/S0006-3495(96)79 |  | CrossRef: hochmuth 1996 — Deformation and flow of membrane into tethers extracted from neuronal  |
| OK | Hosseini2020_AdvSci | https://doi.org/10.1002/advs.202001276 |  | CrossRef: hosseini 2020 — EMT‐Induced Cell‐Mechanical Changes Enhance Mitotic Rounding Strength |
| OK | Hosseini2021_BiophysJ | https://doi.org/10.1016/j.bpj.2021.05.00 |  | CrossRef: hosseini 2021 — EMT changes actin cortex rheology in a cell-cycle-dependent manner |
| OK | Hynes2002_Cell | https://doi.org/10.1016/s0092-8674(02)00 |  | CrossRef: hynes 2002 — Integrins |
| OK | Isambert1995_JBC | https://doi.org/10.1074/jbc.270.19.11437 |  | CrossRef: isambert 1995 — Flexibility of Actin Filaments Derived from Thermal Fluctuations |
| OK | Janmey2019_PhysiolRev | https://doi.org/10.1152/physrev.00013.20 |  | CrossRef: janmey 2020 — Stiffness Sensing by Cells |
| OK | Jetta2023_FCDB | https://doi.org/10.3389/fcell.2023.11981 |  | CrossRef: jetta 2023 — Epithelial cells sense local stiffness via Piezo1 mediated cytoskeleta |
| OK | Kage2017_NatCommun | https://doi.org/10.1038/ncomms14832 |  | CrossRef: kage 2017 — FMNL formins boost lamellipodial force generation |
| OK | Kanchanawong2010_Nature | https://doi.org/10.1038/nature09621 |  | CrossRef: kanchanawong 2010 — Nanoscale architecture of integrin-based cell adhesions |
| OK | Kong2009_Nature | https://doi.org/10.1083/jcb.200810002 |  | CrossRef: kong 2009 — Demonstration of catch bonds between an integrin and its ligand |
| OK | KraningRush2012_PLoSOne | https://doi.org/10.1371/journal.pone.003 |  | CrossRef: kraning-rush 2012 — Cellular Traction Stresses Increase with Increasing Metastatic Potenti |
| OK | Krause2014_NRMCB | https://doi.org/10.1038/nrm3861 |  | CrossRef: krause 2014 — Steering cell migration: lamellipodium dynamics and the regulation of  |
| OK | Kunda2008_CurrBiol | https://doi.org/10.1016/j.cub.2007.12.05 |  | CrossRef: kunda 2008 — Moesin Controls Cortical Rigidity, Cell Rounding, and Spindle Morphoge |
| OK | Labernadie2017_NCB | https://doi.org/10.1038/ncb3478 |  | CrossRef: labernadie 2017 — A mechanically active heterotypic E-cadherin/N-cadherin adhesion enabl |
| OK | Lammerding2004_JCI | https://doi.org/10.1172/JCI19670 |  | CrossRef: lammerding 2004 — Lamin A/C deficiency causes defective nuclear mechanics and mechanotra |
| OK | LeDuc2010_JCB | https://doi.org/10.1083/jcb.201001149 |  | CrossRef: le duc 2010 — Vinculin potentiates E-cadherin mechanosensing and is recruited to act |
| OK | LeimkuhlerMatthews2013 | https://doi.org/10.1063/1.4802990 |  | CrossRef: leimkuhler 2013 — Robust and efficient configurational molecular sampling via Langevin d |
| OK | Levental2009_Cell | https://doi.org/10.1016/j.cell.2009.10.0 |  | CrossRef: levental 2009 — Matrix Crosslinking Forces Tumor Progression by Enhancing Integrin Sig |
| OK | LewisGrandl2015_eLife | https://doi.org/10.7554/elife.12088 |  | CrossRef: lewis 2015 — Mechanical sensitivity of Piezo1 ion channels can be tuned by cellular |
| OK | Li2008_BBRC | https://doi.org/10.1016/j.bbrc.2008.07.0 |  | CrossRef: li 2008 — AFM indentation study of breast cancer cells |
| OK | Licup2015_PNAS | https://doi.org/10.1073/pnas.1504258112 |  | CrossRef: licup 2015 — Stress controls the mechanics of collagen networks |
| OK | Liew2024_CMBE | https://doi.org/10.1007/s12195-024-00811 |  | CrossRef: liew 2024 — Cellular Traction Force Holds the Potential as a Drug Testing Readout  |
| OK | Liu2015_Cell | https://doi.org/10.1016/j.cell.2015.01.0 |  | CrossRef: liu 2015 — Confinement and Low Adhesion Induce Fast Amoeboid Migration of Slow Me |
| OK | Lo2000_BiophysJ | https://doi.org/10.1016/S0006-3495(00)76 |  | CrossRef: lo 2000 — Cell Movement Is Guided by the Rigidity of the Substrate |
| OK | Lomakin2020_Science | https://doi.org/10.1126/science.aba2894 |  | CrossRef: lomakin 2020 — The nucleus acts as a ruler tailoring cell responses to spatial constr |
| OK | Malinverno2016_NatMater | https://doi.org/10.1002/2016gl070096 |  | CrossRef: nole 2016 — Short‐range, overpressure‐driven methane migration in coarse‐grained g |
| OK | MalyBorisy2001_PNAS | https://doi.org/10.1073/pnas.181338798 |  | CrossRef: maly 2001 — Self-organization of a propulsive actin network as an evolutionary pro |
| OK | Manibog2014_NatCommun | https://doi.org/10.1038/ncomms4941 |  | CrossRef: manibog 2014 — Resolving the molecular mechanism of cadherin catch bond formation |
| OK | Marchetti2013_RMP | https://doi.org/10.1103/revmodphys.85.11 |  | CrossRef: marchetti 2013 — Hydrodynamics of soft active matter |
| OK | MarkoSiggia1995_Macromolecules | https://doi.org/10.1021/ma00130a008 |  | CrossRef: marko 1995 — Stretching DNA |
| OK | Masud2025_SciRep | https://doi.org/10.1038/s41598-025-04841 |  | CrossRef: masud 2025 — Exploring the heterogeneity in glioblastoma cellular mechanics using i |
| OK | MattilaLappalainen2008_NRMCB | https://doi.org/10.1038/nrm2406 |  | CrossRef: mattila 2008 — Filopodia: molecular architecture and cellular functions |
| OK | McKenzie2018_SciRep | https://doi.org/10.1038/s41598-018-25589 |  | CrossRef: mckenzie 2018 — The mechanical microenvironment regulates ovarian cancer cell morpholo |
| OK | Moeendarbary2013_NatMater | https://doi.org/10.1038/nmat3517 |  | CrossRef: moeendarbary 2013 — The cytoplasm of living cells behaves as a poroelastic material |
| OK | Mogilner2005_BiophysJ | https://doi.org/10.1529/biophysj.104.056 |  | CrossRef: mogilner 2005 — The Physics of Filopodial Protrusion |
| OK | MogilnerOster1996_BiophysJ | https://doi.org/10.1016/S0006-3495(96)79 |  | CrossRef: mogilner 1996 — Cell motility driven by actin polymerization |
| OK | MogilnerOster2003_BJ | https://doi.org/10.1016/s0006-3495(03)74 |  | CrossRef: mogilner 2003 — Force Generation by Actin Polymerization II: The Elastic Ratchet and T |
| OK | Molter2022_FCDB | https://doi.org/10.3389/fcell.2022.93251 |  | CrossRef: molter 2022 — Prostate cancer cells of increasing metastatic potential exhibit diver |
| OK | MotteKaufman2013_Biopolymers | https://doi.org/10.1002/bip.22133 |  | CrossRef: motte 2012 — Strain stiffening in collagen I networks |
| OK | Mueller2017_Cell | https://doi.org/10.1016/j.cell.2017.07.0 |  | CrossRef: mueller 2017 — Load Adaptation of Lamellipodial Actin Networks |
| OK | Mui2016_PNAS | https://doi.org/10.1242/jcs.183699 |  | CrossRef: mui 2016 — The mechanical regulation of integrin–cadherin crosstalk organizes cel |
| OK | Murrell2015_NRMCB | https://doi.org/10.1038/nrm4012 |  | CrossRef: murrell 2015 — Forcing cells into shape: the mechanics of actomyosin contractility |
| OK | Murrell2015_NRMCB | https://doi.org/10.1038/nrm4012 |  | CrossRef: murrell 2015 — Forcing cells into shape: the mechanics of actomyosin contractility |
| OK | Nakamura2007 | https://doi.org/10.1083/jcb.200707073 |  | CrossRef: nakamura 2007 — Structural basis of filamin A functions |
| OK | Nam2016_PNAS | https://doi.org/10.1073/pnas.1523906113 |  | CrossRef: nam 2016 — Strain-enhanced stress relaxation impacts nonlinear elasticity in coll |
| OK | Notbohm2016_BiophysJ | https://doi.org/10.1016/j.bpj.2016.05.01 |  | CrossRef: notbohm 2016 — Cellular Contraction and Polarization Drive Collective Cellular Motion |
| OK | Ofek2009_JBiomech | https://doi.org/10.1016/j.jbiomech.2009. |  | CrossRef: ofek 2009 — In situ mechanical properties of the chondrocyte cytoplasm and nucleus |
| OK | Omidvar2014_JBiomech | https://doi.org/10.1016/j.jbiomech.2014. |  | CrossRef: omidvar 2014 — Atomic force microscope-based single cell force spectroscopy of breast |
| OK | Oria2017_Nature | https://doi.org/10.1038/nature24662 |  | CrossRef: oria 2017 — Force loading explains spatial sensing of ligands by cells |
| OK | Otto2015_NatMethods | https://doi.org/10.1038/nmeth.3281 |  | CrossRef: otto 2015 — Real-time deformability cytometry: on-the-fly cell mechanical phenotyp |
| OK | Pajerowski2007_PNAS | https://doi.org/10.1073/pnas.0702576104 |  | CrossRef: pajerowski 2007 — Physical plasticity of the nucleus in stem cell differentiation |
| OK | Palchesko2012_PLoSOne | https://doi.org/10.1371/journal.pone.005 |  | CrossRef: palchesko 2012 — Development of Polydimethylsiloxane Substrates with Tunable Elastic Mo |
| OK | Panzetta2019_PNAS | https://doi.org/10.1073/pnas.1904660116 |  | CrossRef: panzetta 2019 — Cell mechanosensing is regulated by substrate strain energy rather tha |
| OK | Park2015_NatMater | https://doi.org/10.1038/nmat4357 |  | CrossRef: park 2015 — Unjamming and cell shape in the asthmatic airway epithelium |
| OK | Park2015_NatMater | https://doi.org/10.1038/nmat4357 |  | CrossRef: park 2015 — Unjamming and cell shape in the asthmatic airway epithelium |
| OK | PathakKumar2012_PNAS | https://doi.org/10.1073/pnas.1118073109 |  | CrossRef: pathak 2012 — Independent regulation of tumor cell migration by matrix stiffness and |
| OK | PaulKonstantopoulos2017_NatRevCancer | https://doi.org/10.1038/nrc.2016.123 |  | CrossRef: paul 2016 — Cancer cell motility: lessons from migration in confined spaces |
| OK | Petitjean2010_BiophysJ | https://doi.org/10.1016/j.bpj.2010.01.03 |  | CrossRef: petitjean 2010 — Velocity Fields in a Collectively Migrating Epithelium |
| OK | Pillarisetti2011_CellReprogram | https://doi.org/10.1089/cell.2011.0028 |  | CrossRef: pillarisetti 2011 — Mechanical Phenotyping of Mouse Embryonic Stem Cells: Increase in Stif |
| OK | Plodinec2012_NatNanotechnol | https://doi.org/10.1038/nnano.2012.167 |  | CrossRef: plodinec 2012 — The nanomechanical signature of breast cancer |
| OK | Plotnikov2012_Cell | https://doi.org/10.1016/j.cell.2012.11.0 |  | CrossRef: plotnikov 2012 — Force Fluctuations within Focal Adhesions Mediate ECM-Rigidity Sensing |
| OK | PollardBorisy2003_Cell | https://doi.org/10.1016/s0092-8674(03)00 |  | CrossRef: pollard 2003 — Cellular Motility Driven by Assembly and Disassembly of Actin Filament |
| OK | PontesGauthier2017_SCDB | https://doi.org/10.1016/j.semcdb.2017.08 |  | CrossRef: pontes 2017 — Membrane tension: A challenging but universal physical parameter in ce |
| OK | Proestaki2019_ExpMech | https://doi.org/10.1007/s11340-018-00453 |  | CrossRef: proestaki 2019 — Modulus of Fibrous Collagen at the Length Scale of a Cell |
| OK | Prost2015_NatPhys | https://doi.org/10.1038/nphys3224 |  | CrossRef: prost 2015 — Active gel physics |
| OK | Provenzano2006_BMCMed | https://doi.org/10.1186/1741-7015-4-38 |  | CrossRef: provenzano 2006 — Collagen reorganization at the tumor-stromal interface facilitates loc |
| OK | Rakshit2012_PNAS | https://doi.org/10.1073/pnas.1208349109 |  | CrossRef: rakshit 2012 — Ideal, catch, and slip bonds in cadherin adhesion |
| OK | Rawicz2000_BiophysJ | https://doi.org/10.1016/S0006-3495(00)76 |  | CrossRef: rawicz 2000 — Effect of Chain Length and Unsaturation on Elasticity of Lipid Bilayer |
| OK | Ray2017_NatCommun | https://doi.org/10.1038/ncomms14923 |  | CrossRef: ray 2017 — Anisotropic forces from spatially constrained focal adhesions mediate  |
| OK | Reffay2014_NCB | https://doi.org/10.1038/ncb2917 |  | CrossRef: reffay 2014 — Interplay of RhoA and mechanical forces in collective cell migration d |
| OK | Riahi2015_NatCommun | https://doi.org/10.1038/ncomms7556 |  | CrossRef: riahi 2015 — Notch1–Dll4 signalling and mechanical force regulate leader cell forma |
| OK | Rosenbluth2006_BiophysJ | https://doi.org/10.1529/biophysj.105.067 |  | CrossRef: rosenbluth 2006 — Force Microscopy of Nonadherent Cells: A Comparison of Leukemia Cell D |
| OK | Ruprecht2015_Cell | https://doi.org/10.1016/j.cell.2015.01.0 |  | CrossRef: ruprecht 2015 — Cortical Contractility Triggers a Stochastic Switch to Fast Amoeboid C |
| OK | Ruscone2024_Bioinformatics | https://doi.org/10.46471/gigabyte.136 |  | CrossRef: noël 2024 — PhysiMeSS - a new physiCell addon for extracellular matrix modelling |
| OK | Saarinen2015_PLoSOne | https://doi.org/10.1371/journal.pone.014 |  | CrossRef: saarinen 2015 — Differential Predictive Roles of A- and B-Type Nuclear Lamins in Prost |
| OK | Salbreux2012_TCB | https://doi.org/10.1016/j.tcb.2012.07.00 |  | CrossRef: salbreux 2012 — Actin cortex mechanics and cellular morphogenesis |
| OK | Salbreux2012_TCB | https://doi.org/10.1016/j.tcb.2012.07.00 |  | CrossRef: salbreux 2012 — Actin cortex mechanics and cellular morphogenesis |
| OK | Saraswathibhatla2023_NRMCB | https://doi.org/10.1038/s41580-023-00583 |  | CrossRef: saraswathibhatla 2023 — Cell–extracellular matrix mechanotransduction in 3D |
| OK | SchwarzSafran2013_RMP | https://doi.org/10.1103/revmodphys.85.13 |  | CrossRef: schwarz 2013 — Physics of adherent cells |
| OK | SerraPicamal2012_NatPhys | https://doi.org/10.1038/nphys2355 |  | CrossRef: serra-picamal 2012 — Mechanical waves during tissue expansion |
| OK | Sharma2014_IntegrBiol | https://doi.org/10.1039/c3ib40246k |  | CrossRef: sharma 2014 — The role of Rho GTPase in cell stiffness and cisplatin resistance in o |
| OK | Smelser2015_BMMB | https://doi.org/10.1007/s10237-015-0677- |  | CrossRef: smelser 2015 — Mechanical properties of normal versus cancerous breast cells |
| OK | Stephens2017_MBoC | https://doi.org/10.1091/mbc.e16-09-0653 |  | CrossRef: stephens 2017 — Chromatin and lamin A determine two different mechanical response regi |
| OK | Stewart2011_Nature | https://doi.org/10.1038/nature09642 |  | CrossRef: stewart 2011 — Hydrostatic pressure and the actomyosin cortex drive mitotic cell roun |
| OK | Storm2005_Nature | https://doi.org/10.1038/nature03521 |  | CrossRef: storm 2005 — Nonlinear elasticity in biological gels |
| OK | Stroka2014_Cell | https://doi.org/10.1016/j.cell.2014.02.0 |  | CrossRef: stroka 2014 — Water Permeation Drives Tumor Cell Migration in Confined Microenvironm |
| OK | Stylianopoulos2012_PNAS | https://doi.org/10.1073/pnas.1213353109 |  | CrossRef: stylianopoulos 2012 — Causes, consequences, and remedies for growth-induced solid stress in  |
| OK | Swaminathan2011_CancerRes | https://doi.org/10.1158/0008-5472.CAN-11 |  | CrossRef: swaminathan 2011 — Mechanical Stiffness Grades Metastatic Potential in Patient Tumor Cell |
| OK | Swift2013_Science | https://doi.org/10.1126/science.1240104 |  | CrossRef: swift 2013 — Nuclear Lamin-A Scales with Tissue Stiffness and Enhances Matrix-Direc |
| OK | Tambe2011_NatMater | https://doi.org/10.1038/nmat3025 |  | CrossRef: tambe 2011 — Collective cell guidance by cooperative intercellular forces |
| OK | Thiery2009_Cell | https://doi.org/10.1016/j.cell.2009.11.0 |  | CrossRef: thiery 2009 — Epithelial-Mesenchymal Transitions in Development and Disease |
| OK | Tinevez2009_PNAS | https://doi.org/10.1073/pnas.0903353106 |  | CrossRef: tinevez 2009 — Role of cortical tension in bleb growth |
| OK | TitushkinCho2007_BiophysJ | https://doi.org/10.1529/biophysj.107.107 |  | CrossRef: titushkin 2007 — Modulation of Cellular Mechanics during Osteogenic Differentiation of  |
| OK | Trepat2009_NatPhys | https://doi.org/10.1038/nphys1269 |  | CrossRef: trepat 2009 — Physical forces during collective cell migration |
| OK | TrepatFredberg2011_TCB | https://doi.org/10.1016/j.tcb.2011.06.00 |  | CrossRef: trepat 2011 — Plithotaxis and emergent dynamics in collective cellular migration |
| OK | Tsujita2021_NatCommun | https://doi.org/10.1038/s41467-021-26156 |  | CrossRef: tsujita 2021 — Homeostatic membrane tension constrains cancer cell dissemination by c |
| OK | VincentEngler2013_BiotechnolJ | https://doi.org/10.1002/biot.201200205 |  | CrossRef: vincent 2013 — Mesenchymal stem cell durotaxis depends on substrate stiffness gradien |
| OK | Vishwakarma2018_NatCommun | https://doi.org/10.1038/s41467-018-05927 |  | CrossRef: vishwakarma 2018 — Mechanical interactions among followers determine the emergence of lea |
| OK | Wahlsten2023_ActaBiomater | https://doi.org/10.1016/j.actbio.2023.08 |  | CrossRef: wahlsten 2023 — Multiscale mechanical analysis of the elastic modulus of skin |
| OK | WalcottSun2010_PNAS | https://doi.org/10.1073/pnas.0912739107 |  | CrossRef: walcott 2010 — A mechanical model of actin stress fiber formation and substrate elast |
| OK | Wolf2013_JCB | https://doi.org/10.1083/jcb.201210152 |  | CrossRef: wolf 2013 — Physical limits of cell migration: Control by ECM space and nuclear de |
| OK | Woodcock2023_Cells | https://doi.org/10.3390/cells12192401 |  | CrossRef: woodcock 2023 — Measuring Melanoma Nanomechanical Properties in Relation to Metastatic |
| OK | XuSulchek2012_PLoSOne | https://doi.org/10.1371/journal.pone.004 |  | CrossRef: xu 2012 — Cell Stiffness Is a Biomarker of the Metastatic Potential of Ovarian C |
| OK | Yamada2019_NRMCB | https://doi.org/10.1038/s41580-019-0172- |  | CrossRef: yamada 2019 — Mechanisms of 3D cell migration |
| OK | YangAnseth2014_NatMater | https://doi.org/10.1038/nmat3889 |  | CrossRef: yang 2014 — Mechanical memory and dosing influence stem cell fate |
| OK | YangKaufman2009_BiophysJ | https://doi.org/10.1016/j.bpj.2008.10.06 |  | CrossRef: yang 2009 — Rheology and Confocal Reflectance Microscopy as Probes of Mechanical P |
| OK | Yao2014_SciRep | https://doi.org/10.1038/srep04960 |  | CrossRef: yao 2014 — Photoacoustic computed microscopy |
| OK | Yap2015_DevCell | https://doi.org/10.1016/j.devcel.2015.09 |  | CrossRef: yap 2015 — Adherens Junctions Revisualized: Organizing Cadherins as Nanoassemblie |
| OK | Yen2020_BBRC | https://doi.org/10.1016/j.bbrc.2020.03.1 |  | CrossRef: yen 2020 — Alteration of Young’s modulus in mesenchymal stromal cells during oste |
| OK | Yonemura2010_NCB | https://doi.org/10.1038/ncb2055 |  | CrossRef: yonemura 2010 — α-Catenin as a tension transducer that induces adherens junction devel |
| OK | Zbiral2023_IJMS | https://doi.org/10.3390/ijms241512208 |  | CrossRef: zbiral 2023 — Characterization of Breast Cancer Aggressiveness by Cell Mechanics |