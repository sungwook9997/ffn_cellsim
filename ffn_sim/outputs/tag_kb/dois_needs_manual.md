# SourceEvidence rows still missing a DOI (manual fill needed)

After auto-backfill, **62** NO_DOI_FOUND rows remain — CrossRef's author+year top-hit
was a wrong/unrelated paper (book, obscure, or the citation key itself is off), so no DOI was
written (avoiding injection of a wrong DOI). A few are real but defeated the journal/topic
matcher (e.g. Helfrich1973, MarantanMahadevan2018_AmJPhys, Ray2021_COCellBiol) — confirm by hand.

| citation_key | short_source | (rejected) auto-candidate |
|---|---|---|
| AlonsoMatilla2023_BiophysJ | Alonso-Matilla et al 2023, Biophys J | 10.1101/2023.11.13.566882 |
| Banerjee2015_PRL | Banerjee et al 2015, Phys Rev Lett 114:228101 | 10.1103/physrevlett.114.228101 |
| Bangasser2017_NCB | Bangasser et al 2017, Nat Cell Biol | 10.3791/56219 |
| Beaune2014_PNAS | Beaune et al 2014, PNAS 111:8055 | 10.4000/books.pur.49676 |
| Bobrowska2021_Cells | Bobrowska J, et al. Biomechanical characterization | 10.31338/uw.9788323548799.pp.90-102 |
| Cai2014_Cell | Cai et al 2014, Cell (collective migration leader) | 10.1201/b17354 |
| Cai2022_FCDB | Cai et al 2022, Front Cell Dev Biol | 10.3389/fcell.2021.817104 |
| Carlsson2018_JPhysCondMat | Carlsson 2018, J Phys Condens Matter | 10.2210/pdb5v7f/pdb |
| CavalcantiAdam2007_BiophysJ | Cavalcanti-Adam et al 2007, Biophys J | 10.1039/b614008d |
| Conklin2011 | Conklin et al 2011 | 10.1515/9781575066288 |
| Doss2020_PNAS | Doss et al 2020, PNAS | 10.1201/9780429344206 |
| EloseguiArtola2016_NatMater | Elosegui-Artola et al 2016, Nat Mater | 10.1016/j.bpj.2017.05.020 |
| FarooquiFenteany2005_JCS | Farooqui & Fenteany 2005, J Cell Sci 118:51 | 10.1242/jcs.01577 |
| Fredberg2015_FASEB | Fredberg 2015, FASEB J (tissue fluidization) | 10.1096/fasebj.29.1_supplement.85.5 |
| Fritzsche2013_MBC | Fritzsche et al 2013, Mol Biol Cell | 10.1007/978-3-642-37495-1_2 |
| Fu2024_PNAS | Fu et al 2024, PNAS (junctional viscosity) | 10.20944/preprints202406.1065.v2 |
| Garcia2015_PNAS | Garcia et al 2015, PNAS 112:15314 | 10.1073/pnas.1510973112 |
| GonzalezRodriguez2012_Science | Gonzalez-Rodriguez et al 2012, Science 338:910 | 10.1126/science.1226418 |
| Guo2017_PNAS | Guo et al 2017, PNAS | 10.1145/3145690.3145698 |
| Han2017_PNAS | Han et al 2017, PNAS | 10.2172/1477879 |
| Han2021_eLife | Han et al 2021, eLife | 10.4211/hs.7c5e032bdc7648a4a8c863a2c175e5b6 |
| Hannezo2014_PNAS | Hannezo et al 2014, PNAS 111:27 | 10.1073/pnas.1312076111 |
| HeadLevineMacKintosh2003_PRE | Head, Levine, MacKintosh 2003, Phys Rev E | 10.1103/physrevlett.91.108102 |
| Helfrich1973 | Helfrich 1973, Z Naturforsch C 28:693 | 10.1515/znc-1973-11-1209 |
| Jansen2018_BiophysJ | Jansen et al 2018, Biophys J | 10.1016/j.jas.2018.02.016 |
| Kim2021_PNAS | Kim et al 2021, PNAS | 10.1145/3508259.3508271 |
| Kothari2018_JApplMech | Kothari et al 2018, J Applied Mechanics | 10.2514/6.2018-0530 |
| Krause2014_NRMCB | Krause & Gautreau 2014, Nat Rev Mol Cell Biol 15:5 | 10.1038/nrm3861 |
| Lindstrom2013_SoftMatter | Lindstrom et al 2013, Soft Matter | 10.1093/obo/9780199766567-0108 |
| Liu2019_ProstateMech | Liu et al. 2019 prostate cancer cell-line mechanic | 10.1007/978-981-13-6962-9 |
| Liu2025_NatPhysics | Liu et al 2025, Nature Physics | 10.22541/au.176275808.84444425/v1 |
| Maitre2012_Nature | Maitre et al 2012, Nature / Science 338:253 | 10.3917/lav.daler.2012.01.0096 |
| Maitre2015_NCB | Maitre et al 2015, Nat Cell Biol 17:849 | 10.1051/shsconf/20152000001 |
| Malinova2021_NatCommun | Malinova et al 2021, Nat Commun | 10.1017/nps.2020.87 |
| MarantanMahadevan2018_AmJPhys | Marantan & Mahadevan 2018, Am J Phys | 10.1119/1.5003376 |
| Mogilner2002_2003 | Mogilner et al 2002/2003, Biophys J | 10.1002/3527601503.ch14 |
| Munster2013_PNAS | Munster et al 2013, PNAS | 10.7551/mitpress/8982.001.0001 |
| Notbohm2016_BiophysJ | Notbohm et al 2016, Biophys J 110:2729 | 10.1016/j.bpj.2016.05.019 |
| Odijk1995 | Odijk 1995 | 10.1037/e493132004-001 |
| Ozawa2020_JCB | Ozawa et al 2020, J Cell Biol | 10.1007/978-981-15-4190-2_7 |
| Pereverzev2005_BiophysJ | Pereverzev et al 2005, Biophys J | 10.1007/s11018-005-0098-9 |
| PerezGonzalez2019_NatPhys | Perez-Gonzalez et al 2019, Nat Physics 15:79 | 10.2172/1597056 |
| PollardBorisy2003_Cell | Pollard & Borisy 2003, Cell 112:453 | 10.1016/s0092-8674(03)00120-x |
| Proestaki2019_ExpMech | Proestaki et al 2019, Experimental Mechanics | 10.1007/s11340-018-00453-4 |
| Pronk2008_PNAS | Pronk et al 2008, PNAS / PRL | 10.1057/dev.2008.28 |
| Ray2021_COCellBiol | Ray & Lele 2021, Current Opinion Cell Biology | 10.1016/j.ceb.2021.05.004 |
| SchwarzSafran2013_RMP | Schwarz & Safran 2013, Rev Mod Phys 85:1327 | 10.1103/revmodphys.85.1327 |
| Steinberg1996_CurrTopDevBiol | Steinberg 1996, Curr Top Dev Biol (DAH) | 10.1387/ijdb.8877443 |
| TapiaRojo2019_SciAdv | Tapia-Rojo et al 2019, Science Advances | 10.1063/1.5126071 |
| Trichet2012_PNAS | Trichet et al 2012, PNAS | 10.3917/bupsy.520.0365 |
| Varma2016_BiophysJ | Varma et al 2016, Biophys J | 10.4135/9781526429629 |
| Vassalli2023_Cancers | Vassalli M, Podesta A, et al. AFM mechanics of ova | 10.3917/mem.086.0017 |
| Venturini2020_Science | Venturini et al 2020, Science | 10.32614/cran.package.dmbc |
| Vishwakarma2018_NatCommun | Vishwakarma et al 2018, Nat Commun 9:3469 | 10.1038/s41467-018-05927-6 |
| Wahlsten2023_ActaBiomater | Wahlsten et al 2023, Acta Biomaterialia | 10.1016/j.actbio.2023.08.030 |
| Wang1997_BiophysJ | Wang et al 1997, Biophys J | 10.1115/97-aa-050 |
| Wei2024_ActaBiomater | Wei et al 2024, Acta Biomaterialia | 10.26434/chemrxiv-2023-8c9vh-v2 |
| Wenger2007 | Wenger et al 2007 | 10.1097/mlr.0b013e31815b97bf |
| WirtzFriedl2011 | Wirtz, Konstantopoulos, Searson 2011 / Friedl migr | 10.1038/nrc3080 |
| Yao2014_NatCommun | Yao et al 2014, Nat Commun 5:4525 | 10.1002/chin.201448205 |
| Yao2014_SciRep | Yao et al 2014, Scientific Reports | 10.1038/srep04960 |
| Yen2020_BBRC | Yen et al. Lineage-dependent cell modulus. Biochem | 10.1016/j.bbrc.2020.03.146 |