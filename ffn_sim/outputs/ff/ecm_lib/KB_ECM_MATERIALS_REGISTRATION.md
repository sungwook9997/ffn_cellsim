# KB registration draft — ECM materials (FF ECM library, 2026-07-10)

New KnowledgeClaims + SourceEvidence for the ECM materials the FF ECM library covers that were **absent
from the KB** (collagen-I and PAA-substrate already present). Values are DOI-verified via the research
workflow (adversarial verify stage). To register: add SourceEvidence rows + KnowledgeClaims to the Notion
Contract-Graph (SoT), then `bash outputs/tag_kb/refresh.sh` + `verify_sources.py --check`. Proposed KB IDs
in the Unit-1.V Microenv namespace (`KB-1.V.3.x`); **PI to ratify IDs before commit**.

## KB-1.V.3.1 — Polyacrylamide (PAA) gel: acrylamide/bis → E calibration
- **Value:** linear-elastic NON-fibrillar chemical gel; E tunable **0.1–40 kPa** by total acrylamide %(w/v)
  + bis-acrylamide %(w/v). ν≈0.45–0.48 undrained (0.30–0.37 fully drained). E(bis) NON-monotonic (softens
  past a crosslink inflection, clustering). Anchors (Subramani 2020, AFM+rheo): 8%/0.01%bis→620 Pa,
  5%/0.1%→1970 Pa, 8%/0.2%→7930 Pa, 8%/0.4%→13 kPa. Does NOT strain-stiffen (contrast collagen/fibrin).
- **Subtopic:** substrate-stiffness / continuum gel. Extends KB-1.21/KB-1.5.
- **Sources (verified):** Tse&Engler 2010 CurrProtocCellBiol `10.1002/0471143030.cb1016s47`;
  Pelham&Wang 1997 PNAS `10.1073/pnas.94.25.13661`; Yeung 2005 CellMotilCytoskeleton `10.1002/cm.20041`;
  Engler 2006 Cell `10.1016/j.cell.2006.06.044`; Denisin&Pruitt 2016 ACS-AMI `10.1021/acsami.5b09344`.

## KB-1.V.3.2 — Fibrin gel mechanics
- **Value:** fibrous clot network; **G0 0.1–2000 Pa** over 0.1–8 mg/mL fibrinogen, power law **G~c^2.3**
  (Piechocka 2010). Fiber diameter ~100–130 nm; single-fiber Young's E **1.7±1.3 MPa** (Collet 2005);
  fibers strain to ~180% without permanent lengthening (Liu 2006) → strong strain-stiffening. Mesh ~1–5 µm.
- **Subtopic:** fibrillar ECM / provisional wound matrix.
- **Sources (verified):** Piechocka 2010 BiophysJ `10.1016/j.bpj.2010.01.040`;
  Collet 2005 PNAS `10.1073/pnas.0504120102`; Storm 2005 Nature `10.1038/nature03521`;
  Liu 2006 Science `10.1126/science.1127317`; Ryan 1999 BiophysJ `10.1016/S0006-3495(99)77113-4`.

## KB-1.V.3.3 — Matrigel / reconstituted basement membrane
- **Value:** soft laminin-111 + collagen-IV BM mimic, near-continuum fine mesh; **E ~30–900 Pa (typ 450,
  AFM 37°C**, Soofi 2009). Batch-variable; poroelastic (permeable polymer network in liquid). ν≈0.45–0.5.
- **Subtopic:** basement membrane / continuum gel.
- **Sources (verified):** Soofi 2009 JStructBiol `10.1016/j.jsb.2009.05.005`;
  Reed 2009 Langmuir `10.1021/la8033098`; Hughes 2010 Proteomics `10.1002/pmic.200900758`;
  Li 2021 PNAS `10.1073/pnas.2022422118`; Fabris 2018 BiophysJ `10.1016/j.bpj.2018.09.020`.

## KB-1.V.3.4 — Agarose gel mechanics
- **Value:** polysaccharide thermogel; **E ~5 kPa–1 MPa (typ ~30 kPa @1%w/v)**, G~c^2.3. Structurally a
  fine sub-isostatic athermal bundle network (mesh ~0.24 µm, pore a~C^-γ γ=0.5–0.75) but at cell scale a
  near-incompressible elastic continuum (ν≈0.5). Cartilage/chondrocyte-encapsulation mimic.
- **Subtopic:** encapsulation gel / continuum.
- **Sources (verified):** Normand 2000 Biomacromolecules `10.1021/bm005583j`;
  Pernodet 1997 Electrophoresis `10.1002/elps.1150180111`;
  Martikainen 2020 Macromolecules `10.1021/acs.macromol.0c00601`; Zignego 2014 JBiomech `10.1016/j.jbiomech.2013.10.051`.

## KB-1.V.3.5 — Hyaluronic acid (HA) gel mechanics
- **Value:** soft flexible polyelectrolyte, usually crosslinked (MeHA); **E ~50 Pa–35 kPa** tunable by
  weight fraction (1–8 wt%) AND crosslink density (independently); brain/neural-mimetic ~0.1–1 kPa. ν≈0.5.
  Typically non-fibrillar (mesh ~0.1 µm).
- **Subtopic:** neural-soft ECM / continuum gel.
- **Sources (verified):** Ananthanarayanan/Kim/Kumar 2011 Biomaterials `10.1016/j.biomaterials.2011.07.005`;
  Burdick&Prestwich 2011 AdvMater `10.1002/adma.201003963`; Ondeck&Engler 2016 JBiomechEng `10.1115/1.4032429`;
  Guvendiren 2013 JMBBM `10.1016/j.jmbbm.2013.11.008`.

## RunResult (FF ECM library validation, 2026-07-10)
- Harness validated on known-modulus continuum (PA E=5000 → uniaxial 5000, shear 0.91× isotropic relation,
  indentation 0.68–0.79×). 5/6 materials' reference modulus IN literature band; agarose modeled as continuum
  (athermal Mikado underestimates dense stiff gels ~60×). Collagen G'(c): absolute Pa at reference matches
  (G(1.5mg/mL)≈11 Pa) but scaling n≈1.07 vs literature 2.05 (athermal density-limit, documented). Alignment
  S∈{0,0.3,0.6,0.85} → anisotropy E∥/E⊥ = 1.0 / 3.1 / 10.3 / 63 (dermis→TACS-2→TACS-3→tendon).
- Code: `ffn_sim/ff/ecm_library.py`, `ffn_sim/ff/ecm_mechanics.py`,
  `ffn_sim/scripts/ff_ecm_validate.py`, `ffn_sim/scripts/ff_ecm_viewer.py`.
