# Parameter-provenance audit (disk + citation-grounded)

Audited **42** headline simulation constants from `params_manifest.yaml` (constant -> KU -> citation -> verdict).

## Summary

| verdict | n | meaning |
|---|---|---|
| UNSOURCED | 10 | no KU / citation_key — constant not traceable to a paper |
| SOURCE_UNVERIFIED | 13 | citation verdict CHECK/NO_DOI_FOUND (real paper, not confirmed-OK) |
| VERIFIED | 19 | value matches disk + KU present + citation verdict OK |

**0 DRIFT rows** (disk worse than declared — these FAIL the CI gate).

## All constants (suspicion-ranked)

| verdict | drift | id | location | KU | citation | note |
|---|---|---|---|---|---|---|
| UNSOURCED |  | bridge-substrate-young-modulus | phase1_unit2_1.yaml:bridge.substrate.young_modulus | KU-1.5 |  | value OK (5000); KU KU-1.5 present but KU->source link not committed |
| UNSOURCED |  | cell-L-cortex-fiber | phase1_unit3.yaml:cell.L_cortex_fiber | KU-3.17 |  | value OK (3e-06); KU KU-3.17 present but KU->source link not committed |
| UNSOURCED |  | cell-R-cell | phase1_unit3.yaml:cell.R_cell | KU-3.17 |  | value OK (1e-05); KU KU-3.17 present but KU->source link not committed |
| UNSOURCED |  | cell-target-coordination-ref | phase1_unit3.yaml:cell.target_coordination_ref | KU-3.17 |  | value OK (3.5); KU KU-3.17 present but KU->source link not committed |
| UNSOURCED |  | cell-water-viscosity | phase1_unit3.yaml:cell.water_viscosity | KU-1.26 |  | value OK (0.0006913); KU KU-1.26 present but KU->source link not committed |
| UNSOURCED |  | cell-xl-stiffness | phase1_unit3.yaml:cell.xl_stiffness | KU-1.28 |  | value OK (0.001); KU KU-1.28 present but KU->source link not committed |
| UNSOURCED |  | ecm-L-fiber | phase1_unit1.yaml:ecm.L_fiber | KU-1.22 |  | value OK (1e-05); KU KU-1.22 present but KU->source link not committed |
| UNSOURCED |  | ecm-water-viscosity | phase1_unit1.yaml:ecm.water_viscosity | KU-1.26 |  | value OK (0.0006913); KU KU-1.26 present but KU->source link not committed |
| UNSOURCED |  | ecm-xl-stiffness | phase1_unit1.yaml:ecm.xl_stiffness | KU-1.28 |  | value OK (0.001); KU KU-1.28 present but KU->source link not committed |
| UNSOURCED |  | junction-kT | phase1_unit4_1.yaml:junction.kT | KU-1.26 |  | value OK (4.28e-21); KU KU-1.26 present but KU->source link not committed |
| SOURCE_UNVERIFIED |  | bridge-catch-bond-F-c | phase1_unit2_1.yaml:bridge.catch_bond.F_c | KU-2.18 | Pereverzev2005_BiophysJ | value OK (7e-12); Pereverzev2005_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | bridge-catch-bond-F-s | phase1_unit2_1.yaml:bridge.catch_bond.F_s | KU-2.18 | Pereverzev2005_BiophysJ | value OK (3e-11); Pereverzev2005_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | bridge-catch-bond-k-off-catch | phase1_unit2_1.yaml:bridge.catch_bond.k_off_catch | KU-2.18 | Pereverzev2005_BiophysJ | value OK (0.4); Pereverzev2005_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | bridge-catch-bond-k-off-slip | phase1_unit2_1.yaml:bridge.catch_bond.k_off_slip | KU-2.18 | Pereverzev2005_BiophysJ | value OK (0.5); Pereverzev2005_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | cell-bending-modulus | phase1_unit3.yaml:cell.bending_modulus | KU-1.1 | Lindstrom2010_Biomaterials | value OK (7e-26); Lindstrom2010_Biomaterials citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | cell-gamma-cortex | phase1_unit3.yaml:cell.gamma_cortex | KU-3.5 | Chugh2017_NatCellBiol | value OK (0.0005); Chugh2017_NatCellBiol citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | cell-gamma-cortex-blebbistatin | phase1_unit3.yaml:cell.gamma_cortex_blebbistatin | KU-3.5 | Chugh2017_NatCellBiol | value OK (5e-05); Chugh2017_NatCellBiol citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | cell-persistence-length | phase1_unit3.yaml:cell.persistence_length | KU-1.1 | Lindstrom2010_Biomaterials | value OK (1.7e-05); Lindstrom2010_Biomaterials citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | cell-stretching-modulus | phase1_unit3.yaml:cell.stretching_modulus | KU-1.2 | Jansen2018_BiophysJ | value OK (8.6e-09); Jansen2018_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | ecm-bead-radius | phase1_unit1.yaml:ecm.bead_radius | KU-1.2 | Jansen2018_BiophysJ | value OK (5e-08); Jansen2018_BiophysJ citation verdict=NO_DOI_FOUND |
| SOURCE_UNVERIFIED |  | ecm-bending-modulus | phase1_unit1.yaml:ecm.bending_modulus | KU-1.1 | Lindstrom2010_Biomaterials | value OK (7e-26); Lindstrom2010_Biomaterials citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | ecm-persistence-length | phase1_unit1.yaml:ecm.persistence_length | KU-1.1 | Lindstrom2010_Biomaterials | value OK (1.7e-05); Lindstrom2010_Biomaterials citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | ecm-stretching-modulus | phase1_unit1.yaml:ecm.stretching_modulus | KU-1.2 | Jansen2018_BiophysJ | value OK (8.6e-09); Jansen2018_BiophysJ citation verdict=NO_DOI_FOUND |
| VERIFIED |  | bridge-fa-growth-A-nascent | phase1_unit2_2.yaml:bridge.fa_growth.A_nascent | KU-2.18 | Bangasser2013_BiophysJ | value OK (1e-12); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-fa-growth-n-clutches-nascent | phase1_unit2_2.yaml:bridge.fa_growth.n_clutches_nascent | KU-2.18 | Bangasser2013_BiophysJ | value OK (50); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-F-stall-per-motor | phase1_unit2_1.yaml:bridge.motor_clutch.F_stall_per_motor | KU-2.18 | Bangasser2013_BiophysJ | value OK (2e-12); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-fa-area | phase1_unit2_1.yaml:bridge.motor_clutch.fa_area | KU-2.18 | Bangasser2013_BiophysJ | value OK (1e-12); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-k-int | phase1_unit2_1.yaml:bridge.motor_clutch.k_int | KU-2.18 | Bangasser2013_BiophysJ | value OK (0.001); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-k-on | phase1_unit2_1.yaml:bridge.motor_clutch.k_on | KU-2.18 | Bangasser2013_BiophysJ | value OK (1); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-n-clutches | phase1_unit2_1.yaml:bridge.motor_clutch.n_clutches | KU-2.18 | Bangasser2013_BiophysJ | value OK (50); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-n-motors | phase1_unit2_1.yaml:bridge.motor_clutch.n_motors | KU-2.18 | Bangasser2013_BiophysJ | value OK (50); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-motor-clutch-v-unloaded | phase1_unit2_1.yaml:bridge.motor_clutch.v_unloaded | KU-2.18 | Bangasser2013_BiophysJ | value OK (1e-07); Bangasser2013_BiophysJ citation verdict=OK |
| VERIFIED |  | bridge-substrate-contact-radius | phase1_unit2_1.yaml:bridge.substrate.contact_radius | KU-1.21 | Discher2005_Science | value OK (1e-06); Discher2005_Science citation verdict=OK |
| VERIFIED |  | bridge-substrate-poisson-ratio | phase1_unit2_1.yaml:bridge.substrate.poisson_ratio | KU-1.21 | Discher2005_Science | value OK (0.45); Discher2005_Science citation verdict=OK |
| VERIFIED |  | bridge-talin-dx-star | phase1_unit2_2.yaml:bridge.talin.dx_star | KU-2.6 | DelRio2009_Science | value OK (1.5e-09); DelRio2009_Science citation verdict=OK |
| VERIFIED |  | bridge-talin-k-u0 | phase1_unit2_2.yaml:bridge.talin.k_u0 | KU-2.6 | DelRio2009_Science | value OK (0.01); DelRio2009_Science citation verdict=OK |
| VERIFIED |  | ecm-target-coordination-ref | phase1_unit1.yaml:ecm.target_coordination_ref | KU-1.3 | BroederszMacKintosh2014_RMP | value OK (3.2); BroederszMacKintosh2014_RMP citation verdict=OK |
| VERIFIED |  | junction-F-per-bond-average | phase1_unit4_1.yaml:junction.F_per_bond_average | KU-4.17 | Buckley2014_Science | value OK (3e-11); Buckley2014_Science citation verdict=OK |
| VERIFIED |  | junction-dx-star-ku417 | phase1_unit4_1.yaml:junction.dx_star_ku417 | KU-4.17 | Buckley2014_Science | value OK (4e-09); Buckley2014_Science citation verdict=OK |
| VERIFIED |  | junction-k-off0 | phase1_unit4_1.yaml:junction.k_off0 | KU-4.17 | Buckley2014_Science | value OK (0.5); Buckley2014_Science citation verdict=OK |
| VERIFIED |  | junction-k-on | phase1_unit4_1.yaml:junction.k_on | KU-4.17 | Buckley2014_Science | value OK (1); Buckley2014_Science citation verdict=OK |
| VERIFIED |  | junction-n-bonds-total | phase1_unit4_1.yaml:junction.n_bonds_total | KU-4.17 | Buckley2014_Science | value OK (100); Buckley2014_Science citation verdict=OK |