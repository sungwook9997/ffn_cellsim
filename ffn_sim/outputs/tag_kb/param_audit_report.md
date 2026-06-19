# Parameter-provenance audit (disk + citation-grounded)

Audited **7** headline simulation constants from `params_manifest.yaml` (constant -> KU -> citation -> verdict).

## Summary

| verdict | n | meaning |
|---|---|---|
| UNSOURCED | 2 | no KU / citation_key — constant not traceable to a paper |
| SOURCE_UNVERIFIED | 2 | citation verdict CHECK/NO_DOI_FOUND (real paper, not confirmed-OK) |
| VERIFIED | 3 | value matches disk + KU present + citation verdict OK |

**0 DRIFT rows** (disk worse than declared — these FAIL the CI gate).

## All constants (suspicion-ranked)

| verdict | drift | id | location | KU | citation | note |
|---|---|---|---|---|---|---|
| UNSOURCED |  | P6-myosin-stall | phase1_unit2_1.yaml:bridge.motor_clutch.F_stall_per_motor | KU-2.18 |  | value OK (2e-12); KU KU-2.18 present but KU->source link not committed |
| UNSOURCED |  | P7-ecm-persistence-length | phase1_unit1.yaml:ecm.persistence_length | KU-1.1 |  | value OK (1.7e-05); KU KU-1.1 present but KU->source link not committed |
| SOURCE_UNVERIFIED |  | P4-cortex-gamma | phase1_unit3.yaml:cell.gamma_cortex | KU-3.5 | Chugh2017_NatCellBiol | value OK (0.0005); Chugh2017_NatCellBiol citation verdict=CHECK |
| SOURCE_UNVERIFIED |  | P5-catch-bond-Fc | phase1_unit2_1.yaml:bridge.catch_bond.F_c | KU-2.5 | Pereverzev2005_BiophysJ | value OK (7e-12); Pereverzev2005_BiophysJ citation verdict=NO_DOI_FOUND |
| VERIFIED |  | P1-junction-dx-star-phase1 | phase1_unit4_1.yaml:junction.dx_star_phase1 | KU-4.17 | Bell1978_Science | value OK (1e-10); Bell1978_Science citation verdict=OK |
| VERIFIED |  | P2-junction-dx-star-ku417 | phase1_unit4_1.yaml:junction.dx_star_ku417 | KU-4.17 | Buckley2014_Science | value OK (4e-09); Buckley2014_Science citation verdict=OK |
| VERIFIED |  | P3-talin-unfold-rate | phase1_unit2_2.yaml:bridge.talin.k_u0 | KU-2.6 | DelRio2009_Science | value OK (0.01); DelRio2009_Science citation verdict=OK |