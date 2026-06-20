# Results-integrity audit (disk-grounded)

Audited **17** headline result-claims from `results_manifest.yaml`.

## Summary

| verdict | n | meaning |
|---|---|---|
| RETRACT | 1 | disk contradicts claim (forbidden metric / value drift) |
| NEEDS_REGEN | 3 | claimed artifact ABSENT — regenerate + commit |
| GPU_UNREPRODUCED | 2 | present but GPU-only, no committed build/CI trace |
| VERIFIED | 11 | artifact present, metric sanctioned, value matches |

**0 DRIFT rows** (disk worse than declared — these FAIL the CI gate).

## All claims (suspicion-ranked)

| verdict | drift | claim | metric | note |
|---|---|---|---|---|
| RETRACT |  | C14-dcm-footprint-report | basal-contact-hull | PI-FORBIDDEN metric 'basal-contact-hull' — A/A0 must be top-down silhouette |
| NEEDS_REGEN |  | C12-dcm-compaction | top-down-silhouette | claimed artifact ABSENT on disk: ffn_sim/outputs/h_dcm_two_stage/two_stage_n400_topdown_summary.json |
| NEEDS_REGEN |  | C13-cleanball-187 | top-down-silhouette | claimed artifact ABSENT on disk: ffn_sim/outputs/h_dcm_two_stage/cleanball_n12.pkl |
| NEEDS_REGEN |  | C17-aa0-fit | top-down-silhouette | claimed artifact ABSENT on disk: ffn_sim/outputs/h_dcm_two_stage/cleanball_aa0_fit.json |
| GPU_UNREPRODUCED |  | C2-compartment-speedup | n/a | value OK: speedups_vs_cupy_cpucomp.native_gpucomp=5.863981474746499 ~= 5.864; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | C3-fullcell-speedup | n/a | value OK: end_to_end_speedup_native_over_cupy=2.432515605855367 ~= 2.4325; GPU-only, no committed build/CI trace |
| VERIFIED |  | warp-B1-baoab-kt0 | n/a | value OK: B1_baoab.kt0.max_abs_pos_diff=0.0 ~= 0.0 |
| VERIFIED |  | warp-B1-baoab-ktpos | n/a | value OK: B1_baoab.ktpos.max_abs_pos_diff=0.0 ~= 0.0 |
| VERIFIED |  | warp-B2-radial-membrane | n/a | value OK: B2_radial_shell.membrane.host_force_rel=0.0 ~= 0.0 |
| VERIFIED |  | warp-B2-radial-nucleus | n/a | value OK: B2_radial_shell.nucleus.host_force_rel=0.0 ~= 0.0 |
| VERIFIED |  | warp-B2-radial-turgor | n/a | value OK: B2_radial_shell.turgor.host_force_rel=4.6959403038584427e-14 ~= 0.0 |
| VERIFIED |  | warp-B4-differentiability | n/a | value OK: B4_differentiability.rel_err_vs_analytic=2.596797331012674e-15 ~= 0.0 |
| VERIFIED |  | warp-Fixman-metric-force | n/a | value OK: Fixman_metric_force.force_rel=6.148529216107672e-16 ~= 0.0 |
| VERIFIED |  | warp-MSHAKE-chain-constraint | n/a | value OK: MSHAKE_chain_constraint.pos_rel=7.365694533460985e-12 ~= 0.0 |
| VERIFIED |  | warp-compartment-harmonic-angle | n/a | value OK: compartment_harmonic_angle.force_rel=2.2066914540404608e-13 ~= 0.0 |
| VERIFIED |  | warp-compartment-harmonic-bond | n/a | value OK: compartment_harmonic_bond.force_rel=4.782658393999819e-16 ~= 0.0 |
| VERIFIED |  | warp-compartment-lj-wca | n/a | value OK: compartment_lj_wca.force_rel=1.753578873717259e-15 ~= 0.0 |