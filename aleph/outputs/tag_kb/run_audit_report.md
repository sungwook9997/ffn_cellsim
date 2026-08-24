# Results-integrity audit (disk-grounded)

Audited **35** headline result-claims from `results_manifest.yaml`.

## Summary

| verdict | n | meaning |
|---|---|---|
| RETRACT | 1 | disk contradicts claim (forbidden metric / value drift) |
| NEEDS_REGEN | 3 | claimed artifact ABSENT — regenerate + commit |
| GPU_UNREPRODUCED | 26 | present but GPU-only, no committed build/CI trace |
| VERIFIED | 5 | artifact present, metric sanctioned, value matches |

**0 DRIFT rows** (disk worse than declared — these FAIL the CI gate).

## All claims (suspicion-ranked)

| verdict | drift | claim | metric | note |
|---|---|---|---|---|
| RETRACT |  | C14-dcm-footprint-report | basal-contact-hull | PI-FORBIDDEN metric 'basal-contact-hull' — A/A0 must be top-down silhouette |
| NEEDS_REGEN |  | C12-dcm-compaction | top-down-silhouette | claimed artifact ABSENT on disk: aleph/outputs/_archive/h_dcm_two_stage/two_stage_n400_topdown_summary.json |
| NEEDS_REGEN |  | C13-cleanball-187 | top-down-silhouette | claimed artifact ABSENT on disk: aleph/outputs/_archive/h_dcm_two_stage/cleanball_n12.pkl |
| NEEDS_REGEN |  | C17-aa0-fit | top-down-silhouette | claimed artifact ABSENT on disk: aleph/outputs/_archive/h_dcm_two_stage/cleanball_aa0_fit.json |
| GPU_UNREPRODUCED |  | AC-bleb-cross-build-parity | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-dt-cause-is-outer-integrator | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-dt-independence-refuted | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-dt-order-undefined | n/a | value OK: config.dt_phys=0.0025 ~= 0.0025; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-kxb-1 | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-kxb-10 | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-kxb-100 | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-kxb-1000 | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-kxb-3 | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-bleb-settled-report | n/a | present; 'n/a' not machine-checkable (.md); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-cortex-ownership-T2-double-count-control | n/a | present; 'measurements.omit_mask_double_count_control.fails_by_exactly_two' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-engine-cortex-AB-native | n/a | present; 'config.engine_cortex' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-engine-cortex-AB-report | n/a | present; 'n/a' not machine-checkable (.md); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-engine-cortex-AB-settled | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-native-cortex-spectrum-stability | n/a | value OK: measurements.native_spectrum.explicit_step_verdict.stability_product_base=0.294585935595612 ~= 0.294586; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | AC-stationarity-kxb1000-timescale-separation | n/a | present; 'measurements.stationarity.gamma_total_pn_per_um.verdict' not machine-checkable (.json); GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | C2-compartment-speedup | n/a | value OK: speedups_vs_cupy_cpucomp.native_gpucomp=5.863981474746499 ~= 5.864; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | C3-fullcell-speedup | n/a | value OK: end_to_end_speedup_native_over_cupy=2.432515605855367 ~= 2.4325; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | inner-budget-down-sweep | n/a | value OK: result.lowest_insensitive_outer=40 ~= 40; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | inner-budget-down-sweep-report | n/a | artifact present; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | observe-crossbridge-is-non-conservative | n/a | value OK: measurements.loop_work_bound.convergence.amplitude_exponent=1.9548221780529091 ~= 1.9548; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | observe-minifilament-internal-zero-modes | n/a | value OK: measurements.minifilament_internal_floppiness.internal_zero_modes_per_minifilament=19 ~= 19; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | observe-tangent-multilinearity | n/a | value OK: measurements.operator_bound.multilinearity.relative_residual=9.373463765719657e-17 ~= 0.0; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | sf-motor-balance-gate-is-the-acceptance-predicate | n/a | value OK: measurements.acceptance.balance_gate.gamma_n=1.0025313912366168e-13 ~= 1.0025e-13; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | sf-motor-balance-gate-positive-control | n/a | value OK: measurements.acceptance.positive_control.balance_residual_pN=0.49999999999999917 ~= 0.5; GPU-only, no committed build/CI trace |
| GPU_UNREPRODUCED |  | sf-motor-energy-balance-residual-is-the-integrator | n/a | value OK: measurements.energy_ledger.mobility_scaling.convergence.residual_exponent=1.002509024280737 ~= 1.0025; GPU-only, no committed build/CI trace |
| VERIFIED |  | L0-run1-void-criterion | n/a | artifact present |
| VERIFIED |  | observe-fisher-effective-dimension | n/a | artifact present |
| VERIFIED |  | observe-floppy-mode-count | n/a | artifact present |
| VERIFIED |  | observe-gershgorin-vs-true-lambda-max | n/a | artifact present |
| VERIFIED |  | observe-timescale-separation | n/a | artifact present |