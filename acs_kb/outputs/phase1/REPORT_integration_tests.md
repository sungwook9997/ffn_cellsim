# Phase 1 — Cross-Unit Integration Tests Report

**Date**: 2026-05-18
**Worker**: A (post-Unit 1.2)
**Branch**: `worker-b/bridge` (checked out at run time; commit is Worker A's)
**Scope**: cross-unit integration tests + Taichi-v4 contamination check
**Repo**: 2 commits added — `acs_kb` complete + `docs/v2` archived

---

## Contamination check (Task 3)

Pattern grep over `acs_kb/` for: `docs/v2 | taichi | mls.mpm | material.point | protrusion_coupled | v2_phase1 | v2_63`.

```
$ grep -rni "docs/v2|taichi|mls.mpm|material.point|protrusion_coupled|v2_phase1|v2_63" acs_kb/
acs_kb/__init__.py:4:(Phase 1: ECM only). Kept separate from `acs/` (the Taichi MLS-MPM
```

**Verdict**: **1 hit, intentional boundary marker, not contamination.**

The single hit is in `acs_kb/__init__.py` line 4, inside the top-level package docstring:

> "Kept separate from `acs/` (the Taichi MLS-MPM production code) so the KB-derived discrete-fiber model can develop independently."

This was written at Unit 1.1 setup to declare what `acs_kb/` is NOT (it is not the Taichi `acs/` track). It is a boundary marker that explains the rationale for the parallel directory, not a reference to v2 plans. **Recommended action: keep.** If the PI prefers zero matches, swap "Taichi MLS-MPM" → "the legacy production code" in the same line.

No hits in any other module, config, test, notebook, or output document.

## Integration tests (Tasks 4 + 5a)

```
tests/integration/ — 10 passed in 0.37 s

test_cell_cortex_uses_ecm.py ............ 4 / 4 PASS
test_motor_clutch_substrate_swap.py ..... 4 / 4 PASS
test_acs_kb_public_api_contract.py ...... 2 / 2 PASS
```

Full log: `acs_kb/outputs/phase1/integration_test_run.log`.

### (a) `test_cell_cortex_uses_ecm.py` — 4 tests

Verifies Worker C's cell/cortex stack reuses Worker A's ECM kernels rather than carrying parallel re-implementations.

| Test | Verifies | KU |
|---|---|---|
| `test_force_balance_imports_real_ecm_compute_forces` | `cell.force_balance` module holds a reference whose value `is` `ecm.fiber_mechanics.compute_forces` | KU-1.24 |
| `test_cortex_imports_real_ecm_cross_links_kernel` | `cell.cortex` module holds a reference whose value `is` `ecm.cross_links._segment_intersections` | KU-1.27, KU-3.1 |
| `test_compute_forces_runs_on_cortex_bead_state` | A Cortex's `bead_positions`, `rest_length`, `cross_links` are shape-compatible with the ECM force kernel; Newton 3rd law holds (Σ F = 0) | KU-1.24, KU-3.1 |
| `test_force_balance_signature_documents_ecm_dependency` | At least one public callable in `cell.force_balance` accepts a `cortex` or `forces` parameter (prevents silent module emptying) | KU-3.17 |

### (b) `test_motor_clutch_substrate_swap.py` — 4 tests

Pins the substrate-stub ↔ future-ECM-adapter swap contract so the
Unit 2.2+ upgrade is mechanical, not architectural.

| Test | Verifies | KU |
|---|---|---|
| `test_substrate_stub_conforms_to_swap_protocol` | `LinearElasticSubstrate` satisfies a `runtime_checkable` `Protocol` with `stiffness` + `compute_displacement` | KU-1.21 |
| `test_two_stub_instances_agree_within_tolerance` | Two stubs at identical KU-1.21 params produce identical `stiffness` and `compute_displacement(F)` for F ∈ {0, 1 pN, 0.1 nN, 5 nN, 10 nN} (within 10 % swap target) | KU-1.21 |
| `test_motor_clutch_only_uses_swap_contract` | Source-level grep on `MotorClutchFA` confirms no `self.substrate.<x>` access outside `{stiffness, compute_displacement}` — i.e. the future ECM adapter does not need to expose anything more | KU-2.4 |
| `test_motor_clutch_step_with_two_identical_substrates_agrees` | 200-step `step()` trajectory with identical RNG seeds + identical stubs is bit-identical (`force_total` per step) | KU-2.4, KU-2.12 |

### (c) `test_acs_kb_public_api_contract.py` — 2 tests + JSON baseline

Pins the inspect.signature of every public symbol that Worker B / C
import by name. New baseline written this run at
`tests/integration/api_contract_baseline.json` (13 symbols / 9.5 kB).

| Test | Verifies | KU |
|---|---|---|
| `test_public_api_all_symbols_importable` | Every entry in the CONTRACT list resolves to an actual object | KU-1.24, KU-1.27, KU-1.28 |
| `test_public_api_signature_baseline` | Parameter names, kinds, defaults, and annotations match the frozen JSON baseline; future drift fails loudly with a diff. Refresh with `ACS_KB_API_REBASELINE=1` and commit the new JSON | KU-1.24, KU-1.27, KU-1.28 |

CONTRACT covers:
- `acs_kb.ecm.fiber_network`: `FiberNetwork`, `generate_2d_fiber_network`, `measure_nematic_order`
- `acs_kb.ecm.fiber_mechanics`: `compute_energy`, `compute_forces`, `energy_of`, `forces_of`, `total_force`
- `acs_kb.ecm.cross_links`: `CrossLink`, `generate_cross_links`, `measure_coordination`, `measure_xl_per_fiber`, `compute_xl_energy_and_forces`

## API / interface issues discovered

None blocking. Minor notes for future passes:

1. Worker C's `cell.cortex` imports `_segment_intersections` (the underscore-prefixed private helper). Test (a) accepts this for Phase 1 — it is genuine kernel reuse — but at Phase 2 the helper should either become a public API (`segment_intersections_2d`) or move into `acs_kb.common`. Filed mentally; not in scope for this report's git history.
2. `acs_kb.ecm.fiber_network` and `cross_links` leak `Any`, `brentq`, `dataclass`, `field`, `ive` into their namespace via wildcard-style `from … import …` lines. The contract test uses a curated list, so leaks are invisible to consumers, but `__all__` declarations would tighten this.

## Skipped-test reactivation procedure

All three test files use `pytest.mark.skipif` guards that key on a
real `importlib.import_module()` attempt:

- `skip_if_worker_a_missing` — `acs_kb.ecm.fiber_mechanics` importable
- `skip_if_worker_c_missing` — `acs_kb.cell.{cortex, force_balance}` importable
- `skip_if_bridge_missing` — `acs_kb.bridge.{substrate_stub, motor_clutch, types}` importable

This run: **all three guards passed**, so every test executed. If a
future tear-down or refactor temporarily removes a worker module, the
guard falls back to `skip` rather than `error`. To re-activate after
the missing module returns, no test edit is needed — `pytest` will
pick the test back up on the next run.

## Files added / modified this session

```
A  acs_kb/                                  (whole subpackage — see Worker A Unit 1.1 + 1.2 reports)
R  docs/v2/ → docs/_archive_taichi_v4_DO_NOT_REFERENCE/
A  docs/_archive_taichi_v4_DO_NOT_REFERENCE/_README.txt
A  tests/integration/__init__.py
A  tests/integration/test_cell_cortex_uses_ecm.py
A  tests/integration/test_motor_clutch_substrate_swap.py
A  tests/integration/test_acs_kb_public_api_contract.py
A  tests/integration/api_contract_baseline.json
A  acs_kb/outputs/phase1/contamination_check.txt
A  acs_kb/outputs/phase1/integration_test_run.log
A  acs_kb/outputs/phase1/REPORT_integration_tests.md   (this file)
```

## Stop

Per PI instruction (Task 5): Worker A halts here. **Unit 1.3, 2.x, 3.x,
4.x, Phase 2** — no autonomous progression. Awaiting next direction.
