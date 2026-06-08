# Compartment ENABLED-PATH smoke harnesses (2026-06-09)

Standalone plumbing harnesses that assemble each default-OFF compartment into a
minimal HOOMD simulation, step it with the project L-M BAOAB integrator, and
confirm the **enabled** build/attach/step path runs without crashing. These are
the FIRST actual executions of each compartment's enabled path — the unit suites
only do `sim.run(0)` static force checks — so they catch crash-on-enable bugs a
static read / coexistence test cannot.

## NOT production physics
- **No band pass/fail, no biological claim.** Each verdict is `PLUMBING_OK`
  (assembled + stepped + finite) or `SCALAR_OK` (runnable scalar laws only).
- **SMOKE-ONLY constants.** Where a compartment's activation constant is `None`
  (PI-pending), the harness uses the `COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09`
  candidate value, labelled `SMOKE-ONLY` at the call site. Production still raises
  until PI ratifies (see `docs/v2_audit/PLATFORM_PI_QUEUE.md`).
- **Physiological baseline.** Drag/kT are at the cytoplasm setpoint (η=65.9 Pa·s,
  kT=4.28e-21 J).
- Every figure carries a `SMOKE` watermark. Outputs go to
  `outputs/compartment_smoke/{figs,json}/`.

## Run
```
python ffn_sim/scripts/compartment_smoke/run_all.py        # all 8 + summary/montage
python ffn_sim/scripts/compartment_smoke/<name>.py         # one compartment
```

## Coverage (8/8)
| harness | kind | what it validates |
|---|---|---|
| `microtubules` | full HOOMD step | stiff aster, CFL gate, rod stability |
| `intermediate_filaments` | full HOOMD step | linear cage, per-r0-bin crosslinks on shared Harmonic, force-free |
| `linc` | full HOOMD step | perinuclear-cap bridges (no zero-bonds trap), force-free |
| `osmotic_regulation` | full HOOMD step | live EnclosedVolumePressure updater (object-identity), RVD sign |
| `cadherin_junction` | full HOOMD step | catch-bond binder forms A↔B trans-dimers (no intra-cell) |
| `stress_fibers` | full HOOMD step | FA→FA passive bundles, 4 sf_ bond types, force-free (NMII off) |
| `membrane_reservoir` | static mesh | own mem_node layer (BLOCKER-1), cross-layer tethers (rupture updater PI-blocked) |
| `junctional_actin` | scalar only | catch-slip shape + engaged fraction (HOOMD build is a STUB → blocked) |

## Found + fixed
- `junction/cadherin.py extend_snapshot_with_cadherins` did not carry the
  per-particle `image` array → state-creation crash on any host frame with image.
  Fixed + regression test (`test_snapshot_extender_extends_velocity_and_image_arrays`).

`_smoke_common.py` holds the shared plumbing (BAOAB stepping sim, inert cortex
stub, gsd round-trip, gamma-map, watermark/fig/json helpers).
