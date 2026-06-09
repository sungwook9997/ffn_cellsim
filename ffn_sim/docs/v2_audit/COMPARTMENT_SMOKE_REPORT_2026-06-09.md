# Compartment ENABLED-PATH Smoke Report (2026-06-09)

**Session:** autonomous platform-deepening, branch `h7/compartment-platform`
(sibling Gate-B session on `h7/full-cell-integration` untouched).
**Scope:** safely deepen the breadth of the 8 default-OFF compartments WITHOUT PI
parameter ratification or shared-file edits. **Deliverable:** the first actual
executions of every compartment's *enabled* build/attach/step path.

> These are **PLUMBING** results — crash/finite/assembly checks with `SMOKE-ONLY`
> candidate constants. **No band pass/fail, no biological claim.** Activation
> blockers (PI/shared-file) are in `PLATFORM_PI_QUEUE.md`; the day-log is
> `PLATFORM_AUTORUN_LOG_2026-06-09.md`.

## Why this matters
The 2-pass physics audit (`COMPARTMENT_PHYSICS_AUDIT_2026-06-09.md`) read the code
and ran `sim.run(0)` coexistence checks, catching 3 crash-on-enable bugs. It did
NOT actually *step* each enabled compartment. This session built standalone
harnesses (`scripts/compartment_smoke/`) that assemble each compartment + an inert
cortex stub, step with the project L-M BAOAB integrator at the physiological
cytoplasm setpoint (η=65.9 Pa·s, kT=4.28e-21 J), and measure a coarse observable.

## Results (9/9 green)

| compartment | kind | verdict | key plumbing result |
|---|---|---|---|
| microtubules | full HOOMD step | PLUMBING_OK | stiff aster, CFL gate honoured, backbone rms dev 0.03% |
| intermediate_filaments | full HOOMD step | PLUMBING_OK | if_backbone + 13 per-r0-bin crosslinks on ONE shared Harmonic; construction 54 kT (force-free) |
| linc | full HOOMD step | PLUMBING_OK | perinuclear-cap → 80 bridges (no zero-bonds trap), force-free (0.0 kT) |
| osmotic_regulation | full HOOMD step | PLUMBING_OK | live EnclosedVolumePressure updater (object-identity), 60 ticks, RVD sign correct |
| cadherin_junction | full HOOMD step | PLUMBING_OK | catch-bond binder → 10 trans-dimers ALL A↔B (0 intra-cell) |
| ventral_stress_fibers | full HOOMD step | PLUMBING_OK | 6 FA→FA bundles, 4/4 sf_ bond types, force-free (passive; NMII off) |
| membrane_reservoir | static mesh | PLUMBING_OK | own mem_node layer → 120 cross-layer tethers (rupture updater PI-blocked) |
| junctional_actin | scalar only | SCALAR_OK | catch-slip F*=6.0 pN, engaged f=0.91; HOOMD build BLOCKED (STUB) |
| surface_manifold | geometry only | GEOMETRY_OK | area→4πR², k-ring reach-coverage master gate = 1.0 at all resolutions |

## What the smoke caught (value-add over the audit)

1. **Crash-on-enable bug FIXED — cadherin `image` array.**
   `extend_snapshot_with_cadherins` extended position/typeid/mass/charge/diameter/
   velocity but NOT the per-particle `image` array → `Snapshot.from_gsd_frame`
   crash (broadcast `(n,3)→(N,3)`) on any host frame carrying an image array
   (every real cortex build-time frame does). Fixed in `junction/cadherin.py`
   (carry image int32, mirroring velocity + the nucleus/MT/IF extenders) +
   regression test `test_snapshot_extender_extends_velocity_and_image_arrays`.
   This is the class of bug only an actual state-creation pass reveals.

2. **Surface-manifold grid-invariance master gate confirmed.** The DERIVED
   `kring_for_reach` keeps the geodesic candidate ring covering the physical reach
   ball (coverage = 1.0) across icosphere subdivisions 0-4 (k adapts 4→8) — the
   grid-invariance the registry requires before any manifold force is sanctioned.
   Done geometry-only (zero force/bond/particle, cross-checked vs the registry
   `geometry_only` invariant).

3. **Honest blockers documented, not bypassed.** junctional_actin's HOOMD build
   raises even when anchored (STUB), and membrane_reservoir's bleb-rupture updater
   raises unconditionally — both were left blocked and recorded for PI rather than
   forced. SF NMII was left off (BLOCKER #2, the `sf_myosin_*` prefix).

## Regression guard
`tests/test_compartment_enabled_canary.py` converts the enabled-path plumbing into
a fast CI guard (MT/IF/LINC minimal assemble+step, ~1 s). Full compartment suite +
canary: **269 passed, 1 skipped**.

## Still PI-gated (see PLATFORM_PI_QUEUE.md)
γ-denylist extension in `cortex/cortical_tension.py`; SF `sf_myosin_*` prefix split
in `cortex/myosin.py`; the 5-point loader wiring (`cell.py`/`manifest.py`);
None-gated parameter ratifications (N_filaments, k_linc, σ_crit_bleb, f_excess, IF
nonlinear curve, MT Y_stretch/L_mt, osmotic Lp); junctional_actin build path; the
two-cell doublet builder for cadherin GATE-J.

## Figures
All under `ffn_sim/outputs/compartment_smoke/figs/` (SMOKE-watermarked):
- `smoke_ALL_montage.png` — suite verdict board (9/9 green).
- `smoke_microtubules.png` — backbone-length distribution + aster x-z projection.
- `smoke_intermediate_filaments.png` — IF backbone lengths + perinuclear cage (x-y).
- `smoke_linc.png` — bridge lengths at r0 + nucleus/acceptor shells.
- `smoke_osmotic_regulation.png` — V0(t) RVD slope + live pressure trace.
- `smoke_cadherin_junction.png` — A↔B trans-dimers + per-dimer tension hist.
- `smoke_stress_fibers.png` — SF backbone lengths + FA→FA bundle layout.
- `smoke_membrane_reservoir.png` — tether lengths + cortex/mem_node cross-layer mesh.
- `smoke_junctional_actin.png` — tensile coupling law + catch-slip off-rate (F* marked).
- `smoke_surface_manifold.png` — area convergence, frame convergence, master-gate coverage.

### Full-cell MORPHOLOGY (PI 2026-06-09 "visualize the constructed cell shape")
Single entry point `ffn_sim/scripts/compartment_vis.py` builds the FULL physiological
MCF7 cell with EVERY LIVE compartment ON (cortex+myosin+xlink spine · baseline
cytoplasm/turgor/nucleus/membrane-surface · osmotic · MT aster · IF cage · LINC
bridges; 19,656 particles) and renders the ACTUAL constructed geometry. Under
`ffn_sim/outputs/h7/figs/`:
- `compartment_cell_overview.png` — 4 panels: 3D scatter of the whole cell;
  equatorial (|z|<0.75 µm) cross-section with the explicit MT/IF/LINC bonds drawn;
  meridional (|y|<0.75 µm) cross-section; and the per-compartment radial density
  profile (each shell at its physiological radius: cortex at R_cell=7.5 µm, filled
  nucleus core to R_nuc, IF cage just outside the nucleus, MT aster from the MTOC).
- `compartment_cell_3d.png` — standalone larger 3D view (spherical cortex shell +
  internal organelles).
- `h7_linc_activation_gate.png` — LINC gate: nucleus↔IF bridges, per-bond EXACT-r0
  histogram (force-free), cortical γ_soft identical OFF vs ON (no contamination).
- `h7_membrane_reservoir_activation_gate.png` — membrane gate: own mem_node layer +
  cross-layer mem_tether mesh, force-free r0 histogram, cortical γ_soft identical OFF vs ON.
  The morphology figs now also render the mem_node membrane layer just outside the cortex.
