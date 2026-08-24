# DCM SimuCell3D upgrades — native mesh shell + bilinear-tent contact

**2026-06-11 · branch `h7/compartment-platform` · commit `2bf782b`.**
Reframed "Phase 2" from the SimuCell3D integration plan
(`docs/v2_audit/_historical/SIMUCELL3D_INTEGRATION_2026-06-11.md`): replace the crude node-LJ cell-cell
contact (which over-dispersed cells → A/A₀ blew to **17**, unphysical) and move the shell to
GPU-native HOOMD forces. **Result: A/A₀ 17 → 2.36 (into the physical band), independently
reproduced.**

## Upgrade A — native mesh shell (`cell/dcm_native_shell.py`)
- **Turgor/volume** → `md.mesh.conservation.Volume` (GPU-native), `U=k(V−V₀)²/2V₀` ⇒
  `p=−K(V−V₀)/V₀` with `k=K·V₀`, K=osmotic bulk modulus = **2500 Pa** (mcf7_p0.xml).
  Per-cell mesh triangle type (`cell{c}`) → each cell conserves its own volume (multi-type
  verified). Replaces the custom `DcmTurgorForce` (so this term no longer forces a per-step
  CPU sync — it is GPU-native).
- **Edges** → existing `md.bond.Harmonic` (k_edge=1e-3 N/m). **Bending** →
  `md.mesh.bending.BendingRigidity`, default OFF (SimuCell3D bending_modulus=0).
- **HOOMD winding gotcha (solved):** HOOMD's mesh `Volume` winding convention is OPPOSITE
  the divergence theorem — an outward-wound mesh gives HOOMD a negative V (rejects V₀≤0);
  flipping the winding trips the mesh bond-builder ("same particle twice in a triangle").
  Fix: keep the consistently-ordered `SurfaceManifold.icosphere` triangulation and **reflect
  vertex x → −x** (handedness flips → HOOMD V positive, builder happy; a reflected sphere is
  still a sphere, physics invariant). V₀ = measured triangulated volume (not the sphere
  formula: subdiv-1 is 13% under sphere, subdiv-2 3%).

## Upgrade B — bilinear-tent contact (`cell/dcm_contact.py`, `DcmTentContact`)
- Symmetric traction-separation **tent in separation d, PEAK at d=c_adh/2**: softening
  `|F|=ω·A·d` (0≤d<c/2), hardening `|F|=ω·A·(c−d)` (c/2≤d<c), repulsion `|F|=ξ·A·(r_contact−d)`
  (overlap). ω=ξ=1e8 Pa/m, c_adh=5e-7 m (mcf7). Matches SimuCell3D Eq.2 (the Lead-corrected
  force law — NOT a 1/d falloff).
- **Geometry (honest):** node↔(other-cell node) neighbour search (scipy cKDTree), inter-cell
  only via mutable `cell_of_node`. The junction-switch seam is preserved
  (`√(cad_mult_i·cad_mult_j)`); all forces capped (5e-8 N). The exact node↔face
  closest-point (Ericson) is the documented next upgrade — but the A/A₀=17 explosion was a
  contact-LAW problem, not a geometry problem, so this first cut already fixes it.

## Validation (CPU BAOAB, independently re-run)
| check | result |
|---|---|
| V1 single cell — native Volume hold | **PASS** V/V₀ 0.95→0.94 (finite, stable 5000 steps) |
| V2 two cells — adhesion | **PASS** 14 contact nodes, sep≈2.05R, no fly-apart/interpenetration |
| V3 14-cell 3D spheroid — A/A₀ | **PASS** A/A₀ [2.09→2.36], **FINAL 2.36** (band ~1.5-4; old node-LJ 17) |
Figure: `figs/dcm_native_validations.png`.

## Honest limits / next
- V₁ settles at 0.94 not exactly 1.0 (edge-spring vs turgor balance) — finite/stable; exact
  setpoint needs `turgor_inflate` or k_edge recalibration.
- node↔node contact is a first cut (not exact non-penetration); node↔face Ericson = the
  documented Phase-1 geometry upgrade.
- The native mesh shell is now GPU-native; the remaining custom forces (tent contact,
  substrate, FA) + the BAOAB per-step sync are the residual GPU items for the large-spheroid
  production run (future). dt=3e-10 and 1e-10 both finite.
