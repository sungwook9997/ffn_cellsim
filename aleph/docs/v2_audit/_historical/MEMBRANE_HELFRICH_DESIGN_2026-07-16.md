---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Membrane Helfrich-bending + real-sheet wiring — verified design (2026-07-16)

**PI (2026-07-16):** the cell membrane is not mechanically OK. Audit (code-grounded, 100% confirmed):
the PRODUCTION cell applies the membrane as a LUMPED uniform Young-Laplace tension `dP_mem=2γ_mem/R_mean`
inward on cortex nodes (Template 1) — not its own surface; the real sheet `ff/membrane_surface.py`
(Template 2: own icosphere, area tension, ERM tether, containment) exists but is UN-WIRED in production
AND has NO Helfrich bending kernel; the only surface-bending in the repo is DCM's biharmonic (wrong operator
for a lipid membrane). So the DEFINING membrane mechanic (bending) is genuinely absent. This is the design
(agent-produced, to implement) to fix it. Key files: `ff/membrane_surface.py`, `ff/network_warp.py`.

## 1. Helfrich bending kernel — Seung–Nelson dihedral hinge
Per interior edge shared by two outward-wound triangles with unit normals n̂₁,n̂₂, flat reference (C0=0):
**E_e = κ̃·(1 − n̂₁·n̂₂)**, summed over edges = discrete Helfrich/Willmore. **κ̃ = 2κ_m/√3 = 0.0956 pN·µm**
(κ_m=0.0828 pN·µm=20 kBT, Rawicz 2000 KB-3.B1.2; the √3 factor is Seung–Nelson 1988 — derived, grid-invariant).

Force = −∂E/∂x. With q₁=n̂₂−(n̂₁·n̂₂)n̂₁, q₂=n̂₁−(n̂₁·n̂₂)n̂₂, A₁,A₂ the triangle areas:
- f_i = −κ̃/(2A₁)·(x_k−x_j)×q₁ − κ̃/(2A₂)·(x_j−x_l)×q₂
- f_j = −κ̃/(2A₁)·(x_i−x_k)×q₁ − κ̃/(2A₂)·(x_l−x_i)×q₂
- f_k = −κ̃/(2A₁)·(x_j−x_i)×q₁
- f_l = −κ̃/(2A₂)·(x_i−x_j)×q₂

(i,j = shared-edge endpoints; k = flap of T1=(i,j,k); l = flap of T2=(j,i,l).) Σf=0 (translation-invariant,
verified analytically); flat hinge → q=0 → zero force. Warp kernel `helfrich_bending_kernel(mpos, hinges,
kappa_tilde, force)` (one thread per hinge, atomic_add, float64) — mirror `forces_warp.cytosim_bending_kernel`.

**Hinge adjacency** `build_membrane_hinges(faces)→(Nh,4)[i,j,k,l]` from DIRECTED edges (keeps each face's native
outward winding → consistent normals — the correctness-critical detail). Closed icosphere: Nh = 3·Ntri/2.

**Validation (`tests/test_membrane_bending.py`, MANDATORY before use):** (1) FINITE-DIFFERENCE the energy →
match −∂E/∂x to <1e-6·‖f‖ (the arbiter of the sign); (2) sphere energy = **8πκ_m = 2.081 pN·µm**, R-independent
(catches winding bugs); (3) flat patch → 0 force; (4) a displaced vertex → finite-width DIMPLE (width ~√(κ_m/σ)),
NOT a spike (the physical acceptance the lumped tension can't give).

## 2. Wiring Template 2 into production (additive, `membrane_surface=None` → bit-identical)
New kwarg `membrane_surface=None` (distinct from the lumped `membrane`). Layout appends a membrane block:
`[cortex(Nc); MT; MTOC; nucleus; membrane(Nm)]`, Nmem0=N_old. Build the sheet from the RESTING cortex
(`build_membrane_mesh(R0+gap, subdiv=4)`), hinges/faces/ERM/containment offset to global indices. Step loop
(both loops + final read): replace the lumped `turgor_kernel(dP_mem_area)` with `membrane_area_kernel` +
`helfrich_bending_kernel` + `erm_tether_kernel` (two-way, ruptures→bleb) + `membrane_containment_kernel`.
When `membrane_surface` set → force `dP_mem_area=0` (the sheet REPLACES the lumped tension; no double-count).
Cortex osmotic turgor `dP_area` unchanged. `reshape_kernel` never touches membrane nodes (not fibers). All
integrator/plate launches dim=N→N' (the AFM plate touches the membrane first).

**Variant A (default):** cortex keeps osmotic turgor; membrane replaces only the lumped tension channel.
**Variant B (PI-gated):** membrane becomes the osmotic envelope (moves dP_area onto the sheet) — larger change,
re-validate the ΔP/γ gate end-to-end.

## 3. CFL / KB / magic numbers (honest)
- κ_m ✓ (KB-3.B1.2), κ̃ ✓ derived, γ_mem/γ_MCA/K_A/f_rupt ✓ (KB-3.B1.1/3/4). Bending is SOFT at subdiv-4
  (κ̃/ℓ³≈1.5 pN/µm) ≪ cortex kmax → dt unaffected. Add all membrane k to kmax; k_wall (containment) is a
  numerical penalty. ⚠ **K_A=2.35e5**: if the reservoir exhausts, kmax jumps ~1000× → dt collapses → keep the
  buffered plateau (K_A off) default, or implicit membrane sub-step, or PI. Never soften silently.
- ⚠ **RESERVOIR_STRAIN=0.60** unsourced magic number (already flagged; MCF7 caveolae-deficient) — PI. Inert
  while K_A deferred. ⚠ **k_erm=50 pN/µm** (CFL-convenience, not single-molecule ezrin) — PI registration.

## 4. Native gate (CLAUDE.md HARD)
Regression: `membrane_surface=None` bit-identical. Bending unit tests (§1). The historical landed gate used
`--cortex-fil 38000`; new-engine reuse must run at the **I0-A-ratified 70,686-filament cortical density** with
`--from-resting --microtubules` + sheet on the production GPU: containment (no actin leak,
straggler count→0); bending sets a finite-width bleb (not a spike, width ~√(κ_m/σ)); tether rupture → bleb
EMERGENT (n_bleb>0); areal strain within the plateau. Interactive 3D HTML: cortex + translucent membrane
sheet + nucleus, full-res, with a bleb visible.

## 5. Risks: dihedral force sign (FD gate is arbiter); double-count if both `membrane`+`membrane_surface`
(guard dP_mem_area=0); containment resolution (per-cortex-node nearest-membrane radius, not a scalar envelope,
else leak); ERM rest built from `--from-resting` (fresh sphere → spurious rupture).

## Change log
- 2026-07-16: created from the verified membrane-design agent. Exact Seung–Nelson dihedral Helfrich kernel +
  force + validation + Template-2 wiring. To implement (with the IF cage) per PI "막도 같이 가자".
- 2026-07-16 (BUILT + native-landed): implemented `helfrich_bending_kernel` + `build_membrane_hinges` +
  `calibrate_kappa_tilde` (membrane_surface.py) + `membrane_surface` channel (network_warp.py) + MANDATORY
  gate `tests/test_membrane_bending.py` (5/5 PASS). ⚠ **CALIBRATION CORRECTION**: the §1 claim "κ̃=2κ/√3 →
  sphere 8πκ" is WRONG for the icosphere — that factor is for the ideal regular lattice (Σ=4√3π≈21.77); the
  icosphere Σ(1−cosθ) converges to **7.50** (measured, resolution/R-independent). Fix (grid-invariant, not a
  fit): κ̃ = 8π·κ_m/Σ_ref, calibrated to the resting mesh. FD-gradient is the force-sign arbiter (PASS).
  Containment deferred (scalar-envelope kernel is broken under non-spherical load, design §5) — ERM tethers
  hold the sheet on the cortex. **Native gate** (Nc=494802, `--mt-reach --membrane-surface`): rest bending
  1.067×8πκ ✓; loaded bending rises 4.13× (s0.30) → 6.56× (s0.45); areal strain −0.05/−0.11 (plateau);
  **ERM rupture→bleb EMERGENT at native** (0 blebs @ s0.30 clean ride → 55 @ s0.45). Viewer browser-verified.
  Variant A wired (cortex keeps turgor; sheet replaces the lumped tension). K_A upturn + full containment = TODO.
