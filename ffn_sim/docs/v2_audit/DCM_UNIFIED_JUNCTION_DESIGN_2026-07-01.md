# DCM cell-cell interaction → ONE junction-mediated interface (new PI principle, 2026-07-01)

**New PI principle** ([[feedback-cell-cell-contact-junction-mediated]]): contact / adhesion /
interfacial-tension / compaction are **ALL one junction-mediated interface** — do NOT split them into
separate lumped forces. Stages (aggregation / compaction / maturation) modulate the **junction STATE /
coupling**, not the mechanism. Generalizes the junction-switch-fine-grained rule.

## The gap — the session's faceted foam uses SEPARATE lumped forces

The aggregation/proliferation work this session built the cell-cell interaction as **several independent
kernels** (confirmed by reading `dcm_warp_decohesion.py`):

| current separate force | kernel | what it lumps |
|---|---|---|
| mechanical contact (rep + adh tent) | `contact_grid_conservative_kernel` / `contact_grid_kernel` | non-penetration + adhesion |
| node-node cohesion | `cohesion_grid_kernel` / `_cad_kernel` | adhesion |
| interfacial tension (Maître/DAH) | `differential_surface_tension_kernel` / `surface_tension_kernel` | γ_ij |
| cadherin catch-bonds | `CadherinBondHost` (Rakshit) | force-dependent adhesion / de-cohesion |

These are superposed as independent forces (the confluent-foam production run used conservative-tent +
differential-γ; the de-cohesion run added cadherin on top). **This violates the new principle** — the
faceting, contact, cohesion and de-cohesion should EMERGE from ONE junction, not four coupled kernels.

## Unification target — the cadherin junction IS the interface

One junction-mediated interface = the **cadherin trans-dimer bond** (Rakshit catch-slip, KU-4.2) is the
single cell-cell interface, and everything derives from its state:

- **Adhesion** = the bond's attractive force (already the cadherin force).
- **Interfacial tension γ_ij** = the junction's line/area energy: a densely-bonded interface has LOW
  γ (Young-Dupré wetting), a sparsely-bonded one HIGH — so faceting/compaction emerge from **bond
  density**, not a separate `differential_surface_tension` kernel. γ_ij = γ_cortex − w_adh(bond density).
- **Contact / non-penetration** = the steric floor of the SAME interface (the bond has a rest length
  r0 + a hard repulsive core below it), not a separate contact tent.
- **De-cohesion / compaction / maturation** = the junction STATE evolving (bond on/off-rate under load,
  density ramp), NOT new forces.

So the ~4 kernels collapse to: **one junction force law** f(bond state, separation) + a junction-state
updater (bond formation/rupture, density). Faceting = high-density junctions flattening contacts;
compaction = junction maturation; de-cohesion = load-driven rupture — all one mechanism, stage-modulated.

## Scope + recommendation (PI decision — this is an architecture reframe)

This is a **significant refactor of a WORKING system** (the faceted-foam production runs use the separate
forces and are validated). The clean path:

1. **Design the unified junction force law** (bond-state → adhesion + steric + γ_ij) grounded in KU-4.2
   (Rakshit) + the Young-Dupré / Manning γ_ij=γ−w_adh relation — one kernel replacing the four.
2. **Parity-gate it against the current foam:** the unified junction at the same physiological bond
   density must reproduce the validated confluent-foam faceting (Q~149, V/V0 1, pen clean) BEFORE
   replacing the separate forces (no regression).
3. Then re-run aggregation → proliferation → (spreading) on the unified junction; stages = junction
   state only.

**Recommend PI confirm the scope/direction before the rebuild** — it re-architects the cell-cell model
the session just validated. Under the 8h mandate I will START step 1 (the unified junction force-law
design) which is non-destructive; the destructive swap (steps 2-3, replacing the working kernels) is the
point to confirm. No magic-number; the unified law's constants are the same lit KU-4.2 / w_cs / γ values.

## The unified junction force law (non-destructive design, drafted 2026-07-01)

State that defines the interface = the local **bond areal density** ρ_b(x) [bonds/µm²] on each apposed
cell-cell face, from the existing `CadherinBondHost` (Rakshit catch-slip, KU-4.2). Everything derives
from ρ_b and the node-face separation d:

**(a) Adhesion** — each bond is a trans-dimer spring with the Rakshit catch-slip off-rate. Per-area
adhesive traction σ_adh = ρ_b · f_bond(d, load), f_bond = k_trans·(d − r0) for d > r0 (pull toward
contact). This IS the current cadherin force — unchanged, it just also drives (b) and (c).

**(b) Steric floor (non-penetration)** — the SAME interface's excluded-volume core below the bond rest
length: σ_rep = k_steric·(r0 − d) for d < r0. One law with (a): a single traction curve σ(d) that is
repulsive for d<r0 and adhesive for d>r0 — replacing the separate `contact_grid_conservative` tent (whose
rep/adh split is exactly this, but decoupled from ρ_b). k_steric from the cortex compressibility (derived).

**(c) Interfacial tension γ_ij — EMERGES from ρ_b (this is the key unification)** — the junction's energy
per unit area is γ_ij = γ_cortex − w_adh(ρ_b), with the adhesion energy density w_adh(ρ_b) = ρ_b · ε_bond
(ε_bond = the trans-dimer binding energy, ~w_cs at saturating density). So:
  - **free surface** (ρ_b = 0): γ_ij = γ_cortex (HIGH) → the face rounds (area-minimising).
  - **bonded contact** (ρ_b high): γ_ij = γ_cortex − ρ_b·ε_bond (LOW) → the contact WETS/flattens
    (Young-Dupré / Manning) → a flat facet.
The force is the usual area-gradient f = −γ_ij·∂A/∂r, but γ_ij is now **per-face from the local bond
density**, NOT a separate `differential_surface_tension` field. **Faceting emerges** from the γ_ij
contrast between bonded (low) and free (high) faces — the exact thing the separate differential-tension
kernel imposed by hand is now a CONSEQUENCE of where bonds are.

**(d) Stages = ρ_b dynamics only** — aggregation: bonds form (ρ_b: 0→); compaction/maturation: ρ_b
ramps up (γ_ij drops further, facets sharpen, Q rises); de-cohesion: load ruptures bonds (ρ_b↓, γ_ij↑,
contact opens). No new forces per stage — only the junction STATE ρ_b evolves. This is precisely the PI
principle: one mechanism, stage-modulated.

**One kernel** `junction_force_kernel(pos, faces, cof, bond_density, params) → force` replaces the four.
Inputs: the CadherinBondHost per-face bond density + the node-face geometry. Outputs: σ(d) traction
(rep+adh) scattered barycentrically + the γ_ij(ρ_b) area-gradient. Parity target: at the physiological
saturating ρ_b the confluent foam must reproduce Q~149 / V/V0 1 / pen-clean before the swap.

**Testable prediction that distinguishes it from the separate-force model:** faceting Q should scale
CONTINUOUSLY with bond density (a maturation sweep ρ_b: low→high gives Q: 113→~150), because γ_ij is
ρ_b-derived — a single knob (junction maturation) spans aggregation→compaction. That sweep is the
clean validation once the kernel exists (still PI-scope-gated for the destructive swap).
