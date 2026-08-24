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

# H.1 / H.5 / H.7 — Manifold-Guided Explicit Mechanobiology — Contact-Manifold Architecture Design Note (2026-06-07)

## 1. Executive summary

- **Proposal.** A "manifold-guided explicit mechanobiology" layer for the cell↔ECM stack (H.1 ECM, H.5 lamellipodium, H.7 cortex): HOOMD stays the **sole explicit force-carrier engine**; a mesh/manifold layer supplies **only** geometry (local normal/tangent frames), contact planes, spatial partitioning / broad-phase, and visualization/export handles.
- **One-line principle.** Broad-phase (patch neighborhood / contact manifold) *narrows candidates*; narrow-phase (Bell-Evans / Hill / catch-bond / WLC on **explicit** bonds) makes the physical decision. Mesh = navmesh + spatial substrate; explicit filaments / ECM fibers / motors / clutches = physics actors.
- **Verdict on the full-fidelity rule: SATISFIED.** Every actin filament, myosin head/minifilament, crosslinker, integrin/talin/vinculin/clutch, and ECM fiber/crosslink stays an explicit particle/bond. The manifold abstracts **zero** filaments into edges; it adds **no** new coarse-graining. The only sanctioned coarse-graining remains the ×40 mesoscopic filament scale.
- **It does NOT lift γ.** The emergent active cortical tension floors ~1×10⁻⁴ mN/m vs the MCF7 datum ~0.27 mN/m (Hosseini 2020). The floor is **generation-limited** in two serial gates and is **upstream of measurement**; a spatial index, contact plane, or soft shell changes none of it.
- **It does NOT remove Gate A or Gate B.** Gate A (myosin heads bind+load to ~0.2–0.44 F_stall but don't walk; grip-stretch s_grip ≈ 0) and Gate B (the rigid M-SHAKE backbone forbids r/r0<1, capping at the dipole ceiling) are the only levers that move γ, and they live in myosin kinetics and the backbone constraint — never in geometry. The soft band *permits* the radial buckling excursion that a hard projection would forbid, but **permitting is not causing**: the manifold supplies no constraint relaxation and produces zero buckled filaments.
- **One new force only.** A soft, **normal-only**, one-sided, default-off shell confinement `U_conf` (`md.force.Custom`, sibling of `cortex/erm.py:ERMHarmonic`). It is provably invisible to the γ estimator (it is not a HOOMD bond, so it cannot enter the method-of-planes sum — the same mechanism by which ERM is deliberately absent).
- **One unified new constant.** `k_conf`, anchored as a fraction of the existing literature-anchored cortex bond stiffness, **k_conf ∈ [1×10⁻³, 1×10⁻²]·bond_k = [3×10⁻⁶, 3×10⁻⁵] N/m** (bond_k = μ/ℓ₀ = 3.0×10⁻³ N/m, verified). This bracket also contains the band-width-fluctuation derivation (6.8×10⁻⁶ N/m ⇒ σ_radial ≈ 25 nm).
- **v0 data model is mesh-free.** No stored triangulation/Delaunay/face-area/edge-adjacency in v0 — only (i) per-bead PCA normal/tangent from `cell/membrane_surface.py:discrete_mean_curvature`, (ii) a per-bead **patch label** (Fibonacci/HEALPix index over (θ,φ)), (iii) a patch→bead lookup for broad-phase. Any true triangulation is render-side only, tagged "geometry, not physics," and forbidden from feeding any reported number.
- **Authoritative γ is unchanged.** It stays `cortex/cortical_tension.py:measure_cortical_tension` verbatim (3 channels: soft method-of-planes, IK whole-shell virial, rigid; B3 separation preserved). Any per-patch stress map is a **non-gate diagnostic** that must sum **exactly** to `gamma_ik` by construction (same bonds, same global 8πR² normalization) — there is no tolerated "quadrature gap."
- **Resolution-invariance is the master gate (VG-1).** Every physical observable (γ, traction, footprint, coordination z, min(r/r0)) must be invariant to patch count `N_patch`; any monotone drift = the manifold doing physics = **halt and surface to PI**.
- **Recommended additive v0.** A 9-step, default-off, `None`-guarded prototype (geometry metadata → broad-phase candidate index → patch-stress/traction *reporting* → optional soft confinement, default-OFF). With the manifold absent, every existing run/test is bit-for-bit identical. The integrator (`integrator/baoab*.py`) is never touched.
- **The manifold earns its keep only as infrastructure** — a better spatial home, contact model, conservative broad-phase, and visualization for the *same* explicit mechanisms. If any physical observable moves with resolution, it is reverted at zero cost.

## 2. Relationship to the cortex-as-mesh REJECTION

This note is **subordinate to** the Lead's rejection in `aleph/docs/v2_audit/H7_CORTEX_AS_MESH_2026-06-07.md` (commit `7b19276`). That decision refused to rebuild the cortex as a **triangulated / icosphere / deformable surface whose edges carry force**, for three reasons that this note adopts wholesale:

1. **A force-bearing surface edge is a 2nd unsanctioned coarse-graining** — one edge abstracts many filaments, and its edge stiffness / area modulus have **no per-filament anchor**, migrating magic numbers.
2. **A closed triangulation is over-constrained** → it **suppresses single-filament buckling**, which is the project's ~10× active-stress amplification lever (Gate B).
3. **It does nothing for the γ-floor**, which is **generation-limited**, not topology-limited.

This note does **not** reopen that. It proposes the orthogonal, complementary thing the rejection left on the table: the mesh/manifold as **geometry / broad-phase / contact-plane / visualization infrastructure**, with every force carrier preserved as an explicit HOOMD particle/bond.

### The bright line (forbidden vs proposed)

| | Mesh-as-physics (FORBIDDEN — already rejected by `7b19276`) | Mesh-as-geometry / broad-phase (PROPOSED here) |
|---|---|---|
| Stores | edge springs, area modulus, force-bearing faces | per-bead normal/tangent frame, patch label, patch→bead lookup |
| Carries force? | **yes** (edge/area elastic energy) | **no in-plane force** — at most one optional soft **normal-only** confinement |
| Replaces explicit objects? | replaces filaments/clutches | replaces **nothing**; explicit filaments/myosin/integrin/clutch untouched |
| Sets γ / traction / spreading? | edge stiffness tuned → γ | **cannot**: adds no in-plane tension, is not a bond, invisible to the estimator |
| Coarse-grains? | edge = many filaments (unsanctioned) | metadata over existing beads; **zero** new CG |
| Buckling lever (Gate B) | suppressed (over-constraint) | untouched (no edge resists r/r0<1) |

The bright line in one sentence: **a mesh edge that pulls two nodes is physics (forbidden); a manifold that says "this bead's normal points here, search for binding partners within this patch, and softly keep this bead inside the 200 nm band" is geometry/broad-phase (proposed).**

**The Lead's sanctioned topology option — "Mikado-on-sphere"** (explicit bead-chain fibers woven on the sphere, crosslinks at geodesic intersections at the physical ~60 nm scale, soft buckling-capable edges, reusing `ecm/mikado.py`) — is **complementary** to and **unaffected by** this note. The manifold is the spatial substrate; Mikado-on-sphere is one allowed way to lay explicit fibers onto it. It would consume the manifold's normal/patch services unchanged. A triangulated icosphere-as-physics remains forbidden.

**Plainly: the manifold does NOT and CANNOT lift γ, and it removes neither Gate A nor Gate B.** Any deliverable citing this note must restate that.

## 3. Why pure 3D bead-filament HOOMD is resource/geometrically inefficient at cell scale

The cell is a **thin shell in a large box**, and naive 3D bead-filament dynamics pays for the volume it does not occupy:

- **×40 inter-filament spacing ≫ binder reach.** The ×40 mesoscale puts ~1000 cortical filaments on a R≈7.5 µm sphere → mean inter-filament spacing **√(A/n) ≈ 1880 nm**, while the physical crosslinker/integrin reach is **~60 nm**. A naive dynamic radius query over all beads at the physical reach finds almost nothing per query; the candidate-generation cost is dominated by traversal of empty space.
- **All-pairs vs surface-local.** Today's binders rebuild a **global** `scipy.spatial.cKDTree` over all actin every batch tick (`cortex/crosslinkers.py:XlinkBondUpdater.act`, `cortex/myosin.py:MyosinStepUpdater`), and the ECM already had to choose `md.nlist.Tree` (BVH) over `md.nlist.Cell` because a uniform cell grid at L/r_cut ≈ 1800 pre-allocates ~5×10⁹ cells and **segfaulted** in the v2 trial (`ecm/mikado.py`). A surface-local patch partition narrows each query to O(patch) candidates instead of O(N).
- **Thin-shell vs full-volume gas.** The cortex occupies a 200 nm band at radius 7.5 µm — a vanishing fraction of the bounding box. A manifold restricted to the shell (cortex) or to the apposed contact zone (cell↔ECM) is the right partition for a structure that is geometrically 2D embedded in 3D.

**Honest caveat (the ×40 erodes the cortex win).** The broad-phase neutrality proofs require `patch_radius ≥` the binder reach. For the cortex that reach is the **×40 coarse-graining-compensated** √(A/n) ≈ 1880 nm (not the physical 60 nm). With only ~1000 mesoscale beads and R≈7.5 µm, a partition whose patches are ≥ 1880 nm radius leaves "patch ∪ 1-ring" close to the whole cap — so the cortex candidate-reduction is **marginal at ×40**. The real broad-phase payoff is at **native (38k filament) scale** and for the **ECM (~66k Mikado beads)** and **cell↔ECM contact** searches, where the occupied sub-volume is a small fraction of the box. v0 must report the realized candidate-count reduction; if it is near-unity for the cortex at ×40, the manifold's cortex value is the **local frame + contact model + visualization**, not speed.

## 4. What the mesh/manifold is ALLOWED to do

1. **Local coordinate frames.** A per-bead / per-patch outward **normal n̂** and **tangent basis (t̂₁, t̂₂)**, computed **mesh-free** from live positions by `cell/membrane_surface.py:discrete_mean_curvature` (PCA tangent-plane quadric fit; smallest-eigenvalue axis = normal, oriented outward). Used to decompose **already-computed explicit forces** into normal/tangential (or radial/circumferential) components for measurement and for candidate generation.
2. **Contact planes.** A substrate/ECM contact plane (flat dish `z = z₀`, n̂ = +ẑ) or a per-contact local tangent plane to the nearest ECM fiber segment, derived from explicit bead positions each query, carrying no stored stiffness.
3. **Spatial partitioning / broad-phase.** A patch label per bead (Fibonacci/HEALPix over (θ,φ)) + a patch→bead lookup, and/or a BVH over patch centroids, returning the **candidate set** of nearby partners. Narrow-phase acceptance is unchanged.
4. **Soft confinement (default-off).** One **normal-only**, one-sided, flat-bottomed radial well `U_conf` keeping already-explicit beads inside the physical 200 nm cortex band, soft enough to leave thermal radial fluctuation, tangential slide, and inward buckling free.
5. **Field placement (geometry, not values).** Spatially *placing* the already-required per-condition substrate values (ligand density `ρ_lig`, compliance `k_sub`) as maps on the contact plane — the manifold does not invent values; the explicit `ecm/substrate.py:SubstrateAnchorSpring` still applies the force.
6. **Visualization / export handles.** Render meshes, per-patch field overlays, broad-phase partition glyphs — read-only over HOOMD positions, provably unable to affect dynamics.

## 5. What the mesh/manifold is FORBIDDEN to do (mapped to guardrails)

| # | Forbidden | Guardrail |
|---|---|---|
| 1 | Become the physics replacing explicit objects | G1 — manifold replaces **nothing**; all filaments/myosin/clutch/ECM stay explicit |
| 2 | A mesh edge carrying elastic force (becoming a filament/area modulus) | G2 — manifold carries **zero** edge/in-plane force; the only force is the soft normal `U_conf` |
| 3 | Mesh/confinement stiffness tuned to fit γ | G3 — `U_conf` is normal-only and not a bond → provably invisible to γ; `k_conf` is set by a thermal/softness requirement, never tuned to γ |
| 4 | Mesh contact force replacing integrin/clutch bonds | G4 — all adhesion stays the explicit catch-bond clutch (`bridge/integrin_bonds.py`, `bridge/clutch_spatial.py`); manifold narrows candidates only |
| 5 | Mesh resolution controlling measured γ / traction / spreading | G5 — **master gate VG-1**: any monotone observable-vs-`N_patch` trend = halt; per-patch maps sum exactly to global, no tolerated gap |
| 6 | Crosslink/motor/integrin reach inflated to fit ×40 spacing | G6 — reach stays physical (60 nm) **or** the already-ratified ×40-compensated √(A/n) bridge; manifold introduces no new reach |
| 7 | Any constant without a Magic-Number Block | G7 — every introduced constant carries a block (§19) |
| 8 | Loosening a gate | G8 — no tolerance is loosened to make a manifold build pass; failing identities surface to PI |
| 9 | Breaking the physiological baseline | G9 — confinement OFF for validation; production runs on the turgor-pressurised, anchored-ECM, 65 Pa·s physiological state |
| 10 | Unilaterally changing a validation contract | G10 — ERM displacement, curved-shell γ, contact-model oracle, etc. are Open PI decisions (§22) |

## 6. H.1 ECM application

**Fibers + crosslinks stay EXPLICIT; mechanics stay emergent.** No change to the runtime: backbone stretch (`ecm/mikado.py` `md.bond.Harmonic`, k = μ/ℓ₀), bending (`md.angle.Harmonic`, Head–MacKintosh), crosslinks (`ecm/cross_links.py` `md.bond.Harmonic` at intersections), WCA excluded volume, frozen L-M BAOAB. `G_0`, the strain-stiffening exponent, and the non-affine 1/r² point-dipole decay (KU-1.30 #1/#2/#3) remain emergent outputs measured by `ecm/shear_protocol.py`; the manifold computes **none** of them and is **forbidden from being read in the KU-1.30 measurement path**.

**A continuum/FEM ECM is out of runtime scope** (it would replace explicit Mikado with mesh elasticity, suppressing the non-affine fiber buckling that produces KU-1.30 #2 strain-stiffening — the ECM analogue of Gate B). If ever wanted, it is a separate **validation oracle** to compare against, surfaced as an Open PI decision (§22), never swapped in.

**What the manifold MAY add (additive, default-off, soft, broad-phase-only):**

- **(6a) Substrate contact plane + local frame.** Analytic plane `z = z₀`, n̂ = +ẑ, tangent (x̂, ŷ). Pure geometry — the half-space the existing point-ligands already live on (`bridge/fa.py` places ligands at z=0). `z₀` is **inherited geometry**, not a free constant. `cell/membrane_surface.py` is NOT this primitive (it is the H.8 closed cortex shell); a flat substrate plane is genuinely new, but any confinement it exerts uses the soft `ecm/substrate.py:SubstrateAnchorSpring` / `cortex/erm.py:ERMHarmonic` external-field template — never a hard projection.
- **(6b) Ligand-density + compliance maps.** Scalar fields `ρ_lig(t̂₁,t̂₂)`, `k_sub(t̂₁,t̂₂)` that *place* the already-required per-condition values (`bridge/ligand_species.py`; `ecm/substrate.py:k_sub`, REQUIRED-no-default) as maps so patterned substrates (durotaxis stripes, the PI Lam4 partial-uniformity Lp≈40 µm pattern) become expressible. The manifold invents no values; the explicit `SubstrateAnchorSpring` still applies the force. The only new number is the **sample spacing `Δ_map`** — a discretization knob gated by resolution-invariance (degenerates to today's scalar for a uniform substrate, bit-identical).
- **(6c) Broad-phase BVH (highest value, lowest risk).** Direct precedent: `ecm/mikado.py` already chose `md.nlist.Tree` over `md.nlist.Cell` because the ECM is so sparse. A BVH over the occupied sub-volume returns the candidate ECM beads per cell-surface particle; **narrow-phase** (WCA, and the future integrin→fiber catch-bond clutch) is unchanged. The broad-phase cell size is **performance-only**, verified by an all-pairs superset cross-check.
- **(6d) Gel boundary.** Analytic half-space/box for a bounded gel/dish wall via the soft external-field template; the gel's *bulk* elasticity stays the explicit Mikado network.
- **(6e) Visualization scaffold.** Plane + ligand-density heatmap + partition + fiber network, read-only over positions.

**`ecm/substrate.py` relationship.** It is **already** the per-ligand explicit force the manifold sits *over*: `SubstrateLigandPin` (rigid, infinite-stiffness backing) and `SubstrateAnchorSpring` (`F = −k_sub(r − a)`, finite compliance). The manifold supplies the anchor geometry and the `k_sub` field *value*; the explicit spring still applies the force, `k_sub` stays REQUIRED-no-default, and the CFL gate `dt ≤ α·γ_ligand/k_sub` still fires. The `effective_E_sub` bridge stays an unratified diagnostic (k_sub→E_sub flat-punch vs Hertzian oracle is an Open PI decision, §22).

**The ECM is currently orphaned from the cell** (grep `actin_ecm` in `cell/`+`bridge/` → 0 hits; FA adheres to flat rigid-pinned point ligands). This manifold is the **seam that could couple them later** (integrin→Mikado-bead clutch) without coarse-graining — but the coupling itself is an Open PI decision (§22), not delivered here. ECM crosslinks remain permanent static bonds (no Bell-Evans/turnover); the manifold adds no remodeling.

## 7. H.5 lamellipodium application

The H.7 spherical-integration work already landed most of the local-frame machinery, so this is largely a **consolidation + observable layer over existing code**:

- **Per-WAVE local tangents already exist.** `cell/lamellipodium.py:generate_lamellipodium_layout(geometry=...)` dispatches to `cell/lamellipodium_basal_ring.py` (outward-radial tangent `t̂ = (cosφ, sinφ, 0)`, ring radius derived from the FA south-cap base radius √(h(2R−h))) and `cell/lamellipodium_polarized_patch.py` (forward tangent = polarization p̂ projected onto the local shell tangent plane; `_tangent_basis(normal)` already builds an orthonormal frame).
- **Region masks already exist** as the two geometry branches (basal-ring azimuthally-uniform; polarized-patch bounded (Δφ,Δθ) cap).
- **Arp2/3, capping, elongation are already local.** Branch daughter tangent = mother tangent rotated 72° about a random perpendicular; elongation advances along `tangent_of[be]` (a real bead appended at `r_be + ℓ₀·t̂`).
- **Mesh-free local normal/curvature already exists** (`discrete_mean_curvature`), the operator needed for any shell-bead frame.

**What the manifold adds (metadata/observable only):**

1. **A shared `SurfaceFrameProvider`** returning `(t̂_growth, n̂_surface)` for any WAVE/tip tag, with `n̂_surface` recomputed from live positions (not a stored mesh normal). This is what the capping updater's `cos θ` should consult instead of the hard-coded `n̂ = +ŷ` TODO (`lamellipodium.py:902`), which is wrong for `basal_ring`/`polarized_patch` once membrane-load is ON. **No rate law changes** — only the `sin θ` between t̂ and the *correct* per-WAVE n̂ is substituted; it stays dormant-correct at F=0. **Membrane-load stays OFF for curved geometries until PI approves** wiring it (touches the KU-5.x membrane module — Open PI decision, §22).
2. **Spreading area from the contact FOOTPRINT** (not a bead-cloud convex hull): explicit cortex/integrin beads in the basal contact zone (the same south-cap mask the FA uses) plus advancing barbed-tip beads projected to the basal plane; area = 2-D projected-hull/alpha-shape at probe radius `r_probe = FA capture_radius` (the existing physical contact scale; no new constant). Reported across an `r_probe` sweep and shown plateau-invariant.
3. **Azimuthal FA traction** as a first-class output: per explicit `fa_actin_clutch` bond compute tension `T = k(L−r₀)` and bin by azimuth φ; report `T(φ)` + first Fourier mode |a₁|/a₀ as the polarization order parameter (≈0 for basal_ring, >0 for polarized_patch). The manifold supplies only the azimuth coordinate; the force is explicit. **Measurement caveat (cross-referenced):** `bridge/fa_growth.py:171` currently computes `F = k·|Δr|` assuming r0=0; with the FA overload-fix per-bond r0-bin family the readout must become `F = k·(|Δr|−r0)`. That r0-subtraction belongs to the **FA overload-fix workstream**, is ratified together with the r0-bin family, and must **not** be folded into a geometry/broad-phase manifold PR (§22).

**No reach constant is introduced.** Branch reach `r_branch_eff = 100 nm` stays the physical Arp2/3 reach (PI 2026-05-29); broad-phase must reproduce the **identical** `d < r_branch_eff` candidate set. A broad-phase cell size altering which mothers are eligible would violate G6 and is forbidden.

**Physiological baseline preserved.** Basal-ring geometry is derived from the resting turgor-pressurised, FA-adhered shell. When membrane-load is eventually enabled it must be at the physiological `γ_mem` (H.8 band), not zero.

**The manifold cannot fix the real lamellipodium gap:** KU-5.1 branch density floored ~0/µm² in the pre-γ baseline because of branch-search reach vs ℓ₀ (a kinetics/reach problem the γ-branching coupling addresses) — a better spatial home for WAVE beads does nothing for it. On the **current** build the lamellipodium is essentially unpopulated, so any traction figure produced today is a near-zero / FA-pin-dominated readout, not evidence of working active mechanics; the readout becomes meaningful only after the branch-reach fix and Gate A/B.

## 8. H.7 cortex application

**The shell already exists implicitly** — `cortex/cortex.py:_project_to_shell_band` projects beads onto the band `[R_cell − thickness, R_cell]` (cortex_thickness = 200 nm, KU-3.17), held at runtime by `cortex/erm.py:ERMHarmonic` (radial spring, k_ERM = 0.1 N/m, KU-3.18). The ERM is a **single-radius, stiff** pin: σ_radial = √(kT/k_ERM) = **0.21 nm** at the canonical kT = 4.28×10⁻²¹ J (note: the erm.py docstring's "0.65 nm" is a pre-existing kT-error to fix separately, not inherited here). The finite-thickness manifold is its **band-valued, soft** generalisation.

**Why soft, not hard.** Hard radial projection every step would re-freeze exactly the radial degree of freedom Gate B needs — the inward buckling/condensation excursion (r/r0<1). **Hard projection is forbidden** precisely because it re-commits the over-constraint sin the cortex-as-mesh rejection warns about.

**Soft confinement `U_conf` (additive, default-off, normal-only).** A `md.force.Custom` sibling of `ERMHarmonic`, a flat-bottomed radial well over `r_i = |r_i − r_c|`, `R_out = R_cell`, `R_in = R_cell − h_cortex`:

```
U_conf(r_i) = ½ k_conf (r_i − R_out)²   if r_i > R_out      (outer soft wall)
            = 0                          if R_in ≤ r_i ≤ R_out (free band — NO force)
            = ½ k_conf (r_i − R_in)²     if r_i < R_in       (inner soft wall)
```

**Inside the band there is no force at all** — beads fluctuate, slide, buckle, condense freely; the wall only resists *leaving* the cortex.

**γ is STILL measured from explicit network stress, never from the manifold.** `cortex/cortical_tension.py:measure_cortical_tension` stays the only sanctioned readout: the soft method-of-planes sum of real per-bond `T = k(L−r0)` over cortical bond types (adhesion denylist applied; `|û·n̂|`; 2πR cut circumference) and the IK whole-shell virial `gamma_ik = Σ T·L·(1−(r̂·û)²)/(8πR²)` (validated <0.5% vs MOP on a synthetic shell). `U_conf` is a **radial force, not a bond**, so it cannot appear in `bonds.types` and `_gamma_soft` never sees it — the **same** mechanism by which ERM is "deliberately absent" from the tension sum. The H.8 membrane area-elastic term is a **separate, additive** `γ_membrane` on the Laplace path, never conflated with cortical tension.

**Two-part γ-invariance requirement (not one).** `U_conf` is *structurally* invisible to the bond-based estimator (proven above), but it is a force that *moves beads*, so it could in principle reshape the realized cortical-bond geometry the estimator reads. Therefore the falsifiable check is **stronger than "γ flat"**: sweep `k_conf` across its admissible bracket and assert that **both** the reported γ **and** the realized cortical-bond geometry — the {L/r0 distribution, min(r/r0), radial-position histogram} — are flat. A flat scalar γ alone is **not** advertised as clean invariance, because on a dead floor (s_grip ≈ 0, ~1×10⁻⁴) γ can be flat simply because the floor is numerically dead. If γ is flat only because of the dead floor, say so.

**Dynamic-binding candidate search becomes patch-local (broad-phase only).** Today `crosslinkers.py` and `myosin.py` rebuild a global `cKDTree` over all actin every tick. A static geodesic patch index (built once from bead (θ,φ), reusing the Fibonacci-sphere machinery in `cortical_tension`) gives O(1)-per-head candidate gather; **narrow-phase is byte-identical** — Euclidean reach, bridge-different-filament rule, bipolar sidedness gate, degree cap, Bell-Evans/Pereverzev kinetics. Patch radius `r_patch = max(max_bind_dist, reach)` (reach = √(A/n), the already-derived mesoscale crosslinker reach) is conservative by construction (no false negatives). This makes the σ_a investigation's long contraction run (drive s_grip→ℓ₀) affordable — but **it does not perform the experiment and does not lift γ**. (Per §3, at ×40 the reach is ~1880 nm, so the cortex candidate-reduction may be marginal; the win is at native/ECM scale.)

**Latent bug guard:** the `myosin.py` `fil_idx = bead_choice // beads_per_filament` map is wrong for the variable-length bimodal cortex and must be fixed (already partially addressed via the explicit `cortex_filament_idx` path) before the faithful bimodal build feeds the patch index, else patches and myosin tangents mis-map.

**ERM displacement is an Open PI decision, not a default.** `U_conf` is proposed as an **additive, default-off** band-soft confinement. Whether it **replaces**, **augments**, or **leaves untouched** the KU-3.18 ERM in production is an Open PI decision and a KU-3.18 contract change — `erm.py` itself states that softening k_ERM "requires PI sign-off vs brief literal KU-3.18." It is **not** an "obviously physiological one-line swap": replacing a sub-nm pin with a ~25 nm band changes the shell's resting boundary condition and the bead positions the γ estimator reads, so it is a candidate requiring its own validation (resting-shell thickness/σ observables preserved + before/after γ under ERM vs U_conf) before any default flips. `U_conf` stays **OFF for all validation references**; ERM stays the production default until validated and PI-ratified. A builder-level "at most one radial confinement attached" assertion guards against ERM+U_conf double-counting, but does **not** pre-decide which wins.

**`connected_mesh.py` is an interim ×40 compensation** (its √(A/n) bridge seeding is the crosslinker analogue of mesoscale myosin force-scaling). The long-term replacement is manifold-guided explicit filament placement at the physical 60 nm scale (Mikado-on-sphere). The manifold makes that natural but does not decide to take it (§22).

## 9. Future cell↔ECM traction mechanics application

The cell↔ECM load path is **already explicit and already wired** into `Cell.build`. The manifold sits **beside** the integrin–ligand arrow as a candidate-scope + coordinate-frame annotation, drawn deliberately **off** the force path:

```
explicit cortical actomyosin   (Stam-Hocky minifilaments bridge/motor.py; cortex/myosin.py — active-stress SOURCE)
   │  (gates A/B floor it HERE)
   ▼
explicit fa_actin_clutch bond  (cell/cell.py, force-free per-bond r0)
   ▼
explicit integrin–ligand catch bond  (bridge/integrin_bonds.py, Pereverzev/Bell-Evans)
   │  ── candidate set scoped by ──▶  CONTACT MANIFOLD (geometry/broad-phase ONLY; local n̂/t̂; NO force)
   ▼
explicit ECM resistance        (SubstrateLigandPin rigid / ecm/substrate.py k_sub spring / Mikado-on-sphere fibers)
```

**Inter-mesh contact manifold (broad-phase / candidate generation).** A per-tick, read-only `bridge/contact_manifold.py:ContactManifold` consuming global positions and emitting **patch-pairs** `{cell_patch_id, ecm_patch_id, n̂, t̂₁, t̂₂, contact_point, gap, area_est, cell_bead_tags, ecm_bead_tags}`. It **creates no adhesion, applies no force, moves no particle**. Cell patches = spatial bins of the south-cap bead cloud (reusing the per-bead PCA frame); ECM patches = the local `z≈0` ligand neighbourhood (flat) or `ecm/mikado.py` fiber-bead bins (fibrous; Mikado-on-sphere slots in here). Pairing = BVH over patch centroids within a conservative broad-phase reach `R_bp = capture_radius_R_FA + r_patch`.

**Narrow-phase = the existing machinery, candidates restricted, local frame supplied.** For each patch-pair, candidates are filtered by physical criteria in order: (1) reach `d ≤ capture_radius_R_FA` (the existing test, now over O(patch) ligands), (2) orientation `n̂_cell · (−n̂_ecm) > cos θ_max` (an integrin can only engage a ligand it faces — a candidate restriction, never an acceptance modifier; v0 `θ_max = 90°`, maximally permissive), (3) availability (existing `_engaged` mask, fan-in cap), (4) **unchanged** stochastic `p_on`/`p_break` with Pereverzev/Bell-Evans `k_off`. **This is the only place an adhesion is decided, and it is identical to today.**

On acceptance, **explicit bonds** are created exactly as now: integrin–ligand Pereverzev catch bond (born force-free via the per-bond r0-bin family — the manifold's pre-screening to the in-contact patch means the realised separation clusters near the nominal ~50 nm gap, so the 1500 pN spurious far-binder is geometrically excluded *before* binning, mutually reinforcing the r0-fix); clutch-actin load-path bond; talin/vinculin WLC + directional-catch when S3/S4 land (the manifold normal feeds the directional kinetics — confirm acceptable vs a fully explicit local actin orientation, §22); ECM anchor (rigid pin or k_sub spring).

**Traction is measured, not prescribed.** Per-bond forces are summed per patch-pair and **decomposed in the manifold frame** (F_normal = peel, F_tangent = shear, plus radial/circumferential), with traction stress `t = ΣF/area_est` (the only place geometry touches a reported number — a denominator for a readout, never a force) and strain energy. This extends the read-only `bridge/fa_growth.py:FAGrowthMonitor` and `bridge/clutch_spatial.py:beta1_distribution_metrics`. **No `t` is ever fed back into HOOMD as a body force.** A lumped field-level traction map (Chan-Odde / Bangasser-Odde traction-vs-stiffness) is permitted **only** as a validation/oracle overlay — same status as Bell-Evans/Hill, never a runtime input (§22). The same `F = k·(|Δr|−r0)` per-bond-r0 caveat applies: `fa_growth.py:171` currently assumes r0=0 and the correction belongs to the FA overload-fix workstream.

**Honesty.** The active stress driving this chain originates at the *top* arrow (cortical actomyosin), which is exactly where gates A/B floor it. The manifold neither adds nor removes a Newton. On the **current** build the ECM is orphaned and γ is floored, so a traction map produced today is a near-zero / FA-pin-dominated readout, not evidence of working active mechanics — meaningful only after gates A/B and the cell↔ECM coupling land.

## 10. Data model

**v0 is MESH-FREE.** Per the pressure-test blocker, v0 stores **no triangulation, no Delaunay/convex-hull build, no per-face `tri_area`, no edge adjacency `geo_adj`** — grep confirms no such primitive exists in-tree (the only area element is the scalar sphere split `A_i = S/N`), so a stored face-area mesh would be net-new geometry whose data structure is the **rejected** triangulated surface (minus, for now, only the edge spring). A cached `tri_area`/edge structure is one line away from `K_A·(A−A0)` area elasticity (the rejected area modulus) and makes a quadrature γ depend on face areas — both forbidden.

v0 stores only:

```
CellSurfaceManifold (cached metadata; rebuilt on the binding-updater cadence)
  vertex_tag      (V,)  i64       # HOOMD tag of each shell bead (the back-pointer)
  vert_normal     (V,3) f64       # per-bead outward normal (discrete_mean_curvature)
  vert_tangent    (V,3,2) f64     # per-bead (t1,t2) tangent basis (PCA frame)
  vert_curv       (V,)  f64       # per-bead mean curvature H [1/m]
  patch_id        (V,)  i32       # Fibonacci/HEALPix label over (θ,φ)  — NO triangle
  patch_lookup    dict            # patch_id -> [bead tags]  (broad-phase)
  patch_adj       dict            # patch_id -> [adjacent patch_ids]  (geodesic ring, from (θ,φ))
  shell_thickness float           # h_cortex = 200 nm (KU-3.17 band)
  centroid, R_mean float          # soft-confinement reference (estimate_area)

PerBeadManifoldRef
  bead_tag        (N,) i64
  nearest_patch   (N,) i32        # broad-phase home patch
  normal_offset   (N,) f64        # signed radial distance to shell (for U_conf only)
  # GLOBAL r stays in HOOMD; normal_offset is DERIVED, never integrated.
```

**Global Cartesian positions are the single source of truth.** The integrator (`integrator/baoab.py`, frozen) advances `r ∈ ℝ³` in global Cartesian with per-particle Stokes drag; the manifold metadata is **recomputed from** those positions, never the reverse. Intrinsic / barycentric surface dynamics is **deferred as an Open PI decision** (§22) — it would require a metric-aware constrained Langevin step (a new integrator, forbidden) and would re-over-constrain buckling. If a true triangulation is ever wanted, it is built **render-side only**, tagged "geometry, not physics," and forbidden from feeding any reported γ/traction/area.

## 11. Force model

**The manifold introduces exactly ONE new force: optional, soft, default-off, purely normal.**

```
U_conf(bead) = ½ k_conf (d_normal − d0)²        d_normal = signed distance past the nearest band edge
∇U_conf = k_conf (d_normal − d0) · n̂            (force purely along n̂; in-plane component ≡ 0)
```

**Proof it adds NO in-plane elastic energy (G1, G2, G3):** the force is **purely along n̂**; the tangential (t̂₁, t̂₂) components are identically zero by construction, so it cannot act like a filament or an edge spring (those are in-plane), and it cannot contribute to γ (the method-of-planes sums `T = k(L−r0)` over **HOOMD bonds**; `U_conf` is a `md.force.Custom` external field, not a bond, so it never enters `bonds.types`, exactly as ERM is "deliberately absent"; nor does a normal force project onto the IK hoop term `1−(r̂·û)²`). γ is **mechanically blind to `k_conf`** — *directly*. The configuration-invariance requirement (§8) covers the *indirect* path.

`k_conf` is constrained to be **softer than the cortex bond stiffness** (k_conf ≪ bond_k) so it can never become a load-bearing in-plane element by the back door and never tightens CFL (it inherits the `dt ≤ α·γ_b/k_eff` gate from `membrane_surface.py`/`erm.py`). An over-stiff wall (k_conf approaching bond_k) would near-pin the radius and thereby *indirectly* suppress the inward buckling excursion → so the admissible `k_conf` band must satisfy **both** σ_radial ∈ [12, 38] nm **and** VG-5 (the r/r0<0.98 buckling-tail count preserved).

**Everything else is unchanged explicit HOOMD physics** (cortex backbone/angles, crosslinker bonds, myosin minifilaments, integrin↔ligand catch bonds, clutch bonds, WCA). The optional substrate normal repulsion `U_sub` (one-sided, soft) is **OFF** until its `k_sub` is supplied by PI (§22).

## 12. Contact model

**Purpose:** form cell-patch↔ECM-patch candidate pairs and report a contact-area estimate, **without contributing any adhesion force** (G4 — all adhesion stays the explicit integrin/talin/vinculin clutch).

- **Pairing:** a cell-patch↔ECM-patch pair is *candidate* when their patches are apposed within a contact band `d_contact = integrin_reach + buffer` (integrin_reach = the existing physical catch-bond capture length; buffer ≈ ½·local patch spacing). Patch-pair candidacy is a cheap patch-centroid KD-tree query (O(10²–10³) patches).
- **A candidate pair opens the broad-phase window** in which the *explicit* integrin↔ligand catch bonds are allowed to form; the **decision to bind is the explicit Pereverzev/Bell-Evans acceptance** at the physical integrin reach. The manifold changes *where the clutch updater searches*, never *whether a bond forms or what force it carries*.
- **Contact-area estimate:** `A_contact = Σ` patch-area over apposed cell patches whose normal faces the substrate within `d_contact`. It is a **reported observable** (the platform analog of PI spreading data), normalized by a **solid-angle share of the global shell area** (never a manufactured face area). It never enters any acceptance and is never fit.

## 13. Neighbor-search model

**The codebase already does a physical-radius broad-phase** (`scipy.spatial.cKDTree.query_ball_point` over reach in `crosslinkers.py`/`myosin.py`/integrin loops; `md.nlist.Tree` for the WCA/LJ pair force). The manifold adds a **patch / geodesic-ring broad-phase** that *narrows* the candidate set the existing radius queries operate on (partner search for a bead on patch p = beads in p ∪ adjacent patches, then the **existing** reach test). It **does not touch HOOMD's nlist** — HOOMD's `md.nlist.Tree` continues to serve the pair-potential narrow-phase; the manifold replaces only the Python-side candidate pre-filter, optionally.

**Cadence/cost:** rebuilt on the binding-updater cadence (every `batch_steps`, **not** every BAOAB step); patch-label assignment is O(N) over (θ,φ). Strictly cheaper than, and amortized over, the binding interval it piggy-backs on; CPU-side alongside the existing CPU binding updaters (not on the GPU-resident hot path).

**PROOF of broad-phase neutrality (G5, the core safety claim):** every binding event is decided by broad-phase (candidate set within physical radius r) then narrow-phase (the physical acceptance). The manifold pre-filter produces `C_mesh ⊇ C_true` **iff** patch diameter + 1 ring ≥ r — a **constructional condition asserted at build time** (`min_patch_geodesic_radius ≥ max(all physical binding radii)`), not a tunable. Given containment, the narrow-phase applies the **same** physical test and acceptance to `C_mesh`; any extra candidate is rejected by exactly the physics that would reject it today → the **accepted set is bit-for-bit identical**, so γ/traction/contact-area/spreading are mesh-resolution-invariant. **Fail-safe:** if the containment assertion cannot hold, the prefilter **disables itself and falls back to the global `cKDTree`** (correct, slower) rather than risk a missed candidate. The reach r it must cover is, for the cortex, the ×40-compensated √(A/n) (G6) — it must read `r_phys`/`reach` from the resolved params and **never define its own reach**.

## 14. Measurement model

**Bright line for measurement:** the manifold enters a measurement **only as a coordinate system and a binning key**. If a proposed measurement cannot be written as "a sum/projection of explicit bond or particle quantities, optionally indexed by a patch label," it is out of scope.

- **M.1 — Global γ (REUSE).** Authoritative γ stays `cortical_tension.measure_cortical_tension` verbatim — soft method-of-planes, IK whole-shell virial `gamma_ik`, rigid — with the B3 rule (turgor never folded into `gamma_structural = gamma_soft + gamma_rigid`) preserved. The `|û·n̂|` absolute value and the `ADHESION_BOND_TYPES` denylist (the √N-vs-N floor audit and the ~100× FA-clutch inflation, respectively) are load-bearing and inherited unchanged.
- **M.2 — Local cortical stress (per-patch, DIAGNOSTIC ONLY).** Bin every cortical bond by the **patch label** of its midpoint; report a per-patch surface-stress tensor from the explicit bond tensions projected onto the patch (n̂, t̂₁, t̂₂) frame. This is a **non-gate diagnostic / colormap**, never reported as γ. **It must sum exactly to `gamma_ik` by construction** — same bonds, **same global 8πR² normalization** (or the per-patch solid-angle share of it), **never a manufactured face area**. There is therefore **no quadrature gap**: if the patch sum disagrees with `gamma_ik`, the **binning is wrong**, not "discretization."
- **M.3 — Contact-patch traction.** Sum the explicit `integrin_ligand` + `fa_actin_clutch` bond forces (the denylisted set), decomposed in the patch frame (normal/tangential/radial/circumferential), normalized by area. Read-only. **Caveat (cross-referenced):** uses `F = k·(|Δr|−r0)`; `fa_growth.py:171` currently assumes r0=0, and the per-bond r0 subtraction is a PI-ratified-with-the-overload-fix change in a separate workstream — not yet in tree, not folded into a manifold PR.
- **M.4 — ECM strain energy.** Explicit sum over ECM bonds `U_ECM = Σ ½ k(L−r0)²` + WLC/bending; the manifold provides only the contact surface and a render overlay.
- **M.5 — Consistency identities (the manifold's correctness gates), both dimensionless ratios:**
  - **(I) Contact-conservation:** `Σ_c (A_c · traction_c) == Σ_b F_b` over all explicit clutch+integrin bonds — an **exact partition-of-unity** to round-off (machine epsilon), independent of patch count. A nonzero residual = a bond dropped/double-counted = a **binning bug**, not physics.
  - **(II) Stress-consistency:** `(1/A) Σ_p A_p·γ_local,p == gamma_ik` **by construction** (same bonds, same normalization), to round-off. The patch integral is the same bonds re-indexed, so there is no tolerance band; `gamma_ik` (zero crossing-variance, all bonds) is the anchor. Any genuine `N_patch` dependence in a reported number is a **halt-and-surface**, not a tolerated gap.

  Per G8/G10, if a real run cannot satisfy these, that is surfaced to PI — the tolerance is **not** loosened to make a manifold build pass.

- **Curved-shell γ (v1, Open PI decision).** On a non-spherical (spread/FA-adhered) shell, **both** the MOP cut circumference (`2πR_cell`) **and** the IK denominator (`8πR²`, origin-centered r̂) are only sphere-valid. The curved-shell generalization must fix **MOP and IK together** (a consistent surface-stress-tensor integration), or neither — fixing only the MOP circumference would break the M.5-II identity for a **geometric** reason and must not trigger a tolerance loosening. v0 keeps the sphere (R_cell global, bit-for-bit identical); the joint generalization is one coupled validation-contract change to PI.

## 15. Visualization / export model

**Principle:** viz is a **debugging instrument**, not decoration; each layer makes a *named* failure mode visible (PI 2026-05-21). Format is split by data character so the bright line is preserved in the files themselves:

- **Explicit particles/bonds → GSD** (reuse the existing `common/gsd_traj.py` writer — the Blender unlock; `scripts/npz_to_gsd.py`, `scripts/render_fresnel.py`). One scalar field per mechanical quantity under `frame.log['particles/<name>']`: cortex filament overlay colored by per-bond tension + connectivity component; lamellipodium leading edge; per-bead local stress/γ (diverging colormap centered at 0). **Do not add a second particle exporter.**
- **Manifold + per-patch fields → glTF/PLY (Blender hero), `.vtp` (ParaView quantitative), JSON sidecar (identity audit)** — explicitly tagged "geometry, not physics": curvature H (from `discrete_mean_curvature`, reuse), per-patch γ_local / principal axes, contact-patch + FA-density heatmap (`f_edge`/`angular_cv`), traction-vector glyphs, ECM `U_ECM` density. The JSON sidecar is the canonical input to the resolution-invariance gate and the M.5 identity audit.

**Viz as bug detector:** wrong normals → render n̂ glyphs (must be uniformly outward; curvature H must be uniform +1/R on a sphere); contact-patch gaps/overlaps → directly visible (these break M.5-I); traction sign errors → glyphs must point centripetally (a sign flip is a color flip); FA mislocalization → density heatmap must concentrate at the periphery; cortex-shell leakage → beads colored by radial distance (a thin annulus); **mesh-resolution artifacts → the same frame at two `N_patch` side-by-side** (if the colormap shifts with patch count, the manifold leaked into physics — the visual twin of VG-1).

**Integrity (non-negotiable):** no axis truncation; colorbars span the full range; signed → diverging, non-negative → sequential; SI scale bar baked into GSD `log_constants`; the literature band [0.35, 0.65] mN/m **and** the MCF7 datum 0.27 mN/m overlaid on every γ figure (overlay, never a fit); ensemble plots show per-realisation thin lines + mean. Auto-viz at run end (best-effort subprocess). **γ on a figure is always the `cortical_tension.py` number, never back-computed from the rendered mesh** — if the mesh render and the explicit γ disagree, the explicit γ is correct and the mesh is the suspect.

## 16. Validation gates

| Gate | Metric | Pass band | Failure signature |
|---|---|---|---|
| **VG-1 (MASTER) — mesh-resolution independence** | γ_soft, γ_ik, γ_rigid, ΣF_clutch, footprint, coordination z, min(r/r0) at N_patch ∈ {320,1280,5120} + no-manifold baseline, ≥5 seeds | each within **±2%** of baseline (stochastic-spread CV, not a physics tolerance) **and** no monotone trend in N_patch | any observable drifts monotonically with N_patch → manifold doing physics → **halt; revert** |
| **VG-2 — patch-anisotropy invariance** | same vector under random SO(3) re-mesh (3 rotations × 5 seeds) + de-rotated per-patch traction histogram | scalars within VG-1 band; histogram KS p>0.05 | traction/γ tracks the icosphere 5-fold defects → discretization imprinting |
| **VG-3 — contact-force conservation** | Σ_patch traction vs Σ_bond F over ADHESION_BOND_TYPES | vector-equal to **1e-9 N** (round-off; bookkeeping) | mismatch > round-off → patch layer generating/absorbing force (G4 breach) |
| **VG-4 — global-vs-local stress consistency** | (1/A)Σ_p A_p·γ_local,p vs gamma_ik (same bonds, same 8πR² normalization) | **exact to round-off** (partition-of-unity; no ±band) | any N_patch dependence → binning bug; converge-to-gamma_ik or **halt** |
| **VG-5 — buckling preservation (GATE-B FIREWALL)** | per-bond r/r0 distribution with vs without manifold under active load | min(r/r0) and the r/r0<0.98 count unchanged within VG-1 band | sub-unity tail thins → confinement acting as a structural shell → reduce k_conf; **hard projection forbidden** |
| **VG-6 — broad-phase neutrality** | realized binding/unbinding event stream + final z, γ_soft: full-distance vs patch-neighborhood enumeration | **event streams bit-identical** (patch set ⊇ physical reach) | a binding event in (a) missing in (b) → broad-phase dropped a reachable candidate (cardinal sin) → halt |
| **VG-7 — visualization integrity** | per-patch exported scalars round-trip to the explicit per-bond quantities | exact for sums; SI units; band overlaid; per-realisation lines | a smooth field where explicit data is noisy → interpolation hiding spread |

VG-4's invariance is also checked as a **convergence-slope test**: the (already near-zero) residual must shrink toward 0 as both n_planes and N_patch rise, converging to gamma_ik. A residual that is small but does **not** shrink is a FAIL (residual grid dependence), per G5.

## 17. Minimal v0 prototype (additive 9-step plan)

| Step | Action | New sibling vs metadata extension |
|---|---|---|
| **G1** | Spherical cell-surface manifold metadata: centroid (reuse `membrane_surface.estimate_area`) + per-bead patch label (Fibonacci/HEALPix over (θ,φ)); **no triangulation** | **NEW** `cell/manifold.py` (numpy only) |
| **G2** | Flat ECM/substrate manifold metadata: plane at basal z of `ecm/substrate.py`, tiled for patch indexing | **NEW** `FlatManifold` in `cell/manifold.py` (owns no physics) |
| **G3** | Bead→nearest-patch + per-patch (n̂, t̂₁, t̂₂); for the sphere reuse the radial unit vector / `discrete_mean_curvature` | thin: a `nearest_patch(positions)` pure function |
| **G4** | Soft shell confinement `U_conf` (`md.force.Custom`, sibling of `erm.py`, same CFL gate), flat-bottom band well, **default-OFF** | **NEW** `cell/manifold_confine.py` (sibling of `erm.py`) |
| **G5** | Patch-neighborhood broad-phase index threaded as optional `candidate_index` into `crosslinkers.make_xlink_updater`, `lamellipodium` Arp branching, integrin loop; `None` ⇒ current full-distance path (bit-identical) | thin: optional arg, `None`-guarded |
| **G6** | Narrow-phase stays physical: Bell-Evans/Hill/catch-bond/WLC/bridge-different-filament — **zero edits**; a test asserts the acceptance code paths are untouched by diff | none (firewall by omission) |
| **G7** | Instrument: candidate-count (on/off), wall-time, z, giant-fraction, γ_soft/γ_ik (call `measure_cortical_tension` unchanged), per-patch traction (VG-3), footprint (explicit hull), min(r/r0) | thin: diagnostics dict in the driver |
| **G8** | Resolution sweep G1–G7 at {320,1280,5120} × ≥5 seeds → VG-1/VG-2 tables | driver-only (`scripts/manifold_v0_sweep.py`) |
| **G9** | Export per-patch fields (traction, local γ, bound-fraction, n̂) to GSD + glTF/.vtp/JSON sidecar with the VG-7 round-trip baked in | extend the existing `*_vis.py` entry-point |

**Success criterion:** a real candidate-count and wall-time reduction **with every physical observable inside the VG-1 band and every gate green**. If any observable moves, revert at zero cost (everything `None`-guarded / default-OFF). v0 introduces **no force-law constant except `k_conf`** (default-OFF), guarded by VG-5.

## 18. Migration path

Stages are ordered by risk; **each is a no-op when its handle is `None`/OFF**, leaving prior runs bit-for-bit identical. **The integrator is never touched.**

- **M0 (prerequisite).** Promote `connected_mesh.build_connected_cortex` (full faithful build) to the production default in `scripts/mcf7_fullcell_stage1.py`, replacing the deferred "fast hybrid" (rejection-doc Rec 1a; independent of the manifold). The manifold is built **on top of** the faithful connected cortex.
- **M1 — geometry only (zero physics).** `cell/manifold.py` (G1–G3): patch labels + per-patch (n̂, t̂₁, t̂₂), written to the sidecar. Nothing enters a force. Revert = delete the import. *Unblocks Blender visualization immediately.*
- **M2 — broad-phase index, candidates-only (G5).** Optional `candidate_index` into the three candidate loops; default `None` ⇒ existing path. **VG-6 bit-identical required before ON.** Revert = pass `None`. **No constant introduced.**
- **M3 — patch-stress + traction aggregation (reporting + VG-3/VG-4).** Pure measurement; revert = stop reporting.
- **M4 — soft confinement (G4), default-OFF.** Only after VG-5 passes on a buckling-active run. Introduces `(k_conf, Δ_conf)`. **CONTRACT-TOUCHING:** this is the **one** stage that is not a pure `None`-guarded no-op — whether it **replaces / augments / leaves untouched** the KU-3.18 ERM is an Open PI decision; its revert is "re-attach KU-3.18 ERM + re-gate," not a one-line flag flip. OFF for all validation references; ON in production only at a PI-ratified physiological band.
- **M5 — substrate/FA manifold contact frames (G2 applied).** Per-patch substrate normal/tangent for cleaner *reporting* decomposition and basal broad-phase. Contact *forces* stay in the explicit clutch (VG-3 firewall). Revert = drop the frame.

**Byte-identical across ALL stages (so the γ-floor work proceeds in parallel, unblocked):** `cortex/cortical_tension.py` (all channels, cut-planes, denylist); the myosin grip_walk kinetics (Gate A); the M-SHAKE/buckling question (Gate B); Bell-Evans/Hill/catch-bond/WLC acceptance (VG-6). Whoever attacks the floor (myosin walking; relaxing M-SHAKE under PI integrator-freeze sign-off) works on the **same explicit filaments** with or without the manifold.

## 19. Consolidated Magic-Number Blocks

| Constant | Value | Source / derivation | Dimensional | Grid-invariance | Not gate-chosen |
|---|---|---|---|---|---|
| **`k_conf`** (soft shell-confinement normal stiffness; the only new force constant) | **[3×10⁻⁶, 3×10⁻⁵] N/m** (= [1×10⁻³, 1×10⁻²]·bond_k; default within ~3× of 6.8×10⁻⁶ N/m, σ_radial ≈ 25 nm) | A **fraction of the existing literature-anchored cortex bond stiffness** bond_k = μ/ℓ₀ = 1.5×10⁻⁹/500 nm = **3.0×10⁻³ N/m** (verified). The bracket brackets the band-width derivation σ_radial = ¼·(h_cortex/2) = 25 nm ⇒ kT/(25 nm)² = 6.848×10⁻⁶ N/m (verified). **Single derivation across the whole note**; the "1×10⁻³ N/m / round-up-a-decade" framing is dropped (it is ~0.33·bond_k, load-bearing, σ_radial 2.1 nm — too stiff). Thermal-hold lower bound √(2γ_b kT/Δt)/(d_shell/2) confirmed to fall inside this bracket (a sanity floor, not binding). | k_conf [N/m]·(Δr)² = J; force k_conf·Δr [N] | **N-invariant**: a ratio to bond_k (the ×40-mesoscale extensional stiffness), so it inherits bond_k's grid-invariance; per-bead **local** wall, independent of N/ℓ₀/dt. Does **not** inherit the membrane `effective_bead_stiffness` 1/N² scaling (that is a global-area per-bead share, not a local stiffness). CFL τ_conf = γ_b/k_conf ≫ dt_cfl, never tightens the step. | Provably γ-invisible (normal-only, not a bond → absent from the method-of-planes sum and the IK hoop term). **Falsifiable gate:** sweep k_conf across the bracket → **both** scalar γ **and** the realized {L/r0, min(r/r0), radial histogram} flat (VG-5). Admissible band must satisfy σ_radial ∈ [12, 38] nm **and** VG-5 buckling-tail preserved. Value outside the bracket = halt-and-surface-to-PI. |
| **`Δ_conf`** (confinement band edge / offset) | h_cortex = 200 nm band; far-field outer offset ~3×bead_radius | KU-3.17 cortex thickness (pre-existing `ResolvedH3.cortex_thickness`, used by `_project_to_shell_band`); far-field offset set so no bead reaches the outer wall in equilibrium (a runaway guard, not a tension source) | [m] | Fixed biological length, independent of N_patch / bead count / dt | Literature value; VG-1/VG-5 falsify any value tuned to move γ |
| **`h_cortex`** | 200 nm | KU-3.17 anchor (pre-existing) | [m] | Fixed biological length | Pre-existing; defines geometry, not a force feeding γ |
| **`r_patch`** (broad-phase patch radius, cortex) | max(max_bind_dist, reach), reach = √(A/n) | No new constant: sized off the binder's own physical/×40-compensated reach (`connected_mesh.mesoscale_reach`) | [m] | Grid-aware: scales with ×40 spacing exactly as reach does | Affects wall-time only; conservative (≥ reach ⇒ no false negatives) |
| **`R_bp`** (cell↔ECM patch-pair reach) | capture_radius_R_FA + r_patch | Purely geometric over-inclusion bound (no dropped narrow-phase pair) | [m] | Scales with patch radius; candidate **union** bounded below by physical capture_radius | Chosen for conservativeness; enlarging costs only CPU |
| **`d_contact`** (cell↔ECM apposition band) | integrin_reach + ½·patch spacing | integrin_reach = existing physical catch-bond capture length; buffer = geometric slack | [m] | Physics term grid-independent; buffer scales with patch spacing | Broad-phase window width; any d_contact ≥ true reach yields the same accepted set |
| **`θ_max`** (orientation acceptance cone) | cos θ_max = 0 (90°, facing hemisphere) | An integrin can engage only the half-space it faces; maximally permissive physical bound on the measured PCA normal (density-invariant) | dimensionless | Inherits the PCA-normal density-invariance; independent of tiling | Excludes only the physically impossible; any tighter cone is a PI param-change |
| **`Δ_map`** (substrate field sample spacing) | discretization knob (set by the patterned condition) | Interpolation spacing over the already-required per-condition k_sub/ρ_lig; introduces no physical value | [m] | Gated by 4× refinement invariance; degenerates to today's scalar for uniform substrate | Convergence-only; cannot be tuned to pass a gate |
| **broad-phase cell / BVH leaf size** | performance-only (matches `md.nlist.Tree` precedent) | Acceleration-structure parameter; carries no physics | [m] | Superset contract: candidate set ⊇ all true contacts within r_cut (all-pairs cross-check) | Cannot enter acceptance; too-coarse → caught as a correctness bug |
| **`N_planes`** | 12 (pre-existing) | Fibonacci isotropy sampling count (`cortical_tension`) | dimensionless | γ_MOP agrees with N_planes-independent γ_ik <0.5% | Predates this note; not tuned to a band |
| **`N_patches`** | set by convergence, not a priori | Geometric label count (Fibonacci/HEALPix) | dimensionless | **Decisive G5:** per-patch map sums exactly to γ_ik; doubling changes nothing | Binning index only; bond tensions unchanged |
| **`k_sub`** (substrate compliance, optional repulsion) | **REQUIRED-no-default / OPEN PI** | Per-condition substrate modulus (KU-1.V.1); **value still unknown** for the GATE-B lamellipodium open item → surfaced to PI, **not chosen** | [N/m] | Per-ligand physical stiffness; field is an interpolation over per-condition values | Stays REQUIRED so it cannot be invented to pass a gate; CFL gate still fires |
| **`z₀`** (substrate plane height) | 0 m (inherited) | Construction z of the substrate ligands (`bridge/fa.py`) | [m] | Single analytic plane, tessellation-independent | No gate reads it |
| **`r_branch_eff`** (Arp2/3 reach) — UNCHANGED | 100 nm (physical) | PI-ratified 2026-05-29 (Funk 2021 + Bieling 2023); the manifold does **not** change it | [m] | Physical reach; broad-phase must reproduce the identical d<r_branch_eff set | Pre-existing; broad-phase forbidden from altering eligibility (G6) |

## 20. Failure modes (consolidated)

1. **Manifold mistaken for the γ-floor fix.** The single biggest risk: a clean traction map or a flat k_conf sweep read as "γ solved." γ is generation-limited upstream (gates A/B). Every deliverable must carry the "does NOT lift γ" banner.
2. **Dead-floor masquerading as invariance.** Because γ is floored ~1×10⁻⁴, a k_conf sweep can show γ flat simply because the floor is numerically dead. Distinguish: report the realized bond geometry, not just scalar γ.
3. **Indirect (configuration) γ leak.** `U_conf` is structurally invisible to the bond-based estimator (not a bond), but it moves beads; an over-stiff wall can reshape realized bond extensions or thin the buckling tail and thus indirectly lower the realized γ. Guard: k_conf ≪ bond_k bracket + VG-5 buckling-tail + the two-part invariance check (scalar γ **and** realized {L/r0, min(r/r0), radial histogram}).
4. **Stored-triangulation creep (re-importing the rejected mesh).** A cached `tri_area`/edge structure is the rejected surface minus the spring and one line from area elasticity. Guard: v0 stores **no** triangulation; any render-side mesh is tagged "geometry, not physics" and forbidden from feeding a reported number.
5. **Hard-projection regression (Gate-B killer).** "Optimizing" the soft band into per-step hard projection re-freezes the radial DOF. Forbidden; the band interior must be strictly force-free (static test: zero force in (R_in, R_out); σ ≈ √(kT/k_conf) realized).
6. **Broad-phase under-inclusion (false negatives).** patch_radius < reach (or too-small ring) drops a reachable candidate, silently changing physics. Guard: build-time containment assertion + fail-safe fallback to the global cKDTree + VG-6 event-stream equivalence.
7. **Mesh-resolution dependence of any reported number (umbrella, G5).** Caught by VG-1 (master) and the exact M.5 partition-of-unity; a single monotone trend → revert.
8. **Mesh anisotropy (icosphere 5-fold defects imprint on per-patch fields).** Detected by VG-2; mitigate by using the patch index for binning only and falling back to the orientation-free `md.nlist.Tree`.
9. **ERM double-counting.** Attaching both stiff ERM and U_conf over-confines and the 0.21 nm ERM dominates. Guard: builder-level "at most one radial confinement attached"; ERM-displacement is an Open PI decision, not a silent swap.
10. **Coarse-graining-compensated reach erodes the cortex speedup.** With ×40 reach ~1880 nm on ~1000 beads, "patch ∪ 1-ring" is near-global → marginal cortex candidate-reduction. Guard: report realized candidate-count; the cortex value at ×40 is the frame + contact model + viz, not speed; the speedup payoff is native/ECM scale.
11. **Contact-area / traction observable mistaken for a force input or a fit target.** A_contact and traction are measurement-only (spreading overlay), never entering acceptance, never fit to PI data.
12. **Measurement-protocol drift (`fa_growth.py:171` r0=0).** The per-bond r0 subtraction does not yet exist; it belongs to the FA overload-fix workstream and must not be smuggled into a geometry/broad-phase PR.
13. **Curved-shell γ one-sided fix.** Fixing the MOP circumference without the IK denominator breaks M.5-II for a geometric reason and tempts a tolerance loosening. Guard: fix MOP and IK together or neither; v0 keeps the sphere.
14. **Hidden integrator coupling via intrinsic coords.** Promoting barycentric metadata to integrated DOFs silently requires a metric-aware constrained Langevin step (forbidden) + re-over-constrains buckling. Guard: v0 keeps global Cartesian as the only integrated coordinate; intrinsic dynamics deferred to PI.
15. **Readout meaningful only after generators are alive.** On the current build the upstream generators are floored (γ ~1×10⁻⁴, branch density ~0/µm², ECM orphaned); any traction map today is near-zero / FA-pin-dominated, not evidence of working active mechanics.

## 21. Answers to the 12 design questions

**1. Does this satisfy the full-fidelity rule?** Yes. Every actin filament, myosin head/minifilament, crosslinker, integrin/talin/vinculin/clutch, and ECM fiber/crosslink stays an explicit HOOMD particle/bond; the only sanctioned coarse-graining (×40) is unchanged. The manifold abstracts **zero** filaments into edges — it is metadata recomputed from explicit positions, carries no force-bearing edge, and the one new force (`U_conf`) is a soft normal-only external field that is provably not a load-bearing in-plane element. It is the navmesh/broad-phase around the physics actors, never a replacement for them.

**2. Minimal manifold layer without rewriting HOOMD.** A plain NumPy/`@dataclass` metadata layer (`cell/manifold.py`) computed from `cpu_local_snapshot` positions on the binding-updater cadence: per-bead PCA normal/tangent (reuse `discrete_mean_curvature`), a Fibonacci/HEALPix patch label, and a patch→bead lookup. The one optional force routes through the existing `md.force.Custom` external-field template (sibling of `erm.py`). The integrator (`integrator/baoab*.py`) is **never touched**.

**3. Global-Cartesian + metadata vs barycentric.** **Global Cartesian positions are the single source of truth**; the manifold stores derived per-bead metadata only. Barycentric/intrinsic dynamics is **deferred (Open PI)** because advancing it correctly needs a metric-aware constrained Langevin step — a new integrator (forbidden) — and would constrain beads *to* the surface, re-introducing the Gate-B-killing over-constraint. Global-Cartesian + metadata touches the integrator zero times and keeps the physics buckling-capable.

**4. Shell confinement without suppressing buckling.** A **soft, one-sided, flat-bottomed** radial well `U_conf` that is **force-free inside the 200 nm band** (beads fluctuate/slide/buckle/condense freely) and softly resists only *leaving* the cortex. **Hard projection is forbidden** (it re-freezes the radial DOF buckling rides on). VG-5 (the r/r0<0.98 buckling-tail count preserved) is the firewall, and the admissible k_conf band must satisfy it.

**5. Inter-mesh contact pairs.** A read-only per-tick `ContactManifold` emits cell-patch↔ECM-patch pairs (with n̂/t̂ and a candidate area) via a BVH over patch centroids within `R_bp = capture_radius_R_FA + r_patch`. It **creates no adhesion and applies no force**; it only opens the broad-phase window in which the explicit integrin↔ligand catch bonds may form.

**6. Exposing planes/normals/tangents to FA/clutch/lamellipodium.** A shared frame provider returns (n̂, t̂₁, t̂₂) recomputed from live positions for any tag (the lamellipodium capping/Arp/elongation already consume local frames; the cortex/FA decompose explicit forces into normal/tangential/radial-circumferential). The frame is used for candidate generation and force *decomposition for measurement*, never to apply a force.

**7. Integrin candidate search via patch neighborhoods preserving physical rules.** Patch-restricted candidate generation (`C_mesh ⊇ C_true` by a build-time containment assertion) feeds the **unchanged** narrow-phase: reach test, facing-orientation filter (90°, excludes only the physically impossible), availability/fan-in, and the stochastic Pereverzev/Bell-Evans acceptance. The accepted set is bit-identical to the global-search path (VG-6).

**8. ECM fibers + ECM surface manifold coexistence.** The ECM stays the explicit Mikado fiber network (`ecm/mikado.py` + `ecm/cross_links.py`); the manifold overlays only (i) a broad-phase BVH over fiber beads and (ii) a per-contact local tangent plane derived from explicit bead positions each query. No continuum mesh replaces Mikado mechanics; Mikado-on-sphere fibers slot in as the ECM bead cloud unchanged.

**9. γ + traction measurement.** γ stays `measure_cortical_tension` verbatim (soft MOP, IK virial, rigid; B3 separation). Per-patch γ/stress and traction are **diagnostics** binned by patch label, summing **exactly** to gamma_ik / Σ_b F_b by construction (no quadrature band). Traction = Σ explicit clutch+integrin bond forces decomposed in the patch frame (with the `F = k(|Δr|−r0)` per-bond-r0 caveat flagged as a separate-workstream change).

**10. Blender/ParaView/Three.js export.** Reuse the existing GSD path (`common/gsd_traj.py` → fresnel/Blender) for explicit particles/bonds; export the manifold + per-patch fields separately to glTF/PLY (Blender), `.vtp` (ParaView), and a JSON sidecar (identity audit), each tagged "geometry, not physics." Three.js consumes the same glTF + JSON overlay. γ on a figure is always the explicit `cortical_tension.py` number.

**11. Failure modes introduced.** See §20 — chiefly: being mistaken for a γ fix; a dead-floor flat-sweep mis-read as clean invariance; indirect γ leak via bead repositioning; stored-triangulation creep; hard-projection regression; broad-phase under-inclusion; mesh-resolution dependence; ERM double-counting; the ×40 reach eroding the cortex speedup; and the readout being meaningful only after the generators are alive.

**12. Validation gates before production.** VG-1 (master resolution-independence) → VG-2 (anisotropy) → VG-3 (contact-force conservation, round-off) → VG-4 (global-vs-local stress, exact by construction) → VG-5 (buckling firewall) → VG-6 (broad-phase neutrality, bit-identical event stream) → VG-7 (viz integrity), plus the k_conf two-part invariance sweep (§19). Any monotone observable-vs-resolution trend, or any failing identity, is a halt-and-surface-to-PI, never a tolerance loosening.

## 22. Open PI decisions

1. **ERM displacement (KU-3.18 contract).** Whether `U_conf` **replaces**, **augments**, or **leaves untouched** the KU-3.18 ERM in production. *Recommended default:* keep ERM as the production default; `U_conf` default-OFF for all validation references; flip only after VG-5 passes, resting-shell observables are validated under U_conf, and PI ratifies the KU-3.18 change.
2. **`k_conf` Magic-Number Block sign-off.** Ratify the single bracket [3×10⁻⁶, 3×10⁻⁵] N/m and the in-bracket production value (~6.8×10⁻⁶ N/m, σ_radial ≈ 25 nm). *Recommended default:* adopt the ratio-anchored bracket; any run needing k_conf outside it is a contract-change trigger.
3. **M0 — promote `connected_mesh.build_connected_cortex` to the production default** in `mcf7_fullcell_stage1.py` (rejection-doc Rec 1a; the manifold's substrate). *Recommended:* land first.
4. **Cell↔ECM coupling.** Whether to couple the cell to the H.1 Mikado (integrin→Mikado-bead explicit clutch) vs keeping H.1 a standalone bulk-rheology deliverable with the manifold as a documented seam. *Recommended:* keep H.1 standalone now; the manifold enables, does not decide.
5. **Curved-shell γ generalization (validation-contract change).** Generalizing the MOP cut circumference **and** the IK 8πR²/origin-centered-r̂ denominator **together** for non-spherical shells. *Recommended:* keep v0 sphere-only (bit-identical); surface the joint generalization as one coupled change; do not loosen M.5-II for a one-sided fix.
6. **`fa_growth.py:171` per-bond r0 subtraction.** A measurement-protocol change coupled to the FA overload-fix r0-bin family. *Recommended:* ratify with the overload-fix; **do not** fold into a manifold PR.
7. **Substrate stiffness `k_sub` value** (GATE-B lamellipodium open item; **still unknown**). *Recommended:* surface to PI; substrate-repulsion term stays OFF and substrate contact stays the explicit clutch until set; the value must come from a literature modulus anchor, never a gate.
8. **k_sub→E_sub contact-model oracle** (flat-punch vs Hertzian) before any manifold-placed compliance field is read as a Pa modulus (`effective_E_sub` diagnostic).
9. **Talin/vinculin directional-catch frame source (S3/S4).** Confirm the manifold normal is acceptable for the vinculin–actin *directional* catch kinetics vs requiring a fully explicit local actin orientation.
10. **Field-level traction map (Chan-Odde / Bangasser-Odde) as overlay-only.** Confirm it stays a validation oracle, never a runtime input (same status as Bell-Evans/Hill).
11. **Patch partition standard** (Fibonacci/HEALPix vs icosphere labels) and canonical N_patches for the resolution-invariance reference grid + auto-viz.
12. **Intrinsic/barycentric surface dynamics (deferred).** Whether a future membrane-species surface diffusion warrants a *separate* metric-aware constrained integrator alongside BAOAB — never a rewrite of it.
13. **Adopting Mikado-on-sphere** (geodesic-intersection crosslinks at ~60 nm, buckling-capable) as the long-term replacement for the ×40 `connected_mesh` compensation. *Recommended:* sequence after gates A+B per the rejection doc; the manifold makes it natural but does not decide it.
14. **Relaxing M-SHAKE (the actual Gate-B γ lever)** is an `integrator/`-freeze change requiring PI sign-off. Independent of the manifold (which neither helps nor blocks it); the migration path keeps it unblocked.
15. **ECM crosslink kinetics** (permanent vs Bell-Evans/catch-bond off-rate for matrix yielding/remodeling) — a ratified ECM-scope call; the manifold adds no remodeling.
16. **Adding contact-area / spread-area and the manifold-prefilter equivalence test as reported quantities / gate-suite members** is a validation-contract addition (G10) — proposed, not enacted.

## 23. Pressure-test summary (honesty appendix)

| Lens | Verdict | Key residual risk |
|---|---|---|
| Coarse-graining-fidelity | pass-with-fixes (applied) | The v0 data model must stay mesh-free (no stored triangulation/tri_area/edges); any future render-side mesh must never feed a reported number. |
| Magic-number / grid-invariance | pass-with-fixes (applied) | `k_conf` unified to one bracket [3×10⁻⁶, 3×10⁻⁵] N/m (= [1e-3,1e-2]·bond_k); the σ_radial(k_ERM=0.1)=0.21 nm correction must propagate (erm.py docstring 0.65 nm is a separate pre-existing bug). |
| γ-floor & measurement | pass-with-fixes (applied) | Two-part γ-invariance (structural **and** configuration) must be tested, not just scalar γ; ERM→U_conf is a candidate requiring validation; curved-shell γ must fix MOP+IK together. |
| Physiological-baseline & contracts | pass-with-fixes (applied) | The physiological-baseline rule must not be used to *mandate* turning U_conf ON; the established cortex tether of record is KU-3.18 ERM; ERM displacement and the per-bond-r0 fix stay Open PI / separate-workstream. |
| Subordination to the cortex-as-mesh rejection | pass (no contradiction) | Keep the "triangulation = index structure, zero edge springs" invariant prominent so a future reader does not re-read the metadata layer as the rejected physics mesh; keep the "permitting buckling ≠ causing it" qualifier mandatory. |
