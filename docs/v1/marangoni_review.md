# Marangoni mechanism + long-time behavior — comprehensive review

PI directive 2026-04-29 (parallel with Production Lam4 finding). Review-
only; no code changes. Goal: identify mechanisms missing from the
current Layer 4 (Marangoni) implementation that could account for the
peak-and-decay long-time behavior observed in Production Lam4.

---

## 0. Production Lam4 trajectory — corrected reading

The previously-committed `docs/v1/production_lam4_finding.md` (commit
b506b57) reports "Production Lam4 80 hr A/A₀_topdown(end) = 2.643" and
interprets this as a plateau near 2.65. This **conflates two distinct
metrics in `metrics.csv`**:

| Metric (column) | t*=0 | peak | peak t* (sim hr) | end (t*=4800, 80 hr) |
|---|---|---|---|---|
| `A_over_A0_topdown` (top-down envelope) | 1.000 | **1.570** | 465 (~7.75 hr) | **1.409** |
| `contact_xy_hull` / A₀ (basal patch only) | 1.000 | **2.643** | 165 (~2.75 hr) | 2.142 |

The 2.643 number comes from the substrate-contact gate
(`contact_area_xy_hull`, runner.py:1056-1071), **not** from the top-
down area metric used throughout the Phase 4-v2 narrative. Both metrics
nevertheless show the **same qualitative behavior**: peak followed by
slow retraction.

The corrected long-time picture:

- **Initial spread (t* = 0–165, 0–2.75 hr)**: contact patch grows
  rapidly to 2.64×; top-down envelope still climbing.
- **Top-down peak (t* = 165–465, 2.75–7.75 hr)**: basal patch already
  retracting (2.64 → ~2.4); top-down envelope reaches its maximum 1.57
  as bulk slumps onto the larger contact area.
- **Long-time decay (t* = 465–4800, 7.75–80 hr)**: both metrics decay
  monotonically. Contact 2.64 → 2.14 (−19%); top-down 1.57 → 1.41
  (−10%). PI experimental Lam4 over the same window grows monotonically
  toward 8–33. Sim/PI ratio degrades from O(1) at hr 1 to ~0.05–0.18
  at hr 80.

**Bucket reclassification** (per `docs/v1/production_lam4_outcomes.md`
boundaries on the top-down endpoint):
- P1 [4.0, 33.1]: NO
- P2 [2.5, 4.0]: NO (1.409 < 2.5)
- P3 [< 2.5]: **YES — strict mechanism-missing**
- P4: NO

The framework does not merely under-spread; it **actively retracts**
in the long-time regime. This is the key qualitative finding the
Marangoni review must address.

---

## 1. Current Layer 4 implementation — accuracy assessment

Code: `acs/physics/mlsmpm.py:307–947`, `acs/runner.py:60, 80–83,
158–161, 275`. Specification: `docs/07_internal_flow_dynamics.md`,
`docs/03_adhesion_dynamics.md`.

### What we model
- γ(φ_p) = γ_max·(1 − φ) + γ_min·φ per particle, with γ_max = 0.020,
  γ_min = 0.003 in dimensionless units (Phase 4-v2 Lam4).
- γ_p scattered to grid (`grid_gamma`, mass-weighted P2G, kernel
  `_scatter_gamma_to_grid` lines 878–915).
- ∇γ via central FD on grid (`_compute_gamma_grad`, lines 918–940).
- Tangent projection `(I − n̂⊗n̂)·∇γ` at boundary cells
  (m > min_cell_mass) inside `_grid_op_overdamped`.
- Impulse `dv += dt·∇_s γ / ρ_local` applied as grid momentum increment
  before G2P.
- Pajic-Lijakovic & Milivojevic 2022 sign convention (flow toward high
  γ region; verified at Stage 1d via retraction direction test).

### What we do not model (gap inventory)
1. **No time dependence** in γ. Static map γ(φ); φ evolves on Cho 2020
   timescale (k_+ ≈ 1.4e-3 s⁻¹, k_- ≈ 4.6e-4 s⁻¹), but γ(φ) reacts
   instantaneously. There is no relaxation time for γ itself.
2. **Boundary-only**: Marangoni impulse applied only at cells with
   surface-normal indicator. Interior φ-gradient (Layer 3 spatial Sₚ
   extension) generates ∇γ in bulk that is currently discarded.
3. **Single source for γ**: only φ. Layer 6 ecm_strength does not enter
   γ at all (it modulates `γ_sub_eff = γ_sub_star · ecm_strength` for
   the substrate-CSF Layer 1a+ Option β, but not the cell-cell γ field).
4. **No surfactant-like transport**: γ-bearing field is not advected
   with interface motion; it is recomputed each step from instantaneous
   φ_p. The Stone 1990 / Aris 1962 surfactant transport equation
   (∂Γ/∂t + ∇_s·(Γ u_s) − D_s ∇_s²Γ = source) is not used.
5. **No instability mode**: γ field is computed deterministically from
   smooth φ; no mechanism to generate Bénard-Marangoni-style cellular
   pattern (Pearson 1958 Ma_c ≈ 79.6) once spreading slows.
6. **No osmotic / volume coupling**: Layer 5 ρ_osm and Layer 4 γ are
   independent. Yadav et al. 2022 (PRF) explicitly observed surface
   tension gradient ∝ mechanical surface strain, both volume-dependent.

### What the implementation does correctly
- Sign convention matches Pajic-Lijakovic & Milivojevic 2022: flow from
  low-γ to high-γ surface, which for our γ(φ) drives retraction once
  boundary φ saturates (consistent with our long-time behavior).
- Tangent projection prevents normal-direction artifacts.
- Mass-weighted P2G for γ keeps surface tension consistent with mass
  distribution.
- Phenotype ordering preserved (Bare < Pre < Lam4 spreading), confirmed
  in Phase 4-v2.

The **direction** of the implementation is correct. The **completeness**
is not. The peak-and-decay behavior is, in fact, the *natural*
consequence of an incomplete Marangoni model: φ saturates at the
boundary (φ_boundary ≈ 0.733 ≈ φ_eq = 0.751), γ becomes uniform on the
boundary band, ∇γ → 0 at the surface but persists from interior φ
decay (φ_min → 0.088), and the implementation discards the latter.

---

## 2. Six candidate mechanisms (A–F)

For each: literature anchor, mechanism, implementation feasibility,
expected long-time impact, paper-narrative value.

### A. Time-dependent γ via active stress reaccumulation
**Lit anchor**: Yadav et al. 2022, *Phys Rev Fluids* **7** L031101
("Gradients in solid surface tension drive Marangoni-like motions in
cell aggregates"). arXiv:2208.01090. Direct experimental observation:
laser-ablation-induced γ gradient drives toroidal flow (τ₁ ≈ 10 s) that
**spontaneously reverses** (τ₂ ≈ 40 s, τ₃ ≈ 70 s) as cells reaccumulate
active surface stress and return to original positions. **Faster than
viscoelastic relaxation timescale** — i.e., γ has its own dynamic.

**Mechanism**: γ is not a function of instantaneous φ but a state
variable that recharges on its own timescale. Replace γ(φ_p) with
ODE per particle:
```
dγ_p/dt = (γ_eq(φ_p) − γ_p) / τ_γ + α·|∇u_p|·(γ_max − γ_p)
```
where the first term is relaxation toward the φ-equilibrium γ and the
second is mechanical-strain-driven reaccumulation (Yadav's surface-
strain term).

**Feasibility**: HIGH. One new per-particle f32 field (γ_p), one ODE
integration alongside φ_p and ρ_osm_p. τ_γ ≈ 70 s ≈ τ_relax (1.0 in
star units), no CFL violation. No new grid kernels.

**Long-time impact**: directly addresses the peak-and-decay anti-
behavior. Reaccumulation prevents γ field from saturating; provides a
*driving* term beyond the boundary band and counters retraction. Yadav
observes reversal but in our spreading geometry (vs static aggregate),
the equivalent is *continued* extension when active stress reloads.

**Paper value**: **HIGHEST**. Direct experimental anchor (Yadav 2022)
that maps 1:1 to our long-time observation. The reversal-of-Marangoni
phenomenon is itself paper-worthy as the first continuum reproduction
of Yadav's mechanism in a spreading geometry.

### B. Multi-component Marangoni (Layer 6 ECM × Layer 4 γ_cc)
**Lit anchor**: Karbalaei et al. 2016 *Micromachines* 7:13 / Tan et al.
2017 *Langmuir* 33:5 1788 ("Marangoni Contraction of Evaporating
Sessile Droplets of Binary Mixtures") + Hu & Larson 2006
*J Phys Chem B* 110:7090 (multi-component selective evaporation
generates surface tension gradients that *reverse* coffee-ring
deposition).

**Mechanism**: combine cell-cell γ(φ) with substrate-coupled
γ_sub(ecm_strength). Currently `γ_sub_eff = γ_sub_star ·
ecm_strength` modulates only Layer 1a+ Option β substrate CSF; it does
not feed into the Layer 4 ∇γ. Adding cross-coupling:
```
γ_total(x) = γ_cc(φ(x)) + ν · γ_sub_eff(ecm_strength(x))
```
with ν a coupling weight (literature: same-order via Bertrand 2024
solid-surface-tension framing, see Bruckner/Yadav references).

**Feasibility**: MODERATE. Requires defining where the substrate γ
contributes (boundary band only? contact band only?), and ν calibration
(no obvious literature value — would trip Magic-Number Block test 1
unless derived dimensionally).

**Long-time impact**: at production scale ecm_strength has hit floor
0.10, so γ_sub_eff is uniformly low under the cell — no spatial
gradient. Mechanism would be most active *during* substrate degradation
(early phase) and saturate at the floor. Likely amplifies peak but not
decay.

**Paper value**: MODERATE. Couples two existing layers but the long-
time saturation limits explanatory power. Useful as a Stage 1d.b
sub-mechanism, not as the primary fix.

### C. Bulk Marangoni (interior ∇γ, not just boundary)
**Lit anchor**: Tan et al. 2017 *Langmuir* 33:5 1788 (binary-mixture
Marangoni *contraction*: bulk flow set up by surface-tension gradient
reaches into the droplet interior). Karbalaei et al. 2019 *PNAS*
116:53 26630 ("Marangoni spreading and contracting three-component
droplets") demonstrates bulk-coupled flow driven by composition
gradients.

**Mechanism**: drop the boundary mask in the ∇γ application. Currently
`(I − n̂⊗n̂)·∇γ` is applied only where the surface normal is defined;
extend to all cells with non-zero ∇γ, projecting onto the local
free-surface tangent or onto a φ-gradient direction (since φ is
the bulk source).

**Feasibility**: HIGH. Remove the boundary mask; add a `||∇γ|| >
threshold` filter to avoid applying force in nominally-uniform bulk.
The grid_gamma_grad field already covers the full grid.

**Long-time impact**: probably *worsens* the retraction. Interior φ →
0.088 means γ_interior → γ_max; surface φ ≈ 0.733 means γ_surface ≈
γ_min. Bulk Marangoni drives flow from low-γ (surface) toward high-γ
(interior) — a *contracting* flow. Without compensating spread driver
this would accelerate the long-time collapse, not arrest it.

**Paper value**: LOW–MODERATE. Mechanistically interesting (the
direction of bulk Marangoni is fundamentally different from boundary
Marangoni) but enacting it without compensation makes the asymptote gap
worse, not better. Useful as a *negative* control to demonstrate that
boundary-only Marangoni is the right level of approximation.

### D. Marangoni instability (Bénard-Marangoni pattern formation)
**Lit anchor**: Pearson 1958 *J Fluid Mech* 4:489 ("On convection cells
induced by surface tension"). Critical Marangoni number Ma_c ≈ 79.6
for hexagonal cellular pattern emergence in thin films with γ-T
coupling. Direct cell-biology analog: Yadav 2022 toroidal flow is a
single-cell version of Bénard-Marangoni in an aggregate.

**Mechanism**: when ∇γ is large enough, perturbations in γ field grow
and form cellular patterns at the spreading edge — converting smooth
spread into discrete intermittent extension events (analog to PI's
"leader cell" picking).

**Feasibility**: LOW–MODERATE. Pure continuum spread does not exhibit
Bénard pattern unless we either (i) add stochastic noise to γ (Magic-
Number trap unless noise amplitude is literature-derived) or (ii)
explicitly include a destabilizing thermal/φ analog (we have φ; Layer
4 already applies ∇γ, so the linear mechanism is already there but
in our overdamped Stokes-like regime the instability is suppressed).

**Long-time impact**: in principle could explain why discrete
lamellipodia events are needed (Stage 1a++.b). The continuum framework
*cannot* exhibit the instability above Ma_c in overdamped MPM without
explicit stochastic forcing.

**Paper value**: HIGH **as analysis**, LOW as implementation. The
finding "the continuum framework cannot exhibit Bénard-Marangoni
above Ma_c because of overdamped dissipation" is a clean theoretical
result that justifies Stage 1a++.b without needing to implement it.

### E. Surfactant-like dynamics (Stone 1990 transport equation)
**Lit anchor**: Stone 1990 *Phys Fluids A* 2:111 ("A simple derivation
of the time-dependent convective-diffusion equation for surfactant
transport along a deforming interface"). Aris 1962 derivation of the
interfacial transport theorem. Karbalaei 2019 *Curr Opin Coll
Interface Sci* 39 review of surfactant Marangoni dynamics.

**Mechanism**: γ becomes a function of an interface-bound concentration
Γ that itself satisfies a transport equation:
```
∂Γ/∂t + ∇_s·(Γ u_s) = D_s ∇_s²Γ + S(Γ, φ)
γ = γ_clean − R T Γ   (Frumkin or Langmuir EoS)
```
Currently Γ is implicitly Γ = γ_max − γ(φ_p), reset every step from
φ_p; the convective and diffusive transport terms are missing entirely.

**Feasibility**: MODERATE-HIGH. Requires: (i) per-particle Γ_p as a
new state variable, (ii) advection via standard MPM transport (free —
particles already advect), (iii) surface diffusion D_s ∇_s²Γ on
boundary band cells, (iv) source term S coupling to φ (E-cad ↔ Int-β1).

**Long-time impact**: solves the φ-saturation problem differently from
Mechanism A. Even when boundary φ saturates, Γ continues to be
transported with surface deformation, generating fresh ∇γ from
geometric stretching. This is exactly the missing convective-stretching
term identified in §1 gap (4).

**Paper value**: HIGHEST (tied with A). Stone 1990 is a foundational
microfluidics paper; importing it into cell-aggregate modeling is a
clean transferable contribution. Mechanism A and Mechanism E are not
mutually exclusive — A is the *temporal* relaxation, E is the *spatial
transport*. A complete model has both.

### F. Coupling with osmotic Layer 5 (volume-tension coupling)
**Lit anchor**: Yadav et al. 2022 *PRF* 7 L031101 (explicit observation:
surface tension gradient is *proportional to mechanical strain at the
surface*, and *both depend on aggregate volume*). Guo et al. 2017
*PNAS* 114:41 E8618 (osmotic-stiffness-coupling foundation already
used in Layer 5).

**Mechanism**: γ_p depends on local ρ_osm_p (cell volume / hydration
state) in addition to φ_p:
```
γ_p = γ(φ_p) · g(ρ_osm_p)
```
where g(·) is a weak monotone increasing function (volume loss → cell
stiffer → higher cortex γ, by Yadav). This couples Layer 5 spreading-
induced volume loss directly into the Marangoni driver.

**Feasibility**: HIGH. No new state variables; one extra multiplicand
in `_scatter_gamma_to_grid`. g(·) form via Yadav linear-strain
approximation: `g(ρ) = 1 + α·(ρ − 1)` with α derivable from Yadav's
surface-strain ↔ surface-tension proportionality (estimated α ≈ 0.5
from Yadav Fig. 3 slope; would need verify against the paper's
quantitative values to avoid Magic-Number trap).

**Long-time impact**: at production scale ρ_osm has *desaturated* (max
1.224 vs pilot 1.600) — i.e., spreading stops being driven by the same
ρ_osm gradient. Coupling γ to ρ_osm would damp this desaturation by
generating a γ-gradient response, partially restoring active driving.

**Paper value**: MODERATE-HIGH. The Layer 5 ↔ Layer 4 coupling is a
clean conceptual addition with direct experimental anchor in Yadav.
Less foundational than A or E but cleaner than B.

---

## 3. Summary table — mechanism scoring

| Mech | Lit anchor | IF | Implementation | Long-time impact | Paper value |
|---|---|---|---|---|---|
| A — γ time dependence | Yadav 2022 PRF | 2.6 | HIGH (1 ODE) | **directly addresses decay** | **HIGHEST** |
| B — Layer 6 × Layer 4 | Tan 2017, Hu-Larson 2006 | 4.7, 3.2 | MODERATE | weak (ECM floor) | MODERATE |
| C — bulk Marangoni | Karbalaei 2019 PNAS | 11.1 | HIGH (mask removal) | **negative** (worsens) | LOW–MOD (negative result) |
| D — Bénard instability | Pearson 1958 JFM | 4.0 | LOW (overdamped) | analytically informative | HIGH as analysis |
| E — Stone surfactant transport | Stone 1990 PoF | 4.6 | MODERATE-HIGH | **directly addresses decay** | **HIGHEST** |
| F — Layer 5 ↔ Layer 4 | Yadav 2022 PRF, Guo 2017 PNAS | 2.6, 11.1 | HIGH | moderate (ρ_osm desat) | MODERATE-HIGH |

The two mechanisms with both strong literature anchor *and* direct
relevance to the peak-and-decay are **A** (temporal γ relaxation /
reaccumulation) and **E** (spatial Stone transport). They are
complementary: A handles the time axis, E handles the convective-
stretching axis. **F** is a clean low-cost addition. **B**, **C**, **D**
are best framed as analysis / negative-result content.

---

## 4. Stage 1d.b vs Layer 4 extension — decision options

Two structurally distinct paths forward. Both can address peak-and-
decay; they differ in *what publication narrative* they support.

### Option α — Stage 1d.b (Layer 4 extension, continuum-only)
Add Mechanisms A + E + F to Layer 4. Stay within continuum framework.

- **Pros**: Self-contained, single-stage extension. Continuum framework
  remains intact. Yadav-anchored. Two new state variables (γ_p, Γ_p)
  plus one ρ_osm coupling. Estimated 200–400 LOC + sanity gate doc.
- **Cons**: Does not address Stage 1a++.b stochastic boundary events
  (lamellipodia / filopodia). The continuum framework has a *separate*
  identified gap that this option does not close.
- **Expected outcome**: peak-and-decay → monotone-spreading-toward-
  asymptote. Final A/A₀_topdown likely 2–4 (Bucket P2). Still below PI
  range [8, 33].

### Option β — Stage 1a++.b (stochastic boundary events)
Add discrete lamellipodia / filopodia / leader-cell stochastic events
at the boundary band. Layer 4 unchanged.

- **Pros**: Closes a different identified gap. PI experimental MCF7
  spheroids visibly exhibit these events; presence in framework is
  qualitative match.
- **Cons**: Stochastic events have no clean Magic-Number-Block-passing
  amplitude calibration without a literature anchor that gives
  per-event impulse magnitude. Risk of fitting amplitude to A/A₀ target.
  Does not fix the Marangoni saturation; the underlying retraction
  driver remains.
- **Expected outcome**: bursty extension events on top of the
  peak-and-decay continuum baseline. Mean trajectory may still decay.

### Option γ — Both, sequenced (1d.b → then 1a++.b)
Run Stage 1d.b first (continuum upgrade); evaluate; if asymptote still
< PI range, add Stage 1a++.b stochastic events on top.

- **Pros**: Each stage is independently publishable. Clear
  attribution: "continuum-only ceiling = X; with stochastic events =
  Y". Matches the Stage 1e Sim A vs Sim B comparison narrative.
- **Cons**: Two more production runs. Six weeks of bench time minimum
  if each takes a Stage-1d-equivalent investigation cycle.
- **Expected outcome**: clean two-axis decomposition of the missing
  ~5–20× factor.

### Option δ — Paper-as-is (no further code changes)
Frame the current Production Lam4 result as the *quantified validity
boundary* of the 5+1 layer continuum framework, analogous to Stage 1e
Sim A vs Sim B.

- **Pros**: No more code. Existing 27-simulation narrative is
  publication-complete. The peak-and-decay finding is itself
  publication-worthy as the first quantitative continuum-only
  asymptote in a 3D MPM cell-spreading model.
- **Cons**: Leaves the most interesting mechanism upgrade (Yadav-style
  γ reaccumulation) unimplemented. Reviewers may reasonably ask "what
  would the framework do *with* time-dependent γ?".

---

## 5. Paper-writing framing options

How the long-time finding fits into the manuscript narrative.

### Framing 1 — "Continuum reaches its ceiling" (Stage 1e analog)
The paper's central comparison is Sim A (1D radial ODE) vs Sim B (3D
MPM continuum) vs PI experimental. Production Lam4 finding adds a
third axis: 3D continuum vs PI experimental, with a quantified gap
attributable to stochastic / time-dependent mechanisms. This is the
*minimal* framing — needs no further code.

### Framing 2 — "Marangoni reversal in spreading geometry" (Yadav extension)
Implement Mechanism A (γ time relaxation) → reproduce a Yadav-style
reversal in a *spreading* (not static) geometry. New finding: the
spreading equivalent of toroidal-reversal is monotone extension when
the τ_γ matches the spreading rate, and decay when it does not. This
positions the paper as the first quantitative test of Yadav's
mechanism in a non-equilibrium geometry.

### Framing 3 — "Surfactant transport in cell aggregates" (Stone extension)
Implement Mechanism E (Stone 1990 transport equation) → demonstrate
that Γ-transport (separate from φ dynamics) is necessary for
sustained spreading. Position as the first import of microfluidic
surfactant theory into 3D cell-aggregate MPM.

### Framing 4 — "Bénard-Marangoni in the overdamped limit" (analysis only)
Derive and report that in the overdamped Stokes regime (Re ≪ 1, our
case), the Pearson 1958 instability is suppressed below a *much
higher* effective Ma_c than in inertia-dominated thin films. This
explains why discrete events (Stage 1a++.b) are necessary even when
the continuum has the linear ∇γ mechanism. Pure analysis; no
implementation needed.

The framings are not mutually exclusive. Framing 1 + Framing 4 is the
*minimum publication path with no further code changes*. Framing 1 +
Framing 2 + Framing 4 corresponds to Option α above.

---

## 6. PI decision options surfaced

1. **Option α + Framing 2 + 4** — implement Mechanisms A + F (Yadav
   anchored), defer E to follow-up paper. Single-stage code addition.
   Matches the "publication-grade single physics extension" pattern of
   Stage 1c → 1d. Estimated 1 week implementation + 1 production run.

2. **Option α-full + Framing 2 + 3 + 4** — implement A + E + F.
   Comprehensive Layer 4 upgrade. Estimated 2 weeks implementation + 1
   production run.

3. **Option γ + Framing 1 + 2 + 4** — Stage 1d.b first (Mechanisms A,
   F), evaluate, then optionally Stage 1a++.b. Most thorough; longest
   timeline.

4. **Option δ + Framing 1 + 4** — no further code; consolidate paper
   draft on the as-is 27-simulation narrative + Production Lam4
   finding. Fastest to manuscript; leaves Yadav extension as future
   work.

PI input required. Auto-entry FORBIDDEN per Production Lam4 outcome
auto-STOP.

---

## 7. Notes for the implementing agent (whichever option is chosen)

- Mechanism A: τ_γ ≈ 70 s (≈ τ_relax = 60 s, i.e., τ_γ_star ≈ 1.0).
  Verify against Yadav Fig. 2 timescales — do not Magic-Number-Block
  fail this.
- Mechanism E: D_s surface diffusion magnitude — Stone 1990 derivation
  gives D_s as an independent input. Literature value for cell-membrane-
  bound proteins ≈ 1e-2 μm²/s (Saffman-Delbrück); convert to dimension-
  less D_s* and run a CFL-style stability check on the diffusion term.
- Mechanism F: α coupling constant from Yadav 2022 Fig. 3 slope; needs
  paper read-through, not just abstract — defer until PI selects this
  option.
- Sanity Gate Protocol mandatory before first execution of any option.
- Hard Rule 9 (Magic-Number Block) applies to every new constant.
- The current `production_lam4_finding.md` (commit b506b57) labeling
  error (top-down vs contact_xy_hull) should be revised in a separate
  documentation pass once PI confirms direction.

---

## Literature anchors (consolidated)

- Yadav, P. et al. 2022. "Gradients in solid surface tension drive
  Marangoni-like motions in cell aggregates." *Phys Rev Fluids* **7**
  L031101. arXiv:2208.01090.
- Stone, H. A. 1990. "A simple derivation of the time-dependent
  convective-diffusion equation for surfactant transport along a
  deforming interface." *Phys Fluids A* **2** 111.
- Aris, R. 1962. *Vectors, Tensors, and the Basic Equations of Fluid
  Mechanics*. (interfacial transport theorem.)
- Pearson, J. R. A. 1958. "On convection cells induced by surface
  tension." *J Fluid Mech* **4** 489.
- Pajic-Lijakovic, I. & Milivojevic, M. 2022. "Marangoni effect and
  cell spreading." *Eur Biophys J* (already cited in our 07_internal_
  flow_dynamics.md).
- Tan, H. et al. 2017. "Marangoni Contraction of Evaporating Sessile
  Droplets of Binary Mixtures." *Langmuir* **33** 1788.
- Hu, H. & Larson, R. G. 2006. "Marangoni effect reverses coffee-ring
  depositions." *J Phys Chem B* **110** 7090.
- Karbalaei, A., Kumar, R. & Cho, H. J. 2019. "Marangoni spreading and
  contracting three-component droplets." *PNAS* **116** 26630.
- Guo, M. et al. 2017. "Cell volume change through water efflux impacts
  cell stiffness and stem cell fate." *PNAS* **114** E8618. (Already
  used as Layer 5 anchor.)
- Moeendarbary, E. et al. 2013. "The cytoplasm of living cells behaves
  as a poroelastic material." *Nat Mater* **12** 253. (Already used as
  τ_relax anchor.)
