# Hard Blocker #3 — FA→ECM Scattering Geometry — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-04 KST
**Status**: design-discussion brief, **opening position only — NOT a
lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under the (B-trigger)
conditions of design-discussion `id=1295/1296` then `id=1326/1328`:
design-discussion has been dormant since 14:30 KST (>5h50m), PI
`id=1300` decision-needed has been pending ~1h55m without response,
and Codex impl explicitly agreed to design-only cross-room dispatch
in `id=1328`.
**Source**: closed-loop ECM gate phased plan
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 (Hard
Blocker #3 enumeration), §6 (Sanity Gate matrix Phase D row), §9
(Phase D blocker resolution checklist).

**Hard contract**: this brief is the **opening position** for the
adversarial design round on Hard Blocker #3. It does not lock the
interface, does not commit the implementation, and does not
authorize any Phase D code. The locked phased plan §1 Phase D
remains blocked until the design-discussion round closes with a
lock artifact (target:
`docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`).

---

## 0. Why this brief exists

The locked phased plan §2 enumerates Hard Blocker #3 as the
"FA→ECM scattering geometry" interface decision: 1 FA at position
`(x, y)` → which ECM grid cells receive its traction signal? With
what kernel? The plan §1 Phase D entry is gated on this blocker
(plus #4) being interface-locked — even if the default Phase D
implementation is identity (no-op).

Without this round, Phase D no-op scaffolding cannot land, which
in turn blocks Phase E active closed-loop response, which is the
end-state Phase 1 deliverable in the closed-loop ECM gate.

The brief routes to design-discussion because the kernel choice
is a **physics-relevant interface**: it shapes what the cell can
do to the ECM and constrains what response laws are coherent. A
nearest-neighbor scattering says "FA traction is a point load on
one ECM cell"; a Gaussian kernel says "FA traction is a smeared
load on a region"; the two have different scaling under refinement
and different coupling to subsequent constitutive laws.

---

## 1. What is already locked (do not redebate here)

Per `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`:

- §2 Hard Blocker #3 itself: phrased as a question ("which ECM
  grid cells receive the FA's traction signal?") with four
  candidate kernels listed (nearest / bilinear / Gaussian / other).
  Phase D entry contract:
  - Inputs: `(adhesions, ecm_state, scattering_params)`
  - Outputs: `traction_density_nN_per_um2[grid_shape]` (units locked)
  - Default no-op: returns zero field (identity)
- §3 effective_stiffness helper guard: any helper that consumes
  `fiber_density` / `orientation_tensor` is ECM→FA bias by another
  name — that is Hard Blocker #4, not #3, but **the kernel choice
  must not bleed into FA-side reads**. The scattering function is
  one-directional FA → ECM, opaque to ECM-side state.
- §6 Sanity Gate matrix Phase D row: per-function units / boundary
  / conservation / numerical / sign / measurement-protocol items
  must each be satisfied at Phase D no-op scaffolding time.
  Specifically: input/output unit chain (`nN per FA` →
  `nN/μm² per ECM cell`), shape/finite/non-negative boundary
  cases, sum conservation (∑ over grid cells equals total FA
  traction projected through the kernel), no auto-shrink, sign
  convention consistent with 6.3a's inward radial-to-centroid
  default.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four candidate kernels (pros/cons)

### Candidate (a) — Nearest-neighbor

A single FA at world position `(x_FA, y_FA)` deposits its full
traction onto exactly one ECM cell — the one whose center is
closest to `(x_FA, y_FA)` (or whose footprint contains the
position).

**Pros**:
- Simplest implementation; no kernel parameter beyond the
  schema's existing `spacing_um`.
- Trivially passes Sanity Gate §3 conservation (each FA's
  traction lands intact on exactly one grid cell).
- Defensible as a "point load" approximation when FA size ≪
  grid spacing.

**Cons**:
- Discontinuous: an FA crossing a cell boundary causes a
  step-function jump in which cell it scatters into. Under grid
  refinement the discontinuity stays. Item 5 sensitivity sweep
  results will show as artifacts of this.
- Anisotropic at cell boundaries — neighboring FAs near the
  same grid line receive identical scattering even with very
  different positions.
- Hard to extend to "one FA whose footprint is comparable to
  the grid" without revisiting the choice.

### Candidate (b) — Bilinear

A single FA deposits its traction onto the **four ECM cells**
whose corners enclose the FA position, weighted by the bilinear
weights of the FA's fractional position within the enclosing
cell. Each weight is the area of the rectangle opposite the
target cell, divided by the cell area; weights sum to 1.

**Pros**:
- Continuous: an FA crossing a cell boundary smoothly
  redistributes weight between the two adjacent cells; the
  scattered traction field has no step-function jump.
- Trivially passes Sanity Gate §3 conservation (the four
  weights sum to 1 by construction).
- Standard finite-element interpolation pattern — derivable
  from first principles on a structured grid; no fitted
  parameter.
- Well-defined at the ECM grid edge: an FA at a corner cell
  receives zero weight on the absent neighbors (or equivalently,
  the absent neighbors are clamped to zero contribution).

**Cons**:
- Smears the load over four cells even when FA size ≪ grid
  spacing; for a sharp point load this is a slight inaccuracy.
- Edge-cell behavior is technically a boundary condition that
  must be explicitly tested (Sanity Gate §2 boundary cases).

### Candidate (c) — Gaussian kernel

A single FA deposits its traction over **N ECM cells within
radius `r_kernel`** weighted by a Gaussian of width
`σ_kernel`. The kernel parameters `(r_kernel, σ_kernel)` are
caller-supplied per the phased plan's caller-supplied discipline
(no project defaults).

**Pros**:
- Smoothest scattering field; suppresses high-frequency
  artifacts in the response field even more than bilinear.
- Defensible from a "FA cluster has finite size on the
  substrate" argument (each FA already represents some local
  averaging).
- Tunable: caller can set `σ_kernel` to match FA size or to
  control coupling range.

**Cons**:
- Two extra caller-supplied parameters (`r_kernel`,
  `σ_kernel`) — Magic-Number Block requires that each is
  derivable from literature or first principles, grid-invariant
  (or its grid-dependence must be documented), and not chosen
  to fit a target. Bigger Magic-Number Block exposure than (a)
  or (b).
- Conservation requires explicit normalization step (the
  truncated Gaussian on a finite grid does not sum to 1 a
  priori); a mistake here breaks Sanity Gate §3.
- Sanity Gate §6 measurement-protocol consistency: the
  effective FA "footprint" the experimentalist measures (e.g.,
  fluorescent paxillin patch size, ~0.5–2 μm) becomes a
  reference for `σ_kernel` choice — a design-discussion item.

### Candidate (d) — Other (literature-derived)

Reserved for any kernel form that has explicit literature
backing: e.g., Green's-function-based scattering on an elastic
substrate, distance-decay forms `1/r²`, stress-from-traction
inversion via Boussinesq solutions. These are heavier and
typically belong to a Phase 2+ unit; flagging here for
completeness so the design round can explicitly reject them in
this round if appropriate.

**Pros**:
- Physically richer if the eventual goal is comparison against
  traction force microscopy (TFM) reconstructions.

**Cons**:
- Heavier Sanity Gate burden (each form has its own conservation
  and unit-chain proofs).
- Unnecessary for Phase D no-op scaffolding (the no-op default
  returns zero field regardless).

---

## 3. Recommended opening position (NOT a lock)

**Candidate (b) bilinear** as the recommended opening position
for the adversarial round.

Reasoning:
1. **Magic-Number Block discipline**: zero caller-supplied
   numerical parameters beyond `spacing_um` (already in the
   schema). No new tunable lands. Compare to (c) which adds
   `(r_kernel, σ_kernel)` and shifts Magic-Number Block work.
2. **Sanity Gate §3 conservation**: trivial, no normalization
   step required. Compare to (c) which requires explicit
   truncated-Gaussian normalization with its own derivation.
3. **Sanity Gate §6 measurement-protocol**: bilinear is the
   standard FE-style interpolation; no new measurement modality
   is invented. The FA position-to-grid-fraction mapping is the
   same mapping used by `make_default_ecm` and existing ECM-OL
   preflight layouts. Compare to (a) which introduces a
   discontinuous mapping that violates §6 modality matching with
   respect to refinement studies.
4. **Phase 1.5 / Phase 2 extensibility**: nothing in (b)
   precludes a later switch to (c) or (d) for a different
   constitutive law. The Phase D no-op scaffolding interface
   contract is identical across choices; only the kernel function
   body differs.
5. **Item 5 sensitivity sweep**: bilinear scattering is expected
   to converge under grid refinement at first-order accuracy
   (kernel error scales as `dx`); this is a clean baseline for
   the closed-loop sensitivity sweep that Phase E adds. (a) is
   non-convergent at cell boundaries.

**Caveat**: the recommendation is conservative and chosen for
gate-discipline reasons, not biological accuracy. If the
design-discussion round wants to argue (c) is biologically more
faithful (FA paxillin patch ~0.5–2 μm versus typical ECM
spacing 0.5–5 μm), that is a valid round 2 challenge — the
opening position is provisional pending literature and the §6
measurement-protocol decision (TFM-style scattering vs
point-load).

---

## 4. Interface contract (restated, locked elsewhere)

Per `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Phase
D entry contract (NOT changeable in this round):

```
def step_fa_to_ecm_response(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm_state: ECMSubstrateState,
    scattering_params: ScatteringParams,
) -> np.ndarray:  # shape == ecm_state.grid_shape, units nN/μm² per ECM cell
    ...
```

- Inputs: `(adhesions, ecm_state, scattering_params)`. The
  `scattering_params` dataclass / mapping is what the round
  decides — its shape depends on which kernel is chosen.
- Output: `traction_density_nN_per_um2` array of shape
  `ecm_state.grid_shape`, units `nN/μm²` per ECM cell.
- Default no-op: returns `np.zeros(ecm_state.grid_shape,
  dtype=np.float64)` — the round must lock this default exactly
  so Phase D no-op scaffolding has nothing to compute beyond the
  zero-array allocation.

---

## 5. Sanity Gate items 6.4-CL-3 will need (preview, NOT commitment)

Per Phase D row of locked phased plan §6 Sanity Gate matrix
applied to the chosen kernel:

- §1 Dimensional: per-FA `nN` → per-cell `nN/μm²` reduction.
  Inline derivation: `nN / (spacing_um · spacing_um)` =
  `nN/μm²`. Hard Rule 10 protected by inline derivation in the
  module docstring.
- §2 Boundary: FA at the exact grid corner / edge / outside the
  ECM extent (out-of-grid policy: clamp to grid? raise? — round
  decides). Empty FA list: returns zero field. Each input
  validation has an explicit failure_kind.
- §3 Conservation: `np.sum(traction_density_nN_per_um2) ·
  spacing_um²` equals the total FA traction magnitude projected
  along the scattering direction; bilinear satisfies this by
  weight-sum-equals-1; nearest-neighbor satisfies trivially;
  Gaussian requires explicit normalization step.
- §4 Numerical: float64 throughout; per-FA scattering is a
  bounded number of cell writes (2×2 for bilinear, 1 for
  nearest, ~`(r_kernel/spacing_um)²` for Gaussian). No dt-rate
  gate at this step (Phase D scattering is a single-step
  algebra, not an explicit-rate update).
- §5 Sign: scattering preserves the sign convention 6.3a uses
  for FA traction (cell-on-substrate inward radial-to-centroid
  default).
- §6 Measurement-protocol: the scattered field's modality must
  match the modality the response law (Phase E) consumes — i.e.,
  if the response law expects `nN/μm² per ECM cell`, the
  scattering returns exactly that. Hard Rule 11 protected by
  a runtime meta-test analogous to Phase B/C precedents
  (`test_scattering_does_not_satisfy_closed_loop_response`).

Magic-Number Block: bilinear is parameter-free (good); Gaussian
adds two parameters (each requires Magic-Number Block proof);
nearest is parameter-free but discontinuous (Sanity Gate §6
modality concern).

---

## 6. What this brief is *not*

- Not a lock. The interface kernel choice + parameters are
  decided by the design-discussion adversarial round; this brief
  opens with a recommendation and expects challenges.
- Not a Sanity Gate document. The Sanity Gate is written by the
  future Phase D 6.4-CL-3 commit, after the lock artifact lands.
- Not a code commitment. No `step_fa_to_ecm_response` function
  is proposed in this brief; the locked plan §2 already
  specifies the function signature, and this round only fills in
  the kernel body's shape.
- Not a closed-loop activation. Phase D no-op scaffolding stays
  default-off until Phase E lands. Activating the scattering does
  not flip on closed-loop ECM remodeling — that's Phase E's
  separate response-law lock.
- Not authoritative for biological parameters. Any kernel
  parameter (e.g., Gaussian `σ_kernel`) requires literature-
  first extraction per Magic-Number Block.

---

## 7. Process expectation (mirrors Hard Blocker #1 / 6.3b precedents)

1. **Cross-room dispatch** to design-discussion with this brief
   + the 4 candidate kernels + the recommended (b) opening
   position. Status: `decision-needed`.
2. **Adversarial round** in design-discussion (Codex/Claude
   adversarial design debate per memory
   `feedback_aggressive_design_debate.md`): challenges, reasoned
   acceptances, candidate compromises.
3. **Unresolved disagreements list** before lock — each unresolved
   item gets a separate sub-round if needed.
4. **Lock artifact** delivered to implementation-work as
   `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`,
   parallel to `docs/v2_63b_protrusion_coupling_locked.md` /
   `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
5. **Cross-room dispatch back to implementation-work** with the
   locked kernel + scattering_params dataclass shape. Phase D
   no-op scaffolding entry permitted.
6. **Hard Blocker #4 round** runs in parallel or sequentially —
   per locked phased plan §2 Phase D requires both #3 + #4
   resolved.

implementation-work stays idle/review-capable while this design
round runs. Per the (B-trigger) ack from Codex impl `id=1326/1328`,
this brief is the impl-work side's "opening dispatch", explicitly
not a lock; the design-discussion round retains all design
authority.

---

## 8. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- 6.3b lock precedent (analogous design-only adversarial round):
  `docs/v2_63b_protrusion_coupling_locked.md`
- 6.3b Sanity Gate (Sanity Gate template post-implementation
  status precedent): `docs/v2_63b_protrusion_coupled_fa_sanity_gate.md`
- ECM-OL preflight (FA traction algebra layer that scattering
  feeds into): `acs/v2/dynamics/ecm_open_loop.py`
- 6.3a static FA traction (the per-FA traction the scattering
  consumes): `acs/v2/dynamics/focal_adhesion.py`
- Forward roadmap: `docs/v2_phase1_forward_roadmap.md`
  ("Closed-Loop ECM Gate" section)
- Memory rules:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md`
  - `feedback_aggressive_design_debate.md`

---

## 9. Cross-room dispatch instruction (impl-work → design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex` (design-discussion-pane Codex receives + Claude
  pane reads)
- `room`: `design-discussion`
- `topic`: `v2-layer-2-hard-blocker-3-fa-to-ecm-scattering`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendation +
  reference to this file
- `refs`: implementation-work `id=1300, 1326, 1328` + locked
  phased plan reference

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context (mirror of `docs/v2_63b_protrusion_coupling_design_brief.md`
post-lock).
