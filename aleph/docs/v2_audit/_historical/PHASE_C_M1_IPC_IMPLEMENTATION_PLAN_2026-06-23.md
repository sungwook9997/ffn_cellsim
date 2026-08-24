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

# Phase C — M1 IPC contact: implementation plan (Warp DCM) — 2026-06-23

The night's 5 prototypes proved no penalty/barrier KERNEL solves the strong-bundle contact failure
(`PHASE_C_M1_CONTACT_FIX_DESIGN`). This is the actionable scope for the real fix, for PI to greenlight.
It is a genuine contact-mechanics project (multi-day), GPU-portable in Warp. NOT started (PI-scoped).

## Why IPC (Incremental Potential Contact, Li et al. SIGGRAPH 2020) over the penalty
The penalty contact equilibrates penetration where contact-force = applied-force; under the lit ×40/×167
bundles that is deep (pen 3.1). The ONLY family that GUARANTEES non-penetration is a barrier method with
**continuous collision detection (CCD)**: the barrier energy →∞ as gap→0, and a CCD-filtered line search
caps each step so no node ever reaches gap≤0. Then the bundle can pull arbitrarily hard and cells still
cannot interpenetrate — they flatten/resist at gap→0+ instead. This is exactly the regime-II contact.

## Components + mapping to the existing Warp engine
1. **Distance / closest-feature** (HAVE): `closest_bary` + the face hash-grid already give node↔nearest-
   face closest point + signed gap. Reuse for the gap `d`.
2. **Barrier energy + force + Hessian** (NEW kernel): IPC barrier `b(d) = -(d-d̂)² ln(d/d̂)` for
   `0<d<d̂` (d̂ = the activation gap, ~c_rep — derived). Force `−b'(d)` along the contact normal; the
   diagonal Hessian `b''(d)` feeds the implicit operator. Both analytic (no FD). Per nearest face only.
3. **CCD — continuous collision detection** (NEW, the hard part): given xₙ and the trial step Δx, find
   the earliest time-of-impact `t∈(0,1]` at which any node first reaches gap=0 against a moving triangle
   (node-vs-moving-triangle TOI = a cubic root, Ericson §5; or conservative-advancement). A Warp kernel
   per node over its grid-local other-cell faces → a per-node `t_node`; reduce to a global `t*`.
4. **Filtered line search** (NEW, in the implicit step): scale the accepted step by `α·t*` (α<1 safety)
   so the configuration stays strictly penetration-free. Replaces the current D8 displacement cap (which
   is a crude scalar version) with the true per-step TOI filter.
5. **Implicit Newton with the barrier** (EXTEND `device_cg`): the barrier Hessian is stiff but PSD and
   ANALYTIC — it goes into the matrix-free operator cleanly (unlike the noisy FD JVP). The adaptive
   barrier stiffness κ (IPC's `κ` update) keeps the barrier ~conditioned; this is where the device_cg's
   pAp-floor guard must be relaxed (the barrier Hessian is genuinely large near contact, not FD noise).

## Hard parts / risks (where the effort is)
- **CCD correctness + cost**: node-vs-moving-triangle TOI is the bulk of the work; must be conservative
  (never miss a collision) yet not over-restrict the step (which would crawl). GPU CCD is well-studied
  but non-trivial; budget the most time here.
- **κ adaptivity**: too small → penetration leaks; too large → CG ill-conditioned. IPC's κ-update +
  the analytic Hessian (not FD) should make this tractable, but it touches the just-stabilized device_cg.
- **Interaction with the frozen-grid implicit operator**: the grid is built on xₙ; CCD must use the same
  frozen connectivity for the step, consistent with the Hessian. Already the engine's pattern.
- **Adhesion coexistence**: keep the multi-face adhesion (it is genuinely multi-contact); only the
  excluded-volume becomes the barrier+CCD. De-cohesion stays cadherin-catch-bond-governed.

## Cheaper alternative to evaluate first (lower risk)
**SimuCell3D's own constraint contact**: relocate a coupled apposed node-pair to their average + share
forces (a hard non-penetration constraint, no barrier/CCD). Memory `project-dcm-physio-nodeface-contact`
notes this is what SimuCell3D does; the deleted `dcm_coupling_host.py` was a (node-NODE) attempt. A
node-FACE constraint version (project the node onto the face + a Lagrange/stiff-spring constraint in the
implicit operator) may reach non-penetration with less machinery than full IPC. Worth a prototype before
committing to CCD. (Earlier coupling attempts were node-NODE + not constraint-correct — this would be the
node-FACE constraint done properly.)

## Recommendation
1. First prototype the **SimuCell3D node-FACE constraint** (relocation/projection + implicit stiff
   constraint) — if it reaches pen<0.3 on the cadherin bundle, it is the cheaper fix.
2. If not, implement **full IPC** (CCD + barrier + filtered Newton) — the guaranteed-correct path.
Either is a scoped project needing PI sign-off (it re-baselines all contact equilibria + touches the
integrator). The night's prototypes have ruled out every cheaper kernel-only option.
