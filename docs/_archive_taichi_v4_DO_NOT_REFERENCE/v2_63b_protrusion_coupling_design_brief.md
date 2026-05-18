# 6.3b Protrusion-Coupled FA Dynamics — Design-Discussion Brief

Date: 2026-05-04 KST
Status: **decision-brief only — no model commitment**. Drafted in
implementation-work after a read-only design pass over
`docs/v2/v2_phase1_forward_roadmap.md`, `acs/v2/dynamics/focal_adhesion.py`,
and `acs/v2/protrusion.py`. Routed to **design-discussion** for an
adversarial round before any 6.3b code is written.

## Why this brief exists

The forward roadmap names 6.3b as **"protrusion-coupled FA maturation"**
(line 39) but does not define the coupling law, the coupling
direction, the scope boundary against the next milestone (closed-loop
Tier 1), or the deterministic-vs-stochastic policy. Locking those
decisions in implementation-work would be design-by-fiat and risks
re-entering the ECM-OL-1 collision pattern (`id=1086`) where two
agents simultaneously coded against incompatible mental models.

## What is already locked

- **6.3a static FA preflight** (`52622b4`, `f43c619`): deterministic
  per-FA `step_focal_adhesions_static`. State table §5 from
  `docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md`. `dt·rate ≤ 0.5`
  safety margin reused from active-contour. No closed-loop, no ECM
  update, no stochastic events.
- **`FocalAdhesionState.linked_protrusion_id`** field exists on the
  schema but **6.3a does not consume it** — this is the explicit hook
  6.3b will inhabit.
- **`ProtrusionEvent` schema** (acs/v2/protrusion.py): schema-only,
  carries `state` (growing/stalled/retracting/ended), `event_type`
  (lamellipodium/filopodium), `start_time_s`, `end_time_s`,
  `length_um`, `boundary_angle_rad`, `force_candidate_nN` (unused),
  and `associated_adhesion_ids: tuple[str, ...]` (FA↔protrusion
  bidirectional link). The schema docstring explicitly notes:
  *"dynamics (event hazard, growth/retraction rules, RNG
  reproducibility) are explicitly out of P0 scope and must write
  their own Sanity Gate before execution."*
- **Phase 1 scope split** (`docs/v2/v2_phase1_forward_roadmap.md`):
  closed-loop Tier 1 (FA → ECM, ECM bias) is the **next** milestone
  after "Separated dynamics modes". 6.3b stays inside Separated
  dynamics modes — protrusion-coupled FA only.

## The 4 open questions

### Q1 — coupling law shape

What does "protrusion couples to FA maturation/binding" actually
mean as a per-step update?

- **A.** Rate multiplier. `k_maturity_eff = k_maturity · (1 +
  α · protrusion_active_indicator)` and similarly for `k_bind_eff`.
  Caller supplies `α` and the per-FA indicator; absent indicator ⇒
  no multiplier (back to 6.3a).
- **B.** State-conditional gate. FA whose `linked_protrusion_id`
  resolves to a protrusion in `state="growing"` follows the
  6.3a-mature rate; FA whose link is `retracting`/`ended` is forced
  toward `slipping`/`unbound`.
- **C.** Force scale. `traction_scale_eff = traction_scale_nN · f(
  protrusion.force_candidate_nN)` per linked FA.
- **D.** State trigger. A protrusion `state` transition (growing →
  stalled, etc.) triggers a deterministic FA state transition
  (nascent → mature on stall, mature → slipping on retract).

A and B are the cleanest per the §5 state-table spirit. C reaches
into 6.3a's traction algebra. D introduces hidden state-transition
machinery.

### Q2 — coupling direction

- **One-way protrusion → FA**: caller supplies a protrusion state /
  timeline; FA reads but never writes back. Sanity Gate burden
  identical to 6.3a's (deterministic explicit-rate preflight).
- **One-way FA → protrusion**: FA `bound_fraction` collapse triggers
  protrusion end. New write-path on the protrusion side.
- **Bidirectional (closed-loop)**: both writes simultaneously.
  Saturation / oscillation / clipping risk per the roadmap "Top
  Risks" §2. Belongs in Tier 1 milestone, not 6.3b.

### Q3 — scope boundary against Tier 1

The next milestone is "Closed-loop Tier 1 single-cell: protrusion
force, FA-protrusion coupling, FA to ECM, ECM bias". 6.3b must
**stay inside Separated dynamics modes**:

- 6.3b: protrusion → FA maturation only.
- 6.3b: **NO** FA → ECM. The ECM update path is owned by Tier 1.
- 6.3b: **NO** ECM bias on protrusion direction. Owned by Tier 1.
- 6.3b: **NO** protrusion-derived force computation that feeds back
  onto active-contour vertices. Owned by Tier 1.

If 6.3b includes any of the three above, it has silently absorbed
the Tier 1 boundary and the "Separated dynamics modes" milestone is
no longer a meaningful separation.

### Q4 — deterministic vs stochastic

6.3a is deterministic (caller-supplied rates, deterministic explicit
update). The protrusion lifecycle is biologically stochastic
(nucleation hazard, branching, retraction). For 6.3b:

- **Deterministic preflight**: caller supplies a fully-prescribed
  protrusion timeline (similar to ECM-OL-1's caller-supplied
  prescribed traction). No RNG. Sanity Gate matches 6.3a / ECM-OL.
- **Stochastic event sampler**: 6.3b owns RNG seeding,
  reproducibility contract, and event-hazard parameters. Sanity
  Gate adds §3 reproducibility, §4 RNG-CFL bound.

The schema-only protrusion module's docstring already commits the
project to writing a *separate* Sanity Gate for stochastic dynamics.
A deterministic 6.3b preflight defers that commitment cleanly.

## Recommended opening position (Codex `id=1183`)

- **Q1**: A + B blend — rate-multiplier coupling (A) gated on
  `linked_protrusion_id` resolving to a non-`ended` protrusion (B
  cheapest version). C and D rejected as out of scope for a
  preflight.
- **Q2**: **one-way** protrusion → FA. No FA → protrusion writeback
  in 6.3b.
- **Q3**: hard boundary. **NO** FA → ECM, **NO** ECM bias, **NO**
  protrusion force on contour in 6.3b. Tier 1 owns those.
- **Q4**: **deterministic** preflight. Caller-supplied protrusion
  state/timeline as test inputs. No RNG. Stochastic sampler
  deferred to its own Sanity Gate.

This default keeps 6.3b on the same Sanity Gate template as 6.3a /
ECM-OL — **deterministic, one-way, caller-supplied, no closed-loop,
no RNG** — and reserves every closed-loop affordance for the Tier 1
milestone where the roadmap already plans the appropriate gate
work.

## Sanity Gate items 6.3b will need (preview, not commitment)

- §1 dimensional: caller-supplied indicator/timeline must have a
  defined dimensionless or 1/s unit chain through the multiplied
  rate.
- §2 boundary: missing `linked_protrusion_id` ⇒ 6.3a behavior;
  protrusion not in supplied registry ⇒ explicit failure_kind.
- §3 conservation: if Q2 is one-way, FA aggregate Newton-3 still
  holds because traction algebra is unchanged; only the rate path
  is altered.
- §4 numerical: `dt_fa · max(rate_eff)` must respect the same
  `_DT_RATE_SAFETY_MARGIN = 0.5` as 6.3a (the multiplier inflates
  the effective rate ceiling).
- §5 sign: rate multiplier must be ≥ 0 and finite; never inverts
  the sign of `k_maturity` / `k_bind`.
- §6 measurement-protocol: per-FA `linked_protrusion_id` round-trip
  through any frame_dump must preserve the multiplier audit chain
  on read-back.

## What this brief is *not*

- Not a model commitment. The opening position above is the
  starting point for a design-discussion adversarial round, not the
  locked contract.
- Not a Sanity Gate document. The Sanity Gate items above are a
  preview to make the gate burden of each Q1–Q4 choice concrete.
- Not an implementation plan. No file paths, no API signatures, no
  test list. Those follow only after design-discussion locks.

## Process expectation

1. **Cross-room dispatch** to `design-discussion` with this brief +
   the four questions + the recommended default.
2. **Adversarial round** in `design-discussion` per the precedent
   set by the Phase-1 forward-roadmap discussion (PI `id=809` /
   `id=814`). Open with strong defaults, expect challenges,
   negotiate concessions, list unresolved disagreements before lock.
3. **Lock artifact** delivered back to implementation-work as a
   `docs/v2/v2_63b_protrusion_coupling_locked.md` note (mirroring
   `docs/v2/v2_p1_derivation_locked.md` for P1 alpha).
4. **Implementation entry** only after the lock artifact exists
   and is mutually acked.
