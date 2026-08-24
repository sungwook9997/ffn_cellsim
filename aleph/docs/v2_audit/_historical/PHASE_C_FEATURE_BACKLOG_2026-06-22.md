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

# Phase C — 8-feature implementation backlog + session state (2026-06-22)

**Boot artifact for context-compression recovery.** PI directed (2026-06-22): implement ALL
of A1, B2, B3, C4, C5, D7, D8, E "지금". This file is the durable plan + live status.

## HARD rules in force (do NOT violate)
- **No parameter-tuning to outcome** (PI 2026-06-22): never sweep/lower a DERIVED param
  (cad-bundle, adh-strength, …) to make spreading/a gate pass = magic-number sin. Legit levers
  ONLY: (1) time/longer integration, (2) add a missing REAL force at its derived value, (3) fix
  a real bug, (4) else surface the limit to PI. R-sweep IS legit (R=controlled variable).
  See memory `feedback-no-param-tuning-to-outcome`.
- Mechanistic de-cohesion must stay catch-bond-governed (cadherin Rakshit catch-slip), not a
  binary/geometric switch. See memory `feedback-junction-switch-fine-grained`.
- Physiological baseline; visualize at every closeout/finding (PI keeps asking for figures).
- All work is on the **Warp DCM engine** (aleph/warp_port/), GPU-resident. HOOMD = parity ref.

## Root finding driving all of this
The proxy-free spheroid does NOT properly aggregate: cells stick + compact slightly but stay
ROUND (Ψ≈0.99 "bag of marbles"), NOT a flattened regime-II tissue. Converged (not timescale).
Cause = sparse cadherin node-NODE point-bonds (~6/junction, mesh-limited) → point contacts.
Fix = SimuCell3D continuous node-FACE adhesion. (You can't spread a bag of marbles like tissue.)

## The 8 features — mechanism change + status

### A1. cell-cell node-FACE coupling ⭐ (IN PROGRESS)
- **Mechanism**: sparse cadherin node-NODE bonds (point contacts, round cells) → **re-enable the
  continuous node-vs-FACE bilinear adhesion** (`contact_grid_kernel`'s `adh` branch = SimuCell3D
  Model-1 traction-separation, mesh-independent) so EVERY apposed node↔face adheres → flat shared
  interface → regime-II flattening. Cadherin catch-bonds stay ON for mechanistic de-cohesion.
- **Status**: IMPLEMENTED as `--coupling` flag + `--adh-strength` (driver). The node-FACE adhesion
  ALREADY existed in the kernel but was disabled (`coh_adh=0`) in cadherin mode. Re-enabling works:
  N=13 free-space aggregation Ψ_final: no-coupling 0.989 → coupling(adh1e7) **0.980** → testing
  lit-anchored adh 4e7/5e7 (memory: DCM-aggregation-resolved gave clean spheroids at 4-5e7; NOT
  tuning — documented value). Awaiting that result to confirm regime-II (Ψ→~0.9) then COMMIT.
- **Open**: de-cohesion reconciliation — layered (node-face geometric release + cadherin catch-slip)
  for v1; the proper unified design is to MODULATE node-face ω by the cadherin catch-bond state
  (node-face geometry + catch-slip kinetics). Surface to PI.

### B2. lamellipodium → explicit substrate clutch (node-to-plane)  [PENDING]
- **Mechanism**: leading node tethered to an abstract advancing actin-bead anchor → the
  lamellipodium's nascent adhesions = integrin clutches gripping the dish (node-to-plane, SAME
  mechanism as the ECM clutch). Removes the abstract tether; traction emergent from substrate clutch.

### B3. filopodia node-FACE tip contact (NEW module)  [PENDING]
- **Mechanism**: NEW. Explicit filopodial protrusions; tips probe + contact cell FACES (node-face)
  and substrate (node-plane), forming adhesions. (Unified-actin-architecture, memory.)

### C4. ligand-density (Bare/Pre/Lam4) → ECM clutch  [PENDING]
- **Mechanism**: clutch k_on/engage-density FIXED → scale by ligand coating condition → reproduce
  the experimental Bare/Pre/Lam4 adhesion-density differences. Values must be lit/experiment-anchored.

### C5. ECM mechano-feedback (strain-stiffen + substrate-E mechanosensing)  [SUBAGENT DONE — review/commit]
- **Mechanism**: rigid plane, fixed clutch stiffness → (1) Chan-Odde substrate-E mechanosensing
  (traction factor E/(E+E_opt), E_opt~5kPa KB-2.8; k_sub=2E*a Hertzian), (2) Storm-MacKintosh
  strain-stiffening hook (K~σ, cap ×100). ADDITIVE, DEFAULT-OFF → rigid-dish runs byte-identical.
- **Status**: subagent implemented in `dcm_substrate_warp.py`, self-test PASS, parity 3/3, KB-anchored
  (KB-2.4/2.8/1.5/1.4). NOT committed. **TODO: review + commit.**

### D7. gravity/buoyancy sedimentation body force  [✅ DONE, committed b6b0a8c]
- **Mechanism**: absent → `f_z=-Δρ·g·v_node` per live node, Δρ≈55 kg/m³ (derived). ~1pN/cell. RHS
  only. `--gravity` (default off). `gravity_body_force_kernel` in dcm_neighbor_warp.

### D8. pen interpenetration fix (G2 FAIL with strong forces)  [PENDING]
- **Mechanism**: strong bundle forces → persistent interpenetration (pen~1.3 > 0.3). Fix: contact
  repulsion stiffness / adaptive dt / active-set so cells don't overlap.

### E. Young-Dupré triplet gate + R-sweep  [SUBAGENT DONE — review/commit]
- triplet gate (cos(φ/2)=η/2, symmetric→120°) added to young_dupre_doublet_gate.py; R-sweep script
  spread_rsweep.py (A/A0=a+b/R+c/R², params FIXED, R=N varied). Subagent done. **TODO: review/commit.**

## Already committed this session (Phase C)
917e51d implicit CG robust (3 guards) · 599a52b decisive (i)/(ii) spread test · 7dad9f2 traj plot ·
b0344b9 bundle FORCE fix (cad40/ecm167, KB-anchored, peeling→cohesive) · 4712693+4a73695 Young-Dupré
oracle+gate (PASS) · 567658e spread_eval tooling · b6b0a8c D7 gravity. (A1 coupling staged, not committed.)

## ✅ ALL 8 FEATURES COMMITTED (2026-06-22)
- D7 gravity/buoyancy `b6b0a8c` · A1 node-FACE coupling `1e39a04` (+actual-cells viz `b0dfbc9`)
- C5 ECM mechano-feedback `fef5cb8` · E triplet+R-sweep `68a22a6` · C4 ligand-density `76f4817`
- D8 pen displacement cap `1e69e2b` · B3 filopodia `d5e593d` · B2 lamellipodium substrate-clutch `a1b92b7`
All GPU-resident Warp, additive/default-off (back-compat), lit/KB-anchored (no outcome-tuning).
New CLI flags: --coupling --adh-strength --gravity --delta-rho --ligand-density --pen-cap-frac
--no-pen-cap --filopodia --lamel-clutch. New modules: dcm_coupling_host.py (node-NODE, SUPERSEDED
by A1's node-FACE re-enable — can delete), dcm_filopodia_{host,warp}.py.

### Documented FOLLOW-UPS (not blocking; surface to PI)
- A1 de-cohesion reconciliation: layered (node-face geometric + cadherin catch-slip) → proper =
  modulate node-face ω by the cadherin catch-bond state (node-face geometry + catch-slip kinetics).
- A1 flattening is REAL but MODEST at lit adh (Ψ 0.989→0.956 at N=13; N=48 interior-cell test
  running to judge proper regime-II — PI flagged N=13 too small vs SimuCell3D).
- D8: validate pen 1.3→low on the N=100 strong-force production config (small tests have pen=0).
- B2: full Pereverzev catch-slip on the lamellipodial clutch (v1 = anchor on dish plane only).
- B3: tip-adhesion stiffness/cap/detach PROVISIONAL (no KB datum).
- C5: needs the deformable/3D substrate (C6 explicit Mikado ECM Warp port) to be meaningful.
- Implicit accelerator ~2-3× net (vs the old native 44× BAOAB); host-hybrid binders (cadherin/
  ecm/lamel KDTree at cadence) are the remaining GPU-main port targets.

## NEXT (resume order)
1. Read A1 lit-adh result (baveae8vy task) → if Ψ→~0.9, COMMIT A1 + viz montage. Else report limit.
2. Review + commit C5 (subagent, dcm_substrate_warp.py) and E (subagent, gate+rsweep).
3. Implement B2, C4, D8, B3 (driver-touching, sequential, myself).
4. Visualize EACH (PI keeps asking): aggregation flattening before/after, montages, gate curves.
5. Then the production spreading test with the FULL stack (A1+D7+...) — does the tissue now spread?
6. KB-check + Notion closeout + receipt line.

## Key file pointers
- Driver: aleph/warp_port/dcm_warp_decohesion.py (--coupling --adh-strength --gravity --cad-bundle --ecm-bundle)
- Kernels: dcm_neighbor_warp.py (contact_grid_kernel adh branch, gravity_body_force_kernel)
- Oracle: validation/oracles/young_dupre.py · gate: scripts/young_dupre_doublet_gate.py
- Eval: scripts/spread_eval.py · agg compare: scripts/agg_quality_simucell3d_compare.py
- gbook: ssh gbook, ~/ffn_phase_c/, A5000, runs via setsid + poller
