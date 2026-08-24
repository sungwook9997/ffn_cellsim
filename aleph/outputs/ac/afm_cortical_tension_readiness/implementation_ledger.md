# Cortical-tension AFM implementation ledger

This ledger maps the seven requested preparation items to reviewable commits. `DONE` means the engine contract
or numerical instrumentation exists and has passed its stated gate; it does not promote a mechanism demo to a
quantitative biological result. `PI GAP` means the runtime deliberately has no default.

| Item | Status | Implementation / evidence | Commit |
|---|---|---|---|
| 1. `Pi_0` provenance | PI GAP guarded | `OsmoticSetpoint` requires value, evidence class and source; a cross-protocol proxy can enumerate only a mechanism demo. Chugh does not source engine `Pi_0`. | `4512278e`, `7bc608ce` |
| 2. Round geometry + holder | DONE contract; viscosity PI GAP | The declaration refuses adherent/spread geometry and FA holding, requiring suspended-round geometry and the T10 exterior-Stokes six-mode holder. Exterior viscosity is required without a default. | `4512278e` |
| 3. Sixth contact edge | DONE runtime/hook; calibration PI GAP | The spherical apparatus calls the existing gap-gated `_unilateral_contact_force_kernel` on every exact inner-solver trial; radius, surface gap and penalty stiffness are sourced inputs. RTX 3090 job 56 passed sign, global-block isolation, adjoint and transaction checks. | `b4589a61`, `7093050b`, `f7550895` |
| 4. Reaction readout | DONE | Device `reaction[0]` is read only between accepted steps; rejected candidates emit no partial curve. | `b4589a61`, `7093050b`, `f7550895` |
| 5. Physical loading speed | DONE mechanism-demo path; apparatus value PI GAP | Accepted depth increments set physical `dt=delta_depth/speed`; the T10 exterior medium shares the transaction. No frame-step speed is accepted. | `4512278e`, `7bc608ce`, `f7550895` |
| 6. Seed scatter | DONE aggregation guard | Enumeration and aggregation require at least three distinct seeds per speed; aggregation uses between-seed sample SD, never within-run SEM. | `4512278e`, `f7550895` |
| 7. C-2 accepted equilibrium | DONE — three-seed native ensemble | The driver regression to membrane subdivision 6 was corrected to the already-established grid-converged subdivision 8. Full-native seeds 0/1/2 accepted one committed physical step each on RTX 3090 at maxPF 0.209269/0.205555/0.208651 pN versus the unchanged ~0.210669 pN predicates. The earlier subdivision-6 plateaus remain rejected controls. | `4d30ed84`; `c2_accepted_ensemble.json` + raw jobs 51/53/54 |

## Experiment split

- **Spherical indenter:** engine perturbation harness. It can run only after C-2 and sourced spherical
  apparatus/contact inputs are supplied.
- **Chugh 2017 reproduction:** a separate tipless flat-cantilever-against-dish apparatus using the published
  geometry inversion. The spherical harness must not be relabelled as this experiment.
- **Result class:** the connected path currently permits passive apparatus mechanism-demo only. Quantitative
  execution additionally requires a post-induction, equilibrated, `dt`-converged active-cortex state handoff;
  direct target-line `Pi_0` is necessary but not sufficient.

## Sweep release condition

The committed three-seed C-2 ensemble satisfies C-2. The host-side
`ac_afm_sweep_preflight.py` validates that evidence and enumerates points only when
`AFMSweepDeclaration.blockers()` is empty for the intended result class. Π₀, spherical apparatus/contact
calibration, exterior viscosity and the accepted active-cortex state handoff remain open, so only numerical
C-2 diagnostics and connector gates may currently be visualized; no quantitative material curve is emitted.
