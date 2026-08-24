# AFM cortical-tension readiness — 2026-08-10

The numerical C-2 prerequisite is now **accepted**, but this directory still contains no AFM material result.
Three independent full-native suspended-cell runs on the shared RTX 3090 accepted one physical step each at
the unchanged projected-force and displacement predicates.  All compartments were present: 70,686 actin
filaments, 494,802 actin nodes, 655,362 membrane vertices and 1,150,806 total nodes.  With membrane subdivision
8 and the production explicit solver, seeds 0/1/2 accepted after 280/300/360 reported inner iterations at
maxPF 0.209269/0.205555/0.208651 pN, respectively, against per-run thresholds near 0.210669 pN.  Each accepted
transaction advanced biological time by 0.01 s.

The decisive change was restoring the membrane subdivision-8 fidelity already established by
`GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §4.  The AFM driver had regressed to subdivision 6, whose membrane
discretization leaves a much larger residual.  The earlier subdivision-6 plateau and its solver variants
remain preserved as rejected controls; they are not evidence of material response and must not be pooled with
the accepted ensemble.  The production default is now 8, and structural guards prevent subdivision 6 from
licensing C-2.

`c2_accepted_ensemble.json` is the compact three-seed evidence record and
`c2_accepted_ensemble.png` is its numerical acceptance visualization.  The raw records and Slurm logs are in
`raw/`.  The older `c2_solver_plateau.png` and `c2_objective_comparison.png` remain historical rejected-state
diagnostics.  None of these figures is an F–δ curve or cortical-tension magnitude.

The spherical-indenter runtime separately passed its CUDA unit gate, then the exact-trial force hook passed on
RTX 3090 (Slurm job 56): separated/touching configurations are force-free, overlap is repulsive,
`reaction[0]` reports the opposite vertical force, only the global membrane block is modified, the
membrane–indenter balance is adjoint to machine precision, and rejected/accepted apparatus positions roll
back/commit correctly. That gate validates connector wiring only; radius, gap and stiffness remain deliberately
unset PI inputs.

Commit `f7550895` connects the apparatus and the T10 exterior-Stokes medium to every exact inner-solver trial
and to the same accepted-step transaction. The loading driver requires the exact C-2 full-native configuration
and a verified build commit, writes no partial curve after a rejected candidate, and reads reaction only between
accepted steps. `ac_afm_sweep_aggregate.py` requires at least three accepted seed paths per speed and uses sample
SD across seeds; `ac_afm_sweep_vis.py` renders only accepted ensembles in the established static style.

`ac_afm_sweep_preflight.py` revalidates the committed C-2 ensemble and writes a fail-closed sweep manifest.
It has no physical defaults: apparatus radius/source, exterior viscosity, speed/depth/seed axes, Π₀ evidence,
and contact distance/stiffness provenance must be declared explicitly. Missing Π₀ or contact calibration
produces `BLOCKED`, an empty point list and exit code 2. A proxy Π₀ can enumerate only when the caller requests
the predeclared mechanism-demo-not-quantitative class. Quantitative enumeration additionally remains blocked
until a post-induction, equilibrated and `dt`-converged active-cortex state handoff closes STATE (c) 3/17; forged
READY manifests and quantitative path records are refused downstream as well.

`protocol_selection.md` records the source audit for Chugh et al. 2017.  Chugh is retained as the rounded-HeLa
cortical-tension magnitude/state oracle, but its apparatus is a tipless flat cantilever against a dish, not the
spherical indenter implemented by this harness.  The spherical harness must not be labelled a Chugh
reproduction, and that paper does not close Π₀, exterior viscosity or numerical-contact calibration.

The remaining release inputs are therefore: (1) a spherical apparatus protocol and radius or authorization
to implement the distinct flat-cantilever Chugh apparatus, (2) target-line Π₀ or an explicit HeLa-proxy
mechanism-demo classification, (3) sourced contact distance and stiffness, (4) exterior-medium viscosity, and
(5) an accepted active-cortex state preparation past induction with a declared time-convergence comparison.
Until they are supplied and connected, the engine must not launch a quantitative force–indentation sweep or
generate a cortical-tension curve. The predeclared protocol uses round suspended geometry, the exterior-Stokes
six-rigid-mode holder, physical speed/depth axes and at least three distinct seeds per point.
