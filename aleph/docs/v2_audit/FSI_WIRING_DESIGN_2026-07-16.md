# FSI two-way wiring — design note (2026-07-16)

**PI runtime contract 2026-07-16:** implement and execute this entire path in Warp on CUDA GPU. No HOOMD or
other simulation runtime is imported/executed for parity, development, fallback, or validation. Python may
orchestrate and inspect outer-step results, but field/state mutation and inner iterations remain device-resident.

**Context.** The existing `aleph/laws/biot_fluid_warp.py` de-risks device-resident field arrays and Peskin-style
spread/interpolation, but it is **not a production cell-domain Biot solver**. Its fixed-cube held-face update
does not implement a conservative moving membrane-minus-nucleus domain and can create fluid content at a
nominal no-flux boundary. Audit reproducer: with `n=9`, one boundary plane initialized to `p=1`, and
`dt=0.1/6`, the nominal no-flux step changed `Σp` from `81.0` to `81.8166667`; moreover, `pore_mass()` is
`∫p dV`, not the general Biot fluid content. The remaining work is therefore a conservative solver rebuild plus wiring, not a
one-line hook into `network_warp.py`. This note fixes that architecture before touching the hot loop.

## The physics the wiring must get right

**The 65.9 Pa·s mapping is retired, not an acceptance target.** Cytoplasmic poroelasticity describes
rate-dependent drainage through a deforming skeleton; it is not equivalent to assigning a single Newtonian
viscosity to every node. The prior `6πηR` mapping was a 0-D clock surrogate, and the audit's dimensional
back-calculation did not reproduce 65.9 Pa·s. The production FSI gate therefore targets storage/flux balance,
poroelastic relaxation, dissipation, and drained↔undrained rate dependence — not an "emergent viscosity."

The current engine represents time and mean pressure through two lumped constructs:

- `gamma_node = 6πηR/Nc` (η=65.9) — used ONLY as the numerical-step→physical-seconds time-scale
  `_dt_real = dt_mu · gamma_node` (`network_warp.py:695,972`). It is NOT a per-node dissipative force in the
  equilibrium loop; it calibrates how many overdamped descent steps equal one physical second.
- the 0-D biphasic block (`network_warp.py:811–825`): `dP_osm` van't Hoff turgor draining on `τ_osm` via
  `V0_eff`, + `dP_solid` Terzaghi effective stress. This is where the current rate-dependence (P4.2) lives.

The new field model separates `p=p_ext+p_bar+p_excess`: `p_bar` is the osmotic/hydrostatic mean set by
membrane water flux and volume closure, while `p_excess` is the spatial poroelastic field. The solid balance
contains `σ_total=σ_eff−αpI`; the fluid-content balance is
`S p_dot+α∇·v_s+∇·q=s_water`, with
`q=φ(v_f−v_s)=−(k/μ)(∇p−ρb)`. Myosin belongs in the solid force balance and must not be added again as a direct
Darcy force. Retire the old `6πηR` clock narrative when this physical scheduler lands, while retaining only a
non-double-counted solid mobility for the chosen inner solve.

## Why this is a time-integration restructure, not a one-line hook

The engine is a **mechanical-equilibrium solver** (ENGINE.md, re-validated T0): `pos += dt_mu·f`,
`dt_mu = 0.1/kmax` is a dimensionless descent step; physical time is *assigned* by `_dt_real = dt_mu·γ`.
`dt_mu` and `_dt_real` are LOCKED through γ. Consequence:

- Retiring 6πηR (η=65.9) → γ_solid (η_solvent≈1e-3, `units.fiber_point_drag_solvent`) shrinks `_dt_real`
  ~65 900×. For a fixed press (physical `t_ramp = R·strain/v_press`, drag-independent), the step count
  `_N_ramp = t_ramp/_dt_real` EXPLODES ~65 900× — computationally infeasible — while the physics is
  unchanged (drainage is driven by physical time, not step count). So a naïve time-scale swap is wrong.

The lock must be broken. Correct architecture = **one physical clock with a nested mechanical solve**:

```
OUTER physical-time step n → n+1       dt_phys obeys the registered field/KMC accuracy limits
  update prescribed loads and membrane water-flux data at physical time t_n
  solve conservative Biot q-p/storage step on the live domain
  reconstruct v_f = v_s + q/φ and advance conservative RAD/KMC with the chosen split
  INNER mechanical solve at frozen outer state
    repeat: assemble solid residual from bending, links, motors, membrane, nucleus, MT, −αpI coupling
            take Warp-GPU explicit descent OR the Warp-GPU implicit step
    until the registered residual/conservation criterion passes
  reject/reduce dt_phys or fail loudly if the inner solve does not converge
```

The exact operator order (staggered, predictor-corrector, or monolithic) is selected by a convergence study,
but every substep advances the same outer physical time. Inner iterations are not seconds, and a fixed inner
budget may not silently substitute for convergence. This is what lets the 0-D `6πηR` clock retire.

## The wiring (production ON; OFF is a regression oracle only)

An additive `biot_fsi=None` path may remain temporarily to prove parity with the frozen Warp-FF arithmetic
path, but it is not
an acceptable production baseline. The new engine's physiological config turns the fluid layer ON from I1.
The ON path contains:

- **field grid/domain:** a conservative finite-volume/cut-cell grid on the live membrane mesh minus the live
  nuclear-envelope mesh. Grid resolution is an I0-B1 accuracy/memory decision, not inherited from P7.F.7.
- **solid→fluid source:** discretize `α∇·v_s` and domain remap conservatively; do not infer volumetric strain
  from cortex radial velocity alone once interior filaments and compartments are live.
- **boundary flux:** outer membrane water flux supplies `s_water`/face flux; the moving nucleus has relative
  no-flux. Impermeable boundaries conserve integrated fluid content.
- **drainage/flow:** solve the storage balance and face flux `q`; reconstruct absolute `v_f` for solute
  advection. The legacy held-face `BiotField.step` is an oracle/parity scaffold, not reused as this solver.
- **fluid→solid transfer:** derive pressure stress/traction from `−αpI` and use an adjoint grid↔node operator;
  gate work sign and closed-system net-force projection. Do not call pressure transfer a solved fluid-momentum
  equation.
- **solid mobility:** use only the mobility required by the I0-A-authorized inner solver; do not convert it
  back into an effective bulk viscosity or a second pore-drag channel.
- **double-count guard:** replace the 0-D transient `p_excess` channel when the spatial field lands; retain
  `p_bar=Π₀` and drained skeleton elasticity only as distinct, sourced equilibrium channels.

## Validation gate (ratified native GPU baseline)

1. **Regression oracle:** FSI-OFF matches the frozen Warp-FF arithmetic path, bit-identical where the same
   operations are retained. No HOOMD/non-Warp executable participates, and this does not authorize an OFF
   production run.
2. **Conservation/numerics before native:** constant-state preservation; Terzaghi/Green oracle; manufactured
   moving-boundary solution; impermeable fluid-content conservation; integrated content change equals outer
   membrane flux; pressure-work sign; adjoint transfer; refinement/convergence study.
3. **Flow kinematics:** Darcy slab flux is linear in `Δp` and `k/μ`; `v_f−v_s=q/φ`; a tracer follows `v_f`,
   not `q`. A Poiseuille profile is required only for a separately authorized Brinkman/free-fluid region.
4. **FSI-ON native reproduces:** drained↔undrained rate dependence and spatial `τ_p` relaxation, with all
   magnitudes reported against registered evidence without retuning `k`, `M`, `μ`, or `α`.
5. **Closed-system diagnostic:** internal pressure/drag coupling does not create cell COM translation; any
   I1b2 relative-velocity drag is equal-and-opposite and dissipative at the discrete level.
6. **Figure:** the CFD-ON native cell HTML (spatial pore-pressure/flux field visible + checkable), per the
   viz-at-closeout rule.
7. **Device residency:** profiler records zero authoritative GPU→CPU state roundtrips inside the outer
   physical-time loop; explicit/implicit mechanics, Biot, RAD, KMC, IBM, and neighbor queries are Warp-CUDA.

## Bounded scope / honest deferrals

- I1a/I1b solve conservative storage plus `q-p` and reconstruct `v_f`. The optional explicit
  relative-velocity solid/fluid drag pair (master I1b2) is a separate PI-KU-gated channel requiring a nodal
  representative volume, equal-and-opposite transfer, and a proof that it does not duplicate pressure work.
- Membrane hydraulic-flux wiring, the deformable nuclear inclusion, cell-type IF replacement, and dynamic
  fiber remodeling are owned by their master increments. Existing membrane Helfrich mechanics are reused;
  FSI wiring completes the **fluid axis**, not "all cell coding."

## Change log
- 2026-07-16: created. Fixes the FSI wiring architecture (inner mechanical equilibrium / outer physical time,
  retirement of the 0-D `6πηR` clock, and double-count guard) before touching the hot loop.
- 2026-07-16 implementation-readiness audit: reclassified the existing field as a component scaffold;
  replaced held-face fixed-box diffusion with a conservative moving-domain contract; corrected `q` vs `v_f`;
  removed the broken emergent-65.9 gate; added storage/flux/work/net-force validation.
