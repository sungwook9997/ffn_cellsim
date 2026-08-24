# ac/fluid — INTEGRATION patch-notes (lead applies in spine order; the fluid track does NOT edit `ff/`)

Per `AC_PARALLEL_SESSIONS_2026-07-16.md` §1.3: where the fluid-spine increments must retire/delete/extend
an `ff/` mechanism, the fluid track does **not** touch `ff/` in its worktree. Each change below is a
documented **patch-note + adapter contract** the lead applies during the serial native-integration step, in
spine order (I1a → I1b → I1c). Every note gives: **file · symbol · before → after · double-count guard ·
same-commit constraint**. Line numbers are from the 2026-07-16 `ff/` state (Explore audit `wf` run); the
lead re-anchors by symbol if they drift.

The new-engine runtime is `ac/fluid/*` (Warp-CUDA source, gated CPU-green via the NumPy oracles). These
notes tell the lead how the old `ff/` mechanisms retire as the `ac/` substrate lands — they are the
"same-commit" guards that stop the *additive-default-off-becomes-production* trap.

---

## PN-1 (I1a) — retire the `6πηR → γ_solid` overdamped CLOCK

- **File · symbol:** `ff/network_warp.py:886` `gamma_node = 6.0*np.pi*eta_bulk_Pa_s*R0/max(Nc,1)`;
  `:887` `_mu_phys = 1/gamma_node`; `:891` `_dt_real = dt_mu/_mu_phys`; secondary OFF-path clock
  `:1426-1432` `_dt_phys = dt_mu*(6πη R0/Nc)`. Retirement target already present:
  `ff/units.py:67-78` `fiber_point_drag_solvent()` (γ_solid at `ETA_SOLVENT=1e-3`, `units.py:64`).
- **Before:** the engine is a dimensionless overdamped descent (`x += dt_mu·F`, `axpy_kernel:325`); physical
  seconds are *assigned* by multiplying the descent step by a lumped whole-cell Stokes drag `6πηR`
  (η = 65.9 Pa·s composite) distributed over `Nc` nodes. `dt_mu` and `_dt_real` are **LOCKED through γ**.
- **After:** physical seconds come from the **outer physical clock** — `ac/fluid/scheduler.py`
  `PhysicalScheduler.outer_step(dt_phys)`, `dt_phys = t_ramp/n_phys` set by the physical process (membrane
  water flux / load ramp), independent of any drag. `γ_solid` (`units.py` `ETA_SOLVENT`) survives ONLY as
  the **inner-solve mobility** of the frozen-state mechanical relaxation, never as the clock. Inner
  iterations converge a registered residual and are never reinterpreted as seconds (I0-A).
- **Double-count guard:** the effective cytoplasm viscosity EMERGES from the Biot storage+mobility (D, M);
  do **not** also scale time by `6πηR` — that double-counts the fluid drag (once in Biot, once in the
  clock). The in-code metric already anticipates this: `:1299` `"6πηR retired; effective η emerges from
  Biot D+M"`.
- **Same-commit:** retire the `6πηR` clock in the SAME commit the conservative Biot substrate (PN-3) lands
  (`ENGINE_ARCHITECTURE_PLAN` hazard: "retire 6πηR the moment ① lands, same commit"). A naïve swap
  `6πηR → γ_solid` as the clock shrinks `_dt_real` ~65,900× and explodes the step count — the fix is the
  outer/inner split, not a drag substitution.

## PN-2 (I1a) — demote the SCALAR TURGOR "balloon" to the mean-pressure channel `p_bar`

- **File · symbol:** `ff/network_warp.py:226-238` `turgor_kernel` (uniform `dP_area` pushes every node
  radially); `:1011-1087` `_refresh_turgor()` closure, esp. `:1064-1071`
  `dP = dP_osm + dP_solid; dP_area = dP*area/Nc` (**one** scalar ΔP from a hull-volume ratio, applied
  identically to all `Nc` nodes).
- **Before:** the cytosol is a 0-D balloon — a single ΔP from `V_cyto = ConvexHull − V_nuc` pushes the shell
  uniformly (`turgor_kernel`), with a 0-D biphasic drainage `τ_osm=(V0_eff−vmin)/(Lp·A·Π_in)` (`:1067`).
- **After:** the pressure splits `p = p_ext + p_bar + p_excess` (`biot_substrate` / `field_grid.p_bar`):
  - `p_bar` = the **mean/osmotic turgor**, seeded at resting `Π₀ = 40 Pa` (I0-B1, `params_i0b1.yaml`), a
    uniform pre-stress — legitimate (the shell IS pre-tensioned; `γ = ΔP·R/2` baseline preserved);
  - `p_excess` = the **spatial** poroelastic field from the conservative Biot solve (`biot_pmass_update_kernel`),
    applied per-node as the LOCAL deviation via `PressureCoupling.accumulate` — NOT a uniform push.
- **Double-count guard:** when the spatial field lands, **REPLACE** the 0-D `p_excess`; keep `p_bar = Π₀` +
  drained elasticity as **distinct channels** (do not sum a 0-D balloon ΔP and a spatial excess — that
  double-counts the osmotic response). Do not force `p_excess` to zero mean during an undrained transient
  (that discards Skempton undrained stiffening).
- **Same-commit:** with PN-3.

## PN-3 (I1a) — replace the HELD-FACE FSI loop with the conservative moving-domain solve

- **File · symbol:** `ff/network_warp.py:1114-1311` the FSI-ON block; esp. `:1182-1184` host-side `num/den`
  mass-lumping spread; `:1187-1190` the reused **non-conservative** `bf.step()` (held-face
  `ff/biot_fluid_warp.py::biot_diffusion_kernel`); `:1195-1203` zero-mean-over-cytosol + clip.
- **Before:** per outer step — `_refresh_turgor` → `fsi_compression_source_kernel` (③→① source) → Peskin
  spread with **host `num/den` normalization** (nonlinear, breaks adjointness) → **held-face
  `BiotField.step`** (creates mass at the "no-flux" box: Σp 81.0 → 81.8167, `FSI_WIRING_DESIGN:9-12`) →
  zero-mean clip → interp back (①→③).
- **After:** the conservative moving-domain FV solve — `ac/fluid/domain.Domain` (live membrane-minus-nucleus
  classification + conservative remap) + `ac/fluid/biot_substrate.biot_pmass_update_kernel` (masked
  conservative stencil) + `ac/fluid/boundary` (membrane hydraulic flux `s_water=L_p(σΔΠ−ΔP)` + nucleus
  relative no-flux). Fluid→solid is the **adjoint** `PressureCoupling.accumulate` (transpose-matched Peskin,
  gated by `ibm_reference`), NOT the nonlinear `num/den`.
- **Double-count guard:** run the conservative path OR the held-face path, never both. The adjoint transfer
  passes the work-sign + closed-system net-force gates (`test_i1a_ibm_gates`). Myosin stays on the SOLID
  balance and is **not** re-added as a Darcy body force (`velocity.py` note).
- **Same-commit:** PN-1 + PN-2 + PN-3 land together (the `①` substrate commit).

## PN-4 (I1a) — reuse `biot_fluid_warp.py` device+IBM pieces; REBUILD its held-face stencil

- **File · symbol:** `ff/biot_fluid_warp.py:29-51` `biot_diffusion_kernel` (freezes `p_new[face]=p[face]` →
  Dirichlet reservoir, not a zero-flux wall) + `:65-66` `biot_gradient_kernel` face-zeroing.
- **Before:** the held/clamped boundary faces silently absorb/emit mass — not a conservative Neumann solver.
- **After:** REUSE only (a) the `_p/_p2` **double-buffer** pattern (`:83-102`) — mirrored by
  `field_grid.p/p_new` + the ping-pong swap in `biot_substrate.step`; (b) the **Peskin IBM** (`_peskin4`
  `:116-128`, `ibm_spread_kernel`/`ibm_interp_kernel` `:131-204`) — mirrored verbatim in
  `ac/fluid/biot_substrate._peskin4` / `ibm_reference` and validated by the adjoint gate. REBUILD (do NOT
  port) the diffusion step as `biot_pmass_update_kernel` (masked conservative FV).
- **Double-count guard:** none (pure replacement). The reused IBM must pass the adjoint gate
  (`Interp = h³·Spreadᵀ`); build the transfer against the RAW spread/interp, NOT the `num/den` wrapper.

## PN-5 (I1a) — `biot_fsi=None` additive channel is a REGRESSION ORACLE, not a production baseline

- **File · symbol:** `ff/network_warp.py:680` kwarg `biot_fsi=None` (default-OFF) on
  `simulate_whole_cell_compression_on_device`; `:1114` `if biot_fsi is not None:` gate.
- **Before:** fluid is an OPTIONAL, default-OFF additive path — so every landed native run was passive
  (memory: "default-off라 native run 다 passive였음"). This is the *additive-default-off becomes de-facto
  production* trap.
- **After:** fluid is **ON from I1** at physiological `Π₀` (the `ac/` engine has no FSI-OFF production
  path). The OFF path survives ONLY as the **FSI-OFF == Warp-FF bit-identical regression oracle** (native
  gate NG-1). `simulate_whole_cell_compression_on_device` itself (`:672-681`, the passive-shell assembly)
  is REBUILT as the fluid-first active assembly (`NEW_ENGINE_BUILD_PLAN §2b ❌`) — it is not the `ac/` entry.
- **Double-count guard:** the physiological-baseline HARD rule — default-off is code hygiene, never a license
  to run physics from an unphysical passive zero.

## PN-6 (I1c) — swap the fixed scalar `G_actin` for the transported monomer FIELD

- **File · symbol:** `ff/polymerization_warp.py:53` `G_actin_uM: float` (scalar field of
  `ResolvedPolymerization`); `:62` `v0 = (k_on*G_actin_uM - k_off)*delta_full_um` (baked ONCE at resolve
  time); `:88` `polymerization_kernel` takes the pre-computed scalar `v0` — **no monomer field, no pool,
  no `c_local`**. Rear half: `ff/motility_warp.py:218` `pointed_end_depoly_kernel`, `:254`
  `fiber_treadmill_kernel` (scalar Δlength budgets, not a field).
- **Before:** `G_actin` is a fixed scalar — spatially uniform, an infinite pool; assembly elongates
  `seg_rest` but decrements no reservoir. Growth is NOT flux-limited.
- **After:** barbed growth reads the **LOCAL** concentration `c_local` IBM-interpolated (Peskin) from
  `ac/fluid/transport.MonomerField` at the barbed-end cell: `v+ = δ(k_on·c_local − k_off)`. Every growth
  event **subtracts** `k_on·c·δ` from that cell and **adds** it to the bound-polymer count; every
  pointed-end/sever event does the reverse (feeding the SAME field, `motility_warp` rear half). Total actin
  `∫φc dV + N_polymer` is conserved to machine precision (`transport_reference` gate).
- **Double-count guard:** the monomer removed from the field == the monomer added to the polymer (exact
  exchange, no free lunch); the pointed-end release feeds the same conserved field, not a second scalar
  budget. Pick ONE turnover channel (do not run the scalar `v0` and the field `c_local` simultaneously).
- **Same-commit:** with the I1c transport landing (needs `v_f` from I1b live).

---

## Frozen interface contracts this track SHIPS (consume read-only)

- **`ac/fluid/domain.py :: Domain`** — the live membrane-minus-nucleus domain. **fluid-spine OWNS it.**
  `set_nucleus_boundary(provider)` (§1.4): I1a ships `StaticSphereNucleusMaskProvider`; the nucleus track
  ships a `DeformableNucleusMaskProvider` conforming to the SAME `NucleusMaskProvider.classify_nucleus(grid,
  mask)` signature (its moving oblate mesh is I1a's inner relative-no-flux boundary). Neither co-edits the
  other's file.
- **`ac/fluid/biot_substrate.py :: PressureCoupling.accumulate(state, out_force)`** — the §1.4 inner
  force-assembly primitive (fluid→solid `-αpI` body force). The lead-owned integrator sums it with
  `MyosinForce` (I3, `ac/motor/`) and `StericForce` (I2b, `ac/solid/`). `state` carries `node_pos`
  (wp.array vec3d) + `node_volume` (wp.array float64).
- **`ac/fluid/transport.py :: MonomerField`** — the I1c conserved G-actin field the I4 weave / lamellipodium
  dendritic nucleation consume as the flux-limit (`c_local`), read-only.
