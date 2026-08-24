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

# FSI two-way coupling to EVERY compartment — verified plan (2026-07-16)

**PI intent:** the reason CFD was specified is that the intracellular fluid did NOT interact with the
existing compartments at all, and it must be done PROPERLY. The current FSI (`biot_fsi`) couples the pore
pressure ONLY to the cortex shell (`fsi_pressure_force_kernel` dim=Nc) → the nucleus + MT are frozen
(INTERNAL_DISPLACEMENT_DIAGNOSIS: RAW |Δ| = 0.0 nm). This plan makes the ONE scalar Biot field a receiver
for every compartment, each with a **dimensionality-matched** back-reaction. Designed + adversarially
verified by a 17-agent workflow (wf_9e422e1c); this note is the reconciled, implementable result.

## The coherent scheme (one field, one source, dimensionality-matched receivers)

ONE scalar Biot field `p(x,t)`. ONE source (cortex compression, dim=Nc, unchanged). Every compartment is
a **receiver** — exactly ONE coupling each (two would double-count):

| Compartment | Dim | Back-reaction | Does the fluid move it? |
|---|---|---|---|
| **Cortex** | 2-D shell (fluid inside) | `+p·(A/Nc)·r̂` **outward** about `centre` (existing, unchanged) | yes (already) |
| **Nucleus** | 2-D shell (fluid outside = inclusion) | `−p·(A_nuc/n_nuc)·r̂` **INWARD** about `centre_nuc` (NEW kernel) | **YES — the real interior-response channel** (oblate flatten, poles-in, V_nuc-conserved) |
| **Microtubule** | 1-D fiber (no enclosed volume) | `−V_seg·∇p` body force (NEW kernel + grid-gradient) | **NO — ~1e-3 pN, negligible** (fiber too slender); wired for dimensional completeness only. MT response = strut (fix C, done) + IF/LINC (B/D) |

**Key simplification (keeps the parity-gated IBM kernels byte-identical):** interpolate the scalar `p` and
the vector `∇p` over **dim=N** (all nodes get their local `p`, `∇p`); each back-reaction kernel indexes its
own slice with an offset. Only the launch `dim` changes (Nc→N); no IBM-kernel signature change.

**Nucleus is a PASSIVE receiver for the first landing** (no forward source): an incompressible nucleus's
true net volumetric source ≈ 0, so a forward source adds sign/self-force/double-count traps for ~0 effect.
It still responds — driven by the **cortex-sourced field's angular (quadrupole) variation across its
surface** (poles see +p → pushed in, equator −p → pushed out ⇒ flatten at constant V_nuc).

## Bugs the verifiers caught in the code (fix first)
- **Div-by-zero:** gating the nucleus centroid/area block on `(_nuc_incompressible or _nuc_fsi)` would run
  `dP_nuc = K_vol·(V0_nuc−V_nuc)/V0_nuc` with `V0_nuc = 0` when the nucleus isn't incompressible → 0/0.
  **Split:** compute `centre_nuc / area_nuc / nuc_faces_d` under the OR; keep `dP_nuc/dP_nuc_area` under
  `_nuc_incompressible` ONLY.
- **Init:** `area_nuc` is a throwaway local; seed a persistent `area_nuc = 4π·R_nuc²` + add to the
  `nonlocal` list; seed `centre_nuc` to the nucleus build centroid (not origin) so the first read isn't stale.
- **Setup:** build `_nuc_faces / nuc_faces_d` + the rehull under `(_nuc_incompressible or _nuc_fsi)`.

## Field-domain fixes (conditional on n_nuc>0 ⇒ cortex-only FSI stays bit-identical)
- **Nucleus-interior no-flux:** mask grid cells inside `R_nuc` (about live `centre_nuc`); use a masked Biot
  diffusion (`step_masked`) so pore pressure does not diffuse through the nucleus (it is a HOLE in the
  cytoplasm) → the nucleus surface samples the genuine cytoplasm pressure. Original kernel untouched (parity).
- **Zero-mean over the cytosol annulus:** subtract the mean over `R_nuc < r < R0` (exclude the nucleus
  interior), re-centred each physical step, not the fixed ball that engulfs the nucleus.
- **Clock:** set `_dt_refresh = dt_phys` before the FSI loop so drainage integrates on the physical clock.
  Do NOT force `_drain_on=False` (REFUTED: the zero-mean field carries the *deviatoric* transient; the
  *uniform* undrained→drained transient rides the 0-D `V0_eff` channel — orthogonal, no double-count).

## Momentum (honest — the PI's constraint met in the weak sense the scalar field allows)
IBM spread/interp adjoint consistency + zero-mean + optional **net-force projection** (subtract `F_net/N`
over immersed nodes → COM drift provably 0). This is NOT strict Newtonian conservation — the scalar Biot
field is massless; source (scalar Δp) and back-reaction (traction) are a *thermodynamic* pair, not a
Newton-third-law pair. **Strict conservation ⇒ the u-p velocity field (roadmap, PI-authored KU).** Do not
overclaim. Adding the interior −∇p forces actually IMPROVES the net-force balance (they are the
divergence-theorem partners the cortex-only path lacked).

## Native validation gate (CLAUDE.md HARD)
`ff_internal_displacement_diag.py`, native full cell, incompressible nucleus + membrane + MT, AFM 0.10.
(a) FSI-OFF == current native (bit-identical); FSI-ON+no-nucleus+no-MT == current cortex-only FSI.
(b) **Nucleus unfreezes THROUGH the fluid:** RAW |Δ| 0.0 → >0, AND oblate flatten (z-half shrinks + equator
bulge) at ~constant V_nuc, AND centroid drift ≈ 0 (a flattening quadrupole, NOT a translation), poles move IN.
(c) MT: RAW |Δ| technically >0 but ~1e-3 nm (still ~frozen) — report as the flagged u-p deficit; do NOT
inflate V_seg (no-param-tuning). (d) net-force ≈ 0 logged. (f) CFD-cell HTML: interior field + flattened nucleus.

## Surface to PI (findings, not knobs)
- `M_fsi` defaults to `K_drained` — conflates the Biot coupling modulus with the drained skeleton stiffness
  (already in `dP_solid`). It pins the whole interior-response magnitude → give it a KB-anchored value; if
  emergent nucleus strain misses the band, surface (do not tune).
- **The MT fluid coupling is negligible; the MT's real fluid response needs the u-p velocity drag
  `f=(η_f/k)(v_s−v_f)` — NO KB claim yet, PI-authored KU required.** (Dominant MT couplings are SOLID:
  strut-tip↔cortex + IF/LINC — out of FSI scope.)
- Design-note correction to `FSI_WIRING_DESIGN_2026-07-16.md`: the field carries the DEVIATORIC transient
  only; the uniform transient stays in the (re-clocked) 0-D channel — "0-D transient REPLACED by the field"
  was subtly wrong.

## Change log
- 2026-07-16: created from the verified 17-agent design+adversarial-verify workflow (wf_9e422e1c).
