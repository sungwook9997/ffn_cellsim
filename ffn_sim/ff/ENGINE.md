# ff — Filament-FEM engine (build spec)

*Stage-6 of the 2026-06-29 two-layer restructure. The going-forward fine-grained engine:
a discrete cross-linked filament network solved as a **mechanical (FEM-like) system** —
energy assembly + implicit nonlinear solve + autodiff — instead of explicit thermostatted MD.*

Physics reference: **Cytosim** (Nédélec & Foethke 2007, *New J. Phys.* **9**:427) — exactly as
the `dcm/` layer references SimuCell3D. Cytosim itself (runnable external C++) is the independent
**parity oracle**, so the active layer is never validated against only itself.

---

## 1. Why MD-free is faithful here (the design rationale, PI-discussed 2026-06-29)

A crosslinked filament network is structurally a **beam/truss finite-element network** (crosslinks =
nodes, filament segments = elastic elements, bending = rotational stiffness). The question is whether
to integrate it as Langevin MD (current `archive/hoomd_legacy`) or solve it as mechanical equilibrium.
Decision rule — **does kT enter the observable through the mean, or through the variance / rare tail?**

| regime | mechanism | treatment | cost |
|---|---|---|---|
| **A — integrable-out** | entropic (WLC) elasticity | closed-form constitutive force (Marko–Siggia); the thermal effect becomes a deterministic term | **free** |
| **B — kinetic layer** | catch/slip bonds, myosin stepping (force-dependent, ATP-driven) | event-driven KMC on the mechanical solve (not per-step thermostat) | cheap; carries its own kinetic timescale |
| **C — must sample** | *pure-thermal* rare-event configurational remodeling | cannot be integrated out; keep stochastic | the only thing that needs sampling |

Key point for an **active** cell: most remodeling is ATP-rate-driven (regime B, clean kinetic layer),
not pure-thermal exploration (regime C). So the MD-free mechanical solve is *more* applicable here than
for a passive thermal gel. Quenched/structural disorder (initial network architecture, bound states)
moves from the time-axis to a **parallel realization ensemble** (Monte-Carlo at init) — exact for slow
disorder, and embarrassingly parallel. See the γ-floor protocol (§4).

> ⚠️ NOT lumping: regimes A and B keep every mechanistic force law (WLC, Bell-Evans, Hill) — only the
> *integration method* changes. This satisfies the project's mechanistic-over-lumped hard rule.

---

## 2. Cytosim physics to implement (Nédélec & Foethke 2007) — LITERATURE-GROUNDING REQUIRED

⚠️ The detailed formulae below must be grounded against the paper before coding (no invented
constants — project hard rule). The paper PDF is **not yet in `ffn_sim/references/`** → fetch +
register a SourceEvidence row first (see §6).

1. **Fiber model.** A fiber = a chain of `n` model points joined by `n-1` segments. Bending elasticity
   from rigidity `κ` (κ = E·I = k_B·T·L_p): `E_bend = (κ/2) ∫ (∂²r/∂s²)² ds`, discretized over
   consecutive-node triples. Optional axial/extensional stiffness per segment (near-inextensible).
   → **literature-ground:** the exact discrete bending operator + κ ↔ L_p mapping.
2. **Integration — the Cytosim move.** Overdamped Langevin `γ ẋ = F(x) + ξ` solved with an **implicit**
   (semi-implicit) scheme so `dt` is bounded by *accuracy*, not the bending CFL → large steps. This is
   exactly the matrix-free linearly-implicit step we already have (§3).
   → **literature-ground:** Cytosim's specific preconditioning / Brownian-term handling.
3. **Hand binding model.** Motors and crosslinkers are **Hands** that attach/detach fiber sites
   stochastically at force-dependent rates (Bell-Evans / catch-slip; Hill for motor stepping). This is
   the regime-B kinetic layer.
   → **literature-ground:** Hand attach/detach rate forms + the motor force-velocity.

---

## 3. Reuse map — ~80% already built in `dcm/` (verified 2026-06-29)

| FF piece | reuse | adaptation |
|---|---|---|
| element force assembly | [`dcm/network_warp.py`](network_warp.py) — harmonic bond + harmonic angle + LJ, one-thread-per-element, atomic accumulate | swap harmonic-bond → **WLC (Marko–Siggia)** constitutive; harmonic-angle → **Cytosim bending** |
| **implicit "FEM solver"** | [`dcm/dcm_warp_implicit.py`](dcm_warp_implicit.py) — `implicit_overdamped_step(pos, force_fn, γ, dt, …)`, matrix-free CG + Newton, **force-agnostic**, precond opt-in | quasi-static γ: `dt→∞` ⇒ Newton on `F(x)=0`; pass the FF force as `force_fn` |
| autodiff (inverse) | [`dcm/dcm_warp_diff.py`](dcm_warp_diff.py) — `wp.Tape` over the whole loop, FD-validated 4e-12 | back-prop γ-loss → fiber κ, prestress |
| geometry | `dcm/geometry.py`, `common/surface_manifold.py` | cortex shell seeding |
| kinetic layer template | `archive/hoomd_legacy` cadherin / junction-switch hosts (force-dependent bond update) | Hand KMC on the mechanical solve |

New code needed: WLC + Cytosim-bending kernels, the Hand KMC layer, the quenched-ensemble harness,
the fiber-network data model ([`fiber_network.py`](fiber_network.py), scaffolded).

---

## 4. First milestone — γ-floor prototype (the decisive experiment)

γ (cortical tension) is **quasi-static** → the load barely changes during measurement → even the
load-coupled myosin has a well-defined steady distribution that can be sampled once (quenched). Protocol:

1. MC-sample N realizations of the cortex fiber network + myosin bound-state at init (regime A/B disorder).
2. Each realization: WLC-beam network + myosin active prestress → **implicit mechanical equilibrium**
   (`implicit_overdamped_step`, `dt→∞`), no explicit BAOAB.
3. Run the N realizations **in parallel** → distribution of γ.
4. **Decisive validation:** compare MD-free γ vs the BAOAB-MD γ from the archived
   `cortical_tension.py` (the γ-floor finding). "Does an MD-free mechanical solve reproduce the same
   γ?" → confirms FF is MD-equivalent for γ, *or* localizes exactly where in-time load-coupling is
   load-bearing.

Ties to the open **γ anchoring** question (Moazzeni MCF7 1e-2 N/m vs SimuCell3D 1e-3 vs emergent;
no clean adherent-MCF7 datum) — γ should be **swept as a controlled variable**, faceting/tension
reported as a function of it (NOT tuned to a target). See `docs/v2_audit/PARAM_AUDIT_SIMUCELL3D_2026-06-25.md`.

---

## 5. Build stages

- **6a** ✅ this spec + fiber-network data model scaffold.
- **6b** WLC + Cytosim-bending Warp kernels (grounded) + wire to `implicit_overdamped_step`; relaxation
  smoke (a bent fiber relaxes to straight at the right rate).
- **6c** Hand KMC layer (myosin/crosslinker) + quenched-ensemble harness.
- **6d** γ-floor prototype run + the MD-free-vs-MD-γ validation + Cytosim parity oracle.

## 6. Open / PI-gated before deep build

- **Cytosim paper grounding:** Nédélec & Foethke 2007 is now in `references/`
  (`nedelec_foethke_2007_cytosim_njp.pdf`, arXiv mirror 0903.5178) — ground the §2 formulae
  (discrete bending operator + implicit scheme + Hand model) from it before any kernel coding
  (hard rule). SE candidate written (`references/SE_REGISTRATION_CANDIDATES_2026-06-29.md`).
- **γ target band:** unresolved (Moazzeni vs SimuCell3D vs emergent) — sweep, don't tune.
- **Active-driver force anchor:** still REFUTED/HALTED to PI (the SF/NMII missing motor-density datum) —
  the spreading magnitude driver remains PI-gated, not auto-tunable.
