# FF ECM library — development roadmap (2026-07-10)

Designed via a 5-lens adversarial design workflow (physics-fidelity · cell-ECM-coupling · active-remodeling ·
validation/KB-rigor · engine-integration), critiqued for hard-rule compliance, and synthesized. Every item
respects the project hard rules (mechanistic full-fidelity · native-scale validation · **NO tuning to an
outcome** · physiological baseline · interactive-HTML viewer per milestone).

## Vision

Advance the library from a **validated single-point-Pa material catalog** into a **living, constitutively-
complete substrate**: a real full-compartment cell **senses, remodels, and is guided** by native-scale ECM
whose nonlinear / anisotropic / time-dependent moduli all EMERGE and match literature master curves, on ONE
grid-invariant stress readout, behind PI-authored gates. Two thrusts: **(1) make the library REACHABLE by a
cell** and prove quasi-static stiffness-sensing NOW (around the crawl-drag blocker); **(2) make the Pa REAL
across the whole deformation regime** (trustworthy virial stress → nonlinear master curve → topological
remodeling → literature-anchored constitutive additions) — converging on a native cell crawling through and
reshaping a native tissue-mimetic matrix once principled Stokes drag unblocks migration.

## NEAR — autonomous, non-gated (next few sessions)

Route the library into the cell scripts, stand up ONE macroscopic stress tensor, and put a resting cell on
real-Pa ECM — zero tuning.

> **Execution status (2026-07-11): NEAR #1, #3, #4, #5 DONE + native/method-verified. #2 & #6 remaining.**

1. ✅ **DONE — Quasi-static stiffness-sensing demo** — `ff_stiffness_sensing.py`: `--from-resting` full-
   compartment cell (Winkler substrate k_sub=2Ea/(1−ν²), Path b). **Native A5000-confirmed** (Nc=266000):
   traction 0→0.090 nN with the engagement threshold (bound 0.02@150Pa → 0.80@2kPa) = the durotaxis basis.
   Monotonic to 40 kPa (interior peak above range). `figs/stiffness_sensing_native.png`. Path (a) — grip the
   ACTUAL library ECM network — is #2 below.
2. **Reachability / Path (a)** — extend `ff_stiffness_sensing.py` (or route the cell scripts) onto
   `build_ecm()` so the cell grips the LIVE library ECM network (traction emergent from its calibrated seg_k),
   `attach_clutches_to_ecm` onto `ecm.net.pos`; guard alignment-S for continuum. *Next autonomous item.*
3. ✅ **DONE — Full macroscopic virial Cauchy stress tensor** — `ecm_material_stress→σ[3,3]` (commit d47b5ae).
   Method-independence PROVEN: PA G_energy=G_virial=G_react=1522; collagen G_energy=16.1≈G_virial=16.4 while
   G_react=1.4 (reaction is the outlier). Resolves the two-readout inconsistency.
4. ✅ **DONE — 6-material native atlas** — `ff_ecm_native_atlas.py`: 6/6 IN BAND at a 40µm native REV, each
   matching its small-REV value (REV↔native consistency), fibrillar via the virial tensor. `figs/native_atlas.png`.
   *Remaining: per-material full-extent annotated HTML viewers (the gallery covers most at moderate box).*
5. ✅ **DONE — Nonlinear σ_xz/N1/K(γ)** — `shear_stress_curve` via the tensor. N1 SIGN independently confirms
   the stretch-dominated regime (model N1>0 vs literature N1<0, Janmey). `figs/normal_stress_N1.png`.
6. **Systematized static traction-remodeling + contact-guidance anisotropy** on validated collagen (extends S6).
   *Next autonomous item (needs Path a / cell-ECM coupling).*

## MID — moderate build or one PI decision

- **ECMRemodeler topological bond-mutation operator** → emergent plasticity + MMP + LOX channels (Kim p07/p39:
  transient Bell crosslinks, covalent = k_off→0 limit, creep-recovery). *Blocked on the matrix Bell-rate datum.*
- **Thermal semiflexible force-extension** (WLC entropic) → real differential-modulus **K~σ^{3/2}** master curve.
  *Changing the production default is PI-gated.*
- **Principled per-segment cylindrical Stokes drag (Kim T6)** — the **migration unblocker** (fixes native crawl
  Σγ∝Nc solver instability). *Blocked on solver stability work.*
- **Long-range stress propagation |σ|(r)~r^-n** via a contractile inclusion (KB-1.10) — needs #3.
- **Kim BD acceptance-oracle cross-validation** for collagen (G, K(γ), N1) — oracle-only, per the rules.
- **Formal ECM ValidationGate family** + disk-grounded checker *(PI-authored gates)*.
- **ν-faithful continuum gel** (break the Cauchy relation with a bending/multibody term) *(PI-gated)*.
- **⟨z⟩(c) LOX crosslink-density law → predicted c²** — literature-anchored, **NEVER an n-fit** *(PI DECISION)*.

## FAR — the living, remodeling, tissue-mimetic matrix (vision, PI-gated)

Durotaxis + stiffness-dependent migration on the gradient/PA substrates · live cell crawling THROUGH the
remodeling matrix (contact-guidance + proteolytic invasion) · 3D confined invasion into collagen · two-phase
poroelasticity (Matrigel/agarose/tissue) · MMP/LOX as secreted reaction-diffusion enzyme fields · tissue-mimetic
ECM atlas (brain/muscle/tumor-stroma/wound) from literature composites · bidirectional matrix-stiffness ↔
cortical-tension feedback + cross-engine DCM⊗FF (DCM cell dividing against the FF library's G'=10/100/400 Pa).

## Sequencing

NEAR #1 (stiffness-sensing) executes now — it is the library's purpose and de-risks the whole migration
program (the same E→k_sub clutch later drives durotaxis). #2/#3/#4 are parallel autonomous enablers. #3 (virial
tensor) unblocks #5 and the MID r^-n. The MID ⟨z⟩(c), ν-faithful, gate-family, and thermal-default items are
the PI decisions; the FAR program follows the Stokes-drag migration unblocker.
