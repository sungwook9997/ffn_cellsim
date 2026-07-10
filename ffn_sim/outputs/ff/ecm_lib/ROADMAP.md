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

1. ⭐ **Quasi-static stiffness-sensing demo (TOP PICK, executing now)** — `ff_stiffness_sensing.py`: a
   `--from-resting` full-compartment cell (cortex+turgor+membrane+nucleus) coupled by basal FA clutches to a
   LIVE library ECM (PA 150 Pa / PA 40 kPa / collagen). Read engaged-clutch **traction + bound-fraction**
   (PRIMARY) and top-down xy-silhouette **spread** (SECONDARY) vs substrate E → emergent biphasic optimum.
   Fresh path (leaves the validated crawl loop byte-identical); coupling de-risked in `ff_ecm_remodel_demo.py`.
   E = literature input, k_sub = DERIVED Winkler/Boussinesq (a=0.05µm), optimum MEASURED. Native A5000 reconfirm.
2. **Reachability wiring** — route `--ecm-material/-alignment-S/-conc` onto `build_ecm()` in the cell scripts;
   `attach_clutches_to_ecm` onto the library net; guard alignment-S for continuum. OFF-path byte-identical.
3. **Full macroscopic virial Cauchy stress tensor** — `ecm_material_stress(ecm,pos)→σ[3,3]` (bond+crosslink+
   bending virial, full tensor). Proves method-independence: virial==energy on dense PA (E=5000), virial vs
   reaction on sparse collagen (energy-route bias quantified). The single grid-invariant stress downstream needs.
4. **Native production for all 6 materials + per-material HTML viewers** — generalize `ff_ecm_native.py`;
   REV↔native consistency; one full-res annotated viewer per material (measured Pa + band + ⟨z⟩/mesh/S).
5. **Nonlinear master curve K(σ) collapse (KB-1.12) + negative normal stress N1** (needs #3) — the thesis-grade
   nonlinear validation; |N1|~|σxz| at ~20% strain is an ungameable magnitude discriminator (Janmey/Kim p13).
6. **Systematized static traction-remodeling + contact-guidance anisotropy** on validated collagen (extends S6).

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
