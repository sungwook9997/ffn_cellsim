# PI-GAP Evidence Cards — one-pass provenance request (2026-07-25)

**What this is.** A single fill-in-the-blanks form so the PI can source, in one sitting, every OPEN PI-GAP
parameter that is currently blocking an `ac/engine` Active-Cell component from going quantitative-native. Each
card names the component(s) it unblocks, the exact symbol + units + code location that consumes it, the current
`null`/provisional value and why it is insufficient, any candidate literature lead already in the KB / the
2026-07-23 sourcing pass (real DOIs only — nothing invented; "NOT FOUND" where none exists), and an empty
**PROVENANCE SLOT** for the PI to complete.

**How many components this unblocks.** **32 parameter slots across 9 components / connectors.** Filling them
closes the magnitude/native gate for: cortex NMII actuator, stress-fiber (SF) NMII, microtubule aster, keratin/
vimentin IF cage, filopodium fascin bundle, lamellipodium Arp2/3 dendritic net, ERM membrane↔cortex clutch,
the osmotic turgor baseline (whole resting setpoint), and the LINC/nucleus tether. The single highest-leverage
group is the NMII force-scale cluster (cards N1–N9): it gates BOTH the cortex γ magnitude (GATE-A/B) AND the SF
traction magnitude, because both read the same per-head force/duty. Turgor (card T1) and the NMII resting
setpoint (card N9) together gate the entire resting-baseline that every other measurement is taken FROM.

**Rules (carried from CLAUDE.md, non-negotiable).** These are **literature-first**: prefer a real in-vivo
sourced value; MCF7/epithelial-specific where it exists, an explicitly-flagged proxy (as HeLa is already used
for cortical tension) where it does not. **No fitting to PI experimental data.** If a quantity is genuinely
unmeasured anywhere, the PI *chooses* a value per the physiological-baseline rule and records it as PI-authored
(owner + expiry_trigger), never tuned to make a gate/band pass. Do **not** leave a production parameter at a
convenient null/zero. For each slot fill:

```
value = ___ | units ___ | source (DOI/PMID) ___ | KU-id ___ | verdict ___
```

`verdict` ∈ {OK-sourced, PROXY-ratified, PI-CHOICE (unmeasured), DEFER}. Confirm any cited KB source's
`source_audit.verdict` is OK before it lands (citation-integrity rule).

---

## Component / gate status at a glance

| # | Component / connector | Ladder state today | Blocked-on cards | Native gate |
|---|---|---|---|---|
| 1 | cortex NMII actuator (`cortex_motor_slice`) | CUDA_UNIT (GATE-B dynamic slice native-validated; tension EMERGES from binding events) | N1–N9 | γ-magnitude INVALID |
| 2 | stress fiber NMII + backbone (`sf_mechanics`/`sf_population`) | KERNEL_BOUND (passive laws launch); NMII magnitude SEAMED | N1–N5, N9, S1 | traction-magnitude INVALID |
| 3 | microtubule aster (`microtubule_rig`) | rig CUDA_UNIT; DI transition layer deferred | M1–M3 | DI dynamics INVALID |
| 4 | keratin/vimentin IF cage (`intermediate_filament_rig`) | linear-tangent reference installable | I1–I2 | WLC strain-stiffening INVALID |
| 5 | filopodium fascin bundle (`params_i0b4`) | analytic shape gates only | F1–F5 | bundle-condensation INVALID |
| 6 | lamellipodium Arp2/3 net (`params_i0b4`) | analytic shape gates only | L1–L4 | protrusion/retrograde-flow INVALID |
| 7 | ERM membrane↔cortex clutch (`erm_cortex_connector`/`_slice`) | CUDA_UNIT (GATE-B ERM native landed) | E1–E5 | shedding-kinetics magnitude INVALID |
| 8 | osmotic turgor baseline (`assemble.PI_0_PA`) | GATE-A **PROVISIONAL** (see Card G-A) at HeLa-proxy 40 Pa | T1, G-A | resting setpoint = proxy |
| 9 | LINC / nucleus tether (`linc_connector`) | KERNEL_BOUND (Hookean+stiffening); magnitude PI-gap | C1–C4 | tether magnitude INVALID |

---

## 1. NMII actomyosin force-scale cluster — cortex γ + SF traction (GATE-A / GATE-B master blocker)

Consumed by the cortex NMII actuator and the SF NMII load-path. Source-of-truth register:
`aleph/components/motor/params_i0b3.yaml`. Native magnitude gate INVALID until N1–N5 close
(`native_run_blocked_on: [F_stall_head, N_side, v0, kappa_hill, k_xb]`). Per PI 2026-07-24, the *representation
of resting myosin* (N9) is additionally an open modeling decision.

### N1 — per-head isometric stall force `F_stall_head` (F_head)
- **Unblocks:** cortex NMII, SF NMII. Units **pN**.
- **Consumes at:** `aleph/components/motor/params_i0b3.yaml:59`; feeds `aleph/components/motor/resting_setpoint.py:250` (`f_head`) and every Hill-step kernel.
- **Current:** `value: null`. Two unaudited claims — **0.5 pN** (AFINES/Freedman 2017; Kovács 2003; Tam 2021) vs **2.0 pN** (Billington structural; `ff/myosin_linear` KB-3.18). Insufficient: neither is an audited single-molecule SourceEvidence row; the 0.5-vs-2 split is partly a units conflation (0.5 is Hill curvature, not a force — see prior pass §1b).
- **Candidate lead:** direct single-molecule NM2A/NM2B isometric stall force **NOT FOUND** (2026-07-23 pass §1b, §7.1). Closest primary is step-size only (Norstrom 2010, `10.1074/jbc.M110.123851`). Muscle-trap history ~1–5 pN makes 0.5–2 pN plausible but not NM2-sourced. **PI to select f_head + NM2A/NM2B isoform.**
- **PROVENANCE SLOT:** `value = ___ | units pN | source (DOI/PMID) ___ | KU-id ___ | verdict ___`  (also record isoform: NM2A / NM2B)

### N2 — explicit heads per anti-parallel half-filament `N_side`
- **Unblocks:** cortex NMII, SF NMII (topology counts, whole-cell unique-ID budget). Units **count**.
- **Consumes at:** `params_i0b3.yaml:78`.
- **Current:** `value: null`. **10** (AFINES simplification) vs **28–30** (Billington 2013 JBC / Niederman-Pollard 1975 EM). Insufficient: fine-grained mandate rejects the AFINES-10 simplification, but 28–30 not yet ratified for I3.
- **Candidate lead:** Billington 2013 (`10.1074/jbc...`) / Niederman-Pollard 1975 EM ~28–30 molecules/minifilament; archived `mcf7_baseline.yaml` production override = 28. Structural count is well-supported → **PI ratify ≈28–30 (N_side≈14–15/half).**
- **PROVENANCE SLOT:** `value = ___ | units count | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N3 — unloaded head stepping velocity `v0`
- **Unblocks:** cortex NMII, SF NMII (Hill FV intercept). Units **µm/s**.
- **Consumes at:** `params_i0b3.yaml:97`.
- **Current:** `value: null`. **0.12 µm/s** (KB-PIV-4; `ff/myosin_linear`) vs **0.2 µm/s** (Kovács 2003). Insufficient: differ by assay temperature / load definition / isoform.
- **Candidate lead:** Kovács 2003 JBC 278:38132 (NM2A single-molecule). **PI bind assay-T + load-definition + isoform.**
- **PROVENANCE SLOT:** `value = ___ | units µm/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N4 — Hill force-velocity curvature `kappa_hill` (a/F₀)
- **Unblocks:** cortex NMII, SF NMII (per-head FV shape). Units **dimensionless**.
- **Consumes at:** `params_i0b3.yaml:114`.
- **Current:** `value: null`. THREE-way: **0.25** (Hill 1938 muscle) / **0.5** (Kovács 2003 NM2A — the value the archived Stam-Hocky port actually used) / **∞** (PI 2026-07-07: non-muscle aggregate FV is linear). Insufficient: unreconciled head-level vs aggregate form.
- **Candidate lead:** Kovács 2003 JBC 278:38132 head-level a/F₀=0.5 is sourced; the linear limit is a ratified aggregate finding. **PI reconcile Kovács-0.5 (head-level) vs Freedman/Tam-linear.**
- **PROVENANCE SLOT:** `value = ___ | units dimensionless | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N5 — per-head crossbridge stiffness `k_xb`  ⚠ MASTER FORCE KNOB
- **Unblocks:** cortex NMII, SF NMII (head strain at load; sets CFL). Units **pN/µm**.
- **Consumes at:** `params_i0b3.yaml:144`.
- **Current:** `value: null`, physical band **100–1000 pN/µm** (~0.1–1 pN/nm). Insufficient: the legacy 1 pN/µm AFINES-soft surrogate collapses stall geometry (working-stroke strain F_stall/k_xb becomes larger than the minifilament).
- **Candidate lead:** band is physically bounded (crossbridge stiffness 0.1–1 pN/nm); no MCF7/isoform single value in KB. **PI source per isoform/assay in 100–1000 pN/µm** (Newton-closure gate certifies the strain is physical — never tune to a γ band).
- **PROVENANCE SLOT:** `value = ___ | units pN/µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N6 — head↔backbone arm rest length `r0_head`
- **Unblocks:** cortex NMII, SF NMII (arm-spring zero-strain length; NOT the crossbridge rest). Units **µm**.
- **Consumes at:** `params_i0b3.yaml:162`.
- **Current:** `value: null`, reference band **0–0.21 µm** (archived head_rest_length = 0.200 µm perpendicular offset). Insufficient: archived value not source-audited; drives the built head offset geometry (`MinifilamentTopology.head_offset_um`).
- **Candidate lead:** archived 0.2 µm perpendicular offset (phase1_h3.yaml) — geometric, not sourced. **PI source myosin head reach / arm rest length.**
- **PROVENANCE SLOT:** `value = ___ | units µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N7 — per-head actin attachment rate `k_on`
- **Unblocks:** cortex NMII, SF NMII (with k_off0 sets EMERGENT engaged fraction φ_b = k_on/(k_on+k_off)). Units **1/s**.
- **Consumes at:** `params_i0b3.yaml:248`; enforced required at `aleph/engine/cortex_motor_slice.py:166,207`.
- **Current:** `value: null`, provisional **50/s** (`ff/hand_kmc` NM2A preset) — config-chosen, NOT an audited head-level SourceEvidence row.
- **Candidate lead:** `ff/hand_kmc` NM2A k_on=50/s (unaudited). **PI audit a head-level NM2A attachment rate.**
- **PROVENANCE SLOT:** `value = ___ | units 1/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### N8 — Pereverzev catch-slip constants `k_catch0` / `x_catch` / `k_slip0` / `x_slip`
- **Unblocks:** cortex NMII, SF NMII — the PHYSIOLOGICAL load-STRENGTHENING detach (duty RISES under load, Kovács 2007), the current default in the dynamic slice. Units: rates **1/s**, lengths **µm/nm**.
- **Consumes at:** `aleph/engine/cortex_motor_slice.py:211-215` (`NMIICatchSlipParams`, mandatory when detach = CATCH_SLIP).
- **Current:** all **None** (GAP). Insufficient: constrained by Kovács 2007 ADP-release slowing (5×/12× for 2A/2B) + NM2B duty 0.2–0.3 but **no direct single-molecule catch-slip fit**. (Prior pass §1 fidelity flag: the pure Bell-slip alternative has the WRONG sign — load must RAISE duty.)
- **Candidate lead:** Kovács, Thirumurugan, Knight, Sellers 2007 PNAS 104:9994 (`10.1073/pnas.0701181104`) — load-dependence magnitudes; Nagy 2013 JBC (`10.1074/jbc.M112.424671`) NM2B duty. No numeric Pereverzev fit exists → **PI author the 4 constants against Kovács-2007 constraints, or ratify a proxy.**
- **PROVENANCE SLOT (×4):**
  - `k_catch0 = ___ | units 1/s | source ___ | KU-id ___ | verdict ___`
  - `x_catch  = ___ | units nm  | source ___ | KU-id ___ | verdict ___`
  - `k_slip0  = ___ | units 1/s | source ___ | KU-id ___ | verdict ___`
  - `x_slip   = ___ | units nm  | source ___ | KU-id ___ | verdict ___`

### N9 — resting bound-head fraction + resting per-head force (resting-γ representation)  ⚠ MODELING DECISION
- **Unblocks:** cortex NMII **static** resting setpoint (how resting myosin appears at t=0). Units: fraction **dimensionless**, force **pN**.
- **Consumes at:** `aleph/components/motor/resting_setpoint.py:62-63` (`fraction`, `per_head_force_pn`).
- **Current:** both PI-GAP (unset → motor left unbound at t0). Per PI 2026-07-24: coarse fiber-quotient does NOT close the native residual — the leftover is discrete-myosin local non-equilibrium, so **how resting myosin is represented in the static baseline is an open PI modeling decision** (continuous prestress vs dynamic-modulator vs percentile), NOT a solver fix.
- **Candidate lead:** resting duty ≈ NM2 unloaded duty rising under load (Kovács 2003/2007; Nagy 2013); per-head force = N1. The *fraction* has no MCF7 resting datum. **PI choose the representation + values.**
- **PROVENANCE SLOT (×2):**
  - `fraction (resting bound-head) = ___ | units dimensionless | source ___ | KU-id ___ | verdict ___`
  - `per_head_force_pn (resting)   = ___ | units pN | source ___ | KU-id ___ | verdict ___`
  - representation decision: `continuous-prestress / dynamic-modulator / percentile → ___`

---

## 2. Microtubule aster — dynamic instability + count

Consumed by `aleph/engine/microtubule_rig.py` + `aleph/components/solid/microtubule.py`. Growth SPEEDS are
per-MT injected arrays; the DI SWITCH is a `DynamicInstabilityRateCard` with NO default.

### M1 — catastrophe / rescue rates `f_catastrophe_per_s` / `f_rescue_per_s`
- **Unblocks:** MT aster DI phase switch. Units **1/s**.
- **Consumes at:** `microtubule_rig.py:590-632` (`DynamicInstabilityRateCard`; `ratified=False` on the proxy factory).
- **Current:** no default (UNSET PI-GAP slot). A caller must supply a card; the Rusan proxy ships `ratified=False` so it cannot silently drive a native run.
- **Candidate lead:** **Rusan 2001 MBoC 12:971** (`10.1091/mbc.12.4.971`), LLCPK-1α interphase plus-end: f_cat **0.026 s⁻¹**, f_rescue **0.175 s⁻¹**. MCF7-specific **NOT FOUND**. In-vivo mammalian-epithelial → recommended **PROXY, PI to ratify** (as HeLa is used for cortical tension). Prior draft KB-DRAFT-3-07 holds in-vitro Walker-1988 values (reject per physiological-baseline rule).
- **PROVENANCE SLOT (×2):**
  - `f_catastrophe_per_s = ___ | units 1/s | source ___ | KU-id ___ | verdict ___`
  - `f_rescue_per_s      = ___ | units 1/s | source ___ | KU-id ___ | verdict ___`
  - 2-state vs 3-state (interphase MTs ~73.5% PAUSED): model PAUSED? `___`

### M2 — plus-end growth / shrink speeds `v_grow` / `v_shrink`
- **Unblocks:** MT aster length dynamics. Units **µm/s**.
- **Consumes at:** `microtubule_rig.py:891-892` (per-MT arrays `v_grow_um_per_s`/`v_shrink_um_per_s`, caller-injected).
- **Current:** caller-injected, no sourced default carried.
- **Candidate lead:** Rusan 2001 plus-end: v_grow **11.5 µm/min ≈ 0.192 µm/s**, v_shrink **13.1 µm/min ≈ 0.218 µm/s** (same proxy dataset as M1 — keep internally consistent). MCF7-specific NOT FOUND.
- **PROVENANCE SLOT (×2):**
  - `v_grow   = ___ | units µm/s | source ___ | KU-id ___ | verdict ___`
  - `v_shrink = ___ | units µm/s | source ___ | KU-id ___ | verdict ___`

### M3 — MT count per cell `n_mt`
- **Unblocks:** MT aster population (whole-cell unique-ID budget). Units **count**.
- **Consumes at:** `aleph/components/solid/microtubule.py:104` (`n_mt: int = 40` — provisional geometry, docstring L29 "not a sourced production count").
- **Current:** provisional **40**. Insufficient: order-of-magnitude placeholder, MCF7-absent.
- **Candidate lead:** **NOT FOUND** — no MCF7 MT number. **PI to source or choose** (interphase mammalian MT number is measurable; surface a target cell-type count).
- **PROVENANCE SLOT:** `value = ___ | units count | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## 3. Keratin / vimentin IF cage — nonlinear WLC card

Consumed by `aleph/engine/intermediate_filament_rig.py`. Linear-tangent reference is installable (k_bb
DERIVED, no magic number); the strain-stiffening WLC card (EA enthalpic wall + x_max α→β crossover) has NO
default. Use AXIAL tensile numbers — bending-derived E (300–900 MPa) is a ~100× anisotropy artifact.

### I1 — keratin enthalpic-wall axial stiffness `IF_KERATIN_WALL_EA_PN` (EA)
- **Unblocks:** keratin (MCF7 primary IF) WLC cage. Units **pN** (EA = axial force scale).
- **Consumes at:** `intermediate_filament_rig.py:75` (`IF_KERATIN_WALL_EA_PN: float | None = None`); card built at `:722-743` (`NonlinearCableCard.axial_stiffness_pn`).
- **Current:** `None` (PI-GAP, no default). Vimentin IS sourced (`IF_VIMENTIN_WALL_EA_PN = 13000.0` pN, Qin 2009; toe EA 0.5–0.9 nN).
- **Candidate lead:** keratin K8/K18 precise EA **NOT FOUND** — Lorenz 2019 PRL 123:188102 (`10.1103/PhysRevLett.123.188102`) / Lorenz 2023 Matter 6:2019 (`10.1016/j.matt.2023.04.014`) are **paywalled** (only qualitative ">2.5-fold, 133% in-cell without rupture"). **PI read full-text/SI or surface.** (Hollow vs solid cross-section A is a PI modeling choice.)
- **PROVENANCE SLOT:** `value = ___ | units pN (EA) | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### I2 — IF entropic→enthalpic crossover `crossover_ratio` (x_max, α→β unfolding)
- **Unblocks:** keratin + vimentin WLC nonlinearity. Units **dimensionless** (x/L_c ∈ (0,1)).
- **Consumes at:** `intermediate_filament_rig.py:722-743` (`NonlinearCableCard.crossover_ratio`, validated in (0,1)); vimentin card factory `:818`.
- **Current:** no default; a caller must supply the card. Vimentin modeling note x_max ~200–300% strain (`IF_VIMENTIN_XMAX_STRAIN=2.5`) is a rupture/plateau NOTE, not the x/L_c crossover.
- **Candidate lead:** vimentin Block 2018 (`10.1126/sciadv.aat1161`) / Forsting 2019 (`10.1021/acs.nanolett.9b02972`); desmin analog 240% (Kreplak 2008 `10.1529/biophysj.107.119826`). Keratin crossover **NOT FOUND** (same paywall as I1). **PI supply crossover per material.**
- **PROVENANCE SLOT (per material):**
  - `crossover_ratio (keratin)  = ___ | units — | source ___ | KU-id ___ | verdict ___`
  - `crossover_ratio (vimentin) = ___ | units — | source ___ | KU-id ___ | verdict ___`

---

## 4. Filopodium fascin bundle

Consumed by `aleph/components/weave/params_i0b4.yaml`. Native protrusion/bundle-condensation gate INVALID until these
close (`native_run_blocked_on`). Values are `null` + `evidence_status: "GAP — PI"`; do NOT script the bundle —
the low-angle tight bundle must CONDENSE from an isotropic seed.

### F1 — fascin crosslink stiffness `k_fascin` (`fascin_bundle_stiffness`)
- **Unblocks:** filopodium bundle crosslinker. Units **pN/µm**.
- **Consumes at:** `params_i0b4.yaml:119`.
- **Current:** `null`. `ff` uses anchored FILAMIN 8.2e5 as a flagged stand-in. Insufficient: filamin ≠ fascin.
- **Candidate lead:** **NO single-molecule fascin stiffness exists anywhere** (KB-DRAFT-3-23; prior pass §4). Do NOT inherit filamin. Honest path = derive from bundle bending stiffness: L_p 10→**150 µm** at fascin:actin 1:2 (Takatsuki 2014 `10.1016/j.bbagen.2014.01.012`), κ_bundle=L_p·k_BT. Coupled-bundle papers (Claessens 2006; Bathe 2008; Vignjevic 2006) are **absent from the corpus — ingest before the derivation is defensible.** **PI approve the bundle-derivation path.**
- **PROVENANCE SLOT:** `value = ___ | units pN/µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### F2 — inter-filament spacing `d_fascin` (`fascin_spacing_um`)
- **Unblocks:** filopodium tight-bundle reach. Units **µm**.
- **Consumes at:** `params_i0b4.yaml:129`.
- **Current:** `null`; claim_ref ~7–8 nm (geometric, architecture_spec) — NOT yet an audited SE row.
- **Candidate lead:** **Courson & Rock 2010 JBC** (`10.1074/jbc.M110.123117`): inter-filament ~8 nm, repeat ~36–38 nm, strictly parallel. Geometric value is well-supported → **PI ratify ~8 nm.**
- **PROVENANCE SLOT:** `value = ___ | units µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### F3 — fascin bundling angle `phi_fascin` (`fascin_bundling_angle_deg`)
- **Unblocks:** filopodium bundle geometry (emergent low angle; a condensation gate target). Units **deg**.
- **Consumes at:** `params_i0b4.yaml:140`.
- **Current:** `null`.
- **Candidate lead:** Courson & Rock 2010 (strictly parallel, low-angle) — no explicit degree in KB. **PI source/choose the bundling angle.**
- **PROVENANCE SLOT:** `value = ___ | units deg | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### F4 — filaments per bundle `N_fil_bundle` (`filaments_per_bundle`)
- **Unblocks:** filopodium population-ledger input (whole-cell unique-ID budget). Units **count**.
- **Consumes at:** `params_i0b4.yaml:150`.
- **Current:** `null`; order candidate 10–30 (architecture_spec); number-of-bundles open.
- **Candidate lead:** order 10–30 (architecture_spec, unaudited); Euler-buckling ~7.7 pN N²-tight-bundle context KB-3.8 (Mogilner & Rubinstein 2005; Pronk 2008). **PI bind per-bundle count + bundle number.**
- **PROVENANCE SLOT:** `value = ___ | units count | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### F5 — tip-nucleator density `rho_tip` (`tip_nucleator_density`)
- **Unblocks:** filopodium/leading-edge tip nucleation loci. Units **1/µm²**.
- **Consumes at:** `params_i0b4.yaml:161`.
- **Current:** `null`.
- **Candidate lead:** per-tip formin count **NOT FOUND** (prior pass §4). **PI source tip-nucleator (formin/Arp2/3) density.**
- **PROVENANCE SLOT:** `value = ___ | units 1/µm² | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## 5. Lamellipodium Arp2/3 dendritic net

Consumed by `aleph/components/weave/params_i0b4.yaml`. Native protrusion/retrograde-flow magnitude gate INVALID until
these close. (Branch-angle knobs θ₀/σ_θ/k_θ are already GROUNDED — Faessler 2020 — and are NOT gaps.)

### L1 — max Arp2/3 branch-nucleation rate `k_arp0`
- **Unblocks:** lamellipodium branch nucleation (enters k_branch = k_arp0·npf·c/(c+K_m)). Units **1/s**.
- **Consumes at:** `params_i0b4.yaml:79`.
- **Current:** `null`. **Candidate lead: NOT FOUND** in KB → **PI source Arp2/3 branch-nucleation kinetics (isoform/assay).**
- **PROVENANCE SLOT:** `value = ___ | units 1/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### L2 — Michaelis monomer half-saturation `K_m` (`k_m`)
- **Unblocks:** lamellipodium monomer flux-limit knee (reads the I1c c field). Units **µM**.
- **Consumes at:** `params_i0b4.yaml:89`.
- **Current:** `null`. **Candidate lead: NOT FOUND** → **PI source Arp2/3 monomer dependence.**
- **PROVENANCE SLOT:** `value = ___ | units µM | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### L3 — capping-protein termination rate `k_cap`
- **Unblocks:** lamellipodium dendritic mesh size (steady length ~ v_grow/k_cap). Units **1/s**.
- **Consumes at:** `params_i0b4.yaml:99`.
- **Current:** `null`. **Candidate lead: NOT FOUND** → **PI source capping kinetics** (do NOT tune to a mesh size).
- **PROVENANCE SLOT:** `value = ___ | units 1/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### L4 — NPF areal density `rho_NPF` (`npf_areal_density`)
- **Unblocks:** lamellipodium leading-edge nucleation gate (npf activity ∈ [0,1]; a declared membrane field). Units **1/µm²**.
- **Consumes at:** `params_i0b4.yaml:109`.
- **Current:** `null`. **Candidate lead: NOT FOUND** (MCF7 leading edge) → **PI source NPF (WASP/WAVE) density.**
- **PROVENANCE SLOT:** `value = ___ | units 1/µm² | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## 6. ERM membrane↔cortex clutch — Bell kinetics + density

Consumed by `aleph/engine/erm_cortex_connector.py` + `erm_cortex_slice.py`. `k_erm` is SOURCED (Braunger
2014 single-bond, 4.6e3 pN/µm, PI-ratified 2026-07-21) and is NOT a gap. The four Bell kinetics are validated at
construction (`erm_cortex_connector.py:174-177` rejects null/non-positive) — no MCF7 ERM on/off datum, so the
GATE-B native slice runs only under an explicit provenance banner.

### E1 — attach rate `k_on`
- **Unblocks:** ERM clutch rebinding. Units **1/s**. **Consumes at:** `erm_cortex_slice.py:126`, `erm_cortex_connector.py:159`.
- **Current:** caller-required, no sourced default. **Candidate lead:** NOT FOUND (MCF7 ERM). **PI source.**
- **PROVENANCE SLOT:** `value = ___ | units 1/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### E2 — zero-force detach rate `k_off0`
- **Unblocks:** ERM Bell-slip shedding prefactor. Units **1/s**. **Consumes at:** `erm_cortex_slice.py:127`, connector `:160`.
- **Current:** caller-required, no sourced default. **Candidate lead:** NOT FOUND (MCF7 ERM). **PI source.**
- **PROVENANCE SLOT:** `value = ___ | units 1/s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### E3 — Bell characteristic force `bell_force_pn` (F₀)
- **Unblocks:** ERM Bell-slip force sensitivity. Units **pN**. **Consumes at:** `erm_cortex_slice.py:128`, connector `:161`.
- **Current:** caller-required, no sourced default. **Candidate lead:** cross-check only — KB single-ERM 11.4 pN / F_unbind ~50 pN (prior pass §5). **PI source F₀.**
- **PROVENANCE SLOT:** `value = ___ | units pN | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### E4 — capture radius `capture_radius_um`
- **Unblocks:** ERM rebind reach. Units **µm**. **Consumes at:** `erm_cortex_slice.py:129`, connector `:162`.
- **Current:** caller-required, no sourced default. **Candidate lead:** NOT FOUND. **PI source/choose (geometric).**
- **PROVENANCE SLOT:** `value = ___ | units µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### E5 — ERM areal density `rho_0`
- **Unblocks:** ERM linker number per node (N_linker/node = ρ₀·A_node, grid-invariant). Units **1/µm²**.
- **Consumes at:** whole-cell ERM population (do NOT hard-code a per-node integer — prior pass §5 wiring).
- **Current:** not yet wired to a sourced ρ₀. **Candidate lead:** **Alert 2015 Biophys J 108:1878** (`10.1016/j.bpj.2015.02.027`) ρ₀ ≈ **100 µm⁻²** baseline (sweep 100–1000); order 10²–10³ robust across 3 methods (Charras 2007 `10.1529/biophysj.107.113605`); active fraction smaller (Bosk 2011 `10.1016/j.bpj.2011.02.039`). MCF7 direct value NOT FOUND. **PI adopt estimate + decide available vs active fraction.**
- **PROVENANCE SLOT:** `value = ___ | units 1/µm² | source (DOI/PMID) ___ | KU-id ___ | verdict ___`  (available / active-scaled?)

---

## 7. Osmotic turgor baseline

### T1 — resting turgor `Π₀` (`PI_0_PA`)
- **Unblocks:** the WHOLE resting cortex/membrane setpoint — every downstream measurement is taken FROM this baseline (physiological-baseline rule). Units **Pa** (1 Pa = 1 pN/µm²).
- **Consumes at:** `aleph/components/incumbent/assemble.py:128` (`PI_0_PA = 40.0`), read by `driver.py:1057`, `dump_state.py:91`, `bleb_growth.py:185`.
- **Current:** **40.0 Pa** — a HeLa-interphase proxy (Fischer-Friedrich 2014 cortical-tension anchor). GATE-A was run at this value and is **PROVISIONAL — PI ratification pending** (Card G-A), and the value is a proxy, not MCF7. Card-3 candidate **72 Pa**.
- **Candidate lead:** Fischer-Friedrich 2014 HeLa 0.2 mN/m → γ≈140 pN/µm at R=7.5 µm (proxy currently used); P2 evidence Card-3 proposes 72 Pa (PI sign-off pending). MCF7-direct turgor NOT FOUND. **PI ratify 40 Pa proxy, adopt 72 Pa, or source MCF7.**
- **PROVENANCE SLOT:** `value = ___ | units Pa | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## 8. Stress-fiber actin axial-stretch stiffness

### S1 — SF actin backbone stretch stiffness `k_axial` (EA_actin)
- **Unblocks:** SF component axial backbone (`sf_mechanics` KERNEL_BOUND passive law). Units **pN/µm**.
- **Consumes at:** `aleph/engine/sf_mechanics.py:215` (`link_k`), required at `:254-271` (REQUIRED, no default — raises if absent), applied `:309`. (The α-actinin crosslink stiffness 4.6e5 pN/µm is SOURCED, Ferrer 2008 — NOT this gap.)
- **Current:** **required, no default** (modeling GAP). Insufficient: NF2007 treats the actin backbone as inextensible (no axial penalty spring), so there is **no sourced EA_actin in `ff.units`**.
- **Candidate lead:** **NOT FOUND** — no sourced actin single-filament EA for an axial Hookean cable in KB. F-actin stretch modulus is measurable (~44 pN·µm² bending / ~4.4e4 pN axial-EA order from Gittes/Kojima single-filament literature is a candidate direction, NOT in corpus). **PI source EA_actin or choose per the inextensible-backbone modeling decision.**
- **PROVENANCE SLOT:** `value = ___ | units pN/µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## 9. LINC / nucleus tether (secondary — from 2026-07-23 pass §6, still OPEN)

Consumed by `aleph/engine/linc_connector.py` (`stiffness_pn_per_um` is a caller-supplied explicit input;
docstring `:771` "LINC stiffness/kinetic magnitudes remain source/PI gaps"; `:1256` surface-to-PI guard). LINC
resting TENSION ~8 pN is sourced (Déjardin 2020, `10.1083/jcb.201908036`) — a tension, not a stiffness.

### C1 — LINC/nesprin tether stiffness `k_linc` (`stiffness_pn_per_um`)
- **Unblocks:** LINC connector (cortex/IF ↔ nucleus). Units **pN/µm**.
- **Consumes at:** `linc_connector.py:763,781` (validated positive-finite) / `:947` (kernel `stiffness`).
- **Current:** caller-required; NO literature stiffness. Déjardin gives a tension only; spectrin-repeat unfolding 25–35 pN (Rief 1999) is a rupture scale, not resting stiffness. **NOT FOUND → PI-authored value required.**
- **PROVENANCE SLOT:** `value = ___ | units pN/µm | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### C2 — LINC areal density
- **Unblocks:** number of nesprin tethers. Units **1/µm²** (or per-nucleus count).
- **Current:** **NOT FOUND** (prior pass §6.1). Zero LINC/nesprin claim rows in KB. **PI source/choose.**
- **PROVENANCE SLOT:** `value = ___ | units 1/µm² | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### C3 — NE rupture threshold
- **Unblocks:** nuclear-envelope rupture criterion (tension/strain/curvature). Units **as chosen**.
- **Current:** WEAK. KB-3.32 asserts ~1% projected-area strain (Zhang-Lele 2018, MCF-10A) but the primary DOI/PMID could NOT be confirmed; Denais/Lammerding 2016 (`10.1126/science.aad7297`) confirms the phenomenon, gives no numeric threshold. **PI resolve primary source or choose.**
- **PROVENANCE SLOT:** `value = ___ | units ___ | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

### C4 — nucleoplasm viscosity `η_nuc`
- **Unblocks:** nucleus interior rheology. Units **Pa·s**.
- **Current:** KB-carried ~52 Pa·s bulk / G′ 18 Pa (Tseng 2004 `10.1242/jcs.01073`; de Vries 2007) but NOT PubMed-confirmed this pass. **PI resolve DOIs before hard-coding.**
- **PROVENANCE SLOT:** `value = ___ | units Pa·s | source (DOI/PMID) ___ | KU-id ___ | verdict ___`

---

## Fill-order recommendation (highest unblock-leverage first)

0. **Card G-A** (§ADDENDUM) — **GATE A closure ratification.** Not a parameter: a gate-contract sign-off that
   is currently blocking the SoT status of the static resting baseline (recorded PROVISIONAL 2026-07-25). Costs
   one decision, unblocks every downstream "GATE A CLOSED" claim.
1. **N1–N5, N7–N9** (NMII cluster) → unblocks cortex NMII **and** SF NMII magnitude at once (2 components, GATE-A/B).
2. **T1** (Π₀) + **N9** (resting representation) → close the resting baseline everything else measures from.
3. **E1–E5** (ERM) → GATE-B ERM native slice already landed; only the Bell/density values gate its magnitude.
4. **M1–M3, I1–I2, F1–F5, L1–L4, S1, C1–C4** → per-component native gates as each component reaches CUDA_UNIT.

No code or constants were changed by this document. It is a sourcing/ratification form for the PI only.

---

## ADDENDUM (2026-07-25, from ECM merge 8b2e7adb + traction-spine scoping)

### Card A1 — α2β1–collagen catch-bond Hand card (NEW — blocks the FA clutch KMC)

- **Unblocks:** `focal_adhesion` + `integrin_collagen_clutch` (the substrate traction spine clutch). Currently SEAMED (series-joint kernel-bound at VS-0, but the catch-bond *rate law* for the collagen ligand is absent).
- **Exact quantity / consumer:** an α2β1–collagen `HandParams` card in `aleph/laws/hand_kmc.py` (mirror `INTEGRIN_A5B1`, hand_kmc.py:150) — `k_on` [1/s], catch-slip `k_catch0`/`x_catch`/`k_slip0`/`x_slip` (Pereverzev/Bell two-pathway), `capture_radius` [µm], `link_k` [pN/µm]. Consumed by the load_path clutch `commit_irreversible` (`aleph/engine/load_path.py:747`).
- **Current state:** only **α5β1–fibronectin** exists (Kong 2009, KB-2.5). The engine clutch chemistry_card is `alpha2beta1_collagen` (contracts.py:281) — no matching rate card. Placeholder rate in the joint commit.
- **Candidate lead:** α2β1–collagen catch-bond — candidate primary literature to source (e.g. Kong/Thomas-lab single-molecule catch-bond studies of α2β1 on collagen). NOT YET SOURCED — PI to source or choose.
- **PROVENANCE SLOT:** `k_on = ___ | k_catch0 ___ | x_catch ___ | k_slip0 ___ | x_slip ___ | capture_radius ___ | link_k ___ | source (DOI) ___ | KU-id ___ | verdict ___`

### Card G-A — GATE A closure ratification (NEW — a gate-contract change, not a parameter)

- **Unblocks:** the SoT status of **GATE A** (converged static resting baseline) and therefore every downstream
  claim that reads "GATE A CLOSED" (R2→R8 native chain, `cell_engine/ROLLING_ROADMAP.md`). Units **n/a** — this
  is a **gate-contract sign-off**, filed here because the queue is where PI-gated items live, not a value GAP.
- **Exact quantity / consumer:** the gate contract itself — observable `max|PF|` (projected per-node net force,
  pN), threshold `< 0.21 pN`, configuration = the cell the gate is measured on. Recorded in
  `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §10–§11; asserted as SoT in `cell_engine/ROLLING_ROADMAP.md:5`.
  Appears in **no runtime code** (the executing criterion is the `sqrt(eps)·ℓ` triple test) — it is a
  *reporting* contract, ~12 hand-copied literals.
- **Current state:** **`PROVISIONAL — PI ratification pending`** (corrected 2026-07-25; previously recorded as
  "✅ GATE A CLOSED"). Insufficient because the closure changed **two things at once** between the FAIL and the
  PASS: (a) the **observable** — raw `max|F|` → projected `max|PF|` (§9, commit `d75e6631`, triggered by an
  **externally relayed** consultation, not internal diagnosis); and (b) the **configuration** — the resting
  myosin was removed from the static baseline (clean turgor + membrane subdiv=8 + ERM + nucleus, no discrete
  heads). The myosin-seeded configuration never got below ≈**1.36** ≈ **6.5×** the gate. `GATE_A_…` §10 ends
  "**SURFACED to PI for ratification**"; no ratification exists in any commit or doc, and GATE A was absent
  from the 2026-07-25 decision queue until this card. Additionally the threshold `0.21 pN` has **no derivation
  anywhere on disk** (first appearance a bare `0.21` in `SESSION_HANDOFF_2026-07-22.md:51`) — a derivation is
  being attempted in parallel and is **not** supplied by this card.
- **Candidate lead / what the PI is actually deciding:** the **physics reframe is not in question** and is not
  proposed for change — the residual is the local non-equilibrium at discrete myosin attachment points
  (`PF ≈ f_head`), ~7× more heads would be needed to get under the gate (unphysical), and resting cortical
  tension is dynamically maintained → GATE B. Numbers as recorded: clean baseline 0.1606 → **0.1398**, seeds
  0/1/2 = 0.1398 / 0.1416 / 0.1586; native coarse A/B OFF 1.48 / A 1.36 / B 1.51. The decision is whether the
  **static gate's observable + configuration change is ratified**, i.e. whether "GATE A" may mean *the clean
  turgor baseline measured in projected force* with resting myosin deferred to GATE B. Related open items that
  travel with it (from `AUDIT_AC_ENGINE_2026-07-25.md:114`): the reported max is **nucleus-block dominated**
  (clean actin max 0.0356) on an un-refined 642-vertex nucleus, the descent used a bespoke fixed-iteration
  harness rather than the driver's accepted-step path, and no machine artifact was archived (one PNG; the three
  ensemble numbers exist only in prose). The physiological-mesh rung is **0.28 > 0.21**.
- **RATIFICATION SLOT:** `verdict = ___ | observable ratified (max|PF|) Y/N ___ | configuration ratified (clean
  turgor baseline, resting myosin → GATE B) Y/N ___ | threshold 0.21 pN accepted pending derivation Y/N ___ |
  artifact required before CLOSED (npz/json + driver-path reproduction) Y/N ___ | nucleus refinement required
  Y/N ___ | date ___ | owner ___`
  `verdict` ∈ {RATIFIED-CLOSED, RATIFIED-PROVISIONAL (stays provisional pending named artifact), REJECTED
  (re-open GATE A), DEFER}.

### Registration note T-collagen — collagen modulus SoT sign-off (formality, not a value GAP)

The collagen modulus "conflict" is ALREADY resolved in-tree (concentration-resolved: 1 mg/mL≈1–5, 1.5 mg/mL=10–16 [Yang-Kaufman 2009], 3 mg/mL=30–100 Pa [Licup]; `ecm_library.py:111`, PI-ratified `AC_DECISION_CARDS_2026-07-22`). What remains is only the **Notion SoT registration + `make kb-check` refresh** per `P2_RATIFICATION_PACKAGE_2026-07-22`. No value decision needed — a sign-off formality.
