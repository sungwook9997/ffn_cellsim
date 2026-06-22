# Phase C wiring audit — is every mechanism actually in the hot loop? (2026-06-23)

**Auditor:** adversarial read-only audit (orchestrated).
**Target:** `ffn_sim/warp_port/dcm_warp_decohesion.py::run_decohesion` (the production hot loop).
**Method:** trace every imported mechanism to the exact `wp.launch` / host-update inside the
device step loop (`step_once`, l.557 + the main `for s in range(1, steps+1)` loop l.870), check
its enabling flag default + whether a production config turns it on at a physiological value, and
look for dead/placeholder/double-count paths. Each "UNWIRED/BUG" finding was refutation-tested
(grep for a launch elsewhere, helper functions, conditional launches) before being reported.

The driver's "full mechanistic production config" is the E R-sweep harness
`ffn_sim/scripts/spread_rsweep.py::run_one` (l.102–112): it turns on
`cadherin, ecm_clutch, lamellipodium, nucleus, surface_tension, bending, edge_edge,
substrate_wetting, use_substrate_well` with `integrator="implicit"`. It does **NOT** set
`coupling, gravity, lamel_clutch, filopodia, ecm_ligand(≠1), remesh_period, junction_switch,
division, necrosis`.

## Table

| Mechanism | Wired into hot loop? file:line | Default | Prod-config turns on? | Issue | Confidence |
|---|---|---|---|---|---|
| **core: turgor** | YES — `dcm_volume_kernel`+`_dp_from_vol`+`dcm_turgor_force_kernel` `dcm_warp_decohesion.py:589–597` | always-on | yes | none | HIGH |
| **core: cohesion+contact** | YES — `cohesion_grid_kernel`/`_cad` l.572–588 + `contact_grid_kernel`/`_cad` l.598–615 | always-on | yes | none (see A1 below for the adh branch) | HIGH |
| **core: cortex edge bonds** | YES — `_bond_accumulate` `:616` | always-on | yes | none | HIGH |
| **core: surface tension γ (B4)** | YES — `surface_tension_kernel` `:638–640` | OFF (`surface_tension=False`) | yes (`spread_rsweep` sets `surface_tension=True`, γ=1e-4) | none | HIGH |
| **core: bending (B5)** | YES — `bending_apply_kernel` (2 umbrella passes) `:628–636` | OFF | yes (`bending=True`) | none | HIGH |
| **core: edge-edge (A2)** | YES — `edge_edge_contact_kernel` `:650–656` | OFF | yes (`edge_edge=True`) | none | HIGH |
| **core: nucleus (E2)** | YES — `nucleus_force_kernel` `:620–627` | OFF | yes (`nucleus=True`) | none | HIGH |
| **core: cadherin (E1)** | YES — host `cad.update/upload` `:897–900`; device `cadherin_bond_force_kernel` `:659–664` | OFF (`cadherin=False`) | yes (`cadherin=True`, bundle 40) | none | HIGH |
| **core: ECM clutch (C6)** | YES — host `ecm.update/upload` `:902–905`; device `ecm_clutch_force_kernel` `:687–690` | OFF | yes (`ecm_clutch=True`, bundle 167) | none | HIGH |
| **core: remesh (A1 remesh)** | YES — `do_remesh()` `:892–895` (and settle `:841`) | **OFF (`remesh_period=0`)** | **NO** — `spread_rsweep` never sets `remesh_period` | not in the production spreading config; backlog says it's the fix for slivers under large spread but the prod sweep runs without it | HIGH |
| **A1 node-FACE coupling** (`--coupling`) | YES (it just sets `coh_adh`, consumed by contact/cohesion kernels l.598–615/572–588) | OFF (`coupling=False`) | **NO** (`spread_rsweep` never sets `coupling`) | **DOUBLE-COUNT + NO-OP-without-cadherin** — see findings | HIGH |
| **B2 lamel substrate clutch** (`--lamel-clutch`) | YES *only as a z-flag* into the SAME tether kernel `lamellipodium_tether_multicell` `:696–701` | OFF | NO | flag only flips anchor `z_basal` (host l.104); needs `--lamellipodium`, else silent no-op; still the floating-tether kernel, not a real node-plane clutch | HIGH |
| **B3 filopodia** (`--filopodia`) | YES — host `filo.update/upload` `:878–881`; device `filopodia_tip_{face,plane}_force_kernel` `:704–715` | OFF | NO | wired+self-tested but `k_tip`/cap PROVISIONAL (no KB datum, backlog l.92) | HIGH |
| **C4 ligand density** (`--ligand-density`) | YES — scales `p_on` in `EcmClutchHost.update` `dcm_ecm_clutch_host.py:120`; passed `ecm_ligand→EcmClutchParams` driver `:454` | 1.0 (baseline) | NO (stays 1.0; only meaningful with `--ecm-clutch`, which IS on) | correct but never set ≠1 in any committed config → Bare/Pre/Lam4 ordering untested in prod | HIGH |
| **C5 ECM mechano-feedback** | **NOT LAUNCHED** — `chan_odde_traction_factor`, `e_sub_to_k_sub`, `run_dcm_substrate_well_mechano_warp`, `dcm_substrate_well_strainstiff_kernel`, `run_dcm_substrate_well_strainstiff_warp` have ZERO callers repo-wide outside their own def file + self-test | OFF (no flag) | **NO** | **fully UNWIRED**; only the rigid-dish well/wetting kernels are imported (`dcm_warp_decohesion.py:34–36`); no `--E-sub`/`--mechano` CLI flag exists | HIGH |
| **D7 gravity** (`--gravity`) | YES — `gravity_body_force_kernel` `:559–561` | OFF | **NO** (`spread_rsweep` never sets `gravity`) | correct + derived (Δρ=55), but not in the production spreading config | HIGH |
| **D8 pen displacement cap** (`--pen-cap-frac`) | YES — `_vaxpy_active_capped` `:725–727` | **ON** (`pen_cap=True`) | **only on the implicit path** (`if implicit:` l.716); `spread_rsweep` uses `integrator="implicit"` so YES | on the explicit BAOAB path the cap is never applied (`_bd_step`, no clamp) — D8 only protects implicit runs | HIGH |
| **E Young-Dupré triplet gate** | N/A to hot loop — standalone gate `scripts/young_dupre_doublet_gate.py` calls `run_decohesion` for a 2/3-cell relax then compares to `young_dupre.triplet_angle` oracle (`cos(φ/2)=η/2`, l.135–143) | gate script | it IS the gate | real gate, not placeholder; correctly external to the loop | HIGH |

## Findings (survived refutation)

### F1 — C5 ECM mechano-feedback is fully UNWIRED [HIGH] (confirms the already-known)
`dcm_warp_decohesion.py:34–36` imports from `dcm_substrate_warp` **only** the rigid-dish kernels
`dcm_substrate_well_accum_kernel, dcm_wetting_scatter_kernel, dcm_wetting_cap_add_kernel,
dcm_wetting_scatter_integrin_kernel`. The C5 mechano API
(`run_dcm_substrate_well_mechano_warp` l.407, `dcm_substrate_well_strainstiff_kernel` l.483,
`run_dcm_substrate_well_strainstiff_warp` l.594, `chan_odde_traction_factor` l.370,
`e_sub_to_k_sub` l.339) has **zero callers anywhere in `ffn_sim/` except its own self-test**
(`grep -rn` repo-wide returned nothing outside `dcm_substrate_warp.py`). No `--E-sub`/`--mechano`
CLI flag exists. **Refutation attempt:** searched for the symbol in the implicit driver, all `.sh`
launchers, `spread_rsweep`, the gate scripts — none import or launch it. → genuinely dormant.
*Fix direction:* add an `E_sub: float|None` arg to `run_decohesion`; when set, replace the
`dcm_substrate_well_accum_kernel` launch (l.667) with `run_dcm_substrate_well_mechano_warp`'s
E-derived `k_sub` + traction factor, and add the strain-stiffen kernel to the basal force; expose
`--E-sub` (physiological E_opt~5 kPa, KB-2.8). Note backlog l.93 already flags C5 needs a
deformable substrate to be meaningful.

### F2 — A1 `--coupling` double-counts adhesion AND is a no-op without `--cadherin` [HIGH]
In **cadherin mode** the design contract (kernel docstring `dcm_neighbor_warp.py:448–450`) is that
node-face adhesion is OFF and the cadherin trans-dimer bond is "the SOLE cell-cell adhesion."
`run_decohesion` sets `coh_adh = 0.0` in cadherin mode (l.410) — **unless `--coupling`**, which
sets `coh_adh = adh_strength`. With `--coupling` + `--cadherin`, apposed cells are then pulled by
**both** (a) the continuous node-face bilinear adhesion tent (`contact_grid_cad_kernel`'s `adh`
branch l.366–374 / `cohesion_grid_cad_kernel` l.303) **and** (b) the cadherin bond
(`cadherin_bond_force_kernel` l.660) → the same interface adhesion is counted twice.
**Refutation attempt:** checked whether the cad kernels exclude the adh branch when bonds exist —
they don't; both run every step on overlapping node pairs. Also: outside cadherin mode `coh_adh`
already defaults to `adh_strength` (l.383), so `--coupling` changes nothing there → it is a no-op
unless `--cadherin` is also set. *Fix direction:* don't re-enable the node-face `adh` branch as a
parallel attraction in cadherin mode; instead modulate the node-face `ω` (adh weight) by the
cadherin catch-bond state (the unified design already flagged in backlog l.86–87), so it is one
adhesion channel, not two. Not currently set in any production config (so harmless today) but the
flag's documented purpose is mis-implemented.

### F3 — D8 pen cap only guards the implicit path [HIGH, lower severity]
`pen_cap`/`pen_cap_frac` are applied via `_vaxpy_active_capped` (l.725–727) **inside `if implicit:`**.
The explicit BAOAB path (`_bd_step` l.731–733) has no displacement clamp. **Refutation attempt:**
searched `_bd_step` and the BAOAB kernel for any cap — none. The production sweep uses
`integrator="implicit"` so D8 is active there; but any `--integrator baoab` run gets no
interpenetration cap despite `pen_cap=True` default. *Fix direction:* either document D8 as
implicit-only or add the same clamp to the explicit step.

### F4 — B2 `--lamel-clutch` is a z-coordinate flag on the old tether, not a real node-plane clutch; silent no-op without `--lamellipodium` [HIGH, lower severity]
`lamel_clutch=True` only changes `LamellipodiumHost.z_basal` from `z0+seed_basal_offset` to `z0`
(`dcm_lamellipodium_host.py:104`); the device force is still the SAME
`lamellipodium_tether_multicell` spring (l.696), not the `ecm_clutch_force_kernel` node-plane
clutch. The backlog itself (l.91) flags "v1 = anchor on dish plane only" — so this matches the
*as-shipped* intent, but the mechanism name "substrate clutch" overstates it. Also `lam` is built
only when `lamellipodium=True` (l.355), so `--lamel-clutch` alone (no `--lamellipodium`) does
nothing and prints no warning. *Fix direction:* either route the lamellipodial anchor through the
real `ecm_clutch_force_kernel` (Pereverzev catch-slip, as backlog l.91 intends), or warn/error when
`--lamel-clutch` is passed without `--lamellipodium`.

### F5 — Production spreading config omits A1-remesh, gravity, coupling, ligand≠1, filopodia [HIGH, scope finding]
`spread_rsweep.py::run_one` (l.102–112) — the committed "full mechanistic" config — leaves
`remesh_period=0` (remesh OFF), `gravity=False`, `coupling=False`, `ecm_ligand=1.0`,
`filopodia=False`. So 5 of the "8 committed features" are NOT exercised by the one production
spreading harness. Of these, remesh (sliver control under large spread) and ligand-density
(the Bare/Pre/Lam4 experimental axis the platform exists to reproduce) are the most consequential
to leave off. *Fix direction:* PI decision — either add them to the prod config at their derived
values (remesh_period, ligand per condition) or document why they're excluded. Note: enabling both
`substrate_wetting=True` and `ecm_clutch=True` is harmless — the wetting proxy is auto-skipped when
`ecm is not None` (l.673), so the clutch wins; but it's a confusing redundant flag in the config.

## Mechanisms confirmed correctly wired (no action)
turgor, cohesion/contact, cortex edges, surface tension, bending, edge-edge, nucleus, cadherin
bonds, ECM clutch, gravity (correct when enabled), D7 derivation, C4 ligand math, the Young-Dupré
gate. All launch inside the device step loop with the host binders refreshed at their cadence in
the main loop (l.870–928).
