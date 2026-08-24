# Blebbing validation plan + FSI verification (2026-07-20)

> Companion to [`_historical/LITERATURE_POSITIONING_2026-07-20.md`](_historical/LITERATURE_POSITIONING_2026-07-20.md).
> PI ratified (2026-07-20) the **blebbing / hydraulic-fracture** regime as the primary scientific target
> (the regime where scale separation breaks down and filament heterogeneity feeds back on the fluid — the
> only defensible ground vs the "why not active gel" critique). This doc records the three follow-ups:
> (1) technical verification of the flagged Biot/FSI consistency issues, (2) bleb-parameter sourcing audit,
> (3) the out-of-band validation-harness design. Verification done by direct code read (agent-assisted).

---

## 1. Technical verification (the §5b FINDINGS + CFL)

| Check | Verdict | Detail |
|---|---|---|
| `M_biot` / storage_S | ✅ **RESOLVED (false alarm for production)** | Production `build_cell` (`assemble.py:434`) + driver + ng4 + ng2_ng3_ng6 all pass `storage_S=1e-4 → M=1e4 Pa`, exactly the c_v-anchored closure (`mobility/S = 5e-3/1e-4 = 50 µm²/s = c_v`). The flagged stale `M=1e3 Pa` survives only in `ng0_parity.py:51`, a param-agnostic Warp-vs-NumPy round-off test that stays self-consistent at c_v=50 via a 10× mobility. **Not a physics defect.** |
| `L_p` | 🔴 **REAL ~62× discrepancy** | Raw "1e-12 vs 1.6e-8" is mostly the m→µm unit factor (×1e6): yaml SI `1e-12 m/(s·Pa)` = `1e-6 µm/(s·Pa)`; code `1.6e-8 µm/(s·Pa)`. Residual ≈ **62×** — the code L_p is ~62× *less permeable* than the Jung 2011 MCF7/AQP5 yaml draft, and is anchored to `gamma_floor` (`assemble.py:114`), NOT to the sourced value. The yaml's "unit reconciliation open" flag is valid. **Directly sets bleb fill rate → must resolve before oracle B.** |
| Biot CFL | ✅ **Correct + respected** | `field_grid.cfl_dt = safety·S·dx²/(2·ndim·mobility)` (`field_grid.py:125-127`), matches the update `p += dt·mobility·acc/(S·dx²)`. Production subcycles: `n_sub = ceil(dt_phys/dt_cfl)`, `dt_sub ≤ dt_cfl` always (`scheduler.py:392-394`, `dt_phys=0.05 s`, safety 0.9). Stability honored regardless of outer dt. |
| **yaml drift (bonus finding)** | ⚠️ **single-source-of-truth hazard** | **No runtime code loads `params_i0b1.yaml`** — it is a ledger referenced only in comments; every runtime value is a hand-duplicated module constant (`BIOT_MOBILITY`, `BIOT_STORAGE_S`, `L_P` in `assemble.py:111-114`). yaml and code can silently diverge. Recommend: load the yaml at construction, or add a test asserting the constants equal the (unit-converted) yaml. |

**Codex handoff (still analytic, not code-value):** Peskin spread↔interp adjoint exactness
(`SolidDilatationCoupling` spread vs `PressureCoupling` interp normalization) and `−α∇p` ↔ `−α·div(v_s)`
energy-conjugacy remain to verify — these were not resolved by the value audit above.

## 2. Bleb-parameter sourcing audit

Nucleation threshold `f_rupt` is **derived, not tuned** — legitimate for validation. Two params need PI action.

| Param | Value | Units | Sourcing | file:line |
|---|---|---|---|---|
| `f_rupt` (ERM rupture threshold) | ≈**11.4 pN** (band 5–40) | pN | ✅ DERIVED `2π√(2κ_m(T_m+γ_MCA))`, KB-3.B1.4 | computed `compartments.py:483`; formula `aleph/laws/membrane_surface.py:223-227` |
| γ_MCA (membrane–cortex adhesion energy) | 10 (band 1–1000) | pN/µm (≈1e-5 J/m²) | ✅ SOURCED KU-3.B1.4 — ⚠️ **3-decade band** | `aleph/laws/membrane_surface.py:31` |
| `gamma_mem` (bilayer in-plane tension) | 10 (band 3–40) | pN/µm | ✅ SOURCED KB-3.B1.1 (Diz-Muñoz 2013) — NOT MCF7 cortical (correctly separated) | `compartments.py:112` |
| `kappa_m` (Helfrich κ) | 0.0828 (=20 kBT) | pN·µm (≈3.4e-19 J) | ✅ SOURCED KB-3.B1.2 (Rawicz 2000) | `compartments.py:113` |
| `k_erm` (ERM spring stiffness) | 100 | pN/µm | 🔴 **GAP / provisional** | `compartments.py:114` |
| ERM rest length | ≈0.10 µm gap (formation length) | µm | geometric (force-free) | `compartments.py:493` |

**PI flags for the blebbing deliverable:**
1. ⚠️ **γ_MCA 3-decade band** (1e-6–1e-4 J/m²) → `f_rupt` inherits that uncertainty. "Sourced but loosely constrained"; report as a band in validation, don't over-claim a sharp threshold.
2. 🔴 **`k_erm` is GAP** — does not enter the threshold `f_rupt` but sets the tension `k_erm·(L−rest)` compared against it, so it governs *when* nucleation fires. Source it before oracle B.
3. **Membrane param block carries no explicit PI-ratified stamp** (unlike the nucleus block). Recommend a PI ratification pass on the membrane block before the blebbing deliverable.
4. `RESERVOIR_STRAIN=0.60` is a flagged magic number but is DEFAULT-OFF for (caveolae-deficient) MCF7 — does not affect the current threshold; leave off.

## 3. Validation harness design (out-of-band; does not change the production runtime)

Two independent oracles, each driving the engine through a scripted protocol + post-processing (+ figures per
the visualize-at-closeout rule). Neither mutates the development loop.

### Oracle A — FSI-on audit (config + profiler gate) — **runnable now (I1a closed)**
Confirms two-way FSI is ON at physiological values, not stubbed.
- **A1** `div_vs` is live: `max|div_vs| > 0` after a contractile perturbation, and the null test
  **rigid translation → `div_vs ≈ 0`** (no spurious source).
- **A2** pressure back-reaction: `−α∇p` body force ≠ 0 on solid nodes under a pressure gradient.
- **A3** physiological setpoint: `storage_S / mobility / α / L_p / Π₀` equal the sourced values (not water
  defaults) — and (per §1 drift finding) equal the unit-converted yaml.
- **A4** GPU residency: zero authoritative GPU→CPU roundtrips inside the physical-time loop (profiler).
- Form: pytest-style gate + a few diagnostic launches. No physics run.

### Oracle B — bleb nucleation + ΔP non-equilibration (acceptance oracle) — **independent of the c_v anchor**
- **Protocol.** From the resting physiological baseline (Π₀=40 Pa, pre-tensioned cortex, ERM intact), apply a
  local perturbation mimicking Charras: (i) local membrane–cortex de-adhesion (drive local ERM tension past
  `f_rupt`), or (ii) a local myosin contraction spike raising local pore pressure.
- **Emergent observables → literature:**
  1. **Nucleation criterion** = pressure vs (adhesion `f_rupt` + tension `gamma_mem`) — must EMERGE from
     `f_rupt / gamma_mem / Π₀`, not be tuned (Charras 2008).
  2. **Pressure non-equilibration** = persistent ΔP gradient over **~10 µm / ~10 s** after nucleation
     (Charras et al. 2005) — THE poroelastic signature, and independent of `c_v` being an input.
  3. **Bleb dynamics** = expansion ~30 s, retraction ~2 min (Charras et al. 2008).
- **Acceptance:** observables 1–3 land in-band **without tuning**. Auto-render figures.

### Prerequisites / gating
| Oracle | Runnable now? | Blockers |
|---|---|---|
| A (FSI-on) | ✅ yes (I1a closed) | none |
| B (bleb) | ⚠️ partially blocked | 🔴 `L_p` ~62× reconciliation (§1) · 🔴 `k_erm` sourcing (§2) · `φ` porosity for cytosol inflow `v_f=v_s+q/φ` (I1b GAP) · membrane-block PI ratification |

## 4. Consolidated PI action list + decisions (2026-07-21)

**Oracle A: DONE** — static wiring/setpoint/drift audit (`test_fsi_on_audit.py`, 6 CPU) + runtime
dilatation-source sign/magnitude (`test_fsi_on_audit_runtime.py`, 2 CUDA); **8/8 verified on the gbook
A5000**. Confirms two-way FSI is wired (not stubbed) and at the physiological setpoint.

| # | Item | PI decision (2026-07-21) |
|---|---|---|
| 1 | `L_p` ~62× gap (code 1.6e-8 vs Jung-2011 draft 1e-6) | **Wait for an MCF7-specific L_p datum** — do NOT adopt the Jung draft now; keep code 1.6e-8 with the tracked invariant. Oracle B's **bleb-fill** part is deferred until the datum lands. |
| 2 | `k_erm` (GAP; KB has NO entry, confirmed) | **Literature-source it (in progress).** New KnowledgeClaim + SourceEvidence needed. Until it lands, oracle B nucleation-*timing* uses the provisional 100 pN/µm and any timing result is reported as provisional (the *threshold* `f_rupt` is already derived/sourced, so nucleation-*criterion* results are firm). |
| 3 | `φ` porosity (GAP; I1b/I1c only) | **Run oracle B pressure-only first** — nucleation criterion + ΔP non-equilibration (I1a, φ-free) now; the bleb-fill velocity `v_f=v_s+q/φ` (needs φ) is deferred to I1b. |
| 4 | Membrane block PI ratification | Pending — stamp together once `k_erm` sources and `L_p` datum land. |
| 5 | yaml drift | **RESOLVED** — `test_reconciled_constants_match_the_params_ledger` pins the reconciled constants; runtime yaml-load is a later hardening, not a blocker. |
| 6 | γ_MCA 3-decade band | Report the nucleation threshold as a **band**, not a sharp value, in oracle B. |

No parameter above may be tuned to make an observable land (Magic-Number Block + no-fitting rules).
**Unblocked now:** oracle B **pressure-only** (nucleation criterion + ΔP non-equilibration). **Still gated:**
oracle B bleb-fill dynamics (needs the `L_p` datum + `φ`).
