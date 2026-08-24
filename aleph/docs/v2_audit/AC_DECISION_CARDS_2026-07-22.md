# AC decision cards — PI-ratified 2026-07-22 (collagen / ERM / γ / sequencing / ac-cell retirement)

PI ratified the following five decisions (final line):
> **(c) conditional collagen gate / ERM production HOLD / re-register Hosseini MCF7 γ direct value /
> run only P0·P1·P2 / `ac/cell` staged strangler retirement.**

These are the SoT for the three parameter/gate decisions + the sequencing + the retirement policy. The
P2 lane (below) executes cards 1–3; P0 executes the physics; P1 the SF component. No P3–P7 coding sessions.

---

## Card 1 — Collagen modulus = (c) CONCENTRATION-RESOLVED conditional gate (contract correction, NOT loosening)

The conflict is a **missing-condition** defect, not a value disagreement: `5–100 Pa`
([ff_ecm_validate.py:37](../../scripts/ff_ecm_validate.py)) is a 1–3 mg/mL exploration envelope; `30–100 Pa`
([ff_ecm_library.py:111](../../ff/ecm_library.py)) is really the **3 mg/mL** range; the ~11 Pa at 1.5 mg/mL
is interpolated from `5 Pa @ 1 mg/mL` × `G′∝c^n`. **Do NOT pin `11 Pa` as a single-number gate.** Split into
conditional gates, each specifying temperature, frequency, strain amplitude, and 3D-bulk vs indentation:

- `VG-ECM-G0(c)` — small-strain storage/shear modulus `G′`, **per-concentration band**
- `VG-ECM-c-scaling` — `G′ ∝ c^n`, **n ≈ 2.0–2.1**
- `VG-ECM-stiffening` — differential modulus increase past a critical strain
- `VG-ECM-normal-stress` — `N1` sign and magnitude

Disposition: `5–100 Pa` → **exploration envelope** (demoted); `30–100 Pa` → **3 mg/mL validation**; 1.5 mg/mL
band → **P2 extracts it directly from the source figure** and registers it as the production gate. Historical
**11–15 Pa results = re-evaluation pending (neither PASS nor FAIL)** until the new contract is PI-approved.

⚠️ **Code inconsistencies P2 must fix (verified this session):** `ecm_library.py` carries `modulus_band_Pa
=(30,100)` AT `ref_conc=1.5 mg/mL` **and** `conc_exponent=0.5` — but the scaling target is `n≈2.0–2.1` (a
second, exponent-level conflict) and 30–100 belongs to 3 mg/mL. The `KB-1.30` note in code is **stale** → use
**KB-1.32 / KB-1.V.2.1**; re-check the Yang–Kaufman DOI in `ecm_library.py`. This is a contract correction
that restores the dropped concentration conditions — surfaced to PI as a gate-contract change, not edited to
pass.

## Card 2 — MCF7 ERM linker density = production HOLD (+ a scientific correction)

**⭐ Correction: `11.4 pN` is NOT the single-ERM rupture force.** `erm_rupture_force = 2π√(2κ_m(γ_mem+γ_MCA))`
is a **continuum membrane tube/tether-extraction force scale** (bilayer bending + membrane–cortex adhesion),
already isolated as diagnostic in code ([assemble.py](../../ac/cell/assemble.py),
[preload_contract.py](../../ac/cell/preload_contract.py)). **This session's resting-baseline diagnosis
mis-labeled it as the per-ERM cap** — see the correction folded into `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md`.
The single ezrin–F-actin measurement (Braunger et al. 2014, doi:10.1074/jbc.M113.530659, purified
reconstitution — NOT MCF7 density) is ~**stiffness 4.6 pN/nm** (this IS the code's `k_erm=4600 pN/µm`),
**unbinding force ~50 pN**, `k_off≈1.3 s⁻¹`, barrier width ~0.7 nm.

Policy:
- `ρ_ERM · F_single ≥ residual pressure` is a **necessary capacity gate ONLY** — do NOT adopt the back-computed
  minimum density as the physiological density.
- `1/node` and zebrafish `600 µm⁻²` stay **diagnostic-only**; production latch stays **false** until an MCF7
  direct value or a PI-approved epithelial proxy is registered.
- Approve as ONE bundle: areal density **+** bound/active fraction **+** Bell `k_on / k_off / F0 / capture`.
- Proxy priority if unavoidable: (1) MCF7 quantitative proteomics × cortical/membrane localization fraction;
  (2) another human epithelial direct areal-density measurement; (3) other-species embryo — diagnostic-only.

## Card 3 — MCF7 γ_cortex = re-register Hosseini 2020 direct MCF-7 measurement (SoT gap, not "no data")

Hosseini et al. 2020 (Adv Sci 7:2001276, doi:10.1002/advs.202001276) IS a direct MCF-7 measurement (suspended
interphase, dynamic AFM parallel-plate confinement, **n=27**): project record central ~**0.27 mN/m**, IQR
**0.18–0.40 mN/m**. **Verified this session:** the KB (KB-6.1.3) only CITES Hosseini 2020/2021 as
"figure-locked/scarce" with **no numeric γ extracted**; the current KU-3.5 target band is **0.35–0.65 mN/m**
(a modeling-contract target, `CORTICAL_TENSION_RECORD_2026-06-30`, NOT a Hosseini measurement) — Hosseini's
0.27 sits below/partly outside it, so registering it may inform/shift the gate.

Policy:
- Slice-1 is suspended/rounded interphase → **Hosseini 2020 is the primary γ candidate** (`0.18–0.40 mN/m`,
  central `0.27`); Hosseini 2021 (~`0.41 mN/m`) is an **upper cross-check** pending protocol reconciliation.
- P2 digitizes both papers' figures into a **dedicated KnowledgeClaim/gate**, comparing cell state, confinement,
  frequency, and effective-tension definition, and reconciles against the existing 0.35–0.65 KU-3.5 target.
- γ is a **physiological initial-state gate the explicit cortex must satisfy**, NOT a runtime lumped
  surface-force substitute.
- **⚠️ Double-count check:** if Hosseini's observable is *effective cortical/surface* tension, re-adding
  `γ_mem` double-counts. Fix to ONE reading:
  - Hosseini γ = total effective tension → `ΔP = 2γ_eff/R`
  - Hosseini γ = cortex-only → `ΔP = 2(γ_cortex + γ_mem)/R`
- First reading at `R=7.5 µm`: `ΔP ≈ 48–107 Pa` (central ~**72 Pa**). The current **40 Pa ⇒ γ_total≈0.15 mN/m**,
  below the Hosseini IQR → **40 Pa stays a HeLa diagnostic proxy only** ([params_i0b1.yaml](../../ac/fluid/params_i0b1.yaml)).

## Card 4 — Sequencing = ONLY P0 / P1 / P2 (no P3–P7 coding sessions now)

- **P0 — critical path (Lead, A5000, solo):** finite-difference-vs-analytic-operator comparison. At one fixed
  state/active-set, compare the true-force directional derivative to the operator action; identify the **first
  mismatching force family**; land the **minimal** fix. **No new solver/preconditioner until this result.**
  Then strict resting GO on a synced A5000 → R2 binding.
- **P1 — SF independent lane (isolated worktree):** SF as a separate component; do NOT touch `ac/cell`, the
  implicit solver, or the world scheduler; may reach CUDA-unit but NO connected/production claim until after P0.
- **P2 — evidence/SoT decision-card lane (expanded):** executes ALL THREE cards — collagen concentration-
  resolved gate, MCF7 ERM density/proxy hierarchy, Hosseini γ figure extraction + measurement mapping.
  **P2 replaces the old P3–P7 as the only non-P0/P1 lane;** P3–P7 (surface/fluid/ECM/NMII/MT-IF binding) are
  **deferred to their R-waves** after P0's native GO, not opened as sessions now.

## Card 5 — `ac/cell` retirement = gate-based 3-stage strangler (NOT a date). Freeze now.

1. **Now — feature freeze.** `ac/cell` is incumbent/reference; **no new biology**. Only P0 blocker fixes,
   instrumentation, and parity work.
2. **R2 + resting Slice-1 PASS →** production default entry switches to `ac/engine`; once SurfaceBody/FluidCore
   own real kernels and pass accepted/rollback/no-DtoH/native gates, `ac/cell` is **reference-only** for that slice.
3. **R4 PASS →** remove `build_cell` from the production CLI/config (after the `SF↔FA↔ECM+NMII` transactional
   vertical slice + 70,686 native gate).
4. **R8/R9 PASS →** move `ac/cell` to **source-only archive/quarantine** (block executable import via a static
   test) after the whole-cell native baseline + core perturbation matrix close.
