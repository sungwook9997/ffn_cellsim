# Parameter-provenance audit — CORTEX + actomyosin (GATE-A / GATE-B)

**Date:** 2026-07-24 · **Scope:** read-only audit of every load-bearing constant on the cortex-geometry,
actomyosin (NMII / ERM), and crosslinker/clutch parameter surfaces. **Trigger:** the PI caught
`density_per_fil=20` (`assemble.py:185` / `architecture_spec.py:92`) as a network-**percolation** spanning
threshold masquerading as biology (`CORTEX_MESH_FIDELITY_2026-07-24.md`). This memo finds the OTHER instances
of the same failure mode against the CLAUDE.md HARD rule: *every parameter + IC is the real in-vivo
physiological value, sourced, from the start — never a convenient null / spanning / placeholder.*

Classification: **SOURCED** (KB-id or a real DOI in code) · **DERIVED** (computed from a sourced quantity) ·
**CONVENIENCE** (percolation / round-number / backward-compat / numerical placeholder silently defaulted) ·
**PI-GAP** (acknowledged unsourced, source-gated, surfaced — an *honest* gap, not a silent violation).

---

## Headline

Besides `density_per_fil=20`, **three more load-bearing physical quantities are set to CONVENIENCE / proxy
defaults of the same kind** (a convenience value baked in as the default, not surfaced as a blocking gap):

1. **Turgor `PI_0_PA = 40` Pa** — the single most load-bearing baseline in the model.
2. **Cortex `seg_um = 0.5` µm** — caps the mesh 5–10× too coarse (the *real* cause behind `density_per_fil`).
3. **Cortex `length_um = 3.0` µm** — 3–30× too long, no Arp2/3 short-filament population at all.

Separately there are **~11 actomyosin/nucleus magnitudes that are unsourced but PROPERLY PI-GAP-gated** (the
code refuses to run without a supplied value + provenance, and the resting baseline leaves them out because the
heads are unbound). Those are honest gaps, *not* violations — but they are still unsourced and load-bearing for
the GATE-B **dynamic** result, so the γ magnitude carries them as a caveat. Of that gated set, **`k_xb` (the
"MASTER FORCE KNOB")** is the one that most directly sets the emergent γ.

**Top 5 to source next (worst-first):** ① turgor Π₀ (MCF7, not HeLa proxy) → ② cortex `seg_um` (mesh) →
③ cortex `length_um` + Arp2/3 population → ④ NMII `k_xb` → ⑤ `density_per_fil` (already flagged; downstream of ②).

The good news: the crosslinker / motor **kinetic** layer (`aleph/laws/hand_kmc.py`) and the fluid Biot closure
(`params_i0b1.yaml`) are genuinely well-sourced, and the ERM/NMII engine slices enforce source-gating rigorously
(`_require_param`, `ERMBellKinetics`, `RestingBoundMyosinSetpoint` all refuse to default a GAP).

---

## Table (sorted worst-first — CONVENIENCE on a load-bearing physical quantity = worst)

| # | Parameter | Location | Value | Class | KB / source | Physiological concern |
|---|---|---|---|---|---|---|
| 1 | **Turgor Π₀** | `assemble.py:128` `PI_0_PA` = 40; `fluid/params_i0b1.yaml` Q3 | **40 Pa** | **CONVENIENCE (proxy)** | Fischer-Friedrich 2014 **HeLa** (KB-6.1.6 / KB-DRAFT-3.B-30) — **no MCF7 datum** | **Worst.** Sets the ENTIRE resting cortical tension γ = ΔP·R/2. MCF7-geometry central estimate ~**72 Pa** (48–107); 40 Pa is ~45% low and **below the Hosseini IQR** (`AC_DECISION_CARDS` §, `AC_ENGINE_INTEGRATION_PLAN` L276). "Sourced" only for the wrong cell — exactly the unphysical-baseline → meaningless-comparison failure the HARD rule names. |
| 2 | **Cortex `seg_um`** | `architecture_spec.py:90`; `assemble.py:183` (`cortex_seg_um`) | **0.5 µm** | **CONVENIENCE (coarse)** | none (0 KB rows) | Caps mesh/pore at ~500 nm vs sourced **50–100 nm** (Morone 2006 / Bovellan 2014 / Chugh 2018). **5–10× too coarse.** The deeper cause of the mesh-fidelity gap — `density_per_fil` is downstream of this. |
| 3 | **Cortex `length_um`** | `architecture_spec.py:90`; `assemble.py:184` | **3.0 µm** | **CONVENIENCE (coarse)** | none | Representative filament vs physiological ~1 µm formin / ~0.1 µm Arp2/3 → **3–30× too long**; `nucleator="formin"` only ⇒ **no Arp2/3 short-filament population** (a missing architecture, not just a number; Bovellan cortex is ~⅓ Arp2/3). |
| 4 | **`density_per_fil`** (ref: already caught) | `architecture_spec.py:92`; `assemble.py:185` | **20.0** | **CONVENIENCE (percolation)** | none (docstring *claims* DERIVED L/δ_filamin; operational note says "keeps the network single-spanning") | Chosen as a spanning threshold, not sourced. ≈1 crosslink / 150 nm contour, on the sparse side. Note: the notes' "DERIVED" story and the "single-spanning" story disagree — a self-inconsistent justification. |
| 5 | **NMII `k_xb`** | `assemble.py:146` `NMII_K_XB_TEST`; `params_i0b3.yaml` k_xb | **1.0e3 pN/µm** | **PI-GAP (master knob)** | GAP — physical band 100–1000 | ⚠ "MASTER FORCE KNOB" — sets head strain F/k_xb and directly scales GATE-B emergent γ. Picked at the **top** of the band. Properly GAP-flagged, but load-bearing + unsourced. |
| 6 | **NMII `f_stall`** | `assemble.py:149`; `params_i0b3.yaml` F_stall_head | **0.5 pN** | **PI-GAP (conflicted)** | GAP — claim_a 0.5 (AFINES/Kovács) vs claim_b 2.0 (Billington) | Per-head stall; unbound at rest so out of GATE-A, enters GATE-B ensemble stall. 4× isoform ambiguity unresolved. |
| 7 | **NMII `kappa` (Hill)** | `assemble.py:151`; `params_i0b3.yaml` kappa_hill | **0.5** | **PI-GAP (3-way)** | GAP — 0.25 (Hill muscle) / 0.5 (Kovács NMIIA) / ∞ (2026-07-07 linear) | Hill curvature; three unreconciled positions. |
| 8 | **NMII `v0`** | `assemble.py:150`; `params_i0b3.yaml` v0 | **0.12 µm/s** | **PI-GAP (conflicted)** | GAP — 0.12 (KB-PIV-4) vs 0.2 (Kovács 2003) | Unloaded stepping velocity; assay/isoform ambiguity. |
| 9 | **NMII `k_on`** | `assemble.py:152`; `params_i0b3.yaml` k_on | **50 /s** | **PI-GAP** | GAP — config-chosen (phase1_h3), **not an audited SE row** | Sets emergent engaged fraction with k_off0. |
| 10 | **NMII `N_side`** | `assemble.py:143` `NMII_N_SIDE` | **10** | **PI-GAP (conflicted)** | GAP — 10 (AFINES) vs 28–30 (Billington 2013 EM) | Heads/side; archived production used 28. Drives head-count → total motor force. |
| 11 | **NMII `r0_head`** | `assemble.py:144` `NMII_HEAD_OFFSET_UM` | **0.200 µm** | **PI-GAP** | GAP — archived 200 nm (band 0–0.21) | Arm rest offset; rests force-free at build so no prestress. |
| 12 | **NMII capture radius** | `assemble.py:177` `NMII_CAPTURE_UM` | **0.210 µm** | **DERIVED** | ⚠ **CORRECTED 2026-08-22 — this row described a defect that had already been repaired ON THE SAME DAY it was written.** It read *"0.05 µm, CONVENIENCE, inconsistent with `hand_kmc` NMIIA 0.210 — two different capture radii for the same head in the same codebase. Reconcile."* `_historical/PARAM_INCONSISTENCY_RECONCILE_2026-07-24.md` reconciled it 0.05 → 0.210 that day, and the constant has carried a provenance comment saying so ever since. This audit is in `STATE.md`'s minimal read set, so the stale row is not inert: it sent a session on 2026-08-22 to the PI reporting a live 4.2x contradiction that did not exist, and to claim PI queue 14's pass/fail line depended on resolving it. Verified by grep before correcting — the whole tree is on one value; `hand_kmc`'s 0.060 (α-actinin, filamin), 0.300 (integrin α5β1) and 0.010 (kinesin) are OTHER MOLECULES, not disagreements. ⚠ Marked rather than rewritten: a dated audit that silently acquires today's answer stops being a record. | 0.210 = r0_head (0.200) + ~10 nm slack; KU-3.5 binding fix. **No open action.** |
| 13 | **NMII `L_bb`** | `assemble.py:143` `NMII_L_BB_UM` | **0.301 µm** | **PI-GAP (well-anchored)** | GAP by isoform, but band = Billington 2013 EM (301±24 nm) | Backbone contour; value itself is on a real EM anchor — closing it is low-risk. |
| 14 | **NMII `k_backbone_factor`** | `assemble.py:148`; `params_i0b3.yaml` | **10× (→1e3 pN/µm)** | **CONVENIENCE (numerical)** | explicitly "numerical (rigid-rod stiffening), NOT a literature constant" | Declared numerical; grid-invariance-blocked. Acceptable but note the per-bond-vs-end-to-end grid-invariance OPEN item (assemble.py:527–534). |
| 15 | **NMII `n_bb`** | `assemble.py:142` | **14** | **CONVENIENCE (numerical)** | archived topology-discretization; count-invariant gate holds for any n_bb≥2 | Declared numerical. |
| 16 | **NMII `k_head_arm`** | `assemble.py:147` `NMII_K_HEAD_ARM_TEST` | **1e2 pN/µm** | **CONVENIENCE (provisional)** | provisional; arm rests force-free (~0 at build) | Out of resting result; provisional. |
| 17 | **ERM `k_erm`** | `compartments.py:123`; `erm_cortex_slice` | **4.6e3 pN/µm** | **SOURCED (nuance)** | Braunger 2014 JBC single-bond, PI-ratified 2026-07-21 (KB-3.B1.6 reg pending) | ⚠ Provenance nuance: Braunger's *most-probable* is **2.2 pN/nm**; 4.6 pN/nm comes from the dissertation event-statistics analysis — a ~2× interpretation choice baked in. |
| 18 | **ERM Bell `k_on/k_off0/F0/capture`** | `ac_gate_b_erm_cortex_native.py:89–92` | 1.0 /s / 0.3 /s / ~4.28 pN / 0.05 µm | **PI-GAP (proxy)** | GAP — no MCF7 ERM on/off datum (2026-07-21 audit) | Loud PROVISIONAL banner; `ERMBellKinetics`/`ERMSliceParams` refuse to default them. Mechanism-shape only, not quantitative. Model-hygiene exemplary; values still unsourced. |
| 19 | **`sigma_EV`** | `assemble.py:136` `SIGMA_EV_UM`; `params_i0b2b.yaml` | **0.007 µm** | **PI-GAP (discretization)** | GAP — physical (7 nm F-actin) known, coarse-graining open | At 0.5 µm node spacing a 7 nm WCA never triggers ⇒ steric effectively inactive at rest. Faithful value, unfaithful discretization. |
| 20 | **`k_EV`** | `assemble.py:137` `K_EV_PROVISIONAL` | **1e3 pN/µm** | **CONVENIENCE (numerical)** | explicitly "NOT a literature value" (Magic-Number Block) | Declared numerical repulsion scale; to be closed at native gate from measured passive F_op, not tuned. |
| 21 | **`steric_force_cap`** | `assemble.py:261` | **1e3 pN** | **CONVENIENCE (numerical)** | declared CFL guard (>> ~0.02 pN physical load) | Declared numerical guard, not the steric verdict. Acceptable. |
| 22 | **NMII `k_off0`** | `assemble.py:153`; `params_i0b3.yaml` | **0.35 /s** | **SOURCED (draft)** | Stam-Hocky 2015 / Tam 2021, KB-3.18 (draft/Medium) | Provisional-sourced; shape gate. |
| 23 | **NMII `f0`** | `assemble.py:154` `NMII_F0` | **≈7.13 pN** | **DERIVED** | = kBT/x_β, x_β=0.6 nm Veigel 2002 (KB-3.18) | Sourced-derived. |
| 24 | **α-actinin `k_off0`** | `hand_kmc.py` ALPHA_ACTININ | **0.066 /s** | **SOURCED (reg pending)** | Ferrer 2008 PNAS (real DOI in code; TAG shows only a filamin-matched draft row — **Contract-Graph registration PENDING**) | Real standard paper; SE row pending (SE_REGISTRATION §5). Not a fabrication. |
| 25 | **α-actinin / filamin `link_k`** | `hand_kmc.py` | **4.6e5 / 8.2e5 pN/µm** | **SOURCED (reg pending)** | Ferrer 2008 PNAS AFM (455 / 820 pN/nm), PI-approved 2026-06-30 | The old 0.1 was the broken AFINES soft surrogate (~4.5e6× too soft) — now fixed. SE rows pending. |
| 26 | **FA clutch (integrin α5β1)** | `hand_kmc.py` INTEGRIN_A5B1 | k_catch0 0.4, Fc 7, k_slip0 0.5, Fs 30 pN | **SOURCED (flag)** | Kong 2009 / KB-2.5 (k_on, k_int PI-gated) | ⚠ KB-consistency flag: recorded params give lifetime-peak F*≈7 pN, not the "F*≈30 pN" the KB-2.5 text claims — used as-recorded, surfaced to PI. Not in the cortex path. |
| 27 | **`R_CELL_UM`** | `assemble.py:117` | **7.5 µm** | **SOURCED** | Wagner 2011 MCF7 Coulter (PMC3147247), reg 2026-07-21 | MEASURED MCF7 suspended radius. Good. |
| 28 | **`CORTEX_MEMBRANE_GAP_UM`** | `assemble.py:126` | **0.10 µm** | **DERIVED** | ½·h_cortex, h_cortex≈200 nm (KB-3.1/3.5) | Half-thickness split not separately resolved; mechanically negligible at rest. |
| 29 | **`n_filaments`** | `assemble.py:173`; `architecture_spec.py:90` | **70686** | **DERIVED (generic)** | = 100 µm⁻²·4π(7.5)²; areal density 100/µm² is KB-3.18 but **generic, non-MCF7** | Order-correct (lands in 30–50 nm mesh regime). Length gap (row 3) is the real issue, not the count. |
| 30 | **Membrane `gamma_mem`** | `compartments.py:121` | **10.0 pN/µm** | **SOURCED** | KB-3.B1.1 (lit-anchored plateau) | In-plane bilayer tension. Good. |
| 31 | **Membrane `kappa_m`** | `compartments.py:122` | **0.0828 pN·µm** | **SOURCED** | 20 kBT, KB-3.B1.2 | Helfrich bending. Good. |
| 32 | **Membrane `f_rupt`** | `compartments.py`; `assemble.py:906` | (derived) | **CONVENIENCE (diagnostic)** | membrane-tube extraction scale — explicitly NOT single-ERM | Correctly labeled `DIAGNOSTIC_...NOT_SINGLE_ERM`; never authorizes production capacity. Good hygiene. |
| 33 | **`membrane_subdivisions`** default | `assemble.py:199` | **3** (prod = 8) | **CONVENIENCE (backward-compat)** | numerical resolution (grid convergence) | Default is the coarse γ-validation sphere; production RESTING = 8 (655k verts). Documented; only a concern if a run forgets to raise it. |
| 34 | **Nucleus `k_lamin_b / k_lamin_ac`** | `compartments.py:107–108` | 15 / 150 pN/µm | **PI-GAP (TEST)** | I0-B2 PI GAP, TEST values | ε=0 ⇒ σ=0 at rest, so OUT of the resting baseline; bite only when the nucleus is loaded. |
| 35 | **Nucleus `k_vol / k_linc / knee_strain`** | `compartments.py:112–113,109` | 1e3 / 1e2 / 0.10 | **PI-GAP (TEST)** | GAP; "8 pN is a TENSION not a stiffness" (k_linc) | Same as row 34 — flag for the loaded-nucleus regime. |
| 36 | **`L_p`** | `assemble.py:133` | **1.6e-8** | **CONVENIENCE (inconsistent)** | "gamma_floor code value; draft/MCF7"; Jung 2011 | ⚠ Code value 1.6e-8 µm/(s·Pa) does **not** match `params_i0b1.yaml` L_p = 1e-12 m/(s·Pa) = 1e-6 µm/(s·Pa). Not load-bearing at rest (uniform p, p_ext balance), but the two surfaces disagree — reconcile. |
| 37 | **Biot `c_v`** | `params_i0b1.yaml` | **50 µm²/s** | **SOURCED** | Moeendarbary 2013 Nat Mater (verified/High) | Primary poroelastic anchor. Good. |
| 38 | **Biot `M / storage_S / mobility`** | `assemble.py:130–131`; yaml | 1e4 Pa / 1e-4 / 5e-3 | **DERIVED** | M = c_v·µ/k (consistency, not independent) | Correctly derived, avoids M=K_drained double-count. Good. |
| 39 | **Biot `alpha`** | `assemble.py:132` | **1.0** | **DERIVED (ratified)** | PI-ratified 2026-07-16 incompressible-constituent limit (Penta 2025 etc.) | Universal modeling limit; documented. Good. |
| 40 | **Straddle windows / grid dims / dx / seed** | `assemble.py:159–160,257–260` | — | **CONVENIENCE (numerical)** | build-time geometric selection / accelerator sizing | Not physics magnitudes; grid-invariant accelerators. Acceptable. |

---

## Reading of the results

**Genuine violations (silent convenience default on a load-bearing physical quantity — the `density_per_fil`
failure mode):** rows **1 (turgor 40 Pa), 2 (seg_um), 3 (length_um), 4 (density_per_fil)**, plus the softer
**12 (NMII capture inconsistency)** and **36 (L_p inconsistency)**. Rows 1–3 are strictly *worse* than
`density_per_fil` because they set the resting-baseline tension and the cortex mesh/length directly, and they
are silently defaulted, not surfaced as blocking gaps.

**Honest, well-gated gaps (unsourced but NOT violations):** rows 5–11, 13, 18, 19, 34–35. The engine enforces
these rigorously — `_require_param` (`cortex_motor_slice`, `erm_cortex_slice`), `ERMBellKinetics`,
`RestingBoundMyosinSetpoint`, and the `CellConfig` ERM/backbone latches all **refuse to construct** without a
supplied value + provenance. This is exactly the discipline the HARD rule wants; they are load-bearing for the
GATE-B *dynamic* γ and should be closed, but they are not hiding.

**Genuinely well-sourced:** the entire `aleph/laws/hand_kmc.py` kinetic layer (Bell/Pereverzev laws, α-actinin/filamin
`k_off0` + `link_k` from Ferrer 2008, NMIIA Bell from Veigel/Kovács/Stam-Hocky), the fluid Biot closure
(Moeendarbary c_v anchor + derived M/S/mobility + ratified α), the membrane γ_mem/κ_m, and R_CELL. Caveat: the
Ferrer 2008 α-actinin/filamin rows and the Braunger k_erm are **DOI-real in code but Contract-Graph registration
is still PENDING** (SE_REGISTRATION §5) — not fabrications, but `tag_query` will not return them yet.

## Caveat on the GATE-A / GATE-B milestones

GATE A (converged resting baseline) and GATE B (emergent dynamic tension) were validated on this **40 Pa /
0.5 µm-mesh / 3 µm-filament** cortex. They correctly establish solver convergence and the dynamic *mechanism*
(tension emerges from binding events), but the **γ magnitude** inherits rows 1–4 as a coarse/proxy baseline —
on top of the already-flagged PI-GAP rate constants. This should be stated in the milestone reports.

## Recommendation (PI decision — not done here, read-only audit)

Source, worst-first: **① MCF7 turgor Π₀** (retire the HeLa 40 Pa proxy; the ~72 Pa MCF7-geometry estimate is
already computed) → **② cortex `seg_um`** (mesh 0.5 → ~0.075 µm, per `CORTEX_MESH_FIDELITY` §7) → **③ cortex
`length_um` + add the Arp2/3 short-filament population** → **④ NMII `k_xb`** (the master force knob, by
isoform/assay) → **⑤ `density_per_fil`** (falls out of ② once the mesh is sourced). Also reconcile the two
NMII capture radii (row 12) and the L_p code-vs-yaml mismatch (row 36), and finish registering the Ferrer 2008
/ Braunger 2014 SE rows so the KB reflects what the code already uses.
