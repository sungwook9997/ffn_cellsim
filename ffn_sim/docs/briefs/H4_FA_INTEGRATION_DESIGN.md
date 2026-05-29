# H.4 → Cell Focal-Adhesion Integration — Design (α restart contract)

> **Status**: DESIGN — PI-ratified path **α** (2026-05-30: keep γ-Phase 1 + integrate FA into `Cell.build`). Implementation contract for the restart. Derived from the 33-agent mechanism audit [`../v2_audit/MECHANISM_AUDIT_2026-05-30.md`](../v2_audit/MECHANISM_AUDIT_2026-05-30.md) §6.
> **Author**: Lead session, 2026-05-30 (autonomous /loop).
> **Architectural principle**: every link is an explicit HOOMD particle/bond with mechanistic kinetics. Chan-Odde, Bell-Evans, Hill, Kong, Riveline/Balaban closed-forms are **acceptance oracles only**, never the runtime mechanism.

## 1. Problem (root cause)

The integrated cell cortex is a **floating shell**: `bridge/` (FA) and `ecm/` (Mikado substrate) are implemented + unit-tested but have **zero call sites in `cell/cell.py`**. `with_fa` (cell.py:584) is a dead no-op; `offsets['fa_integrin']=(cur,cur)` (cell.py:902) hard-zeros the FA tag slot. Consequence:

- **KU-3.5 cortical tension**: cortex has nothing to pull against → myosin stress integrates to ≈0 net tension → γ_total ≈ 1.4×10⁻⁴ mN/m vs target [0.35, 0.65] (~3,000× under).
- **KU-5.1 dendritic density**: lamellipodial barbed ends carry no load → `exp(−Fδ/kT)` capping stays maximal → density ≈ 0.012 vs ~100 /µm² (~10,000× under).

Both failures share one cause: **no mechanical ground + no load path**.

## 2. The force-transmission chain (target architecture)

```
ECM substrate / fixed ligand layer at z=0   (immobile = mechanical ground)
        │  ← integrin–ligand CATCH bond  (Kong 2009 two-state, F*~30 pN)
        ▼
  talin rod   (13-domain serial WLC, per-domain two-state unfold/refold; R3 ~5 pN; +60–80 nm/unfold → exposes VBS)
        │  ← vinculin  (discrete force-gated recruitment to exposed VBS; vinculin–actin DIRECTIONAL catch)
        ▼
  molecular clutch   (cortex / stress-fiber actin → engaged integrin; loaded by myosin retrograde flow)
        ▼
  cortex + lamellipodium actin   ←ERM clutch→   membrane   (Brownian-ratchet barbed-end load; γ_mem ~10–300 pN/µm)
```

Every arrow is an explicit bond population with force-dependent kinetics. The cortical tension and the lamellipodial load both *emerge* from closing this loop.

## 3. Build sequence (wired behind `with_ecm` + `with_fa` in `build_cortex_full_simulation`)

Staged so each stage is independently testable; **do not** jump to KU-3.5 v4 until the full loop is closed + equilibrated. All stages land as **additive, default-off** scaffold first (no regression to current cortex/myosin/lamellipodium behavior), then are turned on for v4.

| Stage | Element | Module(s) | New topology | Oracle |
|---|---|---|---|---|
| S0 | **Substrate anchor at z=0** — fixed/position-restrained ligand layer; integration-group exclusion or stiff pin to initial pos. Density per FA from KU-2.4 (~50 integrins/FA). | new `with_ecm` path in cell.py; reuse `ecm/mikado` or a minimal fixed ligand sheet | `substrate_ligand` particles (immobile) | — |
| S1 | **Integrin–ligand catch bond** engaging cortex actin to anchor. Reuse `IntegrinBondUpdater`. Migrate Pereverzev → Kong 2009 two-state (F*~30 pN). | `bridge/integrin_bonds.py`, `bridge/fa.py::build_h4_state` | `integrin` particles + `integrin_ligand` catch bonds | Kong 2009; Chan-Odde |
| S2 | **`fa_actin_clutch` bond** — NEW topology: cortex/stress-fiber actin → engaged integrin. This is the literal load path; no current module provides it. | new bond type registered in cell.py builder | `fa_actin_clutch` bonds | Chan-Odde traction-vs-stiffness |
| S3 | **Talin multi-domain WLC** between integrin and clutch/actin; per-domain Bell-Evans unfold + refold; R3 first ~5 pN; +60–80 nm contour/unfold. Activate dead `talin` config in `ResolvedH4`. | `bridge/fa.py` (talin block), new updater | talin WLC bonds + unfold-state tag | Yao 2016; Tapia-Rojo 2020 |
| S4 | **Vinculin force-gated reinforcement** — discrete vinculin binds force-exposed VBS (gated by talin unfold state); vinculin–actin directional catch. Replace forbidden scalar `k_eff=k_bare·(1+α·N_vin)`. Activate dead `vinculin` config. | `bridge/fa.py` (vinculin block), new updater | vinculin particles + directional catch bonds | Huang 2017; Elosegui-Artola 2016 |
| S5 | **Close the clutch loop** — myosin (Stam-Hocky, already wired) pulls cortex actin into retrograde flow; engaged integrin/talin clutches resist → network tension = γ. Integration of S0–S4 + existing myosin, not new physics. | cell.py builder | — | Chan-Odde |
| S6 | **Substrate compliance spring** (E_substrate) — anchor each ligand through a harmonic spring whose k encodes E (Hertz/Bangasser-Odde scaling, Magic-Number Block). Add 0.1–100 kPa sweep driver. | new param + sweep driver | `substrate_spring` bonds | Chan-Odde; Bangasser-Odde 2013 |
| S7 | **Explicit membrane object** — deformable surface or tension/area reservoir (κ_m~20 kT, γ_mem~10–300 pN/µm) the barbed ends push against; cortex coupled via ERM. The KU-5.1 fix. | new `cell/membrane.py` | membrane mesh / reservoir | Bieling 2016 f-v |
| S8 | **ERM → Bell-Evans slip clutch population** (not static pin) — per-linker k=1e-4 N/m (Alert 2015), force-dependent slip off-rate (Korkmazhan-Dunn 2022), ρ~100/µm², gated rebind. Reuse `IntegrinBondUpdater`. | `cortex/erm.py` rework | ERM clutch bonds | Alert 2015 |
| S9 | **Integration + equilibration prelude** — register all new bond/WCA params; extend BAOAB `gamma_map`; replace zero-width `fa_integrin` tag slot with real range; run `equilibrate_no_shear` BEFORE the BAOAB production loop. | cell.py builder + `ecm/equilibrate.py` | — | — |

## 4. Co-requisites (must accompany the wiring, not follow it)

The completeness critic established that a turnover-free / pressure-free / membrane-free cell **misses the gates even with FA wired**:

- **Actin turnover** (cofilin severing + pointed-end depoly, τ ~5–30 s; Chugh 2017) — the single dominant determinant of steady-state cortical tension. Add a force/age-aware disassembly Updater to cortex AND lamellipodium. **KU-3.5 unreachable without it.**
- **Enclosed-volume / osmotic pressure** so ΔP = 2γ/R emerges (Stewart 2011; Fischer-Friedrich 2014) — the actual KU-3.1 mechanism. ⚠️ surfaces a gate question: metaphase γ~1.6 mN/m exceeds the KU-3.5 ceiling (1.0) → KU-3.5 may be interphase-biased.
- **Membrane object** (S7) — without it the barbed-end load law is inert (KU-5.1).

## 5. Acceptance (re-run only after the loop is closed + equilibrated)

- **KU-3.5**: γ_total ∈ [0.35, 0.65] mN/m (method-of-planes, γ_soft + γ_rigid). No band change.
- **KU-5.1**: plateau dendritic density → ~100 /µm² (band PI-pending, was placeholder [50,150]).
- **Chan-Odde validation** (S6): traction-vs-E biphasic, optimum ~1 kPa; retrograde flow monotonically ↓ with E. Oracle-only comparison.
- **No magic numbers**: every per-link constant carries a Magic-Number Block from a primary source (see §7 of the audit for the KAIST fetch list).

## 6. PI sign-off items (staged, see `../v2_audit/PI_DECISION_QUEUE_2026-05-30.md`)

This design is additive code; but it interacts with three governed surfaces:
- **frozen-integrator**: the global `dt = min(τ)` reconciliation + `equilibrate_no_shear` placement (S9), and whether the new stiff bonds need constrained-BD.
- **param-change**: integrin Pereverzev → Kong 2009 (F*~30 pN); the per-link FA parameters.
- **gate-contract**: the KU-5.1 band ratification; the KU-3.5 interphase-vs-metaphase question.

Implementation proceeds on the additive scaffold; these three are gated on PI.
