# Ligand-identity FA extension — ready-to-implement SPEC (col-I α2β1 / laminin-111 α6β1)

> **STATUS: SPEC — 2026-06-02. Prep, NOT wired.** Produced in a parallel session
> (`AUTONOMOUS_LOG_2026-06-02.md`); `bridge/fa.py` NOT modified. Implements the Notion
> open item *"sequenced after FA production lands → extend `bridge/fa.py` ligand
> identity"* and PI-exp map §Layer-1. Ships with a standalone, tested data module
> `bridge/ligand_species.py` (import-isolated) the Lead wires in. PI ratified the
> Lam4 = **α6β1** receptor + the affinity+proxy recipe (PI-exp map open items, 2026-06-02).

## Goal

The PI experiment presents **three** ligand conditions on the same surface; the
platform currently has **one** generic integrin↔ligand clutch (FN-like α5β1
catch-slip, KU-2.18). Give the clutch a **ligand identity** so col-I, laminin-111,
and FN are distinct presets, validatable against the per-ligand traction contrast.

| PI condition | Ligand | Integrin | Bond class |
|---|---|---|---|
| Bare / Pre | collagen-I (coated/adsorbed) | **α2β1** (GFOGER) | **slip** (Bell-Evans) |
| Lam4 | laminin-111 (soluble supplement) | **α6β1** (low-affinity LN-111) | **slip (proxy)** |
| (platform baseline / benchmark) | fibronectin | α5β1 (RGD) | **catch-slip** (KU-2.18, exists) |

## The key design fact — ligand identity is a `PereverzevParams` swap, not new mechanism

The runtime off-rate is `IntegrinBondUpdater` → `pereverzev_k_off(F, PereverzevParams(k_s, F_s, k_c, F_c))`
(`bridge/integrin_bonds.py:102,272`). The two-pathway model is
`k_off(F) = k_s·e^{F/F_s} + k_c·e^{−F/F_c}` (slip term + catch term).

- **A pure slip bond is the degenerate case `k_c = 0`.** Bell-Evans slip
  `k_off(F) = k_off0·e^{F·x_β/kT}` maps exactly onto it with
  **`k_s = k_off0`, `F_s = k_BT/x_β`, `k_c = 0`**. (`pereverzev_F_star` correctly
  returns **NaN** for `k_c=0` — a slip bond has no catch peak.)
- **FN α5β1 keeps the full catch-slip** = the platform's existing KU-2.18 default
  (`k_s=0.5 s⁻¹, F_s=30 pN, k_c=0.4 s⁻¹, F_c=7 pN`).

So **no runtime-mechanism change**: ligand identity = selecting the per-species
`PereverzevParams` (+ `k_on`, capture). `bridge/ligand_species.py` provides these
literature-anchored sets + the `k_BT/x_β` conversion + provenance + proxy flags.

## Literature-anchored parameter sets (acceptance values, NOT fitting targets)

`F_s = k_BT/x_β` at **T = 310.15 K** (k_BT = 4.2816×10⁻²¹ J). x_β were measured by AFM
SMFS at RT (~25 °C, k_BT = 4.116×10⁻²¹ J) — a ~4 % systematic on F_s, inside the bands.

| Species | Integrin | model | k_off0 (k_s) | x_β | ⇒ F_s = kT/x_β | k_on | F* / lifetime | Source | Conf | proxy? |
|---|---|---|---|---|---|---|---|---|---|---|
| **col-I (single-cell)** | α2β1 | slip | **1.3 s⁻¹** | **0.23 nm** (2.3 Å) | **18.6 pN** | (generic integrin) | lifetime 0.8±0.7 s; rupture 47±13 pN @500 pN/s | Taubenberger 2007 *MBoC* 18:1634 | high (3-0) | no |
| **col-I (locked-open SMFS)** | α2β1 | slip | **0.44 s⁻¹** | **0.7 nm** | **6.1 pN** | (generic) | (GPP ctrl 11 s⁻¹/0.37 nm) | PMC3588017 (AFM SMFS) | high (3-0) | no |
| **laminin-111** | α6β1 | slip (proxy) | **1.4–2.3 s⁻¹** | **0.28 nm** | **15.3 pN** (≈ reported f_b 14.1±1.3 pN ✓) | (generic) | single-barrier slip | α7β1–**invasin** PMC3882471 (NOT laminin) | med (proxy) | **YES** |
| **FN (benchmark)** | α5β1 | catch-slip | 0.5 s⁻¹ (k_s) | — (F_s=30 pN) | 30 pN | KU-2.4 | catch 10–30 pN; peak 2–10 s @20–25 pN | KU-2.18 / Kong 2009 PMC2712956 | high (3-0) | no |

**Recommended defaults:**
- **col-I → single-cell anchor** (Taubenberger: k_off0=1.3 s⁻¹, x_β=0.23 nm) — it is
  the *cell-context* measurement (whole-cell, in the loading regime the clutch sees),
  more representative than the isolated locked-open I-domain SMFS. Keep the locked-open
  set as a documented variant.
- **laminin → α7β1-invasin shape proxy** (x_β=0.28 nm, k_off0=1.4–2.3 s⁻¹), **scaled
  weaker than FN**: laminin is slip-only (no catch stabilization) AND k_off0 (1.4–2.3)
  ≫ FN k_s (0.5) ⇒ shorter-lived, weaker clutch ⇒ **lower traction** (the lit anchor:
  breast-epithelial traction is lower on laminin-111 than col-I/FN; P=0.016–0.028).
  Every laminin value carries a `proxy=True` flag (direct α6β1–LN-111 force kinetics
  are **genuinely absent** from the literature — PI-exp map §dossier `wf_1494c786-e79`).

**Affinity context (NOT off-rate):** α6β1–LN-111 K_D ≈ 1–20 nM (toward the weak end;
α6β1–LN-511/521 = 0.73 nM). LN-111 is α6β1's **3rd-rank** laminin (LN-511/521 > LN-332
> LN-111). Use the K_D only to set relative on-rate/availability if `k_on` is later
made ligand-specific — not as an off-rate. Cheap follow-up: recover the exact
α6β1–LN-111 K_D from paywalled Nishiuchi 2006 via gbook/KAIST.

## The `bridge/ligand_species.py` data module (ships with this spec — tested, isolated)

Provides, per species: the bond class, `PereverzevParams` (via the `kT/x_β`
conversion), `k_on`, source/conf/proxy. Import-isolated (no existing file imports it
⇒ runtime + suite unaffected). See `tests/test_ligand_species.py` (GREEN in isolation).

```python
from aleph.bridge.ligand_species import (
    LIGAND_REGISTRY, pereverzev_params_for, DEFAULT_LIGAND_FOR_CONDITION,
)
p_col1 = pereverzev_params_for("collagen_I")      # slip-only PereverzevParams (k_c=0)
p_lam  = pereverzev_params_for("laminin_111")     # slip-only proxy
p_fn   = pereverzev_params_for("fibronectin")     # full catch-slip (KU-2.18)
```

## Wiring into `bridge/fa.py` (Lead sequencing step — exact attach points, NOT done here)

`ResolvedH4` already holds a `pereverzev` field built in `resolve_h4` from the
`catch_bond` config block (`bridge/fa.py:254,299`; `pereverzev_F_star/lifetime` at
`:299-302`). The minimal, additive extension:

1. **Config**: add an optional `ligand:` block to the H.4 config —
   `ligand: { species: collagen_I }` (default `fibronectin` ⇒ today's behavior
   bit-for-bit). PI-exp preset (b) sets `species: collagen_I`.
2. **resolve_h4**: if `ligand.species` is set and ≠ fibronectin, replace the
   `pereverzev` field with `pereverzev_params_for(species)` (and the per-species
   `k_on`); else keep the `catch_bond`-derived FN default. One branch, default-off.
3. **No updater change** — `IntegrinBondUpdater` consumes whatever `PereverzevParams`
   it is handed. The slip species (k_c=0) flow through `pereverzev_k_off` unchanged;
   `catch_peak_force`/`_lifetime` become NaN for slip (correct — no catch peak), so
   guard the KU-2.5 catch-peak gate to **skip for slip ligands** (it is a catch-only
   gate). This is the one downstream gate that must learn "slip ligands have no peak."
4. **β1-distribution observable** (PI-exp Layer-1 TODO): a diagnostic over engaged-
   integrin spatial spread (diffuse/peripheral/uniform) — a read-only measurement on
   the existing clutch tags, addable later; not required for the traction contrast.

## Validation acceptance (candidate gates — Layer-1)

| Gate | Criterion | KU / source |
|---|---|---|
| col-I slip kinetics | emergent k_off(F) matches Bell slip `1.3·e^{F·0.23nm/kT}`; single-bond rupture in 38–90 pN over the loading-rate range | Taubenberger 2007 |
| laminin weaker than FN | per-ligand traction(laminin) < traction(col-I/FN) (lower bound-integrin lifetime ⇒ lower transmitted force) | rel-traction anchor (P=0.016–0.028) |
| FN catch peak preserved | FN keeps catch peak F*≈10–30 pN, lifetime 2–10 s (KU-2.5 unchanged) | Kong 2009 |
| single-cell vs collective flag | single-cell laminin traction *lower* (lit) vs Lam4 collective *higher* (poster) — the Layer-1↔Layer-2 split, not a model error | PI-exp map open item |

## Open items for PI / Lead

- [ ] Ratify col-I default = **single-cell anchor** (Taubenberger 1.3 s⁻¹ / 0.23 nm)
      vs the locked-open SMFS variant (0.44 s⁻¹ / 0.7 nm).
- [ ] Confirm `k_on` stays generic (KU-2.4) across ligands for now (no ligand-specific
      on-rate until K_D→k_on is calibrated), so only the off-rate differs by ligand.
- [ ] Guard the KU-2.5 catch-peak gate to **skip slip ligands** (NaN F* is correct).
- [ ] Sequence: wire after the literal-v0 γ run lands (don't perturb the active
      production driver).
