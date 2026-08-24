# ECM Platform — extensibility plan (in vitro → in vivo)

> **STATUS: DRAFT plan — 2026-05-30. NOT a contract.** Produced in a read-only
> research session; no code written. Companion to the Notion KU v2 ECM-V layer
> (KU-1.V.1/2/3) and `CELL_MECHANICS_EXTEND_VS_REBUILD.md`. For PI/Lead review.

## Why this exists (PI 2026-05-30)

PI's experimental reference is **(b) collagen-I COATED on a rigid dish** (glass/TCPS,
backing E ~ GPa) — the standard adhesion-assay condition. The platform must later
**extend to arbitrary ECM conditions spanning in vitro and in vivo** (soft gels,
3D matrices, tumor stroma, tissue, confinement) as **config selections, not code
rewrites**. This plan specifies that extensible platform.

## Current state (code-grounded — what exists today)

The platform already has **two disjoint substrate representations, and NO tunable
stiffness between them**:

| Representation | What it is | Effective backing stiffness |
|---|---|---|
| **H.4 FA substrate-ligand** ([bridge/fa.py](bridge/fa.py): "synthetic substrate plane at z=0"; [cell.py](cell/cell.py) `SubstrateLigandPin`) | ligand particles at z=0, **pinned bit-exactly immobile every step** | **infinitely rigid** (GPa-equivalent) — closest to (b) |
| **H.1 ECM Mikado** ([ecm/mikado.py](ecm/mikado.py)) | 2D soft collagen-I fiber network, free-standing, G_0 ~ few Pa | **soft compliant** (Pa) |

**Gaps (verified by grep — none of these exist):**
1. **No tunable substrate compliance `E_sub`** — substrate is *binary* (rigid pin OR soft network). The whole in-vitro→in-vivo stiffness axis (PA/PEG gels 0.1–40 kPa, tissue, tumor stroma) is unrepresented.
2. **No ligand identity** — the FA "ligand" is generic; not collagen-I vs fibronectin(RGD) vs laminin.
3. **2D only** — Mikado generator is `generate_2d_fiber_network` (z=0 plane); no 3D matrix.
4. The soft Mikado is **not coupled into `Cell.build`** (Lead's "floating shell" finding) — the cell currently feels the rigid FA pin, not the collagen network.

So today, the cell effectively sits on **(b)-like rigid backing** (FA pin) — but without collagen-I ligand identity or any way to soften it.

## Design principle — decompose ECM into 4 orthogonal config axes

Any ECM condition (in vitro or in vivo) = a point in a 4-axis space. Make each axis
a config selection so conditions are presets, not rewrites:

1. **Backing compliance `E_sub`** — rigid coated dish (GPa) | tunable gel (0.1–40 kPa) | tissue (0.1–50 kPa) | none (free network).
2. **Matrix** — thin coating (ligand-only) | 2D fiber gel | 3D fiber gel | dense+crosslinked+aligned tumor stroma. (H.1 Mikado knobs: density ρ, crosslink, alignment S — see KU-1.V.2.)
3. **Ligand presentation** — collagen-I (α2β1/GFOGER) | fibronectin (RGD/α5β1) | laminin | mixed. (FA binding species + density.)
4. **Geometry / confinement** — 2D open | 3D embedded | channel/pore (µm) | dimensionality. (KU-1.V.3.)

## The key unlock — substrate compliance: rigid pin → tunable spring

The single change that opens the entire in-vitro→in-vivo stiffness axis (and is
fully **additive, default = current behavior**):

- **Today:** `SubstrateLigandPin` resets each ligand to z=0 every step ⇒ infinite stiffness.
- **Proposed:** replace the hard pin with an **anchor harmonic spring of tunable stiffness `k_sub`** from each ligand to its z=0 site (a Template-1 custom force over ligand tags, mirroring `erm.py`/`enclosed_volume.py`). 
  - `k_sub → ∞` (or `pin: true`) **exactly recovers today's rigid (b) backing** — so default-rigid is bit-for-bit identical to current FA runs.
  - finite `k_sub` ⇒ **compliant substrate** (PA-gel / tissue). Map `k_sub → effective E_sub` via the ligand areal density (E_sub ≈ k_sub · n_ligand / A_substrate, dimensional bridge; calibrate against KU-1.V.1 bands + a flat-gel oracle).
  - This makes substrate stiffness a **continuous tunable axis** instead of binary, with the rigid coating and the soft gel as the two endpoints — and the cell's FA/clutch (KU-3.5/5.1) then feels the *right* stiffness, which is exactly what the Lead's FA work needs.

## Condition matrix — presets (in vitro → in vivo)

Each row = one `ecm_condition` config preset; columns = the 4 axes + KU anchors.

| Preset | Backing E_sub | Matrix | Ligand | Geometry | KU anchor |
|---|---|---|---|---|---|
| **(b) Col-I-coated rigid dish (PI exp)** | rigid (pin / k_sub→∞) | thin coating (ligand only) | collagen-I | 2D open | KU-1.V.1 stiff cap; KU-3.5 |
| Col-I-coated PA/PEG gel | tunable 0.1–40 kPa | thin coating | collagen-I | 2D | KU-1.V.1 sweep + YAP ~5 kPa |
| 2D reconstituted col-I gel (current H.1) | network self (~few Pa) | 2D Mikado ~1.5 mg/mL | collagen-I | 2D | current baseline; KU-1.30 |
| 3D col-I gel | network self | 3D Mikado (density sweep) | collagen-I | 3D embedded | KU-1.V.2 + V.3 |
| Tumor stroma (in vivo-like) | stiff 5–50 kPa | dense + LOX-crosslinked + TACS-3 aligned | col-I + fibronectin | 3D + aligned | KU-1.V.1/2 tumor |
| Normal tissue (in vivo) | 0.1–10 kPa | physiological 3D | tissue-specific | 3D | KU-1.V.1 normal |
| Confined channel (in vivo migration) | walls + matrix | 3D + pore 3–10 µm | mixed | confinement | KU-1.V.3 (pore 7 µm² → H.9) |

## Platform prep — what to build (additive, default-off / default-rigid)

| Piece | Spec | Reuse / new |
|---|---|---|
| `ecm_condition:` config block | selects {backing, matrix, ligand, geometry} preset → resolves into existing H.1 + H.4 + new substrate params | **new config schema** |
| `ecm/substrate.py` (or `bridge/substrate.py`) | tunable-compliance anchor: `k_sub` spring over ligand tags (Template-1 custom force); `pin:true`/`k_sub→∞` = current rigid | **new module** (the key unlock) |
| Ligand-species abstraction | collagen-I / fibronectin / laminin → FA integrin binding params (k_on, catch params, capture) | **extend bridge/fa.py** |
| 3D Mikado | generalize `generate_2d_fiber_network` → 3D (brief already noted 2D-oracle/3D-runtime) | **extend ecm/mikado.py** |
| Mikado↔Cell coupling | wire the collagen network into `Cell.build` as the compliant matrix (couples to Lead's FA-integration) | **integration** |

Example config sketch (illustrative, not final):
```yaml
ecm_condition:
  preset: col1_coated_rigid        # (b) PI experiment — DEFAULT
  backing:   { mode: rigid }        # rigid | gel | tissue | none ; gel→ E_sub_kPa: <val>
  matrix:    { type: coating }      # coating | mikado_2d | mikado_3d | stroma
  ligand:    { species: collagen_I, density: <val> }
  geometry:  { dimensionality: 2D } # 2D | 3D | channel(width_um)
```

## Sequencing

Rides on **H.1 (ECM) + H.4 (FA)**. The `k_sub` substrate-compliance unlock is the
load-bearing piece and **couples directly to the Lead's FA Cell-integration** (it
sets what stiffness the clutch feels). Order: **after FA production lands** →
add `substrate.py` (tunable compliance, default-rigid = (b)) → ligand identity →
3D Mikado → presets. Each preset validated against its KU-1.V anchor.

## Open items for PI

- [ ] Ratify the 4-axis decomposition + `ecm_condition` config-preset approach.
- [ ] Ratify the **rigid-pin → tunable `k_sub` spring** unlock (default-rigid preserves current FA runs bit-for-bit). Confirm the `k_sub → E_sub` calibration oracle (flat-gel indentation vs KU-1.V.1).
- [ ] Set **(b) col1-coated-rigid as the DEFAULT preset** (PI's experimental baseline) — so out-of-the-box matches the lab.
- [ ] Confirm ligand species set (collagen-I first; fibronectin/laminin later).
- [ ] In-vivo presets (tumor stroma, tissue, confinement) gate on H.9 nucleus (pore limit) — sequence after H.9.
