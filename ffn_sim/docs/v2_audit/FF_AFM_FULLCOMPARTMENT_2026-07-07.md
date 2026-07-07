# Full-compartment virtual-AFM indenter + cortical tension (2026-07-07)

Delivers the PI ask: "인덴터로 누르면서 cortical tension 측정 + 완벽한 풀 컴파트먼트." All four compartments are
now WIRED into `simulate_whole_cell_compression_on_device` (network_warp.py) and indented at native scale
(Nc=494,802 cortex + 480 MT arms + MTOC + 3000-bead nucleus = 498,283 nodes) on the A5000.

## Native strain sweep (cortex + nucleus + membrane + MT, rigid plates)

| strain | V/V0 | dP [Pa] | F_plate [nN] | apparent γ [mN/m] |
|---|---|---|---|---|
| 0.00 | 1.000 | 36 | 0.01 | 0.14 |
| 0.05 | 0.996 | 3404 | 58 | 12.8 |
| 0.10 | 0.983 | 13187 | 471 | 49.6 |
| 0.15 | 0.962 | 29646 | 1581 | 111.5 |
| 0.20 | 0.935 | 52638 | 3920 | 198.1 |
| 0.25 | 0.901 | 84750 | 7883 | 319.4 |

**Cortical tension**: the AFM apparent (Young-Laplace) tension γ=½·ΔP·R rises 0.14 → 319 mN/m over 0–25% strain.
Resting γ≈0.14 mN/m is LOW (MCF7 suspended ~1–10 mN/m) — tied to the known **active-γ floor** (myosin-generated
cortical tension is ~530× under band; the AFM resting tension is turgor-dominated). Under indentation the
response is turgor-driven. (ΔP·R/2 is the effective/apparent tension, per standard AFM cortical-tension protocol;
the active myosin contribution is separately floored — the open scientific gap.)

## Which compartments actually bear load under AFM (coupling verification at 10–25% strain)

| compartment | wired | mechanically significant | evidence |
|---|---|---|---|
| **Cortex** | ✅ | **dominant** (the cell) | the elastic body itself |
| **Nucleus** (bead cloud) | ✅ | **YES, +45% plate force** | removing it: F 471→325 nN, dP 13187→8918 Pa at 10% (V_cyto=V_hull−V_nuc hydrostatic coupling — NOT a passenger, unlike in the crawl driver) |
| **Membrane** (lumped γ_mem) | ✅ | negligible | removing it: ΔF=0.0% (γ_mem≈2.7 Pa ≪ turgor 13–85 kPa) |
| **MT aster** (+tip↔cortex contact) | ✅ | negligible under AFM | ΔF=0.0% at 20%, +0.04% at 25% |

## MT wiring (the new engineering)

Added an additive `microtubule=` param to `simulate_whole_cell_compression_on_device`: the aster merges into
the elastic net ([cortex(Nc);MT_arms;MTOC;nucleus], Ne=Nc+Nmt+1), MT bending + MTOC hub crosslink ride the same
loop, and an **MT-tip↔cortex soft-contact** (excluded-volume) lets the aster prop the cortex under confinement.
`microtubule=None` → **byte-identical to the pre-edit function** (verified on CPU vs git-stashed baseline).

**Why MT is negligible under AFM (physical, not a bug):** the tiny +0.04% at 25% confirms MT IS in the sim
(not silently dropped). But under z-compression against a turgor-pressurized cortex (50–85 kPa), a 40-tube aster
with L_mt=6µm (< R=7.5µm) has almost no arms reaching the compression axis — only ~2 near-±z arms engage, and
40 struts can't compete with 495k turgor-loaded cortex nodes. **L_mt was NOT tuned up to force engagement**
(that would be tuning-to-outcome; `L_mt`/`n_mt` are flagged KB-ungrounded). If MCF7 MTs actually reach/exceed
the cortex radius (real MTs are often 10–25µm and touch the cortex), MT would bear load — needs a KB anchor (PI).

## Honest bottom line

The full-compartment cell is **built and fully wired** (all 4 compartments) and indented at native scale with a
measured force–indentation curve + cortical tension. But **under AFM compression only cortex+nucleus are
mechanically significant**; membrane (lumped) and MT (L_mt<R) are present/coupled but negligible because the
turgor pressure dominates. To make membrane matter → the explicit sheet+ERM (Thread-D); to make MT matter → a
KB-grounded L_mt that reaches the cortex, or a low-turgor / confined-migration regime (not AFM).

## Figures
- `outputs/ff/figs/afm_fullcompartment_cut.png` — cut-view of the full-compartment cell (cortex + MT aster +
  nucleus) at rest; the interactive viewer animates strain 0→0.25 (the AFM flattening the cell), all compartments
  visible via the cut plane. Full-resolution HTML (118 MB, all fibers) regenerable; too large to commit.
