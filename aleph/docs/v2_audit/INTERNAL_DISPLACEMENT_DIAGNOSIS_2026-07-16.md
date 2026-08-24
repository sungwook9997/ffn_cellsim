# Internal-displacement diagnosis — the interior is mechanically DECOUPLED (2026-07-16)

**PI observation (2026-07-16):** "세포 내부 displacement가 이상하다" (the cell's internal displacement looks
odd). Diagnosed on the NATIVE full cell (`scripts/ff_internal_displacement_diag.py`, CLAUDE.md HARD:
native + full + interactive 3D). Method: build the physiological cell (cortex NF=70686 + 3000-bead nucleus
+ reservoir membrane + 40-tube MT aster), relax to the pre-stressed baseline (strain 0), AFM-compress to
strain 0.10, measure per-node displacement Δ = pos(loaded) − pos(rest) by compartment.

## Result (NATIVE, Nc=494,802 — authoritative)

| compartment | n | RAW mean \|Δ\| | RAW std | verdict |
|---|---|---|---|---|
| **cortex shell** | 494,802 | **39.8 nm** | 133.5 nm | deforms (real spread, max 751 nm) |
| **MT aster + MTOC** | 481 | **0.0 nm** | 0.0 nm | **FROZEN — bit-identical rest vs load** |
| **nucleus** | 3,000 | **0.0 nm** | 0.0 nm | **FROZEN — no move, no deform** |

- nucleus flattening under the plate: **0.0 %** (z-half 5.106 → 5.106 µm) — the nucleus does not deform.
- MT tip ↔ cortex minimum gap: **0.900 µm** — the MT tips do not reach the cortex (contact shell 0.5 µm).
- COM drift 0.15 nm (negligible; native compression is symmetric).

**VERDICT: the interior (nucleus + MT aster) is mechanically DECOUPLED from the cortex.** Under
physiological AFM compression the cortex deforms while the nucleus and microtubules do not move or deform
*at all* (displacement identically zero). This is exactly the PI's "odd internal displacement."

## Root cause (grounded, not scale — confirmed identical at coarse + native)

The interior has NO force-transmission path from the deforming cortex at sub-contact strain:
1. **Cytoplasm was 0-D** — no internal stress/flow field to carry cortex compression inward. (The spatial
   Biot FSI added 2026-07-16 pushes only the CORTEX shell, not the nucleus/MT — so even FSI doesn't couple
   the interior yet.)
2. **No cytoplasm→nucleus pressure coupling** — `nucleus_shell_kernel` only restores beads toward R_nuc
   about the cell centroid; the cytoplasm ΔP rise under compression never acts on the nucleus surface.
3. **No IF cage** — vimentin/keratin network + LINC linking nucleus ↔ cortex is genuinely ABSENT (ladder
   finding #4). This is the primary SOLID coupling a real cell uses to transmit cortex strain to the nucleus.
4. **MT tips don't reach the cortex** — 6 µm MTs from a central MTOC end 0.9 µm short of the 7.5 µm cortex,
   so the tip↔cortex contact never engages at physiological strain → MTs can't act as compression struts.

The interior only engages under EXTREME confinement (strain > ~0.23, direct plate contact). At
physiological strain the sub-contact transmission (pressure + cytoskeletal linkage) is missing.

## Fixes (mechanistically-faithful, framework §load-transfer + tensegrity)

- **A. Cytoplasm pressure → nucleus** (quick, hydrostatic): apply the cytoplasm ΔP (and the spatial FSI
  pore pressure) as an inward normal force on the nucleus surface beads (a nucleus-inward `turgor_kernel`).
  → the nucleus compresses/flattens when the cell is squeezed. Reuses existing machinery. Gives PRESSURE
  transmission, not solid-shear coupling.
- **B. IF cage** (proper solid coupling; T5.I build target): vimentin/keratin WLC network with
  strain-stiffening + cytolinkers + LINC, linking nucleus surface ↔ cortex. Transmits cortex DEFORMATION
  to the nucleus, gives internal structure. New compartment (larger build). The framework-faithful fix.
- **C. MT thickness + strut engagement** (tensegrity, PI 2026-07-16): give MTs real 25 nm tube radius +
  shaft excluded-volume, and make them span/engage the cortex (length/MTOC geometry or lateral coupling)
  so they bear COMPRESSION as struts (MT compression vs actin tension, visible under AFM/protrusion).
- **D. LINC** (nucleus ↔ MT/actin anchor): couples the MTOC/cytoskeleton to the nuclear envelope.

## Change log
- 2026-07-16: created. Native diagnosis of the interior-decoupling anomaly + root cause + fix options.
  Harness `scripts/ff_internal_displacement_diag.py`; figure `outputs/mech_hier/figs/internal_disp_native.html`.
