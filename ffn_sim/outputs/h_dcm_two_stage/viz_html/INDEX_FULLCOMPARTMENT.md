# Full-compartment DCM viewers — navigation

GPU-native (A5000) full-compartment DCM spheroid. Each viewer is an **interactive Three.js cell-morphology
viewer** — open the `.html` in a browser, drag to rotate, scroll to zoom. Every cell renders as its
deformable mesh; with nuclei on, a nucleus sphere shows inside each translucent cell (membrane/cortex shell +
cytoplasm interior + nucleus = the full compartment stack). Reconstruction write-up:
`docs/v2_audit/FULLCOMPARTMENT_GPU_RECONSTRUCTION_2026-07-02.md`.

## Viewers (open in a browser)

| File | N | What it shows |
|---|---|---|
| `FULLCOMPARTMENT_n8_membrane_cytoplasm_nucleus.html` | 8 | smallest — clearest single-cell compartment read |
| `FULLCOMPARTMENT_n200_GPUnative_membrane_cytoplasm_nucleus.html` | 200 | GPU-native compartment stack |
| `FULLCOMPARTMENT_n200_FINAL_diffDAH_membrane_cytoplasm_nucleus.html` | 200 | + differential-DAH tension (tightest apposition, gap 0.167 µm) |
| `FULLCOMPARTMENT_n400_GPUnative_membrane_cytoplasm_nucleus.html` | 400 | production-scale spheroid |
| `FULLCOMPARTMENT_n400_FINAL_diffDAH_membrane_cytoplasm_nucleus.html` | 400 | **flagship** — final run + **deep-link support** (below) |

## Deep-link cross-sections (flagship viewer)

The flagship `n400_FINAL` viewer reads URL params so you can jump straight to a specific cross-section
(and it is headless-renderable for proofs):

```
FULLCOMPARTMENT_n400_FINAL_…nucleus.html?section=cut&axis=0&pos=0.5&frame=last
```

| param | values | meaning |
|---|---|---|
| `section` | `cut` \| `peel` \| `slab` \| `off` | cut = true stencil-capped clip (solid interior); peel = whole near-side cells; slab = one-cell layer |
| `axis` | `0` \| `1` \| `2` | clip axis (x/y/z) |
| `pos` | `0`–`1` | clip-plane position along the axis |
| `flip` | `0` \| `1` | flip which side is kept |
| `thick` | `0`–`1` | slab thickness (slab mode) |
| `frame` | `N` \| `last` | trajectory frame (`last` = settled state) |

In the UI: tick **section**, pick the **mode** dropdown, drag **cut pos**. `cut` reveals the interior;
`slab` isolates one cell-layer.

## Honest-state render proofs (PNG)

- `FULLCOMPARTMENT_render_proof_nucleus_visible.png` — nuclei visible inside translucent cells.
- `FULLCOMPARTMENT_n400_render_proof_faceted_nuclei.png` — exterior: faceted tissue + nuclei.
- `FULLCOMPARTMENT_n400_CROSSSECTION_cut_interior_nuclei.png` — mid-plane **cut**: packed faceted cells, a
  nucleus in every cell.
- `FULLCOMPARTMENT_n400_CROSSSECTION_slab_layer_overlap.png` — one-cell **slab**: shows the cell volumes
  crossing toward the dense-nuclei core (the interpenetration, shown honestly not hidden).
- `FULLCOMPARTMENT_n400_CONTROL_explicit_vs_implicit.png` — matched-N=400 control, explicit (baoab) vs implicit
  (IMEX) mid-plane cuts side by side: visually + numerically identical (pen 2.23≈2.11, gap 0.156≈0.158, asph
  0.0543==0.0543) → the integrator is NOT the interpenetration lever (needs IPC, not better integration).

## Honest contact state

- **Floating: FIXED** — inter-cell gap 0.67 → **0.156 µm** (~2 % of R, a physiological cleft) via conservative
  contact + differential interfacial tension.
- **Interpenetration: NOT fixed** — `pen_frac` settles 2.1–2.3 (peak 3.24) = 7–11× the G2 gate (0.3); the run
  is otherwise healthy (V/V0=1.000) so it's a genuine equilibrium overlap, not a blowup. The fix is true
  log-barrier IPC (in progress in the contact-physics core). The exterior renders look faceted-good; the
  cross-section PNGs above show the real interior overlap.
