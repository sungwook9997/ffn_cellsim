# Full-compartment DCM viewers — navigation

GPU-native (A5000) full-compartment DCM spheroid. Each viewer is an **interactive Three.js cell-morphology
viewer** — open the `.html` in a browser, drag to rotate, scroll to zoom. Every cell renders as its
deformable mesh; with nuclei on, a nucleus sphere shows inside each translucent cell (membrane/cortex shell +
cytoplasm interior + nucleus = the full compartment stack). Reconstruction write-up:
`docs/v2_audit/FULLCOMPARTMENT_GPU_RECONSTRUCTION_2026-07-02.md`.

**What is mechanically real vs a render aid (real-vs-drawn honesty):** the cell mesh (membrane/cortex), turgor
(Π₀=133 Pa), cytoplasm drag (η=65.9 Pa·s), membrane surface tension (γ=5e-4), contact, and the nucleus are all
*mechanically active force kernels launched every step* (verified: `nucleus_force_kernel`, `surface_tension_kernel`,
`dcm_turgor_force_kernel` all `wp.launch` in the step loop). The nucleus specifically is a **radial elastic force
field** (E_nuc=4700 Pa, chromatin→lamin bilinear) that pushes the cell-mesh nodes to sustain a core of radius
R_nuc=0.25R — it is *not* a separately-meshed deforming nuclear envelope. The **sphere you see in the viewer is a
render aid** marking that compartment's centre + extent, not a simulated nuclear surface. So the nucleus is real
mechanics, drawn as an idealized sphere.

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
- `FULLCOMPARTMENT_n400_CROSSSECTION_slab_layer_overlap.png` — one-cell **slab**: cells apposed + space-filling
  toward the dense-nuclei core (the confluent Voronoi seam — cells share interfaces; healthy, not pathological).
- `FULLCOMPARTMENT_n400_CONTROL_explicit_vs_implicit.png` — matched-N=400 control, explicit (baoab) vs implicit
  (IMEX) mid-plane cuts side by side: visually + numerically identical (pen 2.23≈2.11, gap 0.156≈0.158, asph
  0.0543==0.0543) → the residual pen is robust to the integrator = it is geometry.
- `FULLCOMPARTMENT_n400_CONTACT_penalty_vs_ipc.png` — penalty tent (pen 2.23, gap 0.156 µm) vs `--ipc` Li-2020
  log-barrier (pen 1.51, gap 0.419 µm) on the SAME confluent full-compartment: `--ipc` lowers the pen ~30 % but
  widens the gap — neither cleans it (both FAIL G2, both V/V0=1.0). The pen is a geometric floor (the confluent
  seam), independently verified on our own runner (not just the concurrent session's number).

## Honest contact state (reframed by audit#7)

- **Floating: FIXED** — inter-cell gap 0.67 → **0.156 µm** (~2 % of R, a physiological cleft) via conservative
  contact + differential interfacial tension.
- **Pressing: YES, and the `pen_frac`≈2.1–2.6 is the confluent geometry, NOT a contact failure.** It is present
  from t≈0 (confluent Voronoi cells share interfaces), robust to the contact method (penalty→2.1, `--ipc`
  log-barrier→2.6 — a *better* contact doesn't lower it) and to the integrator (explicit 2.23 ≈ implicit 2.11),
  at a healthy V/V0=1.000, porosity ~0.10 (space-filling), faceted. The G2 gate (pen<0.3) is an aggregation-regime
  gate for *separate* cells; confluent space-filling tissue inherently exceeds it. `--ipc` (a genuine Li-2020
  log-barrier) already exists and gives the same pen, confirming the overlap is geometric. *(Earlier this file
  called it "interpenetration NOT fixed, needs IPC" — that was over-pessimistic and is corrected here.)*
