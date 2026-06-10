# Helper-session prompt — adopt SimuCell3D cell-cell + spheroid architecture into the DCM

> Paste this whole file into the helper session. Goal: read the SimuCell3D paper + SI
> word-for-word, deep-research the IBER (ETH) tissue/organoid models, inspect the existing
> SimuCell3D assets in this repo, and produce a CONCRETE plan to upgrade our Deformable
> Cell Model's cell-cell contact/adhesion + spheroid construction with SimuCell3D's
> architecture. Read-heavy + research; minimal code (a design + a small proof-of-concept
> only if quick). Report a precise integration plan back.

## Environment
Working dir `/Users/sw1/ffn_cellsim-platform` (HOOMD 7.0.1; branch `h7/compartment-platform`).
Env: `source ~/miniconda3/etc/profile.d/conda.sh && conda activate ffn_sim` (or
`~/miniconda3/envs/ffn_sim/bin/python`). PDFs are read with the Read tool's `pages=` param
(max 20 pages/call; for a long PDF, page through ALL of it — do not stop early).

## What we have now (the thing to upgrade)
Our DCM spheroid (read these): `ffn_sim/cell/dcm.py` — each cell is an **icosphere shell**
(harmonic edge springs = cortical elasticity) + **exact triangulated turgor** (DcmTurgorForce:
divergence-theorem volume + face-normal pressure) + an **adhesive substrate well**
(DcmSubstrateForce) + **cell-cell interaction as a node-level LJ / soft-core+linear-well**
(per-cell-type LJ in dcm.py, or the custom `DcmCellCellAdhesion` by `cell_of_node` in
`ffn_sim/cell/dcm_prolif.py`). `dcm_prolif.py` adds rim-biased contact-inhibited DIVISION
(pre-allocated cell pool, shared membrane type). VERIFIED: division 7→18 cells → footprint
A/A0 1.0→2.38. Master design + calibration: `ffn_sim/docs/v2_audit/DCM_SPHEROID_MASTER_DESIGN_2026-06-11.md`
(+ `references/analysis/_necrosis/`, `_junction/`). KNOWN WEAKNESS: our cell-cell contact is
crude node-pair LJ (porous shells can interpenetrate; not a real face-based contact); the
spheroid morphology (ball→compact→spread) is poor at small N.

## Task 1 — read SimuCell3D paper + SI COMPLETELY (every word)
- `ffn_sim/references/s43588-024-00620-9.pdf` — SimuCell3D (Runser, Vetter, Iber, Nat Comput
  Sci 2024). Read ALL pages.
- `ffn_sim/references/43588_2024_620_MOESM1_ESM.pdf` — Supplementary Information. Read ALL
  pages (the model equations, numerics, parameters live here).
Extract, with equations + parameter values + units:
  1. **Cell representation** — the triangulated deformable surface mesh per cell; the
     membrane/surface energy terms (surface tension, bending/Helfrich, area & VOLUME
     constraints, in-plane elasticity), and how forces are computed on the mesh.
  2. **Cell-cell CONTACT + ADHESION** — how two cell surfaces interact: face-face /
     node-face contact detection, the contact repulsion (non-penetration) + the ADHESION
     energy (this is the part our node-LJ does crudely). The adhesion strength parameters.
  3. **Rheology / integration** — overdamped? the time-stepping; remeshing; numerical
     stability; how they avoid mesh interpenetration.
  4. **Growth & DIVISION** — how a cell grows (volume/area target) and divides (mesh
     splitting, plane selection); proliferation control.
  5. **Spheroid / tissue / organoid construction** — how many cells, how the aggregate is
     built/initialized, boundary/medium, and any spheroid or organoid example.
  6. Performance/scale (cells, runtime), and whether GPU.

## Task 2 — deep-research the IBER models
Deep-research `https://git.bsse.ethz.ch/iber/Publications` (Dagmar Iber group, D-BSSE ETH):
find SimuCell3D and the group's other tissue/organoid morphogenesis models + their
cell-cell / spheroid architectures. Use WebSearch + the PubMed/bioRxiv MCP tools
(ToolSearch). Cite (author, year, DOI). Summarize which model components are reusable for a
fine-grained DCM spheroid.

## Task 3 — inspect EXISTING SimuCell3D assets in the project
There is prior SimuCell3D work — inspect (read-only):
- `~/simucell3d_ref/` — a cloned SimuCell3D source tree + a `build/simucell3d` binary. Map
  its source layout, the core data structures (cell mesh, contact), and the input/config
  format.
- in the sibling worktree `~/ffn_cellsim/ffn_sim/`: `scripts/simucell3d_morphology_vis.py`,
  `scripts/simucell3d_closure_baker.py`, `outputs/simucell3d/` — what was already done with
  SimuCell3D here (what they import, produce, and how it was used).
Report what already exists vs what SimuCell3D natively offers.

## Task 4 — integration plan (the deliverable)
Produce a CONCRETE plan for upgrading our DCM spheroid with SimuCell3D's architecture,
answering:
  - **Cell-cell contact/adhesion:** adopt SimuCell3D's face-based contact + adhesion energy
    to replace our crude node-LJ (fixes shell interpenetration + gives real cell-cell
    mechanics). Exactly which energy terms + parameters, and how to wire into the HOOMD
    BAOAB DCM (as a Custom force) OR whether to call SimuCell3D's engine directly.
  - **Spheroid construction:** how SimuCell3D builds/relaxes an aggregate → adopt for the
    ball→compact→spread morphology + many cells.
  - **Division + growth:** map SimuCell3D's mesh-split division onto our `dcm_prolif`
    pre-allocated-pool division (or replace it).
  - **Reuse vs rebuild:** can we drive SimuCell3D's `build/simucell3d` binary directly for
    the spheroid (and post-process with our necrosis/pressure/junction layers + the
    A/A0=a+b/R+c/R² analysis), OR port specific terms into `cell/dcm.py`? Recommend, with
    effort + the trade-off vs our HOOMD/BAOAB stack and the GPU plan.
  - How it composes with our calibrated necrosis (depth>150µm), bulk-pressure/junction-switch
    (q*≈3.81/5.4, ~0.5–5 kPa), and active traction (lamellipodium/contraction belt) layers.
Write the findings + plan to `ffn_sim/docs/v2_audit/SIMUCELL3D_INTEGRATION_2026-06-11.md`
and report a tight summary (key SimuCell3D equations/params + the recommended integration
path) back to the lead.
