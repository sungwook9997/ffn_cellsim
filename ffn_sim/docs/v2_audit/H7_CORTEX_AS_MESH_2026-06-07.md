# H.7 — cortex-as-mesh investigation (PI proposal, 2026-06-07)

PI proposed building the cortex as a connected MESH (built-in connectivity / guaranteed mesh-mesh
adhesion), worried the independent-filament-lines + emergent-crosslinker design fragments at the
×40 mesoscale; later refined to a DEFORMABLE / game-style variable mesh, and asked whether that is
realistic for ACTIVE SPREADING. A 4-agent Workflow (literature + HOOMD design + critique → synth).

## Verdict: do NOT adopt a triangulated-surface mesh. The instinct is right; the FORM is wrong.

### 1. Mesh/network IS standard — but structure-specific, and MCF7 is on the wrong side
- **RBC membrane-skeleton** (a true quasi-2D spectrin LATTICE): triangulated mesh = gold standard,
  validated vs optical-tweezers/micropipette (Boey-Boal-Discher 1998; Li-Dao-Lim-Suresh 2005;
  Fedosov 2010). BUT the best models use an **IRREGULAR random graph (degree 2-9) with DYNAMICALLY
  BREAKABLE bonds** (Fai 2017; Li-Lykotrafitis 2007) — NOT a rigid fixed triangulation.
- **Bulk crosslinked actin**: the MIKADO / random-fiber model (explicit fibers + crosslinks at
  intersections) is canonical (Wilhelm-Frey 2003; Head-MacKintosh 2003; Broedersz-MacKintosh 2014
  RMP) — i.e. EXACTLY the codebase's own `ecm/mikado.py` and what the current cortex already IS.
- **Bulk cortex as a sheet**: continuum active-surface FE on a triangulated mesh (Salbreux-Jülicher;
  Torres-Sánchez-Arroyo 2019) — but that is a LUMPED continuum, not explicit filaments.
- ⇒ The validated fixed-triangulation success is the PASSIVE RBC membrane-skeleton. The MCF7 bulk
  actomyosin cortex's standard fine-grained representation is a DYNAMIC crosslinked fiber network
  (Mikado/AFINES/MEDYAN/Cytosim) — the current explicit-filament approach. Fragmentation is a ×40
  SCALE-MISMATCH artifact (~1880 nm spacing vs 60 nm reach), not wrong physics.

### 2. Connectivity is already solved (and isn't "more physical" via a mesh)
`connected_mesh.py build_connected_cortex` gives z=3.3, giant 99% with EXPLICIT filaments, seeded at
construction. A mesh would compensate the same ×40 artifact by hard-coding edge length instead of
the √(A/n) bridge reach — neither is more physical at the filament level. **Gap found: the PRODUCTION
script (`mcf7_fullcell_stage1.py`) runs the deferred "FAST HYBRID" (seed onto bare lines), NOT the
full faithful build** — so the real "built-in connectivity" the PI wants is one promotion away.

### 3. THE CRUX — a mesh does NOT touch the γ-floor, and a triangulated one WORSENS it
The floor is GENERATION (s_grip≈0: heads bind+load but don't WALK → tension ~1e-4, ~1000× under the
dipole ceiling — myosin kinetics, orthogonal to topology) gated behind BUCKLING (M-SHAKE forbids
r/r0<1; Ronceray-Broedersz-Lenz 2016 ~10× amplification needs buckling). A closed triangulation is
OVER-CONSTRAINED (each edge in 2 triangles, area-incompressible, 5-6 coordination) → it SUPPRESSES
the single-filament buckling that is the project's amplification lever. The design agent's escape
("soft edges + drop M-SHAKE enable buckling") works on the EXISTING filament lines too — **relaxing
the constraint is the lever, not the mesh.**

### 4. Active-spreading realism (PI's sharp question) + the game-deformable-mesh idea
A FIXED mesh is realistic for a passive stable skeleton (RBC) but WRONG for active spreading, which
is dynamic remodeling (Arp2/3 branched polymerization, treadmilling, flow — adding/removing
filaments). A game-style DEFORMABLE mesh (mass-spring/PBD) handles shape-change, and HOOMD is
already a mass-spring engine — but spreading is mostly GROWTH (polymerization ADDS area), so it needs
dynamic RE-MESHING (add nodes at the front), not just stretching a fixed topology. The literature
agrees the best cortex/skeleton models are IRREGULAR + DYNAMIC, not fixed. **The deformable/dynamic
vision is best realized by the EXPLICIT-FIBER network (dynamic crosslinkers + turnover + lamellipodium
polymerization growth), NOT a triangulated surface.**

### 5. Fidelity (CLAUDE.md)
Triangulated surface mesh = a 2nd unsanctioned coarse-graining (edge abstracts many filaments) →
VIOLATES the explicit-filament hard rule; migrates magic numbers (edge stiffness/area modulus have
no per-filament anchor); forfeits the Chugh-2017 regulated filament-length physics; needs synthesized
myosin walk-paths. **Mikado-on-sphere does NOT violate it** (explicit bead-chain fibers, same class
as ecm/mikado.py; preserves bead order myosin's grip_walk needs).

## RECOMMENDATION (priority order)
1. **~free, unambiguous:** (a) promote the full faithful `connected_mesh.py` build to the production
   default (the real "connectivity built-in at construction" the PI wants; production currently runs
   the fast hybrid). (b) fix the latent bimodal myosin bug `myosin.py:544`.
2. **Attack the wall (the only thing that lifts γ):** gate A (myosin walking + construction-time
   binding pre-equilibration) + gate B (relax M-SHAKE so filaments buckle — ~10× lever; PI
   integrator-freeze sign-off; a mesh would BLOCK this).
3. **ONLY IF built-in topology still wanted after A+B:** Mikado-on-sphere reusing `ecm/mikado.py`
   (NEW `cortex/mikado_sphere.py`; crosslinks at geodesic intersections at the physical 60 nm scale,
   soft edges → buckling-capable; γ measure + myosin need ZERO edits — topology-agnostic). ~3-4 days.
   **Reject the triangulated icosphere outright.**

## Citation caution
Ronceray-Broedersz-Lenz 2016 PNAS (10.1073/pnas.1514208113) + Chugh 2017 were cited from memory by
the design/critique agents — VERIFY before any Notion SourceEvidence use. RBC/active-surface DOIs
were PubMed-checked but spot-check full text.
