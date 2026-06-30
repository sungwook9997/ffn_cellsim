# DCM faceting → production: automated execution plan (PI-greenlit 2026-07-01)

PI greenlit ("계획 잡아두고 자동화 가능 끝까지?") automating the confluent-init faceting work from the
validated prototype through to a production faceted spheroid, and on toward the original /goal
(aggregation properly built → proliferation separation visualized, mesh intact → dynamic spreading).

**Starting state (validated + committed this session):** confluent-init hypothesis confirmed end-to-end —
geometry (icosphere→Voronoi warp, N=8–400, watertight, no-penetration, START FACETED asph 0.04–0.08
Q 150–170) → existing energy MAINTAINS faceting (N=400 capstone asph 0.047 Q157 contact 0.67) → the one
flaw was turgor over-inflation (V/V0 1.37 + interpenetration) from V0=free-sphere. Commits 81a12ba
(viewer), 7b2b556 (prototype), 2c87a80 (init-npz+maintain), 8e9a34c (scale), dfa5969 (capstone).

## Automation policy

- Phases 1–4 (foam → clean production) auto-advance; **report at each phase gate**.
- **Halt + surface to PI** on: any **magic-number trigger** (a derived constant tuned to pass a gate), a
  **gate failure that would need a contract change**, or a **genuine design decision**. (Hard rules.)
- **Checkpoint halt at the 4→5 boundary** — proliferation + spreading (Phases 5–6) are open-ended with
  known hard issues (division stability, prior spreading structural limits); get explicit PI go before them.
- Every parameter stays at its physiological value; faceting must EMERGE from correct init + correct
  physics, never be tuned toward a target Q. If emergent Q < SimuCell3D, that is a finding to surface.

## Phases (each: code → gbook run → measure → gate)

### Phase 1 — V0 = confluent rest-volume setpoint  ✅ GATE 1 PASSED (2026-07-01)

**Result (N=400, `--v0-from-init`):** V/V0 **1.000** (was 1.37), pen-ratio **0.928** (was 0.192), asph
0.048→**0.048** (faceting perfectly held), Q 150→**150**, contact **0.857**. The V0=Voronoi-cell-volume
setpoint was exactly the fix — over-inflation + interpenetration eliminated; the energy maintains a clean
N=400 faceted foam at V/V0=1. Render `viz_html/8_V0FIX_n400_render.png`, `8_V0FIX_n400_clean.html`.

- **Change:** `--v0-from-init` derives each cell's osmotic V0 from its ACTUAL volume in the loaded
  confluent mesh (divergence-theorem volume), not the free-sphere V0. Physiological confluent rest volume;
  derived from geometry, not tuned. (Implemented; mean V0 ratio 0.893 of free-sphere at ε=0.02.)
- **Run:** N=400 confluent init + V0-fix, existing energy (turgor + conservative tent), 15k steps.
- **Gate:** V/V0 ∈ [0.95, 1.10], penetration-ratio > 0.9 (no significant interpenetration), faceting held
  (asph ≥ ~init, Q ≥ ~150), contact high. → removes the capstone's only flaw.

### Phase 2 — Path-B `--builder confluent` production wiring
- **Change:** promote `confluent_init_prototype.py` to a first-class `--builder confluent` in the driver
  (build the confluent foam directly, no npz round-trip). Speed up the Lloyd MC (vectorize assignment /
  scale nsamp) so N≥400 builds in seconds, not the current timeout.
- **Gate:** `--builder confluent` reproduces the prototype geometry (manifold, no-pen, contact, asph/Q
  within tolerance) and runs end-to-end into the maintain energy with the Phase-1 V0 setpoint.

### Phase 3 — sharpen faceting  ✅ CONCLUDED (2026-07-01): Q~150 is a structural ceiling
Swept the legitimate levers at N=200 `--builder confluent`: baseline subdiv2 **Q 150** (asph 0.049);
S2 = REAL differential-γ faceting energy (`--diff-tension --contact-tension-frac 0` at derived w_cs)
**Q 148** (asph 0.030, V/V0 1.000, pen 0.990 — holds faceting clean but does NOT sharpen); S1 = finer
subdiv 3 init **Q 157→161** (modest +7%, 4× heavier). **None reaches SimuCell3D's Q≈250.** Each cell is
a smooth RADIAL deformation of an icosphere filling its Voronoi region — mildly faceted, not razor-sharp
flat-faced polyhedra. Reaching Q≈250 needs a flat-faced cell representation (clip the icosphere to the
Voronoi half-planes, or a different mesh) = a structural change, NOT parameter tuning. **Conclusion (to
PI):** the baseline foam (Q~150, clean, V/V0 1, no-penetration, faceting-held, space-filling) IS a genuine
compact faceted tissue = the production result; razor-sharp Q≈250 is an optional later mesh upgrade.
Production uses subdiv 2 (Q150, clean, fast). No magic-number tuning was used.

### Phase 4 — clean N=400 production aggregation + validation + viz
- **Run:** full N=400 confluent foam, V0-fix, `--builder confluent`, physiological energy, long settle.
- **Validate:** manifold/watertight, no-penetration, V/V0≈1, faceting (asph, Q), contact fraction,
  isoperimetric Q distribution vs SimuCell3D bands. Cut-HTML + render for PI.
- **Outcome:** "aggregation properly built" — a real compact faceted tissue (the original goal's part 1).
- **→ CHECKPOINT: report to PI, request go-ahead for Phases 5–6.**

### Phase 5 — proliferation separation (original /goal part 2) — PI go-ahead first
- Enable division from the confluent foam; cells divide + jam; remesh keeps the mesh intact; visualize
  cells separating cleanly. Known risk: division stability needs remesh-every-step + clean cleavage.

### Phase 6 — dynamic spreading (original /goal part 3) — PI go-ahead first
- From the proliferated/foam state, dynamic spreading; A/A0 top-down silhouette validation. Known prior
  structural limits — likely multiple PI-decision points.

## Execution

Driven by the autonomous `/loop` (this session): each phase = launch gbook run → wakeup → measure →
gate-check → advance or halt-to-PI. gbook GPU (A5000), disjoint from the FF session (ff/ only). Commit
each phase's code + results + viz. Notion closeout + receipt at the eventual session wrap.
