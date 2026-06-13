# Mechanistic GPU lamellipodium graft — build + findings (2026-06-14)

Branch `h7/compartment-platform`. Commit `187da2e`. Follows
`LAMELLIPODIUM_GRAFT_PLAN_2026-06-14.md` (PI chose "build it — full mechanistic graft").

## 1. Why (diagnosis recap, all PI-confirmed)

The DCM cohesive spheroid would not spread with the body-force rim-traction proxy
(`DcmActiveRimTraction*`): every lever was exhausted and ALL contract —
- traction sweep 8×/50×, belt on/off → contract;
- lit-anchored ξ/ω (4e7/5e7), weak cohesion (ω=1e7) → contract;
- strong substrate wetting (W_cs=2.85e-3 lit) + weak cohesion + strong traction →
  **A/A0 0.72** (combined-lever, 2026-06-14).
A body force cannot extend a cell's basal area (flatten it onto the dish); it only
translates nodes, and against cohesion + turgor the net is contraction.

## 2. The graft (`cell/dcm_lamellipodium_gpu.py`, wired into `build_gpu_dcm_simulation`)

GPU-resident, leak-free (dormant-pool, no `set_snapshot`). Three explicit pieces —
NOT a body force:
- **dormant actin-bead pool** in the snapshot (`cell_of_node=−1` → invisible to
  node-face/tent contact; `actin_lamel` typeid → settling skips it; parked → force-
  free). The front ADVANCES by ACTIVATING a parked bead at a force-free site (a
  local-snapshot position write) — same pattern as the remesh/prolif pools, no
  rebuild leak, GPU-resident.
- **FA molecular clutch = per-bead anchor spring** (`ActinClutchAnchorGPU`):
  gripping pins the bead at its birth xy + z0. Folds the explicit ligand grid +
  LigandPin + grip-bonds of the legacy CPU engine into one anchor array (the dense-
  ligand limit), so there is NO bond topology to mutate. The anchor is the FA
  ENSEMBLE (~N_int·KU-2.4 ≈ 5e-2 N/m, ~rigid dish) — it MUST be ≫ k_tether (4e-3)
  or the tether's Newton-3 reaction drags the anchored actin inward (both ends at
  equal drag meet in the middle → no net protrusion). Verified necessary.
- **traction tether** (`LamellipodialTractionTetherGPU`): each rim cell's leading-
  basal membrane node pulled (capped at 5·60Pa·area_per_node = lit MCF7 traction)
  toward the nearest outward anchored actin of its own cell — Newton-3 on the actin.
  Real traction carried by actin/clutch.
- **ratchet** (`GpuLamellipodiumAdvance`, low cadence): a lane advances (activates a
  bead one ℓ0 outward) only when the membrane has CAUGHT UP to its front bead
  (catch-gate), with prob `p = v_front·dt_tick·S/ℓ0`. Front velocity (lit 3–12
  µm/min) is the one constitutive input; footprint/traction/flattening emergent.

Simplifications surfaced to PI (mesoscale construction choices, not a body force):
(1) "born gripped" — a bead anchors at activation (it is born basal, in the contact
band, where it would grip next tick anyway); (2) per-leading-node radial lanes vs
free Arp2/3 branching (membrane-tracked ratchet, as the legacy nucleation updater).

Key build fixes during validation:
- rim cells = the BASAL contact disk (centroid within `contact_band·R` of z0); the
  legacy "few-neighbours periphery" gate wrongly excluded the packed bottom-contact
  cells that actually crawl.
- basal band = `0.3·R` (the substrate contact zone) — NOT `2·mean_edge` (≈9µm),
  which put "leading basal" nodes beyond the tether's reach of a z0-seeded bead →
  the lane re-seeded every tick → pool exhausted in ~1000 steps.

## 3. Validation (the graft engages + spreads)

- **CPU (numpy path), small N** — ratchet engages, seeding controlled, A/A0 RISES
  (1.00 → 1.05) at a moderate clock; the body-force proxy CONTRACTS. Confirms the
  mechanism, not the proxy, spreads.
- **GPU (cupy path), N=400** from the aggregated spheroid — builds + runs stable on
  the A5000 (manifold held, V/V0 healthy 1.125→1.084 at dt=1e-4, S=1), tether
  engages (n_teth=24), no crash. The cupy force paths (scatter_add,
  gpu_local_snapshot, node-face NVRTC) all work.

## 4. The kinetic-budget / quasi-static WALL (the open PI question)

Reaching the PI experimental A/A0 7–10 means the spheroid silhouette growing
73→~200 µm radius = ~240 ℓ0 of collective front advance. This is intrinsically a
~hours real-time process. Two regimes observed:
- **stable but slow**: dt=1e-4, S=1 → V/V0 held but A/A0 ≈ 1.006 in 2000 steps; A/A0
  7–10 would need ~1e8 steps (physiological clock) = wall-clock-infeasible.
- **fast but compresses**: dt=1e-3 (S=5) → A/A0 CONTRACTS to 0.754, V/V0 0.876 in
  2000 steps (volume not conserved at the larger dt); v_front=120 (S≈20, dt=1e-4)
  on small N → A/A0 peaks 1.05 then falls to 0.99 while V/V0 → 0.85. Accelerating
  the clock past the cell's quasi-static spreading rate SQUASHES the cell (the
  tether out-races turgor relaxation) instead of flattening it at constant volume.

This is the same class of result as the Layer-2 magnitude structural limit and the
active-γ floor: the FORM/mechanism is reproduced (cells crawl + flatten), but the
fine-grained mechanistic MAGNITUDE to 7–10 collides with the quasi-static-rate ×
wall-clock budget.

### 4a. N=400 characterisation (dt=1e-4, S=10, band 3.0 → 27 basal rim cells)

The 400-cell cohesive spheroid does NOT spread — it slowly COMPACTS:

| frame | step  | A/A0  | maxZ µm | V/V0  | n_teth |
|-------|-------|-------|---------|-------|--------|
| 0     | 0     | 1.000 | 146.6   | 1.125 | —      |
| 1     | 2105  | 1.006 | 146.8   | 1.083 | 20     |
| 2     | 4210  | 0.992 | 146.3   | 1.054 | 22     |
| 3     | 6315  | 0.959 | 145.1   | 1.032 | 24     |
| 4     | 8420  | 0.917 | 143.7   | 1.007 | 23     |
| 5     | 10525 | 0.877 | 142.0   | 0.983 | 16     |
| 6     | 12630 | 0.843 | 140.6   | 0.959 | 13     |

Two superposed effects, both swamping the lamellipodium: (1) the aggregate starts
OVER-DISTENDED (V/V0=1.125, adhesion-distended from stage 1) and RELAXES toward 1.0
under the soft spread cohesion (ω=1e7) — a contraction the basal traction can't
oppose; (2) the cohesive bulk pulls the periphery inward faster than the ~13–24
crawling leading nodes (of only 27 basal rim cells, vs 400 cohesive cells) pull it
out. NOT a force bug: the SAME engine spreads in isolation (CPU N=4/7 → A/A0 1.05–
1.22). The top-down silhouette is set by the spheroid EQUATOR, so A/A0 cannot grow
until the WHOLE ball slumps into a disk wider than its equator — basal-rim crawling
alone cannot flatten a 400-cell ball against turgor + cohesion in feasible steps.
S=10 vs S=1 barely changes the rate (the catch-gate ties advance to the membrane's
force-velocity), confirming the spread rate is force-balance-limited, not clock-
limited.

## 5. PI decision points

1. **Kinetic acceleration S**: physiological S=1 is wall-clock-infeasible for A/A0
   7–10; a documented S>1 (the old foundation used 6e5) trades quasi-staticity
   (V/V0 drift). What S / V/V0-tolerance is acceptable, or is the magnitude a
   declared structural/compute limit (FORM-only, as Layer-2)?
2. **Collective vs single-cell**: A/A0 7–10 is COLLECTIVE (contact-disk expansion +
   ball slumping into a monolayer + dynamic rim as cells descend), not single-cell
   flattening (the engine's single-cell limit is ~2). Dynamic rim (descended cells
   start crawling) + de-cohesion are needed for the full magnitude — build them?
3. **Force balance** to flatten at V/V0≈1: turgor (resists flattening) vs tether
   (cap = lit traction) vs settle/substrate — tune within lit bands or surface.
