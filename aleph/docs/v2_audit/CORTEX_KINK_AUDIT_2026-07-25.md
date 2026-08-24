# Cortex fiber-kink audit — source of the fine-75nm GATE-A bending-force plateau

**Date:** 2026-07-25
**Scope:** DIAGNOSTIC / read-only. No physics or runtime code changed. This doc is the only committed artifact.
**Trigger:** The fine-75 nm cortex GATE-A residual plateau was force-decomposed (completeness 2e-15) to the
transverse **bending** force — ~5 pN at the worst node, dominant over crosslink/steric/pressure/Jᵀλ. 119/200 of
the worst nodes sit at **crosslink loci** and carry ~**5° local turn angles** vs the ~0.57° expected for a fiber
smoothly following the R = 7.5 µm sphere. The question: are those 5° kinks an **artifact** of construction/
relaxation, or physically intended?

## Verdict (one line)

**Artifact.** The cortex fiber PATHS are built smooth (great-circle arcs, 0.57° natural turn). Crosslink
placement does **not** move nodes. The kinks are created by the **`overlap_free` WCA build-relaxation**, which
pushes individual crossing nodes ~few-nm transverse off their smooth arc *at the crossing/crosslink loci* — a
nm-scale displacement that is negligible at the 0.5 µm coarse mesh but becomes a ~5° kink (and a diverging
bending force) at the 75 nm fine mesh. This is exactly the audit's ~1/L bending divergence with refinement.

---

## (a) How cortex fiber paths are constructed — SMOOTH, not kink-prone

Each cortex fiber is laid as a **centred great-circle arc** on the sphere:

- `aleph/laws/cortex_assembly.py:241-254` — `_great_circle_arc(com_dir, tangent, R, n_beads, seg)`:
  `p(s) = R[cos(s/R)·ĉ + sin(s/R)·t̂]`, beads at arc length `s ∈ [−(n−1)seg/2, +(n−1)seg/2]`. Every node sits
  exactly on radius R; consecutive segments turn by the constant geodesic angle `seg/R`.
- `aleph/laws/cortex_assembly.py:313` — the fiber list is built purely from these arcs.
- Natural per-node turn at the fine mesh = `seg/R = 0.075/7.5 = 0.010 rad = 0.573°` — **this is precisely the
  audit's "~0.57° expected".** At the coarse mesh it is 3.82°/node (still smooth: uniform, no kink).

The FINE 75 nm mesh is built **natively** at `cortex_seg_um = 0.075` (≈41 beads/3 µm fiber), not by resampling a
coarse mesh — `aleph/scripts/ac_cortex_mesh_ab.py:92`, `aleph/components/incumbent/assemble.py:186-203, 359-401`. So the
un-relaxed fine fiber is a smooth arc with 0.57° turns. **The construction of the path itself introduces no
kinks.** The length-distribution option (`length_dist="exponential"`) only changes bead count per fiber; each
fiber is still a smooth arc.

## (b) Crosslink placement and the relaxation — where the kinks actually come from

### Crosslinks do NOT snap/weld/move nodes (ruled out)

Crosslinks are formed by **pairing two near nodes on different fibers** and installing a spring that is
**force-free at the existing separation** — no node is moved to its partner:

- `aleph/laws/gamma_floor.py:159-187` — `_cross_fiber_pairs`: KDTree `query_pairs`, keeps cross-fiber pairs,
  returns index pairs **only** (positions untouched).
- `aleph/laws/weave.py:182-191` and `aleph/components/weave/regions.py:207-211` — `xl_i, xl_j` are those indices;
  `xl_rest = ‖pos[xl_j] − pos[xl_i]‖`, i.e. the **actual current distance** → the crosslink is a zero-force
  tether, not a weld.
- `aleph/components/weave/woven_cell.py:437-438` — after the relax the crosslink rest lengths are **recomputed from
  the relaxed positions**, so the crosslink stays force-free even post-relax.

Therefore a crosslink cannot bend the fiber. The one genuine node-snap in the codebase is the Arp2/3 branch
weld `p[d0] = p[bn]` (`aleph/components/weave/regions.py:442-443`), but that snaps a daughter fiber's **base/end** node
onto its mother's branch vertex (intended θ₀=70° branch topology), not a mid-fiber node at the formin-cortex
crosslink loci the audit flagged.

### The `overlap_free` WCA relaxation deflects crossing nodes → the kink source (confirmed)

After the smooth arcs are placed, the overlap-free build runs a **soft-sphere WCA push**: every different-fiber
node pair closer than the contact distance `r_c = 2^(1/6)·σ_EV` is moved **symmetrically apart by half the
overlap**, iterated to convergence:

- `aleph/laws/cortex_assembly.py:154-186` — `_resolve_cross_fiber_overlaps` (in-builder path).
- `aleph/laws/cortex_assembly.py:189-216` — `overlap_free_shell` (post-placement path).
- `aleph/components/weave/regions.py:388-445` — `_relax_arp23_overlaps` (Arp2/3 leaf, same push).
- `aleph/components/weave/woven_cell.py:433-438` — the **unified cross-region WCA relax** over the combined cortex.
- **Enabled by default** and in the audited run: `aleph/components/incumbent/assemble.py:207`
  (`overlap_free_cortex: bool = True`), `aleph/components/incumbent/driver.py:1293`, and explicitly
  `overlap_free_cortex=True` in `aleph/scripts/ac_gate_a_fq_coarse_test.py:584`.

The push target is a **fixed physical distance**: `r_c = 2^(1/6)·0.007 = 7.86 nm` (F-actin excluded volume),
**independent of mesh**. It shoves a single crossing node transversely off its otherwise-smooth arc while its
in-fiber neighbours stay put → a one-node **kink**. Because two fibers that are close enough to be WCA-pushed
(≤7.86 nm) are exactly the fibers close enough to be **crosslinked**, the deflected nodes coincide with the
crosslink loci — explaining the audit's **119/200 worst nodes at crosslink loci**.

### Why it only bites at the fine mesh (matches the measured ~1/L divergence)

The kink ANGLE from a fixed transverse deflection δ scales as `≈ 2·arctan(δ/seg)` — it **grows as the mesh
refines**, while the deflection itself (set by the mesh-independent 7.86 nm WCA scale) does not:

| quantity | coarse 0.5 µm | fine 0.075 µm |
|---|---|---|
| natural turn / node (`seg/R`) | 3.82° | **0.57°** (= audit "expected") |
| δ needed for a 5° kink | 21.8 nm | **3.3 nm** (well inside the 7.86 nm WCA envelope) |
| same 3.3 nm WCA deflection → kink | 0.76° (negligible) | **5.0°** (= audit "worst node") |
| bending stiffness `α = κ/seg³` (κ=0.073) | 0.6 pN/µm | **173 pN/µm** (×288) |

A ~3 nm WCA nudge that is invisible at the coarse mesh (0.76°, below the 3.8° natural turn) becomes a 5° kink at
75 nm, and the Cytosim bending force `F = α·d` on it is amplified ~288× by `α = κ/seg³`. The crosslink spring is
force-free (rest recomputed post-relax), so the entire residual at that node is intra-fiber **bending** — exactly
the term the force decomposition isolated.

## (c) Assessment: ARTIFACT

The ~5° kinks at crosslink loci are an **artifact of the WCA build-relaxation at fine mesh**, not intended
physics. Physically a stiff actin filament (Lp = 17 µm) crosslinked to a neighbour is **tethered** at the
crossing and bends only gently over its persistence length; it does not acquire a sharp 5° corner at a single
node. The relaxation's job — remove sub-σ interpenetrations so the steric force does not saturate at t0 — is
legitimate, but it currently pays for a clean *steric* start with an artificial *bending* residual that the
coarse mesh hid and the fine mesh exposes. Two design facts make it clearly an artifact:

1. The deflection magnitude is tied to the EV contact scale (7.86 nm), **not** to any fiber-bending physics, so
   the induced kink angle is a pure mesh-resolution accident (`δ/seg`), diverging under refinement.
2. The un-relaxed fine mesh is smooth (0.57° turns); the kinks appear **only** after the relaxation moves nodes.

## (d) Recommended fix direction (RECOMMENDATION ONLY — for Lead/PI)

Ordered by preference; all keep the physiological "start overlap-free" intent while removing the fine-mesh
bending artifact. Any of these is a construction/relaxation change, subject to the Gate-1 bit-identical-cortex
contract (the formin-only cell must stay byte-identical → gate behind the fine/Arp2/3 path, not the default
coarse cortex).

1. **Resolve crossings by radial offset, not transverse in-shell push.** The build already disperses fibers
   across the ~0.2 µm cortex thickness (`cortex_assembly.py:304-312`, `overlap_free_shell`); lean on that to
   separate crossing fibers in RADIUS (a whole-fiber rigid shift preserves arc smoothness) and **suppress /
   cap the in-plane per-node WCA push** so no node is deflected transverse to its own arc. A stiff filament
   should step over its neighbour by sitting at a slightly different radius, not by kinking.
2. **Cap the per-node relaxation deflection to a mesh-fraction of `seg`.** Bound each node's transverse move to
   e.g. `≤ seg·tan(θ_max)` with θ_max at the natural geodesic turn (~0.6°), so the relaxation can never inject a
   kink larger than the fiber's own smooth curvature. Residual sub-cap overlaps are left to the runtime steric
   force (they are tiny and physical), instead of being fully removed at build time at the cost of a kink.
3. **Re-smooth each fiber after the relax.** Post-relaxation, project each fiber's nodes back onto a
   low-curvature arc/spline through its endpoints (a bending-energy-minimising smoothing that preserves fiber
   length and the crosslink attachment nodes to within the tether slack). This directly removes the single-node
   corners the WCA push introduced while keeping the crossings resolved.
4. **(Complementary) let the crosslink absorb the offset as tether slack.** Since crosslink rest length is
   already set from the relaxed geometry, allowing a small finite crosslink rest length (rather than driving the
   two crossing nodes to the EV contact and welding the fiber's shape around it) keeps the tether force-free
   without demanding the fiber bend to meet it.

Do **not** simply disable `overlap_free` — that reintroduces the ~64k t0 interpenetrations and the steric
saturation it was built to fix (`cortex_assembly.py:56-62`). The goal is overlap-free **without** transverse
per-node kinks.

---

### File:line index

- Smooth arc path: `aleph/laws/cortex_assembly.py:241-254, 313`
- Fine mesh built natively at 75 nm: `aleph/components/incumbent/assemble.py:186-203, 359-401`; `aleph/scripts/ac_cortex_mesh_ab.py:90-92`
- Crosslink pairing (no node move): `aleph/laws/gamma_floor.py:159-187`; `aleph/laws/weave.py:182-191`; `aleph/components/weave/regions.py:207-211`
- Crosslink rest recomputed force-free post-relax: `aleph/components/weave/woven_cell.py:437-438`
- Arp2/3 branch weld (end-node snap, not the flagged loci): `aleph/components/weave/regions.py:442-443`
- WCA overlap relaxation (kink source): `aleph/laws/cortex_assembly.py:154-216`; `aleph/components/weave/regions.py:388-445`; `aleph/components/weave/woven_cell.py:433-438`
- Enabled by default / in the audited run: `aleph/components/incumbent/assemble.py:207`; `aleph/components/incumbent/driver.py:1293`; `aleph/scripts/ac_gate_a_fq_coarse_test.py:584`
