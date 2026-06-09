# H.7 active-γ ARCHITECTURE thrust — autorun log (2026-06-09)

PI fork: active-γ closed (generation/architecture-bound). PI direction: "아키텍처
갈아엎어야 되는 거 아닌가?" (scope-pivot rejected = other session's lane). Thrust =
cortex CONSTRUCTION architecture as the transmission lever. Loop: hypothesis → measure →
record → commit.

## LOOP 1: load-path vs connectivity STATIC sweep (cheap, no contraction run)
Built the connected-mesh seed (no sim, ms each) across z_struct / L_long_mean / bundle_mult.
Scale-corrected: percolation needs full ×40 n=1000 (giant 99%; n=300 fragments). RESULT: load-path
is geometrically PINNED ~1.3–1.5 seg at production knobs. `bundle_mult` (not filament length) is the
dominant knob — de-bundling 2→1 lengthens spans but collapses L/lc 5.1→2.5 (violates L/lc≥5.9 gate)
and drops bundling stiffness. A connected long-load-path region exists (z=1.8, L=10µm, bundle=1 →
giant 95%, mean span 2.1 seg ≈1.0µm = ~2× baseline) but is shallow; z=1.0 reaches 2.5 seg but
fragments (giant 78%). ⇒ load-path & connectivity/stiffness are COUPLED in the point-crosslink+bundle
construction; the achievable ~2× is far short of the 30–80× gap and of Chugh/TruongQuang filament-
length-scale (several-µm) overlap. PI's instinct SUPPORTED: real lever = connectivity-primitive
rebuild (carry connectivity on extended parallel bundles/overlap, crosslink bundle ENDS), not a knob.
Tool h7_loadpath_architecture_sweep.py; doc H7_ARCHITECTURE_LOADPATH_2026-06-09.md. NEXT: bound the
transmission payoff — loading-phase WALL-A (g_actin/g_myo) at baseline vs long-load-path knobs.

## LOOP 2: uniform-vs-faithful mesh + WALL-A loading sensitivity (cheap, n=1000 loading)
Wiring the WALL-A test surfaced that production active-γ/Gate-A/B/force-budget ALL run the UNIFORM
cortex (faithful_connected_mesh unset → 3µm uniform fil, no Arp2/3 branches); the bimodal "faithful"
mesh (branches+bimodal lengths) exists in code but is UNUSED for γ. Step-1 sweep was the faithful path.
RESULT: (5) Arp2/3 branching DECOUPLES connectivity from crosslink density — uniform FRAGMENTS at low
z/bundle (giant 2.3% at z=1.8/b=1) so can't lengthen load-path at all; faithful stays connected at the
same sparse knobs (branches carry connectivity) → the long load-path is ONLY available on the branched
mesh. (6) LOADING phase is load-path-INSENSITIVE: faithful baseline vs faithful long-path give identical
g_myo/g_IK/gen_force/g_soft — loading g_soft is dominated by PASSIVE prestress (seed+turgor), not myosin,
so it can't isolate transmission. ⇒ WALL-A transmission test REQUIRES a contraction run (binder-bound GPU
~30min/config). Plumbed --cm-z-struct/--cm-bundle-mult/--faithful into h7_active_force_budget (topology
only, no gate/datum). NET: architecture lever is real+partly-coded (faithful mesh + sparse crosslinks)
but achievable load-path ~2× ≪ 30-80× gap, and its payoff needs the GPU contraction run. Cheap levers
exhausted → PI fork: (A) GPU contraction WALL-A test faithful base vs long-path (decisive, gbook deploy=
PI remote-overwrite), (B) commit to faithful-mesh rebuild as physically-correct architecture, (C) accept
architecture-bound conclusion. HALT→PI. doc H7_ARCHITECTURE_LOADPATH_2026-06-09.md §step2.

## LOOP 3 (PI chose B: faithful-mesh rebuild): execution hits a gate-contract wall
Executing (B) within the hard rules: the cortex-rebuild gates are z∈[3.0,3.5] (Kadzik-Munro), giant≥0.9,
L/lc≥5.9 (Flormann). At physiological L_long=5µm the faithful configs passing ALL THREE (z=2.8/bundle=3 →
z_real 3.30, giant 98%, L/lc 6.3) give mean span **1.26 seg — NO longer (shorter) than uniform ~1.5**:
bundle=3 (needed for L/lc≥5.9) clusters anchors. The long load-path (2.1 seg) needs z_real 2.67 (<3.0) AND
L/lc 1.7 (<5.9) = a gate-contract change on TWO physiological gates → PI sign-off (no inline loosening).
⭐DECISIVE: the short transmission load-path is a CONSEQUENCE of physiological cortical connectivity (real
cortex IS densely crosslinked z~3-4, L/lc~6) — a faithful rebuild REPRODUCES it, doesn't remove it. So
Chugh/TruongQuang overlap tension can't come from longer inter-crosslink spans; it needs a MYOSIN-OVERLAP
mechanism (myosin transmitting along the bundle/antiparallel-actin contour, independent of crosslink
spacing) = a deeper myosin-model change than a mesh rebuild. HALT→PI: (i) adopt faithful for fidelity +
GPU-contraction γ within gates, (ii) authorize gate-contract change (un-physiological per this result),
(iii) reframe lever as myosin-overlap mechanism. doc §step3. No gate/param changed.

## LOOP 4 (PI chose iii: myosin-overlap mechanism): DESIGN doc + per-head force diagnosis
Read myosin.py grip_walk in full. DECISIVE force-budget diagnosis of the active-γ floor into 3 compounding
deficits: (D1) per-head force capped at k·r≈0.74pN ≪ F_stall 8.48pN (~11× under) — grip_walk delivers
F=k_head_actin·min(s_grip,r); once s_grip≥r it SATURATES at k·r; the Hill stall uses k_series=k/2 (the soft
perpendicular offset spring k_head_spring=1pN/µm sits IN SERIES with the power stroke), so the head NEVER
reaches F_stall (would need min(s,r)≈17µm) → walks to the 2ℓ0 cap delivering only k·r. Real cross-bridge
stiffness ~0.3-2 pN/NM = 300-2000 pN/µm — the model's 1pN/µm CONFLATES the structural offset spring with
the load-bearing cross-bridge. With stiff k the Hill stall binds (s_grip≈17nm → F_stall) → ~11× more force.
(D2) no overlap accumulation — local shunting to nearest crosslink, heads don't sum along contour. (D3)
engagement ~11% geometric. Net envelope(18-36×)×1/9×1/11→~260× under band. DESIGN doc written
(H7_MYOSIN_OVERLAP_MECHANISM_DESIGN_2026-06-09.md): O2 overlap-accumulation (in-lane, opt-in) + O1
cross-bridge stiffness (magic-number/PI trigger) + sanity gates + stages M1-M4. ⭐O1 = dominant fixable
deficit AND a magic-number trigger → HALT→PI (stiffness datum sign-off) per hard rules.

## LOOP 5: ROOT CAUSE — cortex head stiffness 1000× softer than the canonical shared motor
The H.3 brief says cortex/myosin shares the Stam-Hocky/AFINES motor via bridge/motor.py (H.4). But the
configs DISAGREE 1000×: bridge/motor head_spring_k=1.0e-3 N/m (1 pN/nm, "AFINES motor stiffness") vs cortex
k_head_spring=1.0e-6 N/m (1 pN/µm), and k_head_actin=1e-6 ("= k_head_spring"). Literature cross-bridge
stiffness ~0.3-2 pN/NM (Veigel/Kaya/Finer) = 300-2000 pN/µm — the canonical 1e-3 is right, cortex 1e-6 is a
pN/µm-vs-pN/nm UNIT SLIP in the H.3 brief literal. This IS deficit D1: at k=1e-6 the Hill stall needs
min(s,r)≈17µm (unreachable → head delivers only k·r≈0.74pN); at the canonical k=1e-3 the stall binds at
≈4nm → head delivers F_stall (~11× more, dominant fixable deficit, makes grip_walk physically correct).
It's a brief-literal CONTRACT value → PI sign-off + literature anchor (NOT inline). CFL: at k=1e-3,
k_backbone=1e-2 → τ_backbone≈39ns≈3×dt — the "far above CFL" note breaks; re-derive (cortex
cfl_safety_factor, not frozen integrator/) or use M-SHAKE rigid. HALT→PI for the contract sign-off.
doc §7. No param changed.

## LOOP 6 (PI: "2D 표면 메시 최적화 (a)(b) 적대적 검토"): 5-agent adversarial workflow
PI asked whether the long-ago "2D surface mesh to optimize the cortex" directive was implemented (it was NOT —
cortex network is pure-3D bead-spring; surface_manifold.py is geometry-only + unused in the cortex/γ path; the
2026-06-07 H7_CORTEX_AS_MESH explicitly REJECTED mesh-as-cortex as coarse-graining). Ran workflow wf_c27dc37d
(advocate-A ∥ advocate-B → cross-adversarial audit ∥ → synthesis). RESULT: (A) search/coord optimization
SURVIVES (moderate) but ORTHOGONAL to γ + CPU search 7-30× SLOWER than scipy cKDTree at ×40 (win only GPU-
resident); real value = deformable coordinate frame for spreading cell; NEW hole = stale bead_tri drops 21-56%
binding pairs on deforming cortex (needs refresh-vs-drift gate). (B) 2D-network reorg DOES NOT SURVIVE (major):
0× of the gap, its geodesic-vs-chord claim ILLUSORY (kring_for_reach is broad-phase over-selector → identical
chord set per VG-6), 4-7d rebuild + fidelity-slip risk → A strictly dominates B; decide A-or-nothing not A-vs-B.
⭐CRUX (both audits+synth): NEITHER fixes transmission; the actual lever is the SEPARATE stiffness unit-slip
contract correction (k_head_actin 1e-6→1e-3 verified in-config phase1_h3.yaml:331,333 vs phase1_h4.yaml:98;
crosslink k_intra=k_attach=1e-7 re-anchor; crosslinkers.py:514 documents the softness) + O2 myosin-overlap.
3 independent PI decisions: (1) transmission stiffness fix = urgent/separate/do-regardless; (2) A-or-nothing;
(3) if A: scope as spatiality not speed, GPU-resident, add staleness gate, defer k_conf. doc
H7_2D_MESH_AB_DECISION_2026-06-09.md.

## LOOP 7 (PI chose A): ManifoldIndex Step-1 — built + VG-6 PASS + staleness gate HOLDS
PI: "A" (manifold as search/coord accelerator, scoped: spatiality not speed, GPU-resident eventual, staleness
gate, defer k_conf, NOT a γ fix). Built cortex/manifold_index.py (ManifoldIndex wrapper over surface_manifold;
persistent bead→patch map + geodesic k-ring candidate gather + drift-staleness policy; ZERO mechanics/γ DOF).
Validation h7_manifold_index_validate.py on the REAL connected-cortex cloud (3279 beads, reach 841nm):
VG-6 candidate-set identity = **PASS** (index within-reach pairs byte-identical to global cKDTree, missed=0
extra=0) at subdiv 3(1280 tri,k=6) AND subdiv 4(5120 tri,k=7). ⭐STALENESS GATE: the audit-predicted 21-56%
miss DOES NOT materialize — geometry-derived conservative kring_for_reach(reach+2·circumradius+safety) absorbs
~1-patch drift: miss ~0% even at 3× drift_threshold (worst 3 pairs=0.02% at subdiv4). So the SAME over-coverage
that guarantees VG-6 ALSO gives drift-robustness → correctness SAFE to wire. Cost: that over-coverage = a large
pool = CPU-slow (matches the audit's "slower than cKDTree" → speed needs GPU-residency). 5 unit tests PASS
(VG-6, staleness-below-gate, derived-k, empty/boundary, no-mechanics fidelity guard). NEXT (pending): cupy
GPU-resident gather + wire the 4 binder cKDTree sites (assert identical bind/break counts) OR the deformable
coordinate-frame use (normals for enclosed_volume/membrane_surface). Transmission stiffness fix stays a
SEPARATE PI-gated track (unchanged). doc H7_2D_MESH_AB_DECISION §A.

## LOOP 8 (PI: "부유 세포 말고 부착 세포 기준 — 계속 이렇게 진행"): ADHERENT-VENTRAL pivot DESIGN
PI reframe (from a cross-session talk, flagged important for me): stop modelling a SUSPENDED/floating cell;
model the ADHERENT cell. PI's diagnosis: the two layers were geometrically inconsistent — S1(surface_manifold=
curved icosphere=suspended geometry) vs B1(filaments=flat disk); the "curvature option" was a SPURIOUS artifact
of that S1/B1 mismatch. Real physics: adherent cell's ventral surface is FLAT (FA+ventral SF+traction
filaments live there), only apical rounded. Fix: S1→flat ventral, consistent with B1, curvature gone.
⭐CONNECTION to my γ work: the pivot RESOLVES the load-path half (the suspended sphere FORCED an isotropic
densely-crosslinked cortex → load-path ~1.3 seg; the adherent traction structure is ventral STRESS FIBERS =
aligned + FA-end-anchored = the long load-path + axial coherence the sphere lacked — we measured the wrong
structure at the wrong operating point) but NOT the soft-coupling half (k_head_actin/k_intra unit-slips are
geometry-independent, same myosin/α-actinin in the SF → separate PI-gated track stays). Observable shifts:
suspended spherical cortical-tension γ (Hosseini) → adherent substrate TRACTION (PI platform: A/A0, integrin
β1). Current state: FA wired but seeds a south-cap contact on the SPHERE (still curved); substrate exists; NO
ventral SF (build needed). Wrote design H7_ADHERENT_VENTRAL_PIVOT_2026-06-09.md (consistent flat geometry +
ventral-SF structure + measurement plan + sanity gates + lane: mine=cortex/SF/myosin on flat ventral, consume
FA/substrate interface, don't touch other session). NEXT: Stage 1 consistent flat-ventral geometry (no S1/B1
mismatch), then ventral SF build, then traction measurement. Soft-coupling fix = separate track.

## LOOP 9: ADHERENT pivot Stage-1a — consistent flat-ventral geometry (S1=B1 flat, no mismatch)
Added press_onto_substrate(manifold, z_basal) to cortex/manifold_index.py (set_verts only, no edit to shared
surface_manifold; clamps verts below z_basal onto the substrate plane → flat ventral + rounded apical = the
adherent shape PI specified). Proved CONSISTENCY: pressed flat-ventral S1 + a flat-ventral B1 actin cloud →
ManifoldIndex VG-6 candidate identity HOLDS (test_adherent_flat_ventral_consistency) — the S1-curved/B1-flat
mismatch (and its spurious "curvature option") is gone. 6 manifold_index tests pass. Geometry only, no
mechanics/γ, in-lane (cortex/). NEXT: Stage 1b flat-ventral actin construction (cortex on the flat disk),
Stage 2 ventral SF (aligned bundle + myosin + α-actinin + FA-anchored ends), Stage 3 traction measurement.

## LOOP 10: ADHERENT pivot Stage-2a — ventral stress-fiber LAYOUT (aligned + FA-anchored ends)
Built cortex/ventral_stress_fiber.py: generate_ventral_sf_layout — explicit bead-spring actin bundle on the
flat ventral plane, ALIGNED along x̂, MIXED polarity (Hotulainen-Lappalainen 2006: ventral SF graded/mixed so
bipolar myosin finds antiparallel overlap to contract), thin bundle cross-section, FA-ANCHOR beads = the two
extreme-x ends of each filament (the substrate-pin / traction-reaction set). This is the aligned+end-anchored
architecture the suspended isotropic sphere lacked (pivot §2/§6). Fine-grained (bead-spring, no lumped bundle).
6 sanity-gate tests PASS (dims/topology, ventral-above-substrate, aligned-along-axis, mixed-polarity,
anchors-are-ends, boundary). Geometry only; myosin/α-actinin = Stage 2b, traction (FA-anchor reaction) = 2c.

## LOOP 11: ADHERENT pivot Stage-2b-1 — ventral SF traction SCAFFOLD (build+run+readout) [WIP]
Built scripts/h7_ventral_sf_traction.py: assembles the ventral SF bundle (layout → HOOMD snapshot: sf_actin +
fa_anchor types, sf-bond Harmonic, sf-angle straight, WCA) + BAOAB, FA anchors held by OVERDAMPED HIGH-DRAG
(γ_anchor=1e4·γ_b, rigid-substrate limit — NOT a re-pin: a hard re-pin after BAOAB is energetically
inconsistent + integrator/ is PI-frozen so no filter). Traction readout = Σ|axial backbone tension| at anchors
(from bond geometry T=k·(r−r0), robust vs HOOMD force-access timing). BUILD+RUN+READOUT all work. ⚠️FINDINGS
for 2b-2: (1) a fiber at FULL CONTOUR extension is TAUT → intrinsic thermal/entropic tension at rest (passive
≠0, expected — a taut WLC carries tension); (2) single-snapshot |T| is noisy (139/109/581 pN across times =
fluctuation, not divergence); ⇒ the Stage-2c SCIENCE measurement must be DIFFERENTIAL (myosin_ON−myosin_OFF) +
TIME-AVERAGED, and/or place FA separation with slack. Myosin (actin-aware Stam-Hocky, reusable) + α-actinin NOT
yet added (2b-2). NEXT: 2b-2 add myosin+α-actinin on the bundle; 2c differential time-averaged traction vs
myosin (+ soft-coupling-fix overlay). Lane clean (consumed nothing of FA/substrate internals; integrator/
untouched; no other-session files).

## LOOP 12 (PI: "단위 슬립 교정하고 나서 마이오신"): stiffness fix is COUPLED to grip_walk r0 conv.
Applied the PI-authorized cross-bridge stiffness correction (k_head_spring/k_head_actin 1e-6→1e-3, anchored to
canonical bridge/motor + Veigel/Kaya) to phase1_h3.yaml and ran h7_active_force_budget (n=1000, mesoscale).
RESULT: per-head bond tension EXPLODED to ~322 pN (F_stall 8.48 pN, 38×), g_myo jumped 135×, γ_soft 114×→14×
under band — but UNPHYSICAL. ROOT: grip_walk attach bond r0≈0 → a freshly-bound head at r~76nm carries k·r
(huge at stiff k; tiny at soft k). A real cross-bridge binds force-free + generates ~F_stall via the nm-scale
power stroke, NOT k·r. Transmission channel (g_actin) stayed floored → confirms PI's point that the CROSSLINK
(k_intra=k_attach=1e-7) is the SEPARATE soft transmission link. ⇒ unit-slip fix CANNOT be config-only; correct
fix is ATOMIC: (1) k=1e-3 + (2) grip_walk r0=bound-length (force=k·s_grip stroke, capped at F_stall, not k·r)
= O2 + (3) CFL re-derive (τ_backbone~39ns) + (4) crosslink k re-anchor. Reverted config to 1e-6 (config-only
1e-3 = over-tensioned broken state). Physical end-state: per-head→F_stall, γ→full-stall envelope (18-36× under
= density/overlap gap), transmission gated by crosslink fix. doc H7_MYOSIN_OVERLAP_MECHANISM_DESIGN §8. HALT→PI:
the unit-slip fix IS the myosin grip_walk surgery (they're coupled, not sequential) — proceed atomically?
