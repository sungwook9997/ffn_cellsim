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

## LOOP 13 (atomic unit-slip fix attempt): bin scheme CANNOT deliver the stiff-cross-bridge stroke
Implemented the atomic fix (k 1e-3 + grip_walk r0=max(r_bind−s_grip,0) re-binning, 5 edits) then hit a HARD
resolution wall: the stiff-cross-bridge stroke is ~4nm (s_grip_stall=F_stall/k_series) but the myosin binding
range is head_actin_max_bind_dist=330nm. The attach-bond r0 is bin-quantized over [0,330nm]; force granularity
= k·bin_width, so resolving the 4nm stroke (granularity<F_stall) needs n_bins≳165 (impractical); at n_bins=60
granularity=23pN=2.7×F_stall → one re-bin overshoots the stall. ⇒ the unit-slip fix needs a CONTINUOUS per-head
force, not bins. THREE nested findings define the fix: (1) k 1e-6→1e-3; (2) r0=bound-length (force=k·stroke not
k·r); (3) CONTINUOUS custom force F=min(k·s_grip,F_stall) per head (md.force.Custom; attach bond kept only for
Bell-Evans off-rate) + likely shrink the 330nm binding range (unphysical for a stiff bridge). Deeper myosin
redesign than config+bins. ALL edits reverted (no under-resolved intermediate). Fully specified in
H7_MYOSIN_OVERLAP_MECHANISM_DESIGN §9. CFL (τ_backbone~39ns) + crosslink re-anchor remain in the atomic set.
HALT→PI: the unit-slip fix = a custom-force myosin redesign (host-sync/GPU-main territory) — execute as focused
sanity-gated effort? doc §9.

## LOOP 14 (§9 execution, Increment 1): continuous per-head custom force — HERMETIC, PASS
Built the §9 fix the bin scheme couldn't deliver (loop13): opt-in stepping_mode="continuous_stroke".
MyosinHeadForce(md.force.Custom) delivers per engaged head F=min(k_head_actin·s_grip, F_stall) along
the head→bead unit vector (reaction on head, Newton-3 pair) — REPLACES the harmonic attach bond as the
force-bearing element, so the k·r explosion (~322pN at stiff k=1e-3, loop12) is capped at F_stall by
construction. Attach bond demoted to k=0 (Bell-Evans off-rate + nlist exclusion ONLY, no double-count);
Bell-Evans Step 1 uses the delivered force for this mode. Reuses grip_walk binding+s_grip-walk+bipolar
gate verbatim; grip_walk/binned_r0 BYTE-IDENTICAL. 13 new sanity-gate tests PASS (dimensional, boundary
0@s=0 + cap@F_stall, Newton-3 Σforce=0 + F[bead]=-F[head], contractile sign, delivered |F|=F_stall NOT
k·r, legacy bond reg unchanged). NO config/CFL change yet (k stays 1e-6, existing gates green). Commit
8dc2680. NEXT (Increment 2): cell.py wiring (append MyosinHeadForce to integrator when continuous_stroke)
+ config k 1e-6→1e-3 (Magic-Number Block, canonical bridge/motor + Veigel/Kaya anchor) + CFL re-derive
(k_backbone=1e-2 → τ~9-39ns, shrink dt via cortex cfl_safety_factor; integrator/ untouched) + crosslink
k re-anchor (lit-verify α-actinin OR PI surface) → h7_active_force_budget smoke: per-head delivered T
caps at F_stall + γ rises to full-stall envelope. GPU cupy port = follow-up (host-sync now).

## LOOP 15 (§9 Increment 2): wire continuous_stroke into cell.py + force-budget tool — VERIFIED
cell.py appends MyosinHeadForce to the live integrator when continuous_stroke (mirrors compartment-force
swap); tool gets --stepping-mode + reads delivered force = min(k·s_grip,F_stall) per-bond (not k·r).
Smoke (n_fil=80, soft k): builds, custom force participates, 2000 steps NO crash, delivered T≈0 at
s_grip≈0 = correct force-free-at-bind (§6.2). Commit 007b0b6.

## LOOP 16 (§9 Increment 3): stiff cross-bridge k + CFL re-derivation — CFL STABLE, cap demonstrated
MODE-COUPLED stiff k (the stiff k explodes harmonic k·r in legacy modes = loop12, so applied ONLY in
continuous_stroke whose custom force caps): configs/phase1_h3.yaml new k_cross_bridge_continuous=1e-3
(Magic-Number Block, canonical bridge/motor + Veigel/Kaya/Finer 0.3-2 pN/nm); resolver applies it only
in continuous_stroke (grip_walk keeps 1e-6, byte-identical). CFL via reconcile_dt (forwarded
build_baseline_cell→Cell.build; dt_reconcile folds myosin k_backbone → integrator dt scalar, NO
integrator/ edit). binding-range 330nm KEPT (§9 shrink was bin-resolution-specific; r doesn't enter the
custom force → no explosion → no shrink; reach gated by capture_perp 210nm). VERIFIED smoke (n_fil=80,
stiff k, reconcile_dt): B1 lowered dt 1.30e-8→9.22e-10 s (14.1×, binding=myosin.k_backbone τ=9.22ns meso),
integration STABLE (3300 steps, all finite, NO blowup), per-head delivered T=1.15pN tracking k·s_grip
capped under F_stall=8.48pN — NOT the k·r=930pN explosion (loop12), demonstrated end-to-end. Magnitude
(s_grip→F_stall, γ→envelope 17.8× under = KNOWN density/overlap gap) needs the GPU contraction run (M4,
~1e7 steps). Commit ab376de. NEXT: crosslink k re-anchor (KB α-actinin lit-verify) → GPU run.

## LOOP 17 (§9 Increment 4): crosslink k KB lit-verify → PI SURFACE (gate-contract)
KB query (tag_query.py): KB-1.28 (verified, High) sets the crosslink bond stiffness k_xl≈1e-4–1e-2 N/m,
DEFAULT 1e-3 (1 pN/nm). The cortex dynamic_crosslinkers k_intra=k_attach=1e-7 is ~10⁴× below that default
(1000× below the range floor). The config anchor "0.1 pN/µm, KU-3.19 (Furuike 2001)" does NOT match a
stiffness datum — KU-3.19 is kinetics-only (off-rates, no k) and Furuike 2001 is filamin unfolding
kinetics; same mis-attribution class the prior audits caught. ⇒ 1e-7 is very likely the §8 "separate soft
transmission link" slip; re-anchor to KB-1.28 (→1e-3) is literature-first BUT a PRODUCTION-WIDE
gate-contract change (all cortex builds, the other session, already-run Gate-A/B) + flexible-crosslinker
molecular-k vs sim-bond-k is genuinely ambiguous (KB range spans 100×). CFL not tightened (τ_xl=391ns≫
9.2ns myosin dt). Per mission rule (crosslink k = lit-verify→PI-surface-if-uncertain + gate-contract→PI),
NOT changed inline → SURFACED to PI (design doc §10). Generation fix needs NO crosslink change to verify;
it is the downstream transmission lever (smoke WALL-A already 66.8% in connected mesh).
STATUS: §9 generation half (items 1-4) DONE + verified (cap + CFL stable, end-to-end). Item 5 (crosslink)
+ M4 (GPU magnitude run) are PI-gated. HALT→PI for the crosslink contract decision before the GPU run.

## LOOP 18 (§9 Increment 5): crosslink k re-anchor 1e-7→1e-3 (KB-1.28, PI-approved)
PI approved (production-wide). dynamic_crosslinkers k_intra/k_attach 1e-7→1e-3 N/m (KB-1.28 verified
default; the 1e-7 "Furuike 2001/KU-3.19" anchor was a mis-attribution — KU-3.19=kinetics-only,
Furuike=unfolding). Magic-Number Block added. CFL not tightened (τ_xl=391ns≫9.2ns myosin dt). Smoke
STABLE (no blowup). ⚠️production-wide (re-run Gate-A/B). Commit 236101e. §9 atomic set (items 1-5) COMPLETE.

## LOOP 19 (ventral SF 2b-2): wire continuous_stroke myosin onto the SF bundle
build_sf_sim with_myosin: actin-aware Stam-Hocky minifilaments along SF tangents (±x̂), MyosinHeadForce
(capped per-head), SF actin = binding substrate; CFL dt re-derived locally (myosin k_backbone, NO
integrator/ edit). main() = 2c DIFFERENTIAL (passive vs active time-avg traction, cancels taut-WLC
thermal tension). CPU smoke: build + passive(2389±496pN) + active + differential pipeline RUN; binding
STARTS (geometry correct); 0-engagement at short smoke = timescale artifact (binding ~1/k_on=0.02s ≈
5e5 steps at CFL dt), NOT a bug. Commit 7237ddc.

## GPU RUNS LAUNCHED (gbook, parallel; PI both-in-parallel directive)
Mac→gbook rsync to FRESH dir ~/ffn_cellsim_h7run (PI-approved remote-overwrite; gbook's own tree had
stale orphaned edits → untouched). (1) FORCE-BUDGET contraction PID 342791: n_fil=500, continuous_stroke,
400k contract steps, reconcile_dt. PROGRESS tick 25000: engagement 9.7%, meanT=7.56pN (→F_stall 8.48pN,
CAPPED, not k·r), gen=4.1nN — the §9 per-head→F_stall target developing on GPU. (2) SF 2c PID 343006:
n_fil=12, 8 motors, continuous_stroke, 120k+180k steps ×2 builds (passive+active differential). Both
results pending (~40min). NEXT: collect γ magnitude (force-budget) + SF differential traction, then
Notion closeout.

## LOOP 20 (ventral SF 2c REDESIGN): same-seed paired differential — noise floor 300× lower
The loop19 2c differenced DIFFERENT realizations → ~10nN taut baseline didn't cancel (got -2495pN,
noise ±900pN). REDESIGN: MyosinHeadForce gains force_scale; 2c builds force-OFF (scale 0) + force-ON
(scale 1) with the SAME seed → identical bundle/myosin/binding/thermostat counter (HOOMD counter-based
RNG → thermal kicks bit-identical per (timestep,tag) regardless of position), so ON−OFF cancels noise +
taut baseline EXACTLY. CPU smoke (20 motors, short): force-OFF 3052.8±698 vs force-ON 3058.0±706 pN
(~identical, same seed), DIFFERENTIAL = +5.22 ± 2.89 pN — POSITIVE (contractile), noise floor ±900→±2.9pN
(~300× lower). Resolvable; 4 engaged heads in smoke → longer run develops engagement. Commit 04bb6e0.

## GPU RUN: force-budget RESULT (plateaued, tick 150k-275k stable)
per-head meanT = 8.44-8.47 pN = F_stall 8.48 pN → PER-HEAD FORCE AT STALL, CAPPED (not k·r ~900pN). ✓
engagement 10.8% (=§9 D3 ~11% geometric limit). γ_soft 7.8e-3 mN/m (~23× under band) = the §9 density/
overlap+engagement gap (band LOCKED, NOT tuned). §9 GENERATION FIX VERIFIED: old soft-k 0.74pN(11×-under)
+ γ~260× under → now F_stall + γ~23× under. Residual gap = documented structural/scope limit (§9 gate 7).

## GPU RUN: SF same-seed differential launched (PID 344489, n_fil=24, 20 motors, 150k+150k ×2)
Develops engagement for a significant positive traction signal. Pending (~50min).

## LOOP 21 (SF 2c coherent observable) + GPU RESULT: alignment+anchoring is NOT sufficient
Fixed the 2c observable: anchor_traction now returns the COHERENT signed contractile traction (two
FA-ends pulled toward each other; thermal cancels) vs the old Σ|T| (fluctuation magnitude, which gave
the misleading −701pN). Commit dd2ee1f. GPU same-seed coherent run (n_fil=24, 20 motors, 68 engaged
heads): COHERENT force-OFF +3476±752 pN (passive EV/bundling pre-tension), force-ON +3414±891 pN,
DIFFERENTIAL = −61.45 ± 27.74 pN — SIGNIFICANT (|−61|>2σ=55) and NEGATIVE (slackening). 68 heads ×
~1pN ≈ 61pN, so the §9 myosin IS exerting force but it is NOT rectified into end-ward contraction.
⭐ FINDING: the §9-corrected motor generates F_stall (verified on cortex), but on a RANDOM MIXED-POLARITY
end-anchored bundle that force does NOT produce net contractile traction — it slightly slackens. The
adherent-pivot §2 claim ("aligned + end-anchored suffices") is REFINED: the missing ingredient is
SARCOMERIC POLARITY ORGANIZATION (Hotulainen-Lappalainen graded polarity: + ends at FAs, − ends
overlapping at centre, periodic α-actinin) that rectifies bipolar myosin sliding into directional
contraction. The same-seed COHERENT method (±28pN floor, resolved a significant 61pN) is validated as
the reusable traction probe. NEXT (PI-gateable, new construction): graded-polarity sarcomeric SF layout
→ re-measure (expect +contractile). docs design §11.

## LOOP 22 (myosin backbone bending) + LOOP 23 (sarcomeric SF) + ⭐DECISIVE RESULT
LOOP22: added opt-in myosin minifilament backbone angle term (rigid-rod fidelity; the brief used
stretch-only = freely-jointed). L_p_myo=17µm derived from Billington 2013 EM-straight (σ_perp≤10% →
L_p≳10µm). Default OFF byte-identical. Diagnostic: backbone doesn't fold either way in this regime
(head-web + stretch + timescale) → faithful but NOT the SF-result cause. Commit 2b68880.
LOOP23: built generate_sarcomeric_sf_layout — Z-bands(α-actinin barbed anchors)↔M-bands(myosin)
periodic, graded polarity (+x left half / −x right half), antiparallel pointed-end overlap at M,
α-actinin Z-disc crosslinks, outer Z = FA anchors. Geometry gates PASS. Commit fd484c4.
⭐⭐DECISIVE GPU RESULT (same-seed coherent differential, 16 motors):
  random mixed-polarity (68 heads): −61 ± 28 pN  (slackening, −0.9 pN/head)
  SARCOMERIC          (49 heads): +131 ± 8 pN  (CONTRACTILE, +2.67 pN/head, 16σ)
The SIGN FLIPPED negative→positive with FEWER heads. {§9-corrected motor} × {sarcomeric structure}
= net contractile traction EMERGES. +2.67 pN/head ≈ the native F_stall (~2 pN) RECTIFIED into the
contractile direction (mixed cancels it). ⇒ adherent-pivot hypothesis CONFIRMED: the active-γ "floor"
on the suspended isotropic cortex was largely a STRUCTURE/observable artifact; the corrected motor on
the correct adherent structure (sarcomeric ventral SF on FAs) produces the platform observable
(traction). The same-seed coherent probe (±8 pN floor) resolved it at 16σ. figure
h7_sf_2c_mixed_vs_sarcomeric.png. NEXT: scale to an SF array → aggregate traction stress (Pa) vs the
PI platform; dynamic α-actinin (vs static Z-disc); cupy gpu_local myosin force port.

## LOOP 24 (2026-06-10): SF-ARRAY scale-up → aggregate traction-stress [Pa] + PI-platform overlay
Scaled the DECISIVE single sarcomeric SF (+131±8 pN, loop23) to a cell-scale ARRAY. Built
generate_sf_array_layout (n_sf parallel sarcomeric SFs across the ventral contact patch, lit
spacing, each an independent FA-anchored validated single-SF unit, WCA-coupled → near-additive;
tracks per_sf_anchor_beads + contact_area). 6 geometry sanity gates PASS. Wired --array into
h7_ventral_sf_traction: aggregate coherent traction → stress [Pa] = ΣFA-reaction / contact
footprint; PI-platform overlay baked in (MCF-7 TFM 102 nN / 63 Pa / 1822 µm², Gil-Redondo 2023
DOI 10.1002/jemt.24368 Table 1 control n=37). CPU smoke (n_sf=3, 900 steps) verified plumbing.
LIT ANCHORS (HARD): FA 43/cell KU-2.4; SF ~20 ≈ FA/2 Hotulainen-Lappalainen 2006; sarcomere
period ~1µm native→~2µm meso Peterson 2004; single active SF ~5-6 nN Kassianidou/Kumar 2017.
Commit 9bb97be. GPU run launched (gbook RTX A5000, ~/ffn_cellsim_h7run): n_sf=4 (2304 beads,
64 minifilaments) primary additivity+aggregate; + n_sf=1 CPU control (reproduce decisive +131
via array path). ETA ~2-3h. KEY invariant: +2.67 pN/engaged-head rectified ≈ native F_stall;
aggregate scales w/ total engaged heads → 102 nN needs ~38k engaged heads = per-SF minifilament
DENSITY lever (= the cortical-γ density/overlap gap, now in the traction observable).

## LOOP 24b (2026-06-10): ⛔ STOP — loop23 "decisive +131 pN" DOES NOT REPRODUCE (sign non-robust)
Before scaling the SF array, a control re-ran the EXACT decisive single-SF config (n_fil=36,
16 motors, equilibrate 150k + contract 150k, 30 samples) across seeds. Per-SF coherent traction
differential is NOT sign-stable: s1=+131 (=decisive exactly), s2=−428, s3=−434, s5=−247, s6=+340 pN.
Across-seed mean −127±156 pN (std 348), 2+/3−, indistinguishable from ZERO. CPU s1=+131.1 reproduces
the GPU decisive → device-independent (not a numerical bug). The "16σ" was WITHIN-realization SEM
(8-23 pN), ~15-40× smaller than the across-seed std (348) = the true error bar. ⇒ loop23 DECISIVE
was a single fortunate seed-1 draw from a high-variance sign-unstable distribution. Scale-up premise
(stable per-SF +131) INVALID → aggregate/[Pa]/Layer-2 seam VALUE suspended (interface design intact).
HALT per stop-rule (verification-failed core-physics). NOT refuted: per-head F_stall fix + the DENSITY
hypothesis — at ~50 heads net direction is placement-noise (small-N ± sum); high-density (thousands of
heads) is untested = the genuine joint test of rectification+magnitude. PI decision (doc §7 BLOCKER):
single high-density SF run (~200-300 minifil) — does per-SF differential converge sign-definite as
N_heads grows? + cheap check that myosin placement enforces antiparallel-overlap engagement. Trail:
outputs/h7/production/h7_sarc_cpu_s{1,2,3,5,6}.json, ensemble.sh. (gpu_s1 + cpu_s4 pending, confirm-only.)

## LOOP 24c (2026-06-10): PI-directed high-density test — does per-SF sign converge as N_heads grows?
PI chose the high-density SF test. Construction-physics insight while sizing it: the myosin 401nm
pairwise-spacing cap means minifilaments pack along the FIBER LENGTH in a thin (400nm) bundle, and
SERIES sarcomeres share tension (don't add force) — so neither lengthening nor thickening-at-fixed-
radius adds parallel force. The force-adding axis is a WIDER bundle (more PARALLEL minifilaments per
cross-section). Sweep: bundle_radius 800/1200/1600nm × n_cross 12/24/48 × n_motors 28/55/100 (P2/P3/P4)
× seeds 1-4 → ~84/165/300 engaged heads (vs baseline ~50). Per level: across-seed mean±std of the
coherent differential → does std/|mean| shrink + sign stabilize as parallel N_heads grows? gbook 12
runs parallel, ETA ~2h. Also running a local placement-sanity check (does the sarcomeric construction
enforce antiparallel-overlap engagement, or same-polarity = the sign-noise source).
