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
