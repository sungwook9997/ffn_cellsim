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
