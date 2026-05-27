# KU-3.5 myosin binding-geometry diagnosis (2026-05-28 /loop autonomous)

**Status:** KU-3.5 (cortical tension) is blocked because cortex myosin does not bind actin → no contraction. Root-caused below (measurement + Codex read-only diagnosis agree). The fix is a **단계-4 myosin placement/binding mechanism change affecting KU-3.5/3.1/3.18 — PI decision required** (not a magic-number patch). The constrained-BD foundation that makes these gates feasible is validated + committed (see `outputs/h3/REPORT.md §constrained-BD 역량`).

## Symptom (measured, 150-fil cortex+myosin)
- myosin head ↔ nearest cortex-actin bead: **median 824 nm** (min 92 nm) vs `head_actin_max_bind_dist = 50 nm` → **0 / 2000 heads in binding range**.
- Over 60 000 steps: `n_engaged ≈ 1`, `n_step_advances = 0`, cortex radius 10.10 → 10.12 μm (**no contraction**).
- The constrained-BD machinery itself is fine: warm-up handoff → stable 54× run, M-SHAKE drift ~1e-15.

## Root causes
1. **Radial head offset (primary).** `generate_cortex_myosin_layout` (cortex/myosin.py:413,416) places heads at `cm ± head_rest_length · n`, where `n = center/R_cell` is the **radial (membrane-normal)** direction. This pushes heads 200 nm radially OFF the actin shell. Fix: offset in the **tangent plane** — `w = normalize(cross(n, u))` (u = backbone axis) — so heads stay in the shell where actin is.
2. **Placement far from actin (primary).** Minifilament centers are sampled uniformly on the sphere, independent of actin. Cortical actin is sparse (1000 fil × 3 μm on 4πR² ≈ 1256 μm² → ~1.5 μm filament spacing), so a random head is ~µm from the nearest filament even with an in-plane offset. Fix: place minifilaments **from actin geometry** — bridging nearby actin segments (biology: myosin minifilaments crosslink/slide adjacent cortical actin).
3. **Bead-center binding (secondary).** Binding searches actin **bead centers** (cortex/myosin.py:747, cKDTree). Bead spacing ℓ₀ = 500 nm → a head within 50 nm of a segment MIDDLE is ~250 nm from either endpoint bead. Fix: **segment-projection** binding (capture 50–63 nm), or — as a documented grid-derived approximation if keeping bead-center — search radius √((ℓ₀/2)² + r_capture²) = √(250² + 50²) ≈ **255 nm** (a discretisation-derived endpoint radius, NOT biological motor reach).
4. **Construction LJ overlap (secondary).** myosin-myosin intra-LJ overlaps at construction blow up even the cfl-dt warm-up at 150 fil (120 fil was benign). Fix: overlap-free myosin placement + robust warm-up (energy-min / dt-ramp).
5. **nlist exclusion cap (tertiary).** Widening the capture radius (diagnostic at 1000 nm) hit HOOMD's compile-time exclusion max (7 bonds/particle) — too many head-actin bonds per particle. Constrains how many heads bind one actin bead; a consideration once placement is fixed.

## Why contraction also needs crosslinkers
Even with binding fixed, a myosin minifilament bound to a single filament only slides along it. Net **cortical contraction** requires the network connectivity (KU-3.19 crosslinkers, already implemented) to transmit motor force into shell tension. KU-3.5 is therefore a full-cell assembly: cortex (rigid backbone) + myosin (binding fixed) + crosslinkers + ERM (shape), at the constrained-BD fast dt.

## Recommended PI decision sequence
1. **myosin placement fix** — tangent-plane head offset + actin-geometry-aware minifilament placement (bridging nearby actin) + overlap-free construction.
2. **binding** — segment-projection (capture 50–63 nm) or grid-derived 255 nm bead-center radius (documented).
3. **assemble full cell** — + crosslinkers + ERM (re-derive k_ERM for the fast dt: at dt=0.001·τ_bend=0.7 μs, the ratified k_ERM=1e-4 gives τ_ERM=3.9 μs which re-violates CFL → soften to ~5.6e-5, or use a smaller dt factor).
4. **tension measurement protocol** — Sanity-Gate measurement-protocol decision: method-of-planes (sum bond/motor/constraint tensions crossing a diametral cut / 2πR) vs box-virial surface tension. Note actin "tension" under rigid constraints is a Lagrange multiplier, not a harmonic-bond force. Target γ ≈ 0.5 mN/m ± 30% (KU-3.5).

Constrained-BD makes the run feasible (KU-3.5 ~0.1–0.5 s sim = ~16 min–1.3 h CPU at 54×); only the myosin binding + assembly + protocol above remain.
