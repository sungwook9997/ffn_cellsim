# KU-3.5 myosin grip-walk — IMPLEMENTATION RECORD (2026-05-31)

> Companion to `KU35_GRIP_WALK_DESIGN_2026-05-31.md` (the ratified design) and
> `KU35_FLOOR_ROOT_CAUSE_2026-05-31.md` (the diagnosis). This records what was
> actually built after PI ratification, the one design nuance that had to be
> reconciled at the HOOMD level, and the validation status.

## 1. PI ratification (2026-05-31)

All four flagged decisions ratified on the brief's recommendations:

| # | Decision | Ratified |
|---|---|---|
| 1.4 | Walk representation | **Continuous sub-bead `pos_a_end`** (float `s_grip`) |
| 2 | Cortex polarity | **Option A — minus = bead 0 (label-only, no topology change)** |
| 3 | Bipolar binding rule | **Antiparallel `sign(m̂·û)` gate + different-filament clause** |
| 4 | Force-ceiling param | **No param change; r0=0 idealized baseline in grip-walk mode** |

## 2. What was implemented (additive, opt-in)

A new `stepping_mode` config on `ResolvedCortexMyosin` (default `"binned_r0"` =
the legacy proxy; `"grip_walk"` = the new motor). With the default, the builder
and updater are bit-for-bit the prior behaviour — the full regression suite and
every committed result are untouched. Changes, all in `cortex/myosin.py` +
the `cell/cell.py` call site:

- **`stepping_mode` field + resolver parse/validate** (`binned_r0`|`grip_walk`).
- **Per-head grip-walk state**: `_head_bound_bead_pos`, `_head_bound_filament`
  (the `(filament, pos)` decomposition the walk increments), `_head_grip_s`
  (the continuous commanded sub-bead stretch, the AFINES `pos_a_end`).
- **Helpers**: `_tag_to_fil_pos` / `_bead_tag` (exact fixed-N bead-tag map),
  `_walk_toward_minus` (Option-A: decrement toward bead 0, clamp = AFINES
  minus-end latch), `_bipolar_accepts` (decision 3 gate).
- **Step-2 binding**: in grip-walk mode, apply the bipolar gate before accepting
  a candidate, then init the grip at the bound bead with `s_grip = 0`
  (force-free construction).
- **Step-3 kernel** (`grip_walk` branch): `F = k·min(s_grip, r)` (honest force
  from the commanded stretch, re-derived from the current `r` every tick so it
  does NOT relax — the proxy's fatal flaw), Hill `v(F)`, `s_grip += v·batch_dt`,
  and on `s_grip ≥ ℓ₀` re-target the bond to the next minus-ward bead
  (material transport) with the shared per-bead degree budget re-checked at
  walk time. The legacy `binned_r0` branch is the original code, untouched.
- **Unbind reset** clears the three new arrays.
- **`cell.py`** threads `ell0_cortex = p_cortex.rest_length` and
  `cortex_beads_per_filament` into the updater (harmless in `binned_r0`).

## 3. Design nuance reconciled at the HOOMD level (PI-FLAG)

Decisions 1.4 (continuous sub-bead `r0_eff`) and 4 (single r0=0 attach type)
**interact**, and the literal "single r0=0 type" cannot realise a *continuous*
commanded stretch: HOOMD's `md.bond.Harmonic` carries `r0` **per bond TYPE**,
not per bond, so a genuinely continuous per-bond `r0_eff` is impossible without
multiple types. Decision 1.4 (continuous) is the load-bearing one (PI: "integer
-only cannot show force at diagnostic scale"), so it wins.

**Realisation:** the continuous commanded stretch is carried by **requantising
`r0_eff = max(r − s_grip, 0)` into the EXISTING bin set every tick** (no new
bond types, no `k` change — honouring the actual param-no-change intent of
decision 4). r0=0 remains the idealised floor (bin 0). This is mechanically
distinct from the condemned proxy on two counts: (a) `r0_eff` is re-derived from
the *current* `r` each tick, so `F = k·s_grip` is HELD as the bead approaches
(the proxy let `F` relax to 0); (b) the grip RE-TARGETS to downstream minus-ward
beads on overflow (the proxy never transported material). The bin range
`[0, max_bind_dist = 330 nm]` covers the `r0_eff` regime (≈ `r − s_grip`), with
±16.5 nm (±0.016 pN) quantisation error.

> **For PI:** this is the one deviation from a literal reading of ratified
> decision 4. It is forced by the HOOMD per-type-r0 constraint and serves
> ratified decision 1.4. No new magic numbers; `k` unchanged.

## 4. Validation status

- **Tier 2 (unit + regression):** 6 new grip-walk tests green (helpers, walk
  minus-end clamp, bipolar orientation gate, different-filament clause, the
  requires-geometry guard, default-is-binned_r0). **Full regression suite
  GREEN** (binned_r0 default byte-unchanged — additivity confirmed).
- **Tier 1 (sustained-tension micro-diagnostic):**
  `scripts/h3_ku35_gripwalk_tier1.py` — A/B (binned_r0 vs grip_walk) at a
  DOCUMENTED diagnostic accelerant on v0 (walk rate only; mechanism is
  v0-independent), measuring `γ_total = γ_soft + γ_rigid` (the dominant rigid
  channel), `s_grip`, re-targets, `r/r0`. n_fil=60, 10 motors, v0_accel=300×,
  6×20 000-step samples, drift ~1e-15 (no blow-up). Result
  (`outputs/h3/production/gripwalk_tier1.json`):

  | quantity | binned_r0 (proxy) | grip_walk |
  |---|---|---|
  | commanded stretch `s_grip` | **0 nm (flat)** | **~250–330 nm sustained** |
  | bead re-targets (transport) | n/a (relabel only) | **12 → 81 (engaged)** |
  | per-head force `k·s_grip` | ~0 | **~0.33 pN (→ F_stall 0.5 pN)** |
  | `γ_total` | 1.8–4.8 ×10⁻⁵ mN/m | 2.4–4.2 ×10⁻⁵ mN/m |
  | `r/r0` | 0.999 | 1.001 |

  **Two-part verdict:**
  1. ✅ **Mechanism validated at the unit level.** grip_walk *sustains* a
     non-zero commanded stretch (~0.33 pN/head, approaching F_stall) and *does
     transport material* (re-targets fire), where the proxy holds `s_grip = 0`
     and never transports. The proxy's diagnosed fatal flaw — force relaxes to
     0, no strain — **is fixed**. Regression-green, opt-in, drift-clean.
  2. ❌ **γ_total does NOT lift off the ~4×10⁻⁵ mN/m floor** — grip_walk ≈
     binned_r0 ≈ the thermal backbone baseline. The myosin contractile signal
     is below the floor at this scale.

## 4b. KEY FINDING — a SECOND floor cause: the myosin force budget (PI)

The Tier-1 γ-non-lift is explained by a force-budget gap the proxy diagnosis
under-weighted (it noted it in passing: root-cause §1 "2 orders short before
any cancellation"). Quantified at the literal parameters (R_cell = 10 µm,
great-circle circumference 62.8 µm):

```
KU-3.5 band 0.5 mN/m  ⇒  great-circle tension needed = 31 416 pN
ideal FULL-PRODUCTION myosin force
  = 100 motors × 2000 heads × 0.5 engaged × F_stall(0.5 pN)  = 500 pN
ideal γ (perfectly aligned, no cancellation)  = 0.008 mN/m   ⇒  63× SHORT
```

So **even an ideal, fully-engaged, perfectly-aligned production myosin
population is ~63× below the KU-3.5 band on force budget alone** — before any
isotropic cancellation. The grip-walk fix was *necessary* (the proxy generated
no sustained force at all) but is *not sufficient* to reach [0.35, 0.65]: the
floor has TWO causes — (1) the lumped proxy (now fixed) and (2) this force-
budget gap.

**Candidate resolutions (PI design decision — NOT taken autonomously):**
- **×40 mesoscale force scaling.** The sanctioned ×40 mesoscopic coarse-graining
  (1000 effective filaments vs ~38 000 native) means each effective minifilament
  stands in for ~40 native ones. If the *per-effective-motor force* (or motor
  count / F_stall) scales ×40 accordingly, 63× → ~1.6× short — within reach.
  This is the most likely fix and is squarely a mesoscale-ratification call.
- Revisit motor density / heads-per-minifilament / F_stall vs the KU anchor.
- Re-examine the KU-3.5 band's mapping to the method-of-planes great-circle
  measurement at the mesoscale.

## 5. Next (PI-gated)

- **PI decides the force-budget resolution (§4b)** — most likely the ×40
  mesoscale force scaling. This is the load-bearing next decision; without it,
  no run length or scale reaches the band.
- PI reviews the §3 HOOMD reconciliation deviation (continuous `r0_eff` via
  per-tick bin requantisation).
- The implementation ships **opt-in (default `binned_r0`)** — do NOT flip the
  default until the force budget is resolved and a production re-measure passes.
- Tier 3 (full [0.35,0.65] re-measure at literal v0, run length sized to ℓ₀/v0)
  is gated on the §4b resolution — a 40 h run before then would land ~63× short
  by construction.
