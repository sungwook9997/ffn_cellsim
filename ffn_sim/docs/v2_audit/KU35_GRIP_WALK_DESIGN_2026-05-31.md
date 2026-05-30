# KU-3.5 myosin grip-point walking — IMPLEMENTATION DESIGN (PI-ratification brief)

> **Status: PROPOSAL. No code written. PI sign-off required before implementation.**
> Companion to the root-cause brief `KU35_FLOOR_ROOT_CAUSE_2026-05-31.md` (the *what/why*);
> this is the *how*. It specifies the exact code-change plan to convert the lumped binned-r0
> ratchet into a true AFINES `pos_a_end` grip-point walking motor on HOOMD, flags the design
> decisions that need PI ratification (filament polarity, bipolar binding rule, force-ceiling
> param), and gives a validation plan that does **not** require a 40 h production run.
>
> Scope rule (CLAUDE.md): every change below is **additive + opt-in** (a new `stepping_mode`
> with default `"binned_r0"` = current behaviour), so the proxy stays the A in an A/B and the
> existing regression suite is untouched until PI promotes the new mode.

---

## 0. TL;DR of the change

| | Current (proxy) | Proposed (grip-walk) |
|---|---|---|
| What a "step" does | decrement `bond_bins[h]` → relabel head↔**same** bead to lower-r0 bin | advance per-head **bound-bead index** toward the filament minus end → re-target head↔**downstream** bead, r0≈0 |
| Where force comes from | a shrinking *idealized* rest length on a bead that relaxes away | the **stretch** `k·(r−0)` of the head spring against the SHAKE-rigid actin shell |
| Force ceiling | bin ceiling 330 nm → `k·330nm = 0.33 pN < F_stall 0.5 pN` | grows ~`k·ℓ₀` (0.5 pN) per net bead walked, multi-bead → ≥ F_stall |
| New per-head state | `bound_to_actin` (bead tag), `step_accum` (frac) | **+ `bound_bead_pos`** (position-in-filament), **+ `bound_filament`** (filament id), **+ minus-end sense per filament** |
| Net contraction | none required at binding (no sidedness) | bipolar binding rule: +/− head sets grip antiparallel/different filaments |

The single load-bearing line-edit is in `MyosinStepUpdater.act()` **Step 3** (`cortex/myosin.py:985-1022`):
replace *"shrink the bin"* with *"walk the bound-bead index and re-target the bond"*. Everything
else (binding bookkeeping, the per-head accumulator, the bond-group rewrite via `set_snapshot`,
the tension measurement) is reused as-is or lightly extended.

---

## 1. Grip-point walking on HOOMD

### 1.1 The mechanism HOOMD already gives us for free

`MyosinStepUpdater.act()` already, every batch tick: reads a snapshot, mutates
`bonds.group` / `bonds.typeid` (to add/remove/relabel attach bonds), and writes it back with
`sim.state.set_snapshot()` (`cortex/myosin.py:1024-1067`). **Re-targeting a head-actin bond to a
different actin bead is the *same* operation as binding a new one** — change the second entry of
the bond pair (`attach_bonds[:,1]`) from bead tag `i` to the downstream tag `i−1`. No new HOOMD
primitive is needed; the walking kernel is a different *content* for the bond-group rewrite that
already runs.

The head **particle** is never moved by the updater (it is integrated by BAOAB). Walking moves
only the *grip point* — which actin bead the head's spring pulls on. Keeping the head in place
while moving the grip downstream is exactly what stretches the spring against the rigid shell.

### 1.2 New per-head state (extends the two arrays at `myosin.py:729-738`)

Add two int64 arrays alongside the existing `_head_bound_to_actin` and `_head_step_accum`:

```python
# myosin.py __init__, after line 738:
self._head_bound_bead_pos = np.full(n_heads_total, -1, dtype=np.int64)  # j within filament
self._head_bound_filament  = np.full(n_heads_total, -1, dtype=np.int64)  # filament id
```

`_head_bound_to_actin[h]` already stores the **global bead tag** the head grips; the two new
arrays store its decomposition `(filament_id, position_within_filament)`, which is what the walk
increments. They are set at bind time (Step 2) and on every successful walk (Step 3), reset to
−1 on unbind (the existing reset block at `myosin.py:849-853`, extended to clear all three).

The bead-tag ↔ `(filament, pos)` map is **exact integer arithmetic for fixed-N cortex**
(`bead_tag = filament*N + pos`, verified: tag 41, N=7 → fil 5, pos 6) and a `searchsorted` on
`filament_starts[]` for variable-N cortex (see §2.4). The cortex bead block is contiguous
`[0, n_cortex_actin)` and tag-stable (FA/myosin appended after), so no re-index is needed.

### 1.3 The walking kernel — replaces `myosin.py:985-1022`

Current Step 3 accumulates `d_bin = v·batch_dt/bin_width` and decrements `bond_bins`. The new
Step 3 accumulates `d_bead = v·batch_dt/ℓ₀` (bead spacings, not bin widths) and **walks the
bound-bead index toward the minus end**, re-targeting the bond:

```python
# ---- Step 3 (grip-walk mode): advance grip point along filament ----
if attach_bonds.shape[0] > 0:
    head_tags  = attach_bonds[:, 0]
    actin_tags = attach_bonds[:, 1]
    r = np.linalg.norm(pos[head_tags] - pos[actin_tags], axis=1)
    F_mag = self.p.k_head_actin * r        # r0 = 0 → stretch carries force
    v_step = hill_velocity_clamped(F_mag, v0=..., F_stall=..., a_over_F_stall=...)
    d_bead = (v_step * self.p.batch_dt) / self.p.ell0_cortex   # float, per-bond

    head_locals = <vectorised inverse map, as today at :1006-1012>
    self._head_step_accum[head_locals] += d_bead
    adv = np.floor(self._head_step_accum[head_locals]).astype(np.int64)   # ≥0 integer beads
    self._head_step_accum[head_locals] -= adv.astype(np.float64)

    # Walk each head's bound-bead pos toward the minus end by `adv` beads,
    # clamped at the minus-end bead (pos == minus_end_pos[filament]); then
    # recompute the new global bead tag and rewrite attach_bonds[:,1].
    for row in np.flatnonzero(adv > 0):
        h    = head_locals[row]
        fil  = self._head_bound_filament[h]
        pos_j = self._head_bound_bead_pos[h]
        new_j = self._walk_toward_minus(fil, pos_j, int(adv[row]))   # §2.3
        self._head_bound_bead_pos[h] = new_j
        new_tag = self._bead_tag(fil, new_j)                          # §2.4
        self._head_bound_to_actin[h] = new_tag
        attach_bonds[row, 1] = new_tag
        self._n_step_advances_total += 1
    # bond_bins is no longer the force-carrying state; in grip-walk mode all
    # attach bonds use a SINGLE r0≈0 bin (see §4) so bond_bins stays 0.
```

Key points:

- **r0 ≈ 0 carries the force** (the `bridge/motor.py:209` `motor_head_actin r0=0` precedent —
  which is registered-but-unused today since no runtime stepper exists on the H.4 side). The bin
  scheme collapses to a single r0=0 attach type in grip-walk mode (§4), so `bin_r0[...]` is no
  longer the force-setting knob.
- **The accumulator already exists** (`myosin.py:732-738, 1014-1016`) and already solves the
  dt-coupling: `d_bead` is sub-1 per tick (≈ 2.7×10⁻⁵ bead/tick at the FA-limited dt, §1.4), so
  integer-bead walks fire rarely; the fractional remainder is retained. This is unchanged — only
  the denominator flips from `bin_width` (33 nm) to `ℓ₀` (500 nm), and the action flips from
  "decrement bin" to "walk bead + re-target bond".
- **Minus-end clamp** replaces the bin-0 clamp: a head that reaches its filament's minus-end bead
  stops walking (AFINES `at_barbed_end` latch, `AFINES_ALGORITHM_NOTES.md:481`) but stays bound
  and keeps pulling at whatever stretch it has. Detachment is via the existing Bell-Evans slip
  (Step 1) — biologically the end-dwelling head.

### 1.4 dt coupling and the *sub-bead force* subtlety (important — reshapes the kernel)

A diagnostic computed at the real parameters (v0 = 0.2 µm/s NMIIA, FA-limited
dt ≈ 6.7×10⁻⁷ s, batch 100):

```
v0·batch_dt (unloaded) = 13.4 pm/tick = 2.68e-5 bead/tick
sim-steps to walk ONE integer bead = 3.7e6   (=ℓ₀/v0 = 2.5 s sim time, dt-independent)
```

**Consequence:** a pure *integer-bead-advance* kernel produces force in a staircase that is FLAT
for ~3.7×10⁶ steps between treads. The default production schedule (80×5000 = 4×10⁵ steps) never
reaches even one tread — the accumulator only gets to 0.11 bead. So a naive integer-only kernel
would *still* read γ≈0 on any run shorter than ~4×10⁶ steps, reproducing the floor for a
different reason.

**Design resolution — continuous sub-bead grip (recommended).** Track a *continuous* grip
fraction and re-target to the **bracketing pair** of beads, OR (simpler, preferred) advance an
*idealized minus-end-ward offset distance* `s_grip[h]` (metres, the AFINES `pos_a_end`) and
realise it as: the head-actin bond targets the *current* bound bead, but with a **non-zero r0
equal to `−s_grip`** projected along the filament tangent — i.e. the rest length is *driven
negative-ward* so the spring is pre-stretched by exactly the walked distance, continuously, even
before a full bead is crossed. When `s_grip` exceeds one bead spacing, re-target the bound bead
downstream and subtract `ℓ₀` from `s_grip` (the AFINES `pos_a_end` overflow → next spring,
`AFINES_ALGORITHM_NOTES.md:478-479, 538-544`).

Concretely the force-carrying quantity becomes `F = k·(r − r0_eff)` with
`r0_eff = max(r − s_grip, 0)` so that the *commanded* stretch is `s_grip` regardless of where the
bead currently sits. `s_grip += v(F)·batch_dt` each tick (continuous), and the integer-bead
re-target keeps `s_grip` bounded to `[0, ℓ₀)` so the bond never points at an absurdly distant
bead. **This is the faithful `pos_a_end`** (a continuous walked-distance) and it generates force
from step 1, not only after 3.7×10⁶ steps.

This means the kernel keeps a per-head **`_head_grip_s` (float metres)** in place of (or in
addition to) the integer `step_accum`; the integer bead-pos walk is the `pos_a_end`-overflow
bookkeeping, and `s_grip ∈ [0, ℓ₀)` is the intra-segment offset. The Hill `v(F)` reads the *true*
current stretch (so the force-velocity feedback is honest), exactly AFINES §4.3.

> **PI decision 1.4** — *integer-bead-only* vs *continuous sub-bead `pos_a_end`*. Recommendation:
> continuous sub-bead, because (a) it is the literal AFINES `pos_a_end` and (b) integer-only
> cannot show force on any feasible diagnostic-length run. Integer-only is simpler but defers the
> floor rather than fixing it at diagnostic scale.

### 1.5 Functions / lines to change (summary table)

| File:line | Change |
|---|---|
| `cortex/myosin.py:729-738` (`__init__`) | add `_head_bound_bead_pos`, `_head_bound_filament`, `_head_grip_s` arrays; add `ell0_cortex`, `minus_end_sense`, `stepping_mode` to ctor args |
| `cortex/myosin.py:646-660` (`register_..._bond_params`) | grip-walk mode: register a single attach bond type at r0=0 (k=k_head_actin) instead of the n_bins r0-binned set (§4) |
| `cortex/myosin.py:859-983` (Step 2 binding) | on bind, also set `_head_bound_bead_pos`/`_head_bound_filament`; apply the bipolar sidedness filter (§3); set `_head_grip_s=0` |
| `cortex/myosin.py:849-853` (unbind reset) | also clear the three new per-head arrays |
| `cortex/myosin.py:985-1022` (Step 3 stepping) | **replace** bin-decrement with grip-walk kernel (§1.3 + §1.4) under `if stepping_mode=="grip_walk"`; keep the old branch under `else` |
| `cortex/myosin.py:1024-1055` (bond rewrite) | unchanged — it already writes `attach_bonds` back; grip-walk just changed `attach_bonds[:,1]` |
| `cortex/myosin.py:194-280` (`resolve_cortex_myosin`) | add `stepping_mode` (default `"binned_r0"`), validate |
| `cortex/myosin.py:1092-1109` (`make_cortex_myosin_updater`) | thread `cortex_minus_end_sense`, `ell0_cortex`, `filament_starts` through |
| `cell/cell.py:1078-1082` | pass the new args (computed from `topology` — already in scope) to `make_cortex_myosin_updater` |

No change to `bridge/motor.py` (its `r0=0` registration is the precedent we adopt; the H.4 FA
motor has no runtime stepper — confirmed only `integrin_bonds.py` has a runtime `act()`).

---

## 2. Cortex filament polarity — **PI DESIGN DECISION (the key open question)**

### 2.1 Finding: cortex actin has NO biological polarity today

`cortex/cortex.py:656-686` (`generate_cortex_topology`) lays each filament's beads along a
tangent: `offsets = (arange(N) − (N−1)/2)·ℓ₀`, `positions[f,j] = com + offsets[j]·tangent[f]`,
with bonds `(f·N+j, f·N+j+1)`. **The bead index order `j=0…N−1` is a pure geometric artifact of
the construction loop direction** — it is NOT a plus/minus end. The tangent azimuth `φ` is drawn
`U(0,2π)` (isotropic), so even the *sign* of the tangent is random per filament. There is no
`polarity`, `plus_end`, `minus_end`, or barbed/pointed attribute anywhere in `CortexTopology`,
the GSD frame, or `ResolvedH3`. (`grep` confirms: zero polarity fields in cortex.)

This is mechanistically a real gap, not just a coding one: myosin walks toward the **pointed
(minus) end** of actin, so a contractile network needs a defined per-filament polarity. The
proxy never needed it (it didn't transport material); the grip-walk motor does.

### 2.2 Options

- **Option A — ASSIGN polarity = bead-index order (minus end at `j=0`).** Declare, by convention,
  that bead `j=0` is the minus (pointed) end and `j=N−1` is the plus (barbed) end of every
  filament. "Walk toward minus" = decrement `j` toward 0. Zero topology change; polarity is a
  pure label carried in the updater (`minus_end_sense = "low_index"`). Because the tangent
  azimuth is already isotropic, this assignment is statistically isotropic in lab frame (no
  artificial global polarity bias) — it just fixes a *per-filament* reference so heads have a
  consistent direction to walk.
- **Option B — ASSIGN random polarity sign per filament.** Draw `minus_at_low_index[f] ∈ {True,
  False}` at construction (seeded), store in `CortexTopology` (+ a `polarity` particle/property or
  a side-array threaded to the updater). Walk decrements or increments `j` per that flag.
  Decorrelates polarity from the (already random) construction order — marginally more "physical"
  but adds a stored array and a plumbing path.
- **Option C — derive polarity from a barbed-end field (future).** When the lamellipodium /
  turnover machinery defines real barbed ends (treadmilling), inherit polarity from it. Out of
  scope for KU-3.5; flagged for Phase 2.

### 2.3 Recommendation + the `_walk_toward_minus` helper

**Recommend Option A** for the KU-3.5 fix: minimal, additive, statistically isotropic, and
sufficient for contractility (what matters is that the *two* head sets of a minifilament walk
*different/antiparallel* filaments — §3 — not the lab-frame polarity distribution). Option B is a
trivial later refinement if PI wants polarity decorrelated from build order.

```python
def _walk_toward_minus(self, fil, pos_j, n):
    # Option A: minus end at j=0. Walk decrements, clamp at 0.
    lo = self._fil_min_pos[fil]   # 0 for fixed-N; filament-local 0 for var-N
    return max(lo, pos_j - n)
```

For Option B, `_walk_toward_minus` branches on `self._minus_at_low_index[fil]` (decrement vs
increment, clamp at the corresponding end `_fil_min_pos`/`_fil_max_pos`).

> **PI decision 2** — adopt **Option A** (minus end = bead 0, label-only, no topology change) for
> KU-3.5; defer Option B/C. *This is the load-bearing polarity decision the diagnosis flagged.*

### 2.4 Bead-tag ↔ `(filament, pos)` map (fixed-N and variable-N)

- **Fixed-N (current production path, `n_filaments × beads_per_filament`):**
  `filament = tag // N`, `pos = tag % N`, `tag = fil·N + pos`. Exact, O(1).
- **Variable-N (`VariableLengthCortexLayout`, `cortex.py:1024-1173`):** carry
  `filament_starts[]` and `n_beads_per_filament[]` into the updater; `filament =
  searchsorted(filament_starts, tag, "right")−1`, `pos = tag − filament_starts[filament]`,
  `tag = filament_starts[filament] + pos`. The Cell builder already has both arrays in scope
  when variable-N is used.

The updater takes `ell0_cortex`, and either `(N,)` for fixed-N or `(filament_starts,
n_beads_per_filament)` for variable-N, via `make_cortex_myosin_updater` (threaded from
`cell/cell.py` where `topology` / `layout` is in scope at `:1078-1082`).

---

## 3. Bipolar sidedness — net contractile dipole (Stam–Hocky)

### 3.1 The gap

Binding (`myosin.py:859-983`) attaches **each head independently** to its nearest eligible actin
bead, with no constraint linking the + and − head sets. So even with sustained per-bond force,
the minifilament has no *organised* dipole: its two sides can grip the same filament, or
parallel filaments, giving zero net contraction (diagnosis §1C). Stam–Hocky contractility
requires the two opposite-polarity head sets to pull **antiparallel** material inward across the
rigid rod.

### 3.2 Proposed binding rule (minimal, additive)

A minifilament has `2H` heads: locals `[0,H)` are **+ side**, `[H,2H)` are **− side** (the layout
convention, `myosin.py:310-313`). Enforce:

- **+ side heads bind filaments whose minus end points toward the rod's `+axis` direction;
  − side heads bind filaments whose minus end points toward `−axis`.** Operationally, at a
  candidate bind, compute the actin filament's **minus-end-ward tangent** `m̂` (the unit vector
  from the bound bead toward bead `j−1`, i.e. toward the minus end under Option A) and require
  `sign(m̂ · û) > 0` for + side heads, `< 0` for − side heads, where `û` is the rod axis
  (`layout.axes[motor]`, already stored). This makes the two sides walk filaments of *opposing*
  minus-end orientation → both sides drag actin toward the rod center → net contraction.
- **Equivalent, simpler variant:** require the + and − sides to bind **different filaments**
  (cheap anti-degeneracy: reject a + bind to a filament already gripped by this motor's − side,
  and vice versa) AND apply the `sign(m̂·û)` orientation gate. The "different filament" clause
  alone removes the worst degeneracy (both sides on one filament = zero dipole); the orientation
  gate adds the antiparallel requirement.

Implementation: in Step 2 (`myosin.py:908-952`), after a candidate `best_bead` is found, compute
`m̂` from `best_bead` and its minus neighbor (`_walk_toward_minus(fil, pos, 1)` gives the
neighbor; both already mapped), look up the head's side (`head_within < H`), and `continue`
(reject this candidate) if the orientation/different-filament test fails. This is a few lines
inside the existing candidate loop; it does not change the binding *probability* model (still
Bell-Evans `k_on`), only *which* acceptors are eligible per head side.

> **PI decision 3** — adopt the **`sign(m̂·û)` antiparallel gate + different-filament clause** as
> the bipolar binding rule. Flag: this slightly lowers realized bound fraction (some nearest
> beads get rejected → head stays free that tick, re-tries next tick), which is biophysically
> correct (myosin engages productively-oriented actin) but should be watched in the engagement
> diagnostic (§5).

### 3.3 Why this is enough (scale check)

With the dipole organised, the ~1660 engaged heads (2000 heads × 0.83 zero-load bound fraction,
§5) each delivering up to `k·ℓ₀ ≈ 0.5 pN` of antiparallel pull, the per-motor contractile force
is O(H · F_stall) and the network sum is the contractile stress the method-of-planes integrates.
The old proxy delivered ~0 because (a) no sustained stretch and (b) no dipole; this fixes both.

---

## 4. Force ceiling to reach F_stall — param decision

### 4.1 Diagnosis recap, quantified

- Single head, r0=0, stretched by one bead spacing ℓ₀=500 nm: `F = k_head_actin·ℓ₀ = 1e-6 · 5e-7
  = 0.50 pN = F_stall` exactly (verified). So **one net bead of grip-walk already reaches
  F_stall** at the head-actin spring — the proxy's 0.33 pN ceiling (bin 330 nm) was an artifact
  of the *bin* ceiling, not the physics.
- **But** the head is also tethered to the rigid rod by `k_head_spring = 1e-6`. The two springs in
  series (head↔actin and head↔backbone) give `k_series = 2.5e-7`, so the *relaxed* tension
  delivered at the actin attach point for a 1-bead offset is `k_series·ℓ₀ = 0.25 pN` (need ~2
  beads of net walk for a sustained 0.5 pN through the series path). Transient (right after a
  re-target, before relaxation) is the full 0.5 pN; time-average sits between.

### 4.2 How grip-walk fixes the ceiling

Grip-walk **removes the 330 nm cap entirely**: the commanded stretch `s_grip` (§1.4) grows with
walking and is bounded only by Hill stall (the motor stops adding stretch when `F → F_stall`) and
by Bell-Evans detachment. So the force self-limits at F_stall by the force-velocity law, which is
the correct physics — no artificial ceiling. Multi-bead walking (sustained engagement across many
ticks) builds the offset the single bin never could.

### 4.3 Param flags for PI

- **`k_head_actin` (currently 1e-6 N/m):** with the series halving, sustained F_stall via the
  head wants either a 2-bead net offset (fine — Hill allows it) or `k_head_actin` raised so a
  smaller offset reaches stall. **Recommendation: keep k_head_actin = 1e-6 (brief-literal,
  KU-anchored, do not tune to pass a gate — CLAUDE.md no-magic-number rule).** The series effect
  is real physics (the cross-bridge is compliant), and Hill stepping naturally builds the offset.
  Only revisit if §5 shows the network can't reach the band with the literal value.
- **Single r0=0 attach bond type:** in grip-walk mode, register ONE attach type at r0=0 (replace
  the `n_bins` r0-binned set, `myosin.py:646-660`). This is the `bridge/motor.py:209` precedent.
  The bin scheme was only there to keep the construction force-free for the *relabel* mechanism;
  grip-walk doesn't relabel r0, so the binning is dead weight. (Binding still records `s_grip=0`
  so the first-tick force is 0 → no construction shock.)

> **PI decision 4** — **no param change** (keep k_head_actin, k_head_spring at brief-literal
> values); replace the r0-binned attach set with a single r0=0 attach type in grip-walk mode.
> Re-open only if the §5 short diagnostic shows the literal stiffness cannot lift γ into
> [0.35,0.65].

---

## 5. Validation plan (no 40 h production run)

Three tiers, cheapest first. Each is a *contract written before the run* (CLAUDE.md
gate-protocol).

### 5.1 Tier 1 — sustained-tension micro-diagnostic (minutes, the decisive test)

**Goal:** prove the head-actin stretch + γ_soft now *persist/grow* instead of relaxing to 0.

The diagnosis showed the proxy holds `r/r0 = 1.000` and γ flat at the floor "even after 1396
advances." The new mechanism must show the *opposite*: a monotone rise in mean head-actin bond
stretch and in γ_soft over a short window.

**Setup (a `/tmp` or `scripts/` micro-driver, NOT the full production driver):**
- Small cortex: `n_fil = 30–50`, fixed-N=7, a handful of motors (5–20), constrained backbone
  (rigid actin via SHAKE — the production physics), grip-walk mode ON.
- **Decouple the dt-walk-rate confound (§1.4):** run a *diagnostic v0* (e.g. 5–20 µm/s) OR a
  larger `batch_dt`, so integer-bead walks fire within ~10⁴–10⁵ steps. This is a *diagnostic
  accelerant*, not a physics change — it only compresses the 2.5 s sim-time/bead into observable
  wall time; the *mechanism* (does stretch persist?) is v0-independent. Document it as such.
  (Alternatively, seed each engaged head with `s_grip = 1–2·ℓ₀` and run ~10⁵ steps to show the
  pre-stretched force does NOT relax — a direct refutation of proxy point A.)
- **Measure each k samples (≈10⁵ steps total):** mean `r_head_actin − r0` (should rise, not
  relax to 0), `γ_soft` (should rise off the 1.4e-5 floor), `n_engaged`, `n_step_advances`.

**PASS contract (Tier 1):**
1. Mean head-actin stretch strictly increases over the first ~5 samples (proxy: flat at ~0).
2. γ_soft rises ≥ 10× off its proxy floor within the window (directional proof; not the full band
   yet).
3. `r/r0 < 1.0` (the cortex mean radius contracts — the contraction proxy the driver already
   logs at `h3_ku35_tension.py:419-432`).
4. No constraint-drift blowup, no LJ-CFL `WARN` (the driver's existing guards, `:435`).

This tier alone *settles whether the redesign works*. If stretch still relaxes, the design is
wrong and we stop before any long run.

### 5.2 Tier 2 — unit tests + regression (CI-fast)

Extend `tests/test_myosin.py` (additive; existing tests stay green because grip-walk is opt-in):

- `test_grip_walk_retargets_downstream`: build a tiny cortex+motor, force one head bound at
  `(fil, pos=3)`, run one tick with a `v0` big enough to fire one integer bead; assert the
  head's `bonds.group[:,1]` moved from bead `fil·N+3` to `fil·N+2` (minus-ward), head particle
  position unchanged, and `_head_bound_bead_pos` decremented.
- `test_grip_walk_minus_end_clamp`: head bound at `pos=0`; assert a walk does NOT decrement past 0
  and the head stays bound (AFINES latch).
- `test_grip_walk_stretch_generates_force`: after a re-target, assert
  `|r_head − r_newbead| ≈ ℓ₀` and `F = k·r ≈ 0.5 pN` (the force-ceiling claim, §4.1).
- `test_bipolar_sidedness_gate`: construct a motor over two antiparallel filaments; assert + side
  heads only bind the `m̂·û>0` filament and − side only `m̂·û<0` (§3).
- `test_polarity_assignment_isotropic` (Option A): assert minus-end labels exist and the lab-frame
  tangent distribution is unchanged/isotropic (no global polarity bias introduced).
- `test_binned_r0_mode_unchanged`: run the *old* mode and assert byte-identical behaviour to a
  saved fixture (the A/B safety net — proves additivity).
- Hill sign/sense tests (`test_myosin.py:192-218`) carry over unchanged.

Plus the existing `tests/test_ku35_v4_driver.py` must stay green in default (binned_r0) mode.

### 5.3 Tier 3 — full KU-3.5 [0.35, 0.65] re-measure (the eventual gate)

Only after Tier 1+2 pass and PI ratifies the mode-flip default. Re-run the (now crash-free, per
root-cause §4) production driver `scripts/h3_ku35_tension.py` with `stepping_mode="grip_walk"` at
production scale (n_fil≈150, full schedule). **Budget the run length to the §1.4 walking
timescale:** at the literal v0=0.2 µm/s, ~10⁶–10⁷ steps are needed for the network to walk
multiple beads and reach steady contractile stress — the default 4×10⁵-step schedule is too short
*by construction* (it would even fail the *new* mechanism). Either extend `--n-sample`/`--interval`
to ≳ 5×10⁶ total steps, or (if PI accepts) measure at the steady-state plateau with a longer
`interval`. Gate: `γ_total = γ_soft + γ_rigid ∈ [0.35, 0.65] mN/m` (the driver already computes
both, `:416-418`). Compare myosin-ON vs myosin-OFF (blebbistatin, KU-3.18) as the A/B that γ rise
is myosin-driven.

> **Validation note (run-length is a first-class finding):** the §1.4 diagnostic shows the
> *previous* production runs were ALSO too short to show contraction even if the mechanism had
> been right (0.11 bead in 4×10⁵ steps). Tier 3 must size the run to ≥ a few ℓ₀/v0 = ≥ several×
> 2.5 s sim time. This is independent of the proxy-vs-grip-walk question and should be surfaced to
> PI regardless.

---

## 6. Risk / scope

### 6.1 Additive vs invasive

- **Additive (low risk):** new per-head arrays; new `stepping_mode` config (default
  `"binned_r0"`); the grip-walk kernel under an `if` branch; new unit tests; the bipolar gate and
  polarity label. With the default, the builder/updater are bit-for-bit the current behaviour, so
  the whole regression suite and every committed result are untouched. **This should ship as
  opt-in first**, exactly like the membrane/turnover/FA additive modules already in the tree.
- **Invasive-ish (contained):** the single r0=0 attach-type registration differs from the binned
  set *in grip-walk mode only*; and Tier 3 changes the production run length. Both are gated
  behind the new mode + PI sign-off.

### 6.2 What could go wrong

1. **Walking timescale (§1.4) — the biggest trap.** At literal v0, force builds over millions of
   steps; an integer-only kernel shows nothing on short runs. *Mitigation:* the continuous
   sub-bead `s_grip` design (§1.4) generates force from step 1; the diagnostic uses an explicit
   accelerant labelled as such. **PI decision 1.4 must be made.**
2. **Series-stiffness force loss (§4.1).** The head↔backbone tether halves the delivered stretch;
   sustained F_stall needs ~2-bead net offset. *Mitigation:* Hill naturally builds it; keep
   literal k (no tuning); re-open only if Tier 3 misses the band.
3. **Re-target + nlist exclusion cap.** Re-targeting moves a head's attach bond from bead `i` to
   `i−1`, changing the *degree* of two beads. The existing degree guard (`_MAX_CORTEX_BEAD_DEGREE
   = 6`, `myosin.py:806, 942`) and per-bead head cap (`MAX_HEADS_PER_BEAD = 3`, `:901`) must be
   re-checked at *walk* time, not just bind time, or a downstream bead can exceed HOOMD's 7-bond
   exclusion cap → the crash root-cause §0 already fixed for binding. *Mitigation:* apply the same
   degree check before accepting a downstream re-target; if the target bead is full, the head
   **stays on the current bead this tick** (walk deferred, `s_grip` retained) — a clean stall, not
   a crash. This is a few lines and MUST be in the kernel.
4. **Bipolar gate lowers engagement.** Rejecting wrongly-oriented acceptors lowers bound fraction;
   if too aggressive, contraction weakens. *Mitigation:* the "different-filament" clause alone
   removes the zero-dipole degeneracy; keep the orientation gate as a soft preference (try
   oriented acceptor first, fall back to any) if §5 shows starvation. Watch `n_engaged`.
5. **Variable-N path.** Production may use variable-length filaments; the `searchsorted` map and
   per-filament min/max-pos clamps must be correct there. *Mitigation:* Tier 2 includes a
   variable-N walk test; fixed-N is the default for the diagnostic.
6. **r0=0 construction shock.** A bond at r0=0 between a head and a bead ~ℓ₀ apart is instantly
   stretched. *Mitigation:* bind sets `s_grip=0` so `r0_eff = r` at bind → zero initial force
   (the §1.4 mechanism makes the *commanded* stretch, not the geometric distance, the force);
   matches the force-free-at-construction invariant.

### 6.3 Recommended sequencing

1. PI ratifies decisions **1.4 / 2 / 3 / 4** (this brief).
2. Implement grip-walk as **opt-in** (`stepping_mode="grip_walk"`, default off). Land Tier-2 unit
   tests green; binned_r0 fixture proves additivity.
3. Run **Tier 1** micro-diagnostic. If stretch persists/grows → mechanism validated.
4. PI reviews Tier-1 evidence, ratifies flipping the default (or keeps opt-in for one more cycle).
5. **Tier 3** full re-measure at correctly-sized run length → KU-3.5 [0.35,0.65] decision.

No step commits to a 40 h run until Tier 1 (minutes) has proven the core physics.

---

## 7. Flagged PI decisions (consolidated)

| # | Decision | Recommendation |
|---|---|---|
| **1.4** | integer-bead-only vs continuous sub-bead `pos_a_end` walk | **continuous sub-bead** (literal AFINES; generates force at diagnostic scale) |
| **2** | cortex filament polarity: assign how? | **Option A**: minus end = bead 0, label-only, no topology change |
| **3** | bipolar binding rule | **antiparallel `sign(m̂·û)` gate + different-filament clause** |
| **4** | force-ceiling param | **no param change** (keep literal k); replace r0-bin set with single r0=0 attach type in grip-walk mode |
| **(impl)** | ship opt-in first (default = current binned_r0) for safe A/B | **yes** |
| **(val)** | production run length must be sized to ℓ₀/v0 (≫ current 4×10⁵ steps) | surface regardless of mode |

## 8. Key files (for the implementer)

- `cortex/myosin.py` — `__init__:694-745`, binding `:859-983`, **stepping `:985-1022`** (the
  replace), bond rewrite `:1024-1067`, bin registration `:646-660`, resolver `:197-280`, factory
  `:1092-1109`.
- `cortex/cortex.py` — topology `:590-686` (polarity gap; `bond_groups`/`tangents`/`positions`
  are the inputs), variable-N `:1058-1173` (`filament_starts`).
- `cell/cell.py:1078-1082` — `make_cortex_myosin_updater` call site (thread new args here).
- `bridge/motor.py:209` — `motor_head_actin r0=0` precedent (adopt). No runtime stepper there.
- `scripts/h3_ku35_tension.py` — production driver + γ_soft/γ_rigid method-of-planes (`:51-115`,
  `:192-259`); the Tier-3 re-measure host.
- `docs/AFINES_ALGORITHM_NOTES.md:466-565` — the `pos_a_end` mechanism this implements.
- `docs/briefs/H3_cortex.md:56-65` — the Stam–Hocky/Hill spec the proxy departed from.
