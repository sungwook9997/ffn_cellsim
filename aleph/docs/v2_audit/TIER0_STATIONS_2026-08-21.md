# TIER 0 — the built cell DOES admit bipolar stations, and the thing stopping them is placement, not polarity

**2026-08-21, 06:15 KST.** Zero of the eight crossbridge PI-GAPs used. No stiffness, no rate, no
force, **no magnitude at any reach.** Record: `world_tier0/stations.json`, self-stamped.

**Why this ran before PI queue item 11 is decided, not after.** If the built cell admitted no legal
bipolar station, **no value of any of the eight constants would make this family contractile** —
`STATE.md` (e) 5 records exactly that outcome for SF geometry (`N_stations: 0`, for a model reason,
with the standing instruction *do not force it*). The question was pre-registered by the lane that
is blocked, in the same message that told me the answer might be zero.

## The curve

Reach is a **declared axis**, swept over multipliers of the two length scales the family module names
— the built NMII head offset (0.200 µm) and the cortex thickness — declared before the run and not
narrowed after. 442 minifilaments, H = 30 (⚠ UNRATIFIED, PI queue item 2).

| reach [µm] | stationed | unstationed | ⚠ no cortex in reach | ⚠ no anti-parallel partner | minus-row first choice |
|---:|---:|---:|---:|---:|---:|
| 0.05 | **0** | 442 | **442** | 0 | 0 |
| 0.10 | **0** | 442 | **442** | 0 | 0 |
| 0.15 | **0** | 442 | **442** | 0 | 0 |
| 0.20 | **0** | 442 | **442** | 0 | 0 |
| 0.30 | 89 | 353 | 351 | 2 | 38 |
| 0.40 | 323 | 119 | 106 | 13 | 137 |
| 0.60 | **442** | 0 | 0 | 0 | 192 |
| 0.80 | **442** | 0 | 0 | 0 | 192 |
| 1.20 | **442** | 0 | 0 | 0 | 192 |

## What the numbers say, arithmetically

**1. The answer is not zero.** The cell admits a full set of legal bipolar stations — all 442 — at a
pairing reach of 0.60 µm and above. The blocking outcome that would have made item 11 moot did not
occur.

**2. The two failure reasons do not compete — one dominates everywhere.** ⚠ They are never summed,
because they have **different owners**: `no_cortex_in_reach` is a PLACEMENT problem owned by
`build/nmii.py`; `no_antiparallel_partner` is a POLARITY problem owned by `build/cortex.py`'s
`polarity_mix` (PI decision 4). At every reach where either is nonzero, placement is larger by one to
two orders — 442 vs 0, then 351 vs 2, then 106 vs 13. **Polarity is not what is stopping the stations.**

**3. The threshold sits on top of the built head offset.** At reach = 0.200 µm — *exactly* the
`head_offset_um` the minifilaments were built with — the count is still **0**. It first becomes
nonzero at 0.30 µm, 1.5×. So the minifilaments are standing further from the nearest actin than their
own heads extend. That is a statement about `build/nmii.py`'s radial placement against
`build/cortex.py`'s shell, and it is the arithmetic, not a diagnosis.

**4. The anti-parallel requirement is not free, even where nothing fails.** At saturation
`n_minus_row_first_choice` is **192 of 442**. The minus row takes its nearest filament in 192 cases;
in the other 250 the anti-parallel constraint moved the choice. So polarity is an active constraint at
large reach even though it causes no failures there — which is a different fact from "polarity is
unconstrained", and the two would be indistinguishable without this column.

## What this run does NOT say

⚠ **No magnitude, at any reach.** This counts what *could ever be paired at build time*. The kinetic
capture radius that decides what binds on a tick is a different quantity and is a PI-GAP — conflating
them would let a build detail set a rate.

⚠ **No station count is quotable without its reach.** The count is monotone in reach, so a single
value spans 0 to 442 by choice of one argument. Picking the reach that yields an acceptable count is
precisely what `CLAUDE.md` forbids. The curve is the result; any row of it alone is not.

⚠ **Three PI-GAPs set this answer**: `n_bb`, `head_offset_um`, `backbone_length_um` fix the straddle
geometry, so they fix what a reach can reach. They are recorded with their `PI_GAP` labels beside the
values, not in a footnote.

⚠ **Interpretation is not mine.** The blocked lane specified the readout and pre-registered the one
conclusion available before the numbers existed (`n_stationed == 0` ⟹ item 11 is moot). It did not
occur; which failure reason dominates decides the next question, and that reading belongs to them.
