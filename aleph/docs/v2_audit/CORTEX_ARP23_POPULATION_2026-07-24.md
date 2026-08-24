# Mixed formin + Arp2/3 cortex population — the missing native architecture

**Date:** 2026-07-24
**Status:** FLAGGED / DEFAULT-OFF prototype. Default (`cortex_arp23_fraction = 0.0`) keeps the current
pure-formin cortex **bit-identical**. This is a **"PI A/B/C" modeling option** — sourced, but the
representative Arp2/3 length + branch treatment are choices to ratify, not silent defaults.
**Branch:** `codex/ff-ac-codex`. **Not committed.**

---

## 1. The gap

The model cortex is **formin-only**: `aleph/laws/architecture_spec.py::CORTEX` is
`nucleator="formin", length_um=3.0` for all 70,686 filaments. The native cortex is **not** formin-only —
it is a **mixed** network:

> Cortical F-actin is ~2/3 formin-nucleated (longer, unbranched) + ~1/3 Arp2/3-nucleated
> (short, branched), **by mass** (Bovellan et al. 2014, *Curr Biol* 24:1628; KB-3.18).

KB-3.18 already flags cortical filament length as a **"PI A/B/C decision"**: native cortical filaments are
SHORT (~hundreds of nm) and high-turnover, ~1/3 Arp2/3. The current 70,686 × 3 µm are *coarse representative*
formin filaments (see `CORTEX_COUNT_CROSSCHECK_2026-07-21.md`), with **no Arp2/3 branched sub-population at
all**. This memo adds that sub-population as an opt-in.

---

## 2. Mass fraction ≠ count fraction (the honest arithmetic)

The physiological datum is a **mass** fraction (~1/3 of cortical actin is Arp2/3-nucleated). But we build
discrete filaments, and **mass ∝ n · length**. Because Arp2/3 filaments are ~10–20× **shorter** than the
representative formin filaments, the Arp2/3 **count** fraction is *much larger* than 1/3.

Fix the cortex filament **budget** (`N_total`, so total areal density ~100/µm² is preserved) and solve:

```
f = N_a·L_a / (N_a·L_a + N_f·L_f)    (Arp2/3 mass fraction)
N_a + N_f = N_total                  (fixed budget → areal density preserved)
```

⇒ `M = N_total / (f/L_a + (1−f)/L_f)`, then `N_a = round(f·M/L_a)`, `N_f = N_total − N_a`.

(`aleph/components/weave/regions.py::derive_cortex_arp23_split`.)

### Worked cases

| f (mass) | L_f | L_a | N_formin | N_arp23 | count ratio a/f | Arp2/3 count frac |
|---|---|---|---|---|---|---|
| 1/3 | 1.0 µm | 0.1 µm | 11,781 | 58,905 | **5.00×** | 0.833 |
| 0.33 | 3.0 µm | 0.15 µm | 6,514 | 64,172 | 9.85× | 0.908 |
| 0.33 | 1.0 µm | 0.15 µm | 16,502 | 54,184 | 3.28× | 0.766 |

Row 1 is the task's illustrative check (**~5× more Arp2/3 by count** for 1/3 by mass at a 10× length ratio).
Row 2 is the **CellConfig default** (`cortex_length_um=3.0`, the coarse representative formin length). Row 3
is what the **fine-cortex** config (formin `length_um=1.0`) would give — the count ratio is sensitive to the
formin length, which is itself a coarse representative, so it is reported per-config, never hard-coded.

**Key point:** ~1/3 by mass ⇒ Arp2/3 is the numerical **majority** of cortical filaments. That is physical
(each Arp2/3 daughter is a short branch), and it is why a formin-only census under-represents the branched mesh.

---

## 3. The representative Arp2/3 length + branch treatment — "PI A/B/C" choices

These are **modeling choices**, flagged (not tuned to an outcome, not sourced point values):

- **Arp2/3 length `L_a` (CellConfig `cortex_arp23_length_um = 0.15 µm`).** Native cortical Arp2/3 filaments are
  short (~0.1–0.2 µm; hundreds of nm, KB-3.18 / Bovellan 2014). 0.15 µm is the **midpoint** of that range.
  *Option A* 0.10 µm (short end, matches the illustrative 5× case with L_f=1.0), *B* 0.15 µm (midpoint,
  default), *C* 0.20 µm (long end). The split scales with `L_a`; the memo table shows the sensitivity.
- **Arp2/3 discretization `seg_um` (`cortex_arp23_seg_um = 0.05 µm`).** A filament must have **≥3 nodes** to
  carry an interior branch bead (mother→branch→daughter). At `L_a=0.15`, `seg=0.05` → 4 nodes / 3 segments —
  branchable, ~2 interior beads. A too-coarse seg (< 3 nodes) **raises**, not silently degrades to a rod.
- **Branch topology / mother fraction (`cortex_arp23_mother_fraction = 0.2`).** The Arp2/3 population is split
  into branch **roots** (mothers) + **daughters**; each daughter attaches at one Arp2/3 branch junction. 0.2 ⇒
  ~4 daughters per mother (a dendritic fan-out echoing the lamellipodium). This is a topology choice, not a
  rate — the native flux-limited nucleation KMC (`aleph/components/weave/nucleation.py`) is the dynamic story; here we seed a
  resting census.
- **Branch angle θ₀ = 70° ± 9°.** NOT a choice — **reused** from the Arp2/3 angle-harmonic path (Fäßler 2020
  in-cell cryo-ET; `ARP23_BRANCH_ANGLE_RAD` / `ARP23_K_THETA` via equipartition σ_θ=9°).

---

## 4. What was added (all additive, default-off)

- **`aleph/components/weave/regions.py`**
  - `derive_cortex_arp23_split(n_total, arp23_mass_fraction, formin_length_um, arp23_length_um)` — the mass→count
    map above (pure, testable).
  - `cortex_arp23_region(...)` + `_cortex_arp23_arch(...)` — an Arp2/3 cortical-sub-population `RegionSpec`
    (manifold `"sphere_dendritic"`, `nucleator="arp23"`, θ₀=70°/σ_θ=9°).
  - `_build_cortex_arp23_region(...)` — sphere-native dendritic builder. **REUSES the Arp2/3 branch path**:
    `branch_angle.sample_branch_angles` + `ARP23_THETA0_RAD`/`ARP23_K_THETA`, and emits `branch_triples` /
    `branch_anchors` in the **same layout** as `build_lamellipodium`, so they feed the **same device
    `branch_angle_kernel`** the lamellipodium uses. The only NEW code is the sphere placement (isotropic bases +
    tangent-plane branch rotation) — the 72°/70° angle-harmonic branching is **not reinvented**.
  - `build_region` dispatch for `"sphere_dendritic"`.
- **`aleph/components/incumbent/assemble.py`**
  - `CellConfig`: `cortex_arp23_fraction=0.0` (mass fraction; 0 = OFF), `cortex_arp23_length_um=0.15`,
    `cortex_arp23_seg_um=0.05`, `cortex_arp23_mother_fraction=0.2`.
  - `_cortex_regions(cfg)` — returns the **single** legacy formin region when fraction==0 (byte-identical to the
    prior `_cortex_region(cfg.n_filaments, …)` call → Gate-1 parity), else `[formin(N_f), arp23(N_a)]`.
  - `build_cell` now weaves `_cortex_regions(cfg)` (was the inline single region).
- **`tests/ac/cell/test_cortex_arp23_population.py`** — 13 CPU tests (below); CUDA guarded (no `build_cell`).

Bit-parity is preserved: the formin sub-population is still the delegated `ff.weave.weave` sphere cortex (just
with `N_f` filaments); the Arp2/3 region has no myosin (like the lamellipodium); crosslinks in the Arp2/3
region are the stiff daughter-base↔branch anchors (force-free at the branch geometry).

---

## 5. Test result

`pytest tests/ac/cell/test_cortex_arp23_population.py` → **13 passed** (0.5 s, CPU). Covers:
- `fraction=0.0` ⇒ single region == legacy `_cortex_region` (bit-parity), no branches.
- mass→count split preserves the budget + round-trips the mass fraction; the illustrative **5.00×** case;
  the default-lengths **(6514, 64172)** split; bad-input guards.
- `fraction=0.33` ⇒ two regions, derived counts, total == 70,686 (areal density preserved), same shell radius.
- the built Arp2/3 region carries a valid branch per daughter, junction angles cluster at θ₀ (70°), and a
  too-coarse seg raises.
- `weave_cell([formin, arp23])` (small) concats into one network with the branches + a correct population ledger.

Existing regressions still green: `test_cortex_resolution_config`, `test_cortex_parity`, `test_regions`,
`test_woven_cell`, `test_lamellipodium_oracle`, `test_branch_angle_oracle` and the full `tests/ac/weave/`
(77 passed, 1 CUDA-skip).

---

## 6. The native call (for the Lead)

The prototype builds the **topology** (host NumPy). To close the loop physically:

1. **Build + compare on the gbook A5000** at full native population:
   `build_cell(CellConfig(cortex_arp23_fraction=0.33))` vs the formin-only baseline. Report node/branch/byte
   deltas (Arp2/3 at default = 64,172 filaments × 4 nodes ≈ 257k nodes + ~51k branch junctions, on top of the
   6,514 formin × 7 nodes). Verify the `branch_angle_kernel` is wired to the woven `branch_triples` so the
   Arp2/3 junctions carry real angle-harmonic force (the CPU seed only lays geometry).
2. **PI A/B/C ratification** of `L_a` (0.10 / 0.15 / 0.20 µm), the mother fraction, and whether the formin `L_f`
   used in the mass split should be the coarse 3.0 µm or a shorter native value — the count ratio is sensitive
   to both (§2 table).
3. **Inter-population crosslinking** (Arp2/3↔formin binding into one mesh) is a follow-up: today the Arp2/3
   region only carries its own branch anchors. A cross-region bond (`CrossRegionBond`) or a shared-reach
   crosslink pass would bind the two sub-populations — flagged, not built here.
4. **Per-stage visualization** (isolation + cumulative) per the PI 2026-07-22 ladder gate once it builds native.
