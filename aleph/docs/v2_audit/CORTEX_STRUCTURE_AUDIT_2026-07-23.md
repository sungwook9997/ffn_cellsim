# Cortex structure audit (2026-07-23) — the native cortex mesh is structurally broken (5 problems)

Branch `codex/ff-ac-codex`. PI intuition ("our cortex structure seems to have a lot of problems") — confirmed
hard, on a **FRESH native A5000 build of current code** (not a stale dump). This is the upstream root of the
resting-gate blockers (convergence rollback + no transmission + head-binding cap); parameter sourcing sits below
this. Deferred from the resting/sourcing work per PI priority (first: gbook sync + native baseline).

## Method

Fresh `aleph.components.incumbent.dump_state --from-resting --n-filaments 70686 --n-inner 60` on the A5000 (gbook, current
code, md5-synced), then a pure-NumPy/SciPy structural audit of `pos_post` (494,802 actin nodes), `fiber_offsets`
(70,686 fibers), `xl` (crosslinks), `mf_head_node` (8,840 myosin heads). Reproduced identically to the 2026-07-22
dump ⇒ these are properties of the committed weave, not an artifact.

## The 5 confirmed problems (fresh native)

| # | problem | measured | should be |
|---|---|---|---|
| 1 | **zero cortex thickness** | every actin node at R = 7.400 µm, **std = 0.0000** (perfect math sphere; filaments all tangential) | ~0.2 µm shell (Clark/Salbreux) |
| 2 | **uniform filaments** | every fiber = exactly 7 nodes = 6×0.5 µm = 3 µm | exponential length dist (KB-3.18) |
| 3 | **heads float off actin** | myosin head→nearest actin node median **0.133 µm** (head radial std 0.14); only **15.9%** within capture 0.05 µm | heads ON their actin (bind gate works) |
| 4 | ⭐ **fragmented network** | segments+xl → **11,430 connected components**, largest only **79.6%**; ~20% (≈100k nodes) in isolated fiber-fragments; **1.00 crosslink/fiber** (70,686 xl) | 1 spanning component; dense crosslinks (α-actinin/filamin) |
| 5 | **steric interpenetration** | **64,517** actin nodes with steric overlap > 0 (~13% of nodes) at t0 | ~0 (excluded-volume respected) |

## Why this is THE root of the resting-gate blockers

- **Convergence rollback (native #1/#2, diagnoses 22c–23d):** a network that is 20% disconnected, crosslinked at
  ~1/fiber (right at the percolation threshold — a random graph with ~1 edge/node fragments), zero-thickness, and
  interpenetrating in 64k nodes is **under-constrained with floppy modes** and cannot equilibrate a balanced hoop
  tension. No parameter set fixes a broken mesh.
- **No transmission (23c §2, native inert-seed):** myosin loading a fragmented, sparsely-crosslinked net cannot
  propagate tension across the ~11k fragments → the cortical tension never reaches the membrane.
- **Head-binding cap (Codex finding 7, native realized-duty 0.17):** problem #3 directly — heads at 0.133 µm can't
  bind at capture 0.05 µm.

## Reconstruction targets (the actual fix — "other work", PI-deferred)

1. **Crosslink density → percolation + carrying capacity.** 1/fiber → many/fiber (physiological α-actinin/filamin
   density) so the network is a single spanning component that can carry balanced hoop tension. (Memory
   `cortex-construction-rebuild` 2026-06-04 claimed "fragmented→CONNECTED SPANNING MESH 99%"; this native build is
   79.6% ⇒ that rebuild is NOT wired into the current weave, or regressed — investigate `ac/weave` /
   `connected_mesh.py`.)
2. **Head-on-actin placement.** Seed minifilament heads at/near their actin so capture 0.05 µm binds them (not
   splayed 0.13 µm off).
3. **Steric-clean t0.** Resolve the 64,517 interpenetrations at build (relax/repack) so the solve doesn't start
   inside a steric wall.
4. **(Lower priority) cortex thickness + filament length distribution** — shell of ~0.2 µm, exponential lengths.

## Status

Native-confirmed on current code. Gate stays OPEN. Reconstruction is the next substantive lane once PI greenlights
(it is the "다른 작업" deferred behind the sync + native-baseline task). Artifacts: the fresh audit ran on
`/tmp/assembled_state_fresh.npz` (gbook); the committed 2026-07-22 dump `outputs/ac/cell_assembled/
assembled_state.npz` reproduces every number.
