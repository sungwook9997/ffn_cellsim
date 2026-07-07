# GPU cadherin port — engineering plan (2026-07-07)

## Problem

The DCM cadherin junction (`CadherinBondHost`, `dcm_cadherin_host.py`) manages bonds **on the CPU**:
every `batch_steps` (=50) the loop does

```
cad.update(pos_d.numpy())   # GPU→CPU download of ALL node positions
  ├─ BREAK: numpy force→koff→stochastic rupture over existing bonds
  └─ FORM : scipy cKDTree over free nodes → query_pairs(r_bind) → greedy ≤1-bond accept
cad.upload(device)          # push new bond list back to GPU
```

The **bond FORCE** is already a Warp kernel (`cadherin_bond_force_kernel`, on GPU). Only the
**bond MANAGEMENT** (break/form) is CPU. This violates the GPU-main directive (PI 2026-05-31) and
is the N-scaling bottleneck: the `pos_d.numpy()` download + `cKDTree` build/query scale with node
count, every 50 steps.

**Measured (2026-07-07, N=64 A5000):** implicit + `accel_dt=8e-5` (10× dt) is stable
(V/V0=1.000, cfl=0) but only **~1.5× faster** overall — the implicit CG solve + the CPU cadherin
still dominate. So `implicit` alone does **not** make N=2000 feasible; the CPU cadherin must go GPU.
Projected N=2000 full-compartment: ~2–3.5 h (baoab) / ~1.3–2.3 h (implicit) — both too slow, and
the CPU cadherin fraction grows with bond count (~70k bonds at N=2000).

Already validated this session: implicit+cadherin runs stable; per-node **stress** (`fmag`) +
per-node **junction count** (`nbond`) + per-frame **bond pairs** (`bondpairs`, for drawing the
junctions as lines) now save from the sim (`dcm_warp_decohesion.record()`), verified at N=64.

---

## Option 1 — GPU cadherin port (CHOSEN)

Move break + form + bookkeeping onto Warp kernels; bonds live on the device; **no GPU↔CPU
round-trip per batch**. Reuse the existing `wp.HashGrid` infrastructure (already used for
contact/cohesion neighbour search in `dcm_neighbor_warp.py`).

### Kernels (new: `dcm/dcm_cadherin_gpu.py`)

1. **`cad_break_kernel`** — one thread per existing bond `b=(i,j)`:
   - drop if `cof[i]<0 or cof[j]<0` (remesh-collapsed node);
   - `L=|P[i]-P[j]|`, `F=k_trans·max(0,L-r0_trans)`;
   - `koff` by linear interp of the **catch-slip table on the device** (small `wp.array` of
     `(fs, koff)` uploaded once); `p_break=1-exp(-koff·dt_batch)`;
   - keep if `wp.randf(state) >= p_break`; write survivor into a compacted array via `wp.atomic_add`
     on a device counter (stream compaction).

2. **`cad_partner_kernel`** — one thread per **free** node `i` (cof≥0, currently unbonded):
   - `wp.hash_grid_query(grid, P[i], r_bind)` → nearest node `j` with `cof[j]≠cof[i]` **and** `j` free;
   - stochastic gate `wp.randf(state) < p_on`; else `partner[i] = -1`; store `partner[i]=j`.

3. **`cad_form_kernel`** — one thread per free node `i`:
   - **mutual-nearest**: if `partner[i]==j and partner[j]==i and i<j` → append bond `(i,j)` via
     `wp.atomic_add(count,1)` then write. This **guarantees ≤1 bond per node** (each node has one
     `partner`), needs **no fragile claim-atomics**, and is fully parallel + deterministic.

   > Divergence from the CPU greedy accept: mutual-nearest forms a subset (only reciprocated
   > nearest pairs) → slightly fewer bonds per batch, but the same steady-state class (bonds keep
   > forming each batch until apposed free nodes are matched). Parity target below covers this.

### Bond-array management
Fixed-capacity device arrays (`bonds: wp.array(vec2i)`, `count`), grown ×2 on the host only when
`count` approaches `cap` (rare). Break compacts into a scratch array; form appends. `bonded[]`
(uint8, N) recomputed each batch by a scatter kernel over surviving bonds.

### RNG
`wp.rand_init(seed, batch_index*N + tid)` per thread — deterministic, reproducible, matches the
"one stochastic draw per candidate per batch" semantics.

### Integration
- `CadherinBondHost.update_gpu(pos_d, grid, batch_index)` — new device path; the old `update()`
  (numpy) stays for CPU/parity. `run_decohesion` calls `update_gpu` when `device!='cpu'`.
- The existing `cadherin_bond_force_kernel` consumes the device `bonds`/`count` directly (no change).
- `_dev` grows to hold `partner`, `bonded`, scratch, counter, koff table.

### Parity validation (before any production run)
At **N=64**, run 5000 steps with CPU cadherin and with GPU cadherin from the **same seed/init**;
accept if:
- steady-state **bond count** within ±15% (mutual-nearest forms a subset — a band, not exact);
- **inter-cell gap** median and **cad-churn** (cum_broken slope) match within ±15%;
- cohesion holds (no cluster fly-apart), V/V0=1.000, cfl bounded.
If parity fails, fall back to Option 2 and iterate the GPU rule.

### Then
Re-run N=2000 full-compartment Voronoi+cadherin (implicit + GPU cadherin + stress/junction save),
retrieve, render the FEM-stress + junction-line viewer.

### Risk / effort
~150 lines of Warp + host wiring + parity harness. Core shared files touched
(`dcm_cadherin_host.py`, `dcm_warp_decohesion.py` call-site, new `dcm_cadherin_gpu.py`). Concurrent
FF session does **not** touch these. Medium effort, medium risk (mitigated by the N=64 parity gate).

---

## Option 2 — pragmatic intermediate (NOT chosen; documented)

Keep the CPU cadherin but cut its frequency + use implicit:
- `integrator=implicit`, `accel_dt≈8e-5` (stable, measured) → fewer steps for the same physical time;
- `cad_batch` 50 → **200** → 4× fewer CPU cadherin calls (each still a `cKDTree`, but 4× rarer).
- Net: roughly ~4–6× faster than baoab+cad50. Enough to bring **N≈1000–1500** to a tolerable
  ~30–60 min without touching kernels.

**Limits:** still O(N) CPU cKDTree per batch (won't reach N=2000 comfortably); `cad_batch=200`
coarsens the catch-slip kinetics (bonds updated 4× less often — a modeling change to sanity-check);
does not satisfy the GPU-main directive. A stopgap, not the fix.

---

## Decision
PI (2026-07-07): **write both plans, then do Option 1.** Proceeding with the GPU port + N=64
parity gate before the N=2000 production.
