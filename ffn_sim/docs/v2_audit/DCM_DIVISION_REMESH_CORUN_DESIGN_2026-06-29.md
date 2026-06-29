# DCM division + remesh co-run — design note (2026-06-29)

PI goal: "proliferation을 어떻게든 잘 분리되는 것까지 … 메쉬 깨지지 않도록." Clean proliferation
with the daughters' stretched/cut meshes kept in-band needs **division and remesh to co-run**.
They are currently mutually exclusive — `run_decohesion` force-disables remesh whenever
`--division` is set ([dcm_warp_decohesion.py:280](../../dcm/dcm_warp_decohesion.py#L280)):

```
if division and remesh_period:
    print("  [warn] division + remesh both requested — disabling remesh this run …")
    remesh_period = 0
```

The header comment says the blocker is a one-line cof-sentinel disambiguation (remesh pool
== −1 vs parked cell == −2). **On closer reading it is not that simple** — recording why, so the
gbook session / PI can pick an approach rather than ship a fragile sentinel swap blind.

## Why the one-line sentinel swap breaks cleave

Two division mechanisms exist:

- **Mitotic** (`DivisionHost.update`, default, `--division`): activates a whole parked icosphere
  from a reserved cell-ID slot, indexed by **fixed npc-blocks** (`mother*npc:(mother+1)*npc`).
- **Cleave** (`DivisionHost.cleave_plan` + `cleave_cell`, not CLI-exposed): carves the mother
  shell along the Hertwig plane; **topology-agnostic** (works after the mesh is relabelled).

The dormant-node bookkeeping is entangled:

1. **Parked-cell node blocks double as the dormant NODE supply for cleave.** `cleave_cell` draws
   its cut-ring nodes from `cof < 0` ([dcm_cleave.py:88](../../dcm/dcm_cleave.py#L88)), and the
   build allocates **no separate remesh pool when remesh is off** (`pool_factor=(pool_factor if
   remesh_period else 0.0)`, [dcm_warp_decohesion.py:323](../../dcm/dcm_warp_decohesion.py#L323)).
   So cleave relies on the parked-cell nodes being `cof < 0`. Re-tagging parked cells to `−2`
   and making `cleave_cell` draw only `−1` would starve cleave of ring nodes.

2. **Mitotic is fundamentally incompatible with remesh** regardless of the sentinel: remesh
   SWAP/SPLIT/COLLAPSE **relabels `cof` and resizes faces**, after which `update`'s fixed
   npc-block indexing (`mother*npc:…`) points at the wrong nodes. A remesh COLLAPSE that parks
   a cell's *first* node also flips `_cell_active` (which tests `cof[cell*npc] >= 0`) to dead.

## Two viable approaches (pick with PI / on gbook)

**A. Index-range disambiguation (least invasive, mitotic stays remesh-free).**
Keep all dormant nodes at `cof == −1`. Disambiguate by **node-index range**, which the build
already separates: parked cells occupy `[n_cells·npc : n_total·npc]`, the remesh pool occupies
`[n_total·npc : ]`. Pass `remesh_pool_start` into `split_edge` so remesh SPLIT only activates
nodes `≥ remesh_pool_start`; cleave keeps drawing from the parked-cell range. No sentinel change,
cleave untouched. Still does **not** fix mitotic+remesh — only enables **cleave + remesh**, which
is the correct (topology-agnostic) co-run anyway. Requires exposing `--cleave` on the CLI.

**B. Unified pool manager (correct, larger).**
A single `NodePool` object owns all dormant nodes and hands them out by purpose (remesh-split vs
cleave-ring vs mitotic-cell-activation), tracks frees (COLLAPSE returns, division consumes), and
exposes `is_active(cell)` topology-agnostically (any node with `cof==cell`, not the first block
node). This removes the fixed-npc-block assumption entirely and lets *either* division path
co-run with remesh. This is the SimuCell3D/CellSim3D pattern.

## Known caveat (must be checked on gbook, not assumed away)

A prior note ([dcm_warp_decohesion.py:758](../../dcm/dcm_warp_decohesion.py#L758)) records that a
single `remesh_pass` right after a cleave made mesh **quality WORSE** — the fresh high-valence
cut/cap config defeats the conservative COLLAPSE. So cleave+remesh needs either a remesh cadence
that lets the cut config relax first, or a cleave cap that emits in-band triangles. Validate on
gbook at N≥400 before declaring the co-run "clean proliferation".

## Recommendation

Approach **A** + expose `--cleave`, validated at N≥400 on gbook (the only place the
separation+mesh-integrity claim is real). Deferred from the CPU autonomous loop because (1) it
can't be production-validated offline, (2) the quality caveat is a real risk, (3) it is a
gate-adjacent architecture change worth a PI glance. The cof comment at line 277-279 should be
corrected to reference this note (the "−2" one-liner is insufficient).
