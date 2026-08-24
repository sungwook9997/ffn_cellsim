# INTEGRATION — protrusion (SEAMED → KERNEL_BOUND, lamellipodium branch-angle)

**Owner file touched:** `aleph/engine/protrusion.py` + `aleph/tests/ac/engine/test_protrusion.py` only.
No shared file (dispatch/world/contracts/load_path) was edited. This note lists what the Lead must do to
compose the new binding and to advance it from KERNEL_BOUND to CUDA_UNIT on the A5000.

## What landed

`LamellipodiumBranchAngleMechanics` — a `MechanicsContributor` delegate that binds the **real** Warp kernel
`aleph.components.weave.branch_angle_warp.branch_angle_kernel` (Arp2/3 3-body angle-harmonic) through the existing
mechanics seam. Evidence state for the lamellipodium **branch-angle term** moves SEAMED → **KERNEL_BOUND**:
the kernel runs through the seam with **no private state path** — it reads the actor-owned
`branch_triples_d` / `branch_active_d` and scatters force directly into the actor-owned `force` array the seam
passes in.

- Constants are the sourced **Faessler 2020** anchors only (`theta0 = 1.2217304764 rad`,
  `k_theta = 0.173 pN·µm/rad²`, equipartition-derived), imported from `ac.weave.branch_angle`. Off-source
  values and a non-real kernel are rejected at construction (no magic numbers).
- `LamellipodiumStateOwner.__post_init__` now calls `mechanics.assert_bound_to(...)` when the injected
  mechanics is this delegate, so a delegate carrying a **private copy** of the branch topology is rejected —
  the kernel is proven to launch against the actor's authoritative arrays.
- The `launch` seam defaults to `wp.launch`; tests inject a capturing double (`_LaunchSpy`) so the structural
  gate is CPU-green on the dev Mac (no CUDA here).

## Composition (Lead, A5000)

Construct the delegate with the **same** branch arrays the owner owns, then inject it:

```python
mech = LamellipodiumBranchAngleMechanics(
    n_branch=owner_n_branch,
    branch_triples_d=branch_triples_d,   # the SAME arrays passed to LamellipodiumStateOwner
    branch_active_d=branch_active_d,
    # launch defaults to wp.launch; theta0/k_theta default to the Faessler anchors
)
LamellipodiumStateOwner(..., branch_triples_d=branch_triples_d, branch_active_d=branch_active_d, mechanics=mech, ...)
```

No `dispatch.py` / `world.py` change is required for this term: composition is through the existing
`MechanicsContributor` injection. If/when the Lead registers a dedicated dispatch task for the lamellipodium
`ASSEMBLE_INTERNAL` phase, the delegate's `accumulate(pos, force)` is the canonical single call.

## Remaining gap to CUDA_UNIT (A5000, GPU lane — cannot run on the Mac)

Run the branch-angle kernel through this seam on CUDA and close the gates already specified in
`branch_angle_warp.py` / `branch_angle.py`:

1. **Force/oracle parity:** device force == host oracle `ac.weave.branch_angle.angle_forces` to tolerance.
2. **Newton 3rd / sign:** the three nodal forces of each triple sum to zero; a bent junction restores toward
   `theta0`; dormant (`active==0`) branches contribute zero force.
3. **Precision/dimensional:** float64, pN·µm units, on a non-empty native branch population.
4. **Rollback:** bit-exact restore of `force`/topology under a rejected predicate (transaction seam already
   present; needs the real snapshot/rollback kernel bound, still a `_TransactionSpy` here).

## Still SEAMED (not advanced this increment — honest scope)

- **Filopodium mechanics is deliberately NOT bound.** Its bundle cross-link stiffness `k_fascin`, spacing,
  bundling angle, and per-bundle count are all `value: null / GAP — PI` in `ac/weave/params_i0b4.yaml`.
  Binding a filopodium spring kernel now would require a magic number → blocked on PI source cards.
- **Lamellipodium stretch/bend/excluded-volume + membrane/cortex/cytosol/FA connector loads** remain seam
  protocols/spies; only the branch-angle term is kernel-bound.
- **Reaction kernels** (branching/polymerization/capping/severing commits) and **crosslink KMC**
  (`crosslink_kmc_warp`, which also depends on GAP hand params `k_off0/f0`) remain accepted-step spies.
