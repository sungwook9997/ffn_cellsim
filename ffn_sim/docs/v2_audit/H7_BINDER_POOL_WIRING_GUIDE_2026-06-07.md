# Binder → attachment-pool wiring guide (the ~8.4–16× lever) — 2026-06-07

The dominant remaining Gate-A cost is the binders' per-firing **global `set_snapshot`**
(51.9 ms/firing — `h7_native_gate_a_profile`). The native `FFNAttachmentSpringForce`
pool (3a26565, parity 3.4e-16, toggle 17 µs, end-to-end 5.3× in `h7_binder_pool_demo`)
removes it. This is the **exact minimal change** to wire it — physics-preserving, so it is
low-risk; the binding *decisions* (Steps 1–2) are untouched, only the bond-mutation
mechanism (Step 3) changes. Sub-session enabler is built + validated; the wiring is the
Lead's myosin/xlink lane (their `act()`), hence this guide rather than a blind port.

## XlinkBondUpdater (cortex/crosslinkers.py) — the change

**Build/attach (once):** replace the HOOMD `xlink_attach_b*` bonds + their harmonic bond
force with one pool force, sized to the head count:
```python
from ffn_hoomd_plugin import NativeAttachmentSpringForce
self._attach_force = NativeAttachmentSpringForce(pool_size=self.n_heads)
integrator.forces.append(self._attach_force)         # replaces the attach bond force
# track the bind bin per head (for per-bond r0); -1 = unbound
self._head_bin = np.full(self.n_heads, -1, dtype=np.int64)
```
Set `self._head_bin[head_local] = bin_idx` at each bind (Step 2, line ~1117) and `= -1`
at each break (Step 1, line ~1038), alongside the existing `_head_bound_to_actin` updates.

**act() Step 3 (lines 1129–1172) — REPLACE the whole snapshot rebuild + `set_snapshot`
with:**
```python
bound = self._head_bound_to_actin                    # (n_heads,) actin_tag or -1
active = bound >= 0
head_tags = (np.arange(self.n_heads) + self.n_cortex_actin).astype(np.int32)
actin_tags = np.where(active, bound, -1).astype(np.int32)
bin_r0 = xlink_attach_bin_rest_lengths(self.p.n_bins, self.p.max_bind_dist)
r0 = np.where(active, bin_r0[np.clip(self._head_bin, 0, None)], 0.0)
kv = np.where(active, self.p.k_attach, 0.0)
self._attach_force.set_attachments(head_tags, actin_tags, kv, r0)   # 17 µs, NO set_snapshot
```
The **other bonds (backbone, FA clutch, myosin backbone)** stay HOOMD bonds — they are not
mutated by this updater, so they need no snapshot at all. The `get_snapshot()` READ at the
top stays for now (3 ms, not the wall); a later optional win is `gpu_local_snapshot` for the
positions feeding the off-rates + KDTree.

Note the per-cortex-bead degree budget (`_MAX_CORTEX_BEAD_DEGREE`) was a HOOMD bonded-exclusion-cap
workaround; with the pool the attach springs are NOT HOOMD bonds, so they no longer consume the
7-exclusion budget — the cap can be relaxed for attach (keep it only for real HOOMD bonds).

## MyosinStepUpdater (cortex/myosin.py) — same shape
Same pattern: the head↔actin grip springs become a second `NativeAttachmentSpringForce` pool
(or share one with a head-tag offset); the grip_walk `act()` keeps its s_grip stepping +
bind/unbind decisions but writes the pool via `set_attachments` instead of `set_snapshot`.
The active stepping force (Hill F–v) is the spring with the walking r0 — set `r0` to the
current grip rest length per head.

## Validation (before trusting the swap)
Statistical parity vs the current set_snapshot binder over a warmed cell: same per-tick
bind/unbind counts, same γ_soft + γ_rigid distributions (Welch |z|<3), nonconv=0. Harness
pattern: `h7_native_gamma_parity` / `h7_native_gate_a_verify` (compare an arm with the pooled
binder to one with the set_snapshot binder).

## Projected: removes the 51.9 ms/firing set_snapshot (the 62.6% step slice) → ~8.4–16× full-cell
(the binding decisions, vectorised as today, are cheap once the snapshot rebuild is gone). Pair
with `FFN_GPU_DEVICE_COMPARTMENTS=1` + `FFN_GPU_DEVICE_BAOAB=1` (see H7_GATE_A_NATIVE_SMOKE).
