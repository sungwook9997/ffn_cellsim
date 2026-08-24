---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# MDA-MB-231 Migration Failure — Root-Cause Analysis (2026-07-14)

**Result.** Native MDA-MB-231 crawl reaches steady **~1.0 nm/s** (mda231_diag_v3_long, 300 s), **~6× below**
the 3D-collagen target (0.39 µm/min = 6.5 nm/s). This is a genuine **FAIL** reported honestly (not tuned).
This document is the root-cause decomposition — quantitative force balance + code-level mechanism trace.

---

## 1. The quantitative force balance (from data, both cells)

At steady state, propulsion = drag, so the net forward force is **F_net = γ_COM · v** with
γ_COM = 6πηR (the physical `--com-drag` COM drag):

| | v (steady) | γ_COM | **F_net** | total clutch traction | F_net / F* |
|---|---|---|---|---|---|
| **MDA-231** | 1.08 nm/s | 1513 pN·s/µm | **1.6 pN** | 756 pN | 0.23 clutch |
| **MCF7** | 0.76 nm/s | 9316 pN·s/µm | **7.1 pN** | 596 pN | 1.0 clutch |

→ **The net forward force is ¼–1 clutch's worth (F* = 7 pN).** To reach 6.5 nm/s at MDA drag you need
F_net ≈ **9.8 pN (6× the current 1.6)**.

**Note (honesty):** the "756 pN total" is the sum of clutch-**spring** extension magnitudes — the radial
contractile GRIP, which is *supposed* to be radial. It is a different budget from the propulsion channel
(see §3). Do NOT read "1.6/756 = 0.22% efficient"; read it as "the radial grip vector-cancels laterally,
and the separate propulsion channel under-delivers."

## 2. Where the traction goes (ECM displacement, both cells)

Near-footprint collagen displacement, decomposed vs the footprint centre (radial) and migration axis (polar):

| | radial | polar | slab bulk drift | far-field nodes |
|---|---|---|---|---|
| **MDA-231** | **−11 nm** (centripetal) | +8.5 nm | +4 nm (≈0) | 1.8 nm (pinned) |
| **MCF7** | **−12 nm** (centripetal) | +0.1 nm | — | — |

→ Both cells **contract the collagen inward (centripetal) ~symmetrically**. The slab is **properly
anchored** (bulk drift ≈ 0, far-field pinned, median node disp = 0). So **anchoring is NOT the problem and
drag is NOT the problem** — the cell has a fixed substrate and low physical drag, yet still doesn't move.

## 3. Code-level mechanism — there is **no locomotory force dipole** (radial-by-construction)

Traced in `fa_ecm.py`, `motility_warp.py`, `network_warp.py`, `ff_crawl_on_substrate.py` run(). Exactly two
traction terms; only one carries the polarity axis `phat`:

- **Dominant term = radial contractile grip (the 596–756 pN).** `clutch_ecm_spring_kernel`
  ([fa_ecm.py:38-50](../../ff/fa_ecm.py)) is a geometric Hookean spring pulling actin toward its gripped
  collagen node — **no polarity axis**. It is loaded by `myosin_kernel`
  ([network_warp.py:105-120](../../ff/network_warp.py)) = **constant f_myo, isotropic, no polarity** → the
  cortex contracts every basal clutch's collagen node centripetally → radial. Summed over the ~symmetric
  footprint the lateral vector sum ≈ 0. **This is the radial contraction; it is correctly radial and cancels.**
- **Only forward channel = slip traction (≤ F*·N_bound ≈ 160 pN), and it self-cancels.**
  `clutch_ecm_slip_traction_kernel` ([fa_ecm.py:53-84](../../ff/fa_ecm.py)) adds a coherent **+phat** push,
  but (i) capped at F* = 7 pN/clutch, and (ii) its **−phat reaction lands on the MOVABLE collagen node**,
  which the spring then drags the cell back toward → the forward push is largely returned/cancelled. Net ≈ 1
  clutch.
- **`polarize` / `ecm_regrip` add NO force — they only bias BINDING.** polarize front-bind
  ([ff_crawl_on_substrate.py:731-733]) forms nascent adhesions only at the front; rear de-adhesion
  ([:748-751]) simply **UNBINDS** rear clutches (`bd[_rear]=0`) → the rear generates **zero** forward force;
  ecm_regrip re-maps anchors. Net-forward is therefore only a *residual of biased binding* + the capped,
  self-cancelling slip term.
- **`myo_rear_bias` (rear-biased contraction) is DROPPED** — numerically unstable at native (OverflowError
  in from-resting), all presets set 0.0 ([cell_type.py:66]). **Protrusion injects no net force** — it is
  COM-conserving rest-length growth by design ([motility_warp.py:210-213]); the old leading-edge body-force
  push is retired.

## 4. Root cause (one sentence)

**The FF crawl engine has no front-back locomotory traction DIPOLE: the dominant force is an isotropic
myosin-driven radial contraction (which correctly cancels laterally), and the only polarized channel — the
F*-capped slip traction — self-cancels because its rearward reaction is delivered to the movable matrix and
the rear edge merely unbinds (zero force) instead of retracting the body forward. Net-forward = ¼–1 clutch
(1.6–7 pN), 6× short.**

This is a **coherence / cancellation** deficit, **NOT**: drag (com_drag gives physical low γ), substrate
anchoring (slab pinned), per-clutch magnitude (F* is physiologically correct), or total force (596–756 pN of
grip exists). It is **universal** (MCF7 and MDA are the identical mechanism); MDA is faster only because its
6× lower η lets the same tiny net force move faster — **not** because 231 mechanics improved propulsion.

Deceleration (2.95→1.0 nm/s) = the cell densifies/contracts the local collagen symmetrically (recruit grows
+185→+230 nm) and embeds in a symmetric well, opposing the residual polar force → v decays.

## 5. Where the fix must go (the missing dipole)

Convert ~2 clutches of *uncancelled* dipole (not more magnitude): 

1. **Rear-retraction traction (largest gap).** [ff_crawl_on_substrate.py:748-751] — rear clutches currently
   only UNBIND. Add a +phat retraction traction keyed on the rear cap (`_s < 0`) so the trailing edge snaps
   the body forward as it releases (the back half of the dipole).
2. **Uncancelled front anchorage.** [fa_ecm.py:82-84] — deliver the front slip's −phat reaction to the
   pinned matrix (or hold the front with a non-slipping engaged clutch) so the forward push survives instead
   of being dragged back through the spring.

Edit sites: `clutch_ecm_slip_traction_kernel` (fa_ecm.py:53-84) + dish twin `clutch_slip_traction_kernel`
(motility_warp.py:308-327) + the polarize turnover block. Keep the F* cap (per-clutch magnitude is correct).

**This supersedes the earlier "drag-solver / free-draining" framing** in the G3 gate / RunResult notes: the
drag is not the bottleneck for the COM (com_drag handles it). The bottleneck is the **missing locomotory
force dipole** — a mechanism gap, PI-gated (touches the propulsion force model). Contractility↑ and raw
adhesion↑ were already shown NOT to help (M6) — consistent: more isotropic myosin = more radial contraction,
not more dipole.

## 6. Verification status
- Force balance F_net=γv: computed from mda231_diag_v3_long + emt_flow_clamp npz (steady windows).
- Radial/polar + slab-anchor: computed from ecm_pos0/posf + ecm_frames.
- Mechanism: code trace (read-only) confirmed the two-term structure + binding-only polarize + dropped
  rear-myosin + retired protrusion-force.
- Open (optional): a run with a rear-retraction dipole term would test §5 directly (PI-gated force change).
