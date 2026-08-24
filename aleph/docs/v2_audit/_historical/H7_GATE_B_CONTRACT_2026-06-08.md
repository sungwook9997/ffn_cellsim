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

# H.7 Gate-B — validation contract (written BEFORE the run)

**Status:** DRAFT for PI ratification. Branch `h7/full-cell-integration`. Author: Lead.
This is a gate contract (sanity-gate protocol): the predicate, observable, pass/fail
criterion, reference datum, controls, and magic-number discipline are fixed here,
before the relaxed-constraint mode is built or run. No gate-loosening: the criterion
is not editable to make a run pass — a wrong criterion is surfaced to PI for a
contract change.

Supersedes the (deleted) `cortex/CORTICAL_TENSION_PROTOCOL_REVIEW_2026-06-08.md`,
which was built on two errors corrected here: (a) the 0.391 mN/m "in-band structural"
was a discarded FA-integrin-overload artifact, not a measurement; (b) the relaxed
lever is the **compression/buckling** degree of freedom, not axial stretch.

---

## 1. What Gate-A established (locked — not re-opened)

`H7_GATE_A_RESULT_2026-06-08.md` (REFUTE, firm, 35-sample plateau), full physiological
MCF7 cell, grip_walk myosin genuinely walking, native GPU stack:

- s_grip develops to ℓ₀-scale (0.024 → 0.58, 22×) and heads bead-step (real material
  transport).
- γ_soft (active soft-MOP) rises **16.2×** from the loading floor (1.89e-4) to a stable
  **contraction plateau 3.06e-3 ± 0.39e-3 mN/m** — a REAL sustained contraction.
- It caps there: **1/114 of the proxy band, ~1/88 of the MCF7 datum 0.27**.
- `r/r0 = 1.0000` throughout: the constrained backbone does **not** condense
  (no buckling/densification — the Murrell/Lenz/Miyazaki symmetry-breaking).

**Verdict (locked):** the active-γ wall is NOT generation (heads bind, load, walk, step,
contract) and NOT s_grip delivery (γ_soft responds 16×). It is **transmission/lever**:
the real contraction does not transmit into spanning hoop tension on a non-condensing
backbone.

## 2. Gate-B predicate (the question)

> Does permitting filament **buckling / condensation** (relaxing the M-SHAKE constraint
> on the COMPRESSION side only) transmit Gate-A's real contraction into spanning
> cortical tension — i.e., does the emergent cortical γ rise from 3.06e-3 mN/m toward
> the MCF7 datum?

CONFIRM ⇒ buckling/condensation is the transmission lever (the model reaches the right
order untuned). REFUTE ⇒ the wall is elsewhere (soft long-range transmission, or a
generation-density shortfall) and buckling alone does not close it.

## 3. Mechanism under test + the degree-of-freedom correction

The lever is **compression-side yielding**, NOT axial compliance:

- **Tension side** (bond stretched, |r| > ℓ₀): stays ~inextensible. Physiological —
  actin axial stiffness k_axial = EA/ℓ₀ ≈ 154 pN/nm (E≈2 GPa, r≈3.5 nm, ℓ₀≈500 nm),
  so a ~4 pN myosin load strains a segment ~0.026 nm ⇒ r/r0 ≈ 1.00005 (matches Gate-A's
  1.0000235). Relaxing this would be unphysical and is NOT done.
- **Compression side** (|r| < ℓ₀, or filament under axial compression): permit the
  filament to **buckle** (Euler, bending-governed) / the network to **condense**, the
  symmetry-breaking a fixed-length constraint suppresses. Real slender actin buckles
  under compression.

**Design-investigation item (not a contract term):** the backbone already has soft
bending angles (`connected_mesh.py` cortex-angle harmonic), so buckling is geometrically
possible in principle. Why Gate-A still shows no condensation must be isolated during
relaxed-mode design — candidate suppressors: (i) the M-SHAKE fixed bond length removing
the compression DoF that seeds buckling; (ii) bending stiffness angle_k above the Euler
threshold at the segment scale; (iii) loading geometry keeping the backbone in tension.
The relaxed mode must target whichever genuinely suppresses condensation, and the
controls (§5) must show condensation actually engaged.

## 4. Observable + measurement protocol

Report ALL channels separately (turgor never folded in — B4 discipline). The GATE
criterion is on the **myosin-attributable** cortical tension, isolated by differencing
over the physiological baseline (the physiological-baseline HARD-rule worked example:
measure on the resting cell, add myosin as the modulator):

> **γ_active = γ_cortical |(myosin ON) − γ_cortical |(myosin OFF)**, same build, same
> turgor, at the settled operating point.

where γ_cortical is the constraint-inclusive structural cortical stress (soft-MOP +
Lagrange-MOP), so a contraction transmitted either through soft cortical bonds OR
through condensation is counted, and passive load (turgor/thermal/EV) cancels in the
difference rather than being dropped by channel omission. Raw γ_soft, γ_rigid,
γ_passive, and the membrane component (γ_mem ≈ 0.10 mN/m, currently omitted — add as a
4th reported component) are all reported alongside.

## 5. Controls (mandatory, written before the run)

1. **Rigid-limit parity:** with buckling disabled (compression constraint restored),
   reproduce Gate-A γ_soft = 3.06e-3 mN/m to within seed noise. Proves the relaxed mode
   is a clean superset of rigid M-SHAKE (nonconv = 0 preserved).
2. **Tension-side inextensibility preserved:** bonds under tension keep r/r0 ≈ 1
   (k_axial physiological). The relaxation is compression-only — NOT a soft-everywhere
   backbone.
3. **Condensation actually engaged:** report a buckling/condensation observable
   (filament end-to-end shortening, mean curvature, or pair-density increase). γ moving
   without condensation engaging would mean the wrong DoF was relaxed → invalid.
4. **Myosin-OFF baseline:** the differencing reference; isolates the myosin rise from
   turgor/thermal structural load.

## 6. Operating point (RATIFIED PI 2026-06-08)

**Gate-B is measured on a suspended/rounded interphase cell** (turgor ON at the resting
setpoint, FA off/minimal), observable = cortical-γ. This matches the in-vitro condition
of the only real MCF7 datum (Hosseini 2020, suspended/rounded interphase), avoids the
FA-integrin-overload artifact (the source of the discarded 0.391), and is the cell state
for which cortical-γ is the correct observable. The FA-adhered operating point and the
traction-observable pivot were considered and not chosen for this gate.

Note vs the physiological-baseline HARD rule: suspended/rounded-with-turgor is a fully
physiological cell state (not a null/unphysical baseline) — it is the state of the
experiment being reproduced. The adherent collective operating point remains the
ultimate platform target (PI-exp overlay) and is gated separately (likely on traction).

## 7. Pass / fail criterion + reference datum (RATIFIED PI 2026-06-08)

Reference: **Hosseini 2020 (Adv Sci 7(19):2001276; AFM parallel-plate confinement,
Fischer-Friedrich method), MCF7 control, suspended interphase: γ ≈ 0.27 mN/m, IQR
0.18–0.40, whiskers 0.05–0.62.** The CONFIRM band is the **Hosseini IQR [0.18, 0.40]**.
The [0.35,0.65] band is a non-MCF7 (HeLa/L929) proxy, retained only as a rounded-cell
envelope, not the MCF7 anchor.

| outcome | criterion (on γ_active, §4) | reading |
|---|---|---|
| **CONFIRM** | rises into **[0.18, 0.40] mN/m** (Hosseini IQR) | buckling IS the transmission lever; model reaches the right order untuned |
| **PARTIAL** | rises ≥10× above the Gate-A plateau (≥~0.03) but < 0.18 | buckling contributes but is insufficient; quantify the residual gap + next lever |
| **REFUTE** | < 10× the Gate-A plateau (< ~0.03) | buckling is NOT the lever; wall is soft long-range transmission or generation density |

## 8. Magic-number discipline

The compression-yield / buckling threshold is the **physiological Euler buckling force**
F_crit = π²κ/L_seg² from the actin bending stiffness κ = ℓ_p·k_BT (ℓ_p ≈ 10–17 µm,
KB-anchored) and the segment length L_seg — **derived, grid-invariant, not a free knob**.
There is no parameter tuned to land γ_active in the pass range; if the derived threshold
cannot be set without a free parameter, halt and surface to PI.

## 9. What Gate-B does NOT do

- Does not set an unphysical axial k_axial (tension side stays inextensible).
- Does not re-open the Gate-A REFUTE (γ_soft genuinely plateaus at 3.06e-3 on the rigid
  backbone).
- Does not chase the band; the criterion is the named MCF7 datum, fixed here.
- Does not fold turgor into the active number (differencing keeps active/passive split).

## 10. Gate-contract decisions — RATIFIED PI 2026-06-08

1. **Operating point:** ✅ suspended/rounded + cortical-γ (§6). FA-adhered and the
   traction pivot were not chosen for this gate.
2. **Pass band:** ✅ Hosseini MCF7 interphase IQR **[0.18, 0.40]** (§7).
3. **Adherent-cell observable:** the adherent collective operating point's gate (likely
   **traction**, not cortical γ) is tracked separately as the platform-target gate; not
   part of this Gate-B.

**Contract status: LOCKED.** Next = build the buckling-permitting relaxed-constraint
native mode (§3 design-investigation: isolate the condensation suppressor), then run
Gate-B on the suspended/rounded cell per §4–§7.
