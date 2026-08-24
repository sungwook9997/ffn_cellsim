# GATE A — derivation of the resting-convergence threshold (`max|PF| < 0.21 pN`)

**Date** 2026-07-25 · **Status** DERIVATION + RECOMMENDATION, nothing changed · **Owner** PI decision
**Executable companion** `aleph/scripts/ac_gate_threshold_derivation.py` (9/9 checks, CPU, no Warp import)
**Trigger** `AUDIT_WHOLE_REPO_2026-07-25.md` finding 9 · `TRAJECTORY_HOW_WE_GOT_TO_AC_2026-07-25.md` §"stop
treating it as a solver problem", item (2): *"Nobody currently knows whether 0.63 pN is a failure or a pass,
because the threshold has no derivation on disk."*

> **No gate was edited.** `ac_resting_converge.py`, `ac_gate_a_figure.py` and the other 21 files carrying the
> literal are untouched. Everything below is a derivation plus a proposal that requires PI ratification,
> because a gate-contract change is PI-owned (CLAUDE.md, *no gate-loosening*).

---

## 0. Verdict in one page

**The `0.21` literal is not arbitrary, and it is not a physics threshold. It is the force-equivalent of the
runtime's own float64 round-off displacement tolerance, evaluated once at the CLAUDE.md baseline mesh and then
frozen by hand into 54 sites.**

Exactly:

```
executing predicate (aleph/components/incumbent/inner_mechanics.py:296)   dt_mu * max|PF|  <=  tolerance_um
  tolerance_um = sqrt(eps64) * ell        (driver.py:506, ell = assemble.py seg_mean)
  dt_mu        = 0.1 / kmax               (assemble.py:964)

=> implied per-node force bound
   F_pred = tolerance_um / dt_mu = 10 * sqrt(eps64) * ell * kmax
          = 10 * 1.4901161193847656e-08 * 0.4999048934 um * 2828099.2886817832 pN/um
          = 0.2106697371 pN
```

The frozen literal `0.21` is that number to **0.34 %**. Reproduced to **exact float64 round-off**
(`rel = 0.000e+00`) against the committed native artifact
`outputs/ac/resting_native/resting_native_probe2_capture0.6_2026-07-23.json`
(`inner_tolerance_um = 7.449163398500275e-09`, `inner_dt_mu = 3.535943748517097e-08`,
`ledger.kmax_pN_per_um = 2828099.2886817832`). Someone read `tolerance/dt_mu` off a native run in July, wrote
"gate ≈ 0.21" in `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md:79`, and never wrote down where it came from.

**Consequences, in order of severity:**

| # | Finding | Factor |
|---|---|---|
| 1 | The gate is a **machine-precision** criterion, not a biological one. Recompiled in float32 the same formula gives **4882 pN**. | ×23 170 |
| 2 | It is proportional to `kmax`, which at the native config is set by the **WCA near-core CFL guard** — `max_core_stiffness(f_cap=1e3 pN, σ=7 nm, ε from K_EV_PROVISIONAL=1e3)`, a constant whose own comment says **"PROVISIONAL Magic-Number-Block scale"**. Turn steric off and the same formula gives **0.0611 pN**. | ×3.4 |
| 3 | It is proportional to `ell`. The literal was frozen at `ell = 0.5 µm` and then applied unchanged at 0.2 µm and 0.075 µm, where the runtime's own bound is **0.0843** and **0.0316 pN**. The frozen literal has been **silently loosening the gate as the mesh refined**. | ×2.49 / ×6.64 |
| 4 | `max_i \|PF_i\|` is **not a grid-invariant quantity at all** — a nodal force is a force *density* × `ell`, and different force families scale with opposite exponents. No single per-node force literal can survive a refinement program. | see §5 |
| 5 | At 310 K the **equipartition force floor** `sqrt(kBT·k)` on a myosin-loaded node is **2.07 pN** — ten times *above* the gate. The seeded `f_head = 1.5 pN` is itself only **0.73×** the thermal force noise of its own crossbridge. The 0.5–3.6 pN plateaus two weeks of solver work fought are **inside the thermal band of the bonds that carry them**. | ×9.8 |

**Answer to the trajectory audit's question — is 0.63 pN a failure or a pass?** Both statements are true and
they are different statements: 0.63 pN at `ell = 0.5 µm` is **2.99× over** the runtime's own numerical
convergence bound (the minimiser has *not* converged — a real failure of route (c)), and simultaneously
**0.30×** the thermal force floor of the crossbridge carrying it (so it cannot be claimed as a *physical*
non-equilibrium either). The project has been reporting a numerical-convergence failure as if it were a
physics failure.

**Recommendation (§7):** keep route (c) but **compute it per run instead of freezing it** — that is a
*tightening* at every mesh finer than 0.5 µm and a no-op at the baseline. Add a genuinely grid-invariant
physics criterion (route d1, the signed radial residual bias on the reported cortical tension). Use the
thermal floor (route a2) only as a **ceiling on what may be claimed**, never as a pass criterion — it is
looser than 0.21 and adopting it as a pass would be gate-loosening.

---

## 1. Scope, and what was actually read

Every number below comes from a definition site on disk. Nothing is invented; where a quantity genuinely does
not exist in the code, §4.1 says so instead of guessing.

| quantity | symbol | value | definition site |
|---|---|---|---|
| thermal energy at 310 K | `kBT` | 4.280011900e-3 pN·µm | `aleph/laws/units.py:45` |
| actin bending modulus | `κ` | 7.276020230e-2 pN·µm² | `aleph/laws/units.py:50` (`= kBT·L_p`, `L_p = 17 µm`, Gittes 1993) |
| cortex segment rest length (config) | `seg_um` | 0.5 µm | `aleph/components/incumbent/assemble.py:228` — its own comment says **"CONVENIENCE (coarse): 5–10× above the sourced 50–100 nm mesh"** |
| realized mean segment length | `seg_mean` (= `ell`) | 0.4999048934 µm | `aleph/components/incumbent/assemble.py:741`, back-computed from the committed `inner_tolerance_um` |
| steric diameter | `σ_EV` | 0.007 µm | `aleph/components/incumbent/assemble.py:158` |
| WCA contact stiffness | `k_EV` | 1.0e3 pN/µm | `aleph/components/incumbent/assemble.py:159` — **"PROVISIONAL"** |
| WCA operating force cap | `f_cap` | 1.0e3 pN | `aleph/components/incumbent/assemble.py:329` |
| NMII crossbridge stiffness | `k_xb` | 1.0e3 pN/µm | `aleph/components/incumbent/assemble.py:168` (PI-GAP, "MASTER knob") |
| ERM tether stiffness | `k_erm` | 4.6e3 pN/µm | `aleph/components/incumbent/compartments.py:123` (PI-ratified 2026-07-21) |
| α-actinin crosslink | `k_xl` | 4.6e5 pN/µm | `aleph/laws/hand_kmc.py:132` (Ferrer 2008 PNAS) |
| filamin crosslink | `k_xl` | 8.2e5 pN/µm | `aleph/laws/hand_kmc.py:148` (Ferrer 2008 PNAS) |
| resting turgor | `Π₀` | 40 Pa = 40 pN/µm² | `aleph/components/incumbent/assemble.py:147` via `common/turgor_pi0.py:60` — **CONVENIENCE, HeLa proxy, PI-GATED** (see §6.4 caveat) |
| cell radius | `R` | 7.5 µm | `aleph/components/incumbent/assemble.py:124` |
| CFL numerator | — | 0.1 | `aleph/components/incumbent/assemble.py:1001`, `inner_mechanics.py:267` |
| float64 epsilon | `eps64` | 2.220446049250313e-16 | `np.finfo` |

Line numbers are as of this commit; `aleph/components/incumbent/assemble.py` is being edited concurrently, so the companion
script reads these constants **by symbol** (with an ordered fallback chain) rather than by line, and halts if a
symbol disappears instead of defaulting.

**Units convention** — the project's engine units: length µm, force pN, energy pN·µm, stiffness pN/µm,
pressure pN/µm² (`1 Pa == 1 pN/µm²`). Every derived threshold in this document is dimensionally verified
through an explicit exponent-tuple algebra in the companion script (§8, check `dimensions`).

---

## 2. Two different gates have been conflated

There are **two** convergence criteria in this codebase and only one of them executes.

**(i) The executing criterion** — `ac/cell/inner_mechanics.update_convergence_kernel` (line 296), a
three-part device-side predicate, all three parts in **µm**:

```
max_displacement      <= tolerance_um
max_constraint_error  <= tolerance_um
max|PF| * dt_mu       <= tolerance_um        <-- the force test, expressed as a displacement
```

with `tolerance_um = sqrt(eps64) * ell` (`driver.py:506`, `preload_force_probe.py:341`,
`ng2_ng3_ng6.py:1253`). Its own docstring is explicit and, importantly, **correct** about what it is:

> *"The tolerance is `sqrt(float64 epsilon) * ell_ref` … This is a derived numerical-resolution criterion,
> grid-invariant and not tuned to make a biological gate pass."*

That claim is true **of a displacement tolerance**: `sqrt(eps)·ell` is the smallest position change that is
meaningful in float64 when positions are differenced against a segment of length `ell`. It is a round-off
criterion and it is honest about being one.

**(ii) The reported criterion** — the bare literal `0.21`, in **54 sites across 22 files** (measured, not
quoted: `aleph/**/*.py` excluding `archive/`, minus 2 unrelated false positives). It appears in
`scripts/ac_resting_converge.py:1,40,57`, `scripts/ac_gate_a_figure.py:5,18`,
`scripts/ac_erm_schwarz_probe.py:33` (`GATE = 0.21`), `scripts/ac_myosin_stacking_probe.py:60,61,80,81`,
`aleph/components/incumbent/erm_schwarz.py:6,31`, `aleph/components/incumbent/erm_gauss_seidel.py:7`, `aleph/components/incumbent/assemble.py:242,245`,
`tests/ac/cell/test_resting_balance_oracle.py:15`, and 12 more `scripts/ac_gate_a_*.py` files. **It has no
definition site.** No module exports it, no YAML carries it, no `ValidationGate` row anchors it.

Criterion (ii) is a *reporting contract*, and §3 shows it is a hand-frozen snapshot of criterion (i).

---

## 3. Route (c) — the accepted-step predicate's own force floor. This is where 0.21 came from.

### 3.1 Derivation

The force limb of the executing predicate is a displacement test in disguise. Solve it for the force:

```
dt_mu * max|PF| <= tolerance_um
=>  max|PF| <= tolerance_um / dt_mu
             = [ sqrt(eps64) * ell ] / [ 0.1 / kmax ]
             = 10 * sqrt(eps64) * ell * kmax                    ... (C)
```

**Dimensions:** `[1] × [µm] × [pN/µm] = [pN]` ✓. And `dt_mu = 0.1/kmax` carries `[µm/pN]`, so
`dt_mu × [pN] = [µm]` is a length — the predicate really is comparing lengths, as its docstring says.

### 3.2 What sets `kmax`

`assemble.py:942-964` takes the max over every composed stiffness. Reconstructing all seven terms at the
native config:

| term | value [pN/µm] |
|---|---|
| actin bending `κ/ell³` | 0.582 |
| NMII crossbridge `k_xb` | 1.0e3 |
| WCA contact stiffness `U''(r_min)` | 1.0e3 |
| membrane ERM `k_erm` | 4.6e3 |
| α-actinin crosslink | 4.6e5 |
| filamin crosslink | 8.2e5 |
| **WCA near-core `U''(r_eq(f_cap))`** | **2 828 099.2886817832** ← dominant |

The winner is the **steric near-core stiffness at the operating force cap** —
`wca_analytic.max_core_stiffness(f_cap = 1e3 pN, σ = 0.007 µm, ε = 8.574462700673442e-4 pN·µm)`, where `ε`
itself is back-derived from `K_EV_PROVISIONAL = 1e3 pN/µm`. Reproduced to `rel = 0.0e+00` against the
committed `ledger.kmax_pN_per_um`.

### 3.3 The arithmetic, and the match

```
ell = 0.4999048934 µm   (realized seg_mean)      kmax = 2828099.2886817832 pN/µm
tolerance_um = 1.4901161193847656e-08 * 0.4999048934 = 7.449163398500275e-09 µm   [matches artifact exactly]
dt_mu        = 0.1 / 2828099.2886817832           = 3.535943748517097e-08 µm/pN   [matches artifact exactly]
F_pred       = 7.449163398500275e-09 / 3.535943748517097e-08 = 0.2106697371 pN
```

**`0.21` vs `0.2106697` → 0.34 % low.** At the config default `ell = 0.5 µm` (rather than the realized
`seg_mean`) the formula gives `0.2107098` — the 1.9e-4 relative gap is the beads-per-filament rounding, 0.095 nm
per segment. The literal is route (c) at the baseline, to well inside a rounding.

**This is the provenance.** It is not thermal, not a displacement-physics choice, and not fitted to any
result. It is float64 round-off × baseline mesh × a provisional steric guard.

### 3.4 What the gate therefore secretly depends on

| dependency | nature | if it changed |
|---|---|---|
| `sqrt(eps64)` | machine precision | float32 → **4882 pN** (×23 170 looser) |
| `ell = seg_mean` | mesh | 0.2 µm → 0.0843; 0.075 µm → 0.0316 (§5) |
| `kmax` = WCA core at `f_cap` | **provisional numerical guard** | steric off → **0.0611** (×3.4 tighter); `f_cap` or `K_EV_PROVISIONAL` moves it proportionally-ish |
| `0.1` CFL numerator | duplicated bare literal (host `assemble.py:964` + device `inner_mechanics.py:267`) | linear |

A "physics gate" that moves by four orders of magnitude when you change floating-point precision, and by 3.4×
when you toggle excluded volume, is a numerical-convergence statement wearing a physics label. **That is the
core finding of this document**, and it does not depend on which alternative route one prefers.

One further compounding: the `ell` the literal was frozen at is itself now labelled **"CONVENIENCE (coarse):
5–10× above the sourced 50–100 nm mesh"** at its own definition site (`assemble.py:228`). So the gate is
`round-off × a convenience mesh × a provisional steric guard` — three separately-flagged non-physical inputs,
multiplied.

---

## 4. Route (a) — thermal noise floor

### 4.1 (a1) Langevin `sqrt(2·γ·kBT/dt)` — **NOT EVALUABLE, and would be wrong if it were**

Two independent reasons, both from the code:

1. **There is no per-node drag `γ` in `ac/`.** `driver.py:7` states it explicitly: *"the physical clock comes
   from the process, NOT a drag-scaled descent step (the retired 6πηR clock)"*, and `aleph/components/fluid/__init__.py:14,28`
   records the `6πηR → gamma_solid` **retirement**. The resting solve is a quasi-static projected descent
   `x += dt_mu·P·F`; `dt_mu` has units **µm/pN**, not seconds. There is no mobility and no physical time step
   to substitute. The one viscosity that does exist, `mu_pore = 1.0e-3 Pa·s`
   (`ac/fluid/params_i0b1.yaml:61`), is the **Darcy pore-fluid** viscosity, and its own policy field forbids
   this use: *"MUST stay separate from BOTH the RETIRED 65.9 Pa·s bulk drag AND the gamma_solid inner-solve
   mobility."* Synthesising a per-node `γ = 6πηa` from it would be inventing a coupling the model does not
   have — refused per the *never invent a number* rule.
2. **Even with a `γ`, the expression is not a threshold.** `sqrt(2γkBT/dt) → ∞` as `dt → 0`: it characterises
   the *integrator's injected noise*, not the *state*. A quasi-static equilibrium criterion cannot depend on a
   step size that is being driven to zero.

**Route (a1) is rejected on both grounds.** This is a PI-GAP only in the trivial sense that a per-node drag
would have to be created; it should not be created *for this purpose*.

### 4.2 (a2) The equipartition force floor `sqrt(kBT·k)` — the correct thermal statement

For a degree of freedom held by stiffness `k` at 310 K, equipartition gives `⟨x²⟩ = kBT/k`, hence an
instantaneous net force fluctuating with rms

```
F_th = sqrt(kBT * k)          [pN·µm · pN/µm]^(1/2) = pN  ✓          ... (A2)
x_th = sqrt(kBT / k)          the companion position uncertainty
```

Both are **γ-free and dt-free** — the right form for a quasi-static state. Evaluated at the model's own
stiffnesses (`ell = 0.5 µm`):

| mode | `k` [pN/µm] | `F_th` [pN] | `x_th` [nm] |
|---|---|---|---|
| actin bending soft mode `κ/ell³` | 0.582 | **0.0499** | 85.7 |
| NMII crossbridge `k_xb` | 1.0e3 | **2.069** | 2.07 |
| membrane ERM tether `k_erm` | 4.6e3 | **4.437** | 0.965 |
| α-actinin crosslink | 4.6e5 | **44.37** | 0.0965 |
| filamin crosslink | 8.2e5 | **59.24** | 0.0723 |
| steric WCA core at `f_cap` | 2.83e6 | **110.0** | 0.0389 |

**Independent limiting-case verification.** `x_th = sqrt(kBT/k)` at `k = 0.1 N/m = 1e5 pN/µm` gives
**0.2069 nm**; `_historical/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md:120` computed **0.207 nm** for exactly
this quantity, independently, seven weeks earlier (and used it to catch a pre-existing "0.65 nm" kT error in
`erm.py`). The machinery reproduces a number this project derived by another hand. ✓

**What this says about the incumbent gate — plainly, and it is uncomfortable:**

- A gate on `max_i |PF_i|` at **0.21 pN** is **below the thermal force floor of every stiff bond in the
  model** — 9.8× below a bound crossbridge, 21× below an ERM tether, 210× below an α-actinin crosslink.
- The seeded `f_head = 1.5 pN` is **0.73×** `sqrt(kBT·k_xb)`. The physiological load the resting cortex is
  asked to balance is *smaller than the thermal force noise of the spring that applies it*.
- Therefore the 0.5–3.6 pN plateaus (`0.51–0.63` GATE-B native queue; `~3.6` solver consult; `0.776` ERM
  floor) all sit **inside** the thermal band of the bonds that carry them: 0.63 pN is `0.30 × sqrt(kBT·k_xb)`.

**Route (a2) is LOOSER than 0.21 — by ×9.8 at a crossbridge node and ×210 at a crosslinked node.**
🚩 **It must NOT be adopted as a pass criterion.** Doing so would loosen the gate by two orders of magnitude,
and that is a PI-only contract change (CLAUDE.md, *no gate-loosening*). Its correct use is §7 G3: a **ceiling
on what may be claimed** — a residual below `sqrt(kBT·k_local)` may not be reported as "physically
non-equilibrium", only as "the athermal minimiser has not converged".

### 4.3 `kBT / segment-length` — arithmetically available, physically vacuous

The task asked to consider it. `kBT/ell = 4.28e-3 / 0.5 = 8.56e-3 pN` (and `1.585 pN` at 2.7 nm, `0.611 pN` at
σ_EV = 7 nm). It is dimensionally a force, but it is not the thermal force *on anything*: it is the force
scale of a thermal energy dropped across an arbitrary length, and it inherits whatever length you pick, over a
2-decade range, with no mechanical justification. `sqrt(kBT·k)` is the same idea done correctly, with the
length supplied by the mechanics (`x_th`) rather than by choice. **Rejected as arbitrary.**

---

## 5. Route (b) — physically-meaningless-displacement, and why it collapses into (a2)

The form is `F = k · δ` (dimensions `pN/µm × µm = pN` ✓), using the **actual** per-node stiffness rather than
a guess. The trouble is `δ`:

| `δ` | bending 0.582 | `k_xb` 1e3 | `k_erm` 4.6e3 | α-actinin 4.6e5 | WCA core 2.83e6 |
|---|---|---|---|---|---|
| actin monomer rise **2.7 nm** | 0.00157 | 2.70 | 12.4 | 1242 | 7636 |
| steric **σ_EV = 7 nm** | 0.00408 | 7.00 | 32.2 | 3220 | 19 800 |

A single fixed `δ` spans **more than five decades** of implied threshold across stiffnesses that all coexist
on the same mesh. Route (b) with a fixed `δ` is therefore **ill-posed for a multi-stiffness network**: it is
absurdly tight on the bending mode and absurdly loose on the crosslinks.

`δ` must be **mode-matched**, and the only mode-matched length the physics supplies is the node's own thermal
position uncertainty `x_th = sqrt(kBT/k)` — below which a displacement is not observable even in principle.
Substituting it:

```
F = k * sqrt(kBT/k) = sqrt(kBT*k)          == route (a2)
```

**Routes (a2) and (b) are the same criterion.** That coincidence is a good sign for both, and it means there
is no independent third answer hiding in (b). The "energy" phrasing of the same statement — *the residual must
not be able to release more than ~kBT of work as it relaxes*, `|PF|²/(2k) ≤ kBT` — differs only by a factor
`sqrt(2)` and is the form to quote in a gate contract, because the *criterion* (residual stored energy per DOF
< kBT) is grid-invariant even though the per-node force it implies is not.

---

## 6. Route (d) — grid invariance. This is the part that matters most.

### 6.1 How route (c) scales

`kmax` is set by the WCA near-core term, a function of `(f_cap, σ_EV, ε)` **only** — no `ell`, no `n_nodes`
(`wca_analytic.max_core_stiffness` docstring asserts exactly this grid-invariance, and it holds). So

```
F_pred(ell) = 10 * sqrt(eps64) * ell * kmax     =>   F_pred ∝ ell^(+1),   independent of n_nodes
```

Verified numerically: `F_pred/ell` is constant across `ell ∈ {0.5, 0.2, 0.075}` to `1.3e-16`
(check `route_c_scales_linearly_in_seg`). While steric dominates `kmax` the scaling is *exactly* linear; if
`f_cap`/`K_EV_PROVISIONAL` were lowered enough for a crosslink to take over, `kmax` would become
`ell`-independent in a different way but still `n_nodes`-independent, and if bending ever dominated
(`κ/ell³`, which needs `ell < 4.46 nm` to beat filamin 8.2e5 and `ell < 2.95 nm` to beat the WCA core — i.e.
below the actin monomer, so **unreachable in practice**) the scaling would flip to `ell^(-2)`. All three
regimes are stated so a future refinement does not have to rediscover them.

### 6.2 The incumbent literal under refinement

| `ell` | runtime's own bound `F_pred` [pN] | frozen `0.21` is | reported `max\|PF\|` | verdict vs the runtime's own bound |
|---|---|---|---|---|
| 0.500 µm | **0.2107** | 1.00× | 0.14 | PASS |
| 0.200 µm | **0.0843** | **2.49× looser** | 0.21 | **FAIL by 2.5×** |
| 0.075 µm | **0.0316** | **6.64× looser** | 3.42 | **FAIL by 108×** |

(reported values: `FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md:429`.)

**The frozen literal loosened the gate by inaction.** Two of the three rungs in the fine-mesh series — including
the 200 nm rung the multigrid design calls *"the known-good rung (weak solver already reaches 0.21)"*
(`FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md:293`) — are **re-opened by arithmetic**, using nothing but the
runtime's own criterion evaluated at the correct `ell`. No new physics, no new run.

### 6.3 A per-node force gate is the worst possible metric for a refinement program

A nodal force is a force **density** times a length: `PF_i ≈ f(s_i) · ell`. Different families therefore scale
with opposite exponents at fixed physical state:

| residual source | nodal force scaling | why |
|---|---|---|
| extensive load (turgor, tension, motor load) | `∝ ell^(+1)` | density × node spacing |
| bending against a **fixed geometric** imperfection δ | `∝ ell^(-3)` | `κ·δ/ell³` |
| bending against a **resolution-limited** imperfection (`δ ∝ ell²`) | `∝ ell^(-1)` | `κ·ell²/ell³` |

Measured exponents from the committed series (companion script prints these):

```
seg 0.075 -> 0.200 :  max|PF| ∝ ell^(-2.845)      route-(c) bound ∝ ell^(+1.000)
seg 0.200 -> 0.500 :  max|PF| ∝ ell^(-0.443)      => gap widens as ell^(-3.845) / ell^(-1.443)
```

The `-2.845` exponent over the fine rungs is close to the pure `ell^(-3)` **fixed-kink bending** signature,
and it corroborates the fine-mesh doc's own node-level forensics — *"worst node 2169521: bending = 5.00 pN (of
the 3.42 pN |PF|)"* (`FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md:497`). The 3.42 plateau is a **geometric kink
being resolved**, not "more accurate physics"; and the measurement diverges from the gate as `ell^(-3.8)`.
**Chasing a per-node force literal under mesh refinement is arithmetically hopeless.**

(Caveat, stated because it weakens the series and the whole-repo audit already flagged it: the 0.075 µm rung
also changed `cortex_density_per_fil` 20→40, so the three rungs are not a clean one-parameter refinement.
That does not affect the *gate* arithmetic here, which uses only `ell`.)

### 6.4 The grid-invariant alternatives

Grid-invariance requires a **dimensionless** ratio, or a **total** whose continuum limit exists.

**(d1) — signed radial residual bias on the reported observable. RECOMMENDED.**
Only the *coherent radial* part of a residual field can contaminate the measured cortical tension; randomly
signed local residuals cancel. Define

```
S       = Σ_i ( PF_i · r̂_i )                                    [pN]   (signed, not |·|)
ΔP_spur = S / (4πR²)                                             [pN/µm²]
γ_spur  = ΔP_spur · R / 2 = S / (8πR)                            [pN/µm]
gate    : γ_spur  ≤  f · γ_target ,   γ_target = Π₀R/2 = 150 pN/µm
```

Dimensions: `pN/(µm) → pN/µm` ✓, and dividing by `γ_target` gives a dimensionless ratio containing **neither
`n_nodes` nor `ell`**. Budgets at `R = 7.5 µm`, `Π₀ = 40 Pa`:

| `f` (fraction of `γ_target`) | `|S|` budget | as a fraction of `F_load = Π₀πR² = 7068.6 pN` |
|---|---|---|
| 0.1 % | 28.3 pN | 0.004 |
| **1 %** | **282.7 pN** | 0.04 |
| 5 % | 1413.7 pN | 0.2 |

⚠️ **Caveat that travels with G2's denominator.** `γ_target = Π₀R/2` inherits Π₀, which as of 2026-07-25 is
provenance-gated in `aleph/components/incumbent/assemble.py:135-147` and classified **CONVENIENCE — Fischer-Friedrich 2014 HeLa
interphase, wrong cell line**. Three values are live in the tree (40 Pa here, 133 Pa band-implied ⇒ circular,
~72 Pa MCF7-geometry estimate) — a **3.3× spread**. So a *relative* G2 gate (`γ_spur/γ_target ≤ f`) is robust
to that spread (both sides scale with Π₀ and it cancels), but the **absolute `|S|` budgets in the table below
move by 3.3×** if Π₀ is re-anchored. State G2 in the relative form; treat the pN budgets as derived-at-40-Pa.

This is the criterion that answers *the question the gate exists to answer*: "would starting from this state
bias the tension I am about to measure?" **It cannot be evaluated retrospectively** — `S` is not computed in
any committed artifact. Adding it is a signed reduction alongside the existing per-family `max` telemetry in
`driver.py`'s force-family block (the machinery already loops the families and already has `r̂` available from
the membrane/nucleus radial masks).

**(d2) — L1 total `Σ|PF_i|` against the load. NOT recommended.** It converges (`Σ|f_i|ell → ∫|f|ds`) so it is
formally grid-invariant, but it **double-counts cancelling residuals** and so is not a measure of measurement
contamination. Indicatively: `mean|PF| × N` is `0.027 × 494 802 ≈ 1.34e4 pN` at the coarse mesh and
`0.129 × 2 898 126 ≈ 3.74e5 pN` at the fine mesh — 1.9× and 53× the *entire* physiological load, which is
mostly cancellation, not error. Reported only to show why L1 is the wrong norm.

**(d3) — grid-convergence diagnostic, not a gate.** Because `PF_i` is a density × `ell`, a *converging*
sequence must show `mean|PF| ∝ ell^(+1)`, i.e. `mean|PF|/ell` constant. Indicatively the two available means
give `0.027/0.5 = 0.054` vs `0.129/0.075 = 1.72` — a **32× divergence**, i.e. the fine mesh is far from grid
convergence. ⚠️ The two means come from *different states* (coarse = seeded t0; fine = post-multigrid
plateau), so this number is indicative only and must be re-measured on comparable states before it is quoted.

---

## 7. Recommendation

### The one change that is unambiguously a tightening — recommend adopting

**G1 — stop freezing route (c); compute it.** Replace the 54 hand-copied `0.21` literals with the run's own

```
F_pred = tolerance_um / dt_mu          (both already in every report: `inner_tolerance_um`, `inner_dt_mu`)
```

and report `max|PF| / F_pred` as a dimensionless ratio. Effects: at `ell = 0.5 µm` it is a **no-op**
(0.2107 vs 0.21 — the literal is 0.3 % *tighter*, negligible); at every finer mesh it **tightens** the gate to
the runtime's own criterion; and it makes the gate self-consistent under any future change to `ell`, `f_cap`,
`K_EV_PROVISIONAL`, or the CFL numerator. This is still a gate-contract change and so PI-owned, but it
**cannot loosen anything** and it eliminates the magic number outright. The two numbers needed are already in
every committed report JSON.

Re-scored against G1, for the record (all at `ell = 0.5 µm`, `F_pred = 0.2107`):

| result | `max\|PF\|` | ratio to G1 | reading |
|---|---|---|---|
| clean cortex, no myosin (`consult §3`) | 0.0356 | 0.17 | PASS with 5.9× margin |
| GATE-A static closure (`GATE_A_RESTING_CONVERGENCE:300`) | 0.1398 | 0.66 | PASS with 1.5× margin |
| GATE-A seed-2 (same doc) | 0.1586 | 0.75 | PASS |
| physiological-mesh rung (audit §6) | 0.28 | 1.33 | **FAIL** |
| GATE-B native queue (`NATIVE_QUEUE_RESULTS:6`) | 0.51 – 0.63 | 2.42 – 2.99 | **FAIL** |
| solver-consult plateau | ~3.6 | ~17 | **FAIL** |

### The genuinely new physics gate — PI decision required

**G2 — route (d1), the signed radial bias.** Grid-invariant, tied to the observable, and it is the only
candidate here that would still be the right gate at 2.9 M nodes. **The form is derived; the fraction `f` is
not** — `f` is a reporting-precision contract and only PI can set it. Three real options:

| option | `f` | `\|S\|` budget | trade-off |
|---|---|---|---|
| **A** conservative | 0.1 % | 28.3 pN | Below any plausible measurement precision; may be unreachable for the same reason the per-node gate is — untested. |
| **B** measurement-matched | **1 %** | **282.7 pN** | 1 % of `γ_target` = 1.5 pN/µm, far inside the 10–30 % spread of the cortical-tension literature the project compares against. My recommendation. |
| **C** permissive | 5 % | 1413.7 pN | Comfortably reachable, but 7.5 pN/µm of tension bias — larger than the entire `γ_network = 3.70 pN/µm` currently reported by GATE B, so it would make the flagship number meaningless. Not recommended. |

**Prerequisite:** `S` is not computed anywhere. G2 cannot be ratified before the reduction exists and has been
run natively at least once, so PI should ratify the *form* now and the *fraction* after the first measurement.

### The interpretive floor — must NOT become a pass criterion

**G3 — route (a2) as a claim ceiling.** Report, per run, `max_i |PF_i| / sqrt(kBT·k_i)` using each node's own
local stiffness. Rule: a residual with ratio **< 1** may not be described as a physical non-equilibrium — only
as an unconverged athermal minimisation. 🚩 G3 is **looser than 0.21 by ×9.8 to ×210** and is being surfaced
to PI precisely so it cannot be quietly swapped in as a pass criterion. It changes no gate; it changes what
sentences are allowed in a closeout.

### Also recommended: fix the framing, not just the number

The resting solve is an **athermal (T = 0) energy minimisation**. That is a legitimate and probably correct
modelling choice — but it means the resting-baseline claim is *"we located the mechanical equilibrium about
which thermal fluctuations occur"*, **not** *"the cell is at rest to within 0.21 pN"*. The second sentence is
not supportable at 310 K by §4.2 arithmetic, at any threshold. Every GATE-A closeout should carry the first
sentence.

---

## 8. Verification — dimensional analysis and independent limiting cases

Run: `python aleph/scripts/ac_gate_threshold_derivation.py` → **9/9 checks pass, exit 0**.

| check | what it verifies |
|---|---|
| `dimensions` | all four route formulas reduce to `pN^1` through an explicit (µm, pN, s) exponent algebra — not a comment, a computation that raises on a non-integral `sqrt` exponent |
| `dt_mu_units` | `0.1/kmax → µm¹·pN⁻¹`, so the force limb of the predicate really is a length comparison |
| `reproduce_native_kmax` | derived `kmax` == committed `ledger.kmax_pN_per_um`, `rel = 0.0e+00` (tol 1e-12) |
| `reproduce_native_predicate_force_exact` | route (c) at the realized `ell` == committed `tolerance/dt_mu`, `rel = 0.0e+00` (tol 1e-12) |
| `reproduce_native_predicate_force_from_config` | route (c) from the *config default* `ell` matches to `1.9e-4` (tol 1e-3); the gap is quantified as the 0.095 nm beads-per-filament rounding |
| `frozen_literal_is_route_c` | the literal read out of `ac_resting_converge.py` matches route (c) to 0.34 % (tol 1 %) |
| **`limiting_case_sigma_radial`** | **independent limiting case:** `sqrt(kBT/k)` at `k = 0.1 N/m` → 0.2069 nm vs the 0.207 nm computed independently in `_historical/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md:120` |
| **`limiting_case_steric_off`** | **independent limiting case:** steric disabled → `kmax` = filamin 8.2e5 → route (c) = 0.0611 pN vs the ~0.06 pN stated independently in `AUDIT_WHOLE_REPO_2026-07-25.md` finding 9 |
| `route_c_scales_linearly_in_seg` | `F_pred/ell` constant to `1.3e-16` across the three rungs — the `ell^(+1)` claim of §6.1 |

### Negative control — the checks were watched failing

```
$ python aleph/scripts/ac_gate_threshold_derivation.py --assume-kmax 8.2e5
!! NEGATIVE CONTROL: a deliberately wrong input was injected (kmax=820000.0, seg=None) — checks MUST fail.
   [FAIL] reproduce_native_kmax: derived kmax=820000 (ASSUMED (negative control)) vs committed ledger
          2828099.289  rel=7.101e-01 (tol 1e-12)
   [FAIL] reproduce_native_predicate_force_exact: ... 0.06108313987 pN vs committed ... 0.2106697371 pN
          rel=7.101e-01 (tol 1e-12)
   [FAIL] reproduce_native_predicate_force_from_config: ... rel=7.100e-01 (tol 1e-3 ...)
   [FAIL] frozen_literal_is_route_c: frozen literal 0.21 pN vs derived route-(c) 0.0610948 pN -> rel=70.9073%
RESULT: 5/9 checks pass  -> 4 FAILED           exit 1
```

`8.2e5` is not a random wrong number: it is the filamin crosslink stiffness that
`AUDIT_WHOLE_REPO_2026-07-25.md` finding 9 assumed when it estimated the implied bound at "~0.06 pN". That
estimate is **arithmetically right for a steric-free build and wrong for the native one** — the native `kmax`
is the WCA near-core term, 3.45× larger, which is exactly why the implied bound *is* 0.21 rather than being
"tighter than 0.21". The negative control therefore also corrects a finding in the audit that commissioned it.

A second negative control (`--assume-seg 0.075`, exit 1, 3 FAILED) confirms the `ell` limb is live too.

---

## 9. What this derivation does NOT settle

1. **It does not tell you whether the resting state is physically correct.** Route (c) certifies a minimiser
   converged; route (d1) would certify the measurement is uncontaminated. Neither certifies that
   `f_head = 1.5 pN`, `fraction = 0.5`, `k_xb = 1e3` are the right physiology — all three remain PI-GAP.
2. **G2 has never been measured.** `S = Σ PF·r̂` does not exist in any artifact. The fraction cannot be
   ratified from data that has not been taken.
3. **The `ell^(-2.845)` exponent rests on a 3-point series with a confounded variable** (density 20→40 at the
   finest rung). The direction is robust; the exponent is indicative.
4. **`d3`'s 32× grid-divergence compares two different states** (seeded t0 vs post-multigrid plateau). Must be
   re-measured on comparable states.
5. **The route (a2) table uses a single stiffness per mode**, not the true per-node Hessian diagonal after the
   NF2007 inextensibility projection. The projector removes the stiff *stretch* directions, so the true
   per-node `F_th` is somewhere between the bending value (0.05 pN) and the crosslink value (44 pN) depending
   on local topology. Computing it properly means a per-node Hessian-diagonal reduction — worth doing, and it
   is the same reduction G3 needs. **The qualitative conclusion is unaffected:** even the *softest* candidate
   stiffness a myosin-loaded node has is `k_xb = 1e3`, giving 2.07 pN ≫ 0.21.
6. **It does not address the `0.1` CFL numerator**, a bare literal duplicated in host (`assemble.py:964`) and
   device (`inner_mechanics.py:267`) code, which multiplies straight into the gate. Small, separate, real.
7. **G2's absolute budgets inherit the Π₀ CONVENIENCE claim** (40 Pa HeLa proxy; 3.3× spread across the three
   live values). The relative form is Π₀-invariant; the pN numbers are not. Stated in §6.4.
8. **Nothing here was run natively.** No CUDA on this host; everything is arithmetic over committed artifacts
   and constants, which is the correct scope for a threshold derivation — but the G1 re-scoring table inherits
   whatever provenance the source documents have (the audit's Tier-(d) warnings apply to `~3.6` and `0.28`).

---

## 10. Surface to PI

| # | item | why it is PI's, not mine |
|---|---|---|
| 1 | **G1** — replace the 54 frozen `0.21` literals with the per-run `tolerance_um/dt_mu`. Tightens at every mesh < 0.5 µm, no-op at baseline. | gate-contract change |
| 2 | **The 200 nm and 75 nm rungs are re-opened by arithmetic** (2.5× and 108× over the runtime's own bound). The 200 nm rung is currently documented as "known-good". | a closed result becomes open |
| 3 | **G2** — ratify the *form* of the signed-radial-bias gate now, the *fraction* (A/B/C in §7) after the first native measurement of `S`. | new gate contract |
| 4 | **G3 is LOOSER than 0.21 by ×9.8–×210.** Surfaced explicitly so it is never swapped in as a pass criterion. Recommended use: a ceiling on claims. | no-gate-loosening rule |
| 5 | **`f_head = 1.5 pN` is 0.73× the thermal force noise of its own crossbridge** (`sqrt(kBT·k_xb) = 2.07 pN`). Either `k_xb` (PI-GAP "MASTER knob") is too stiff, or a *single* resting crossbridge is not a meaningful mechanical unit and the resting prestress must be represented as a coarse-grained mean field. This is the same modelling question `project-ac-gate-a-p02-p03-native` already escalated, now with a number attached. | modelling decision |
| 6 | **The gate is proportional to `K_EV_PROVISIONAL` and `steric_force_cap`**, both provisional numerical guards. Any Magic-Number-Block closure on those two silently moves the convergence gate. | parameter closure has an undeclared side effect |
| 7 | **`AUDIT_WHOLE_REPO_2026-07-25.md` finding 9 needs a correction**: the implied bound at the native config is 0.21, not "~0.06 pN, tighter than 0.21" — that figure is the steric-off build. | audit correction |
| 8 | **Reframe every GATE-A closeout** from "at rest to within X pN" to "athermal mechanical equilibrium located; thermal force noise at the loaded bonds is Y pN". | wording of a scientific claim |

---

## Figures

None. This deliverable is arithmetic over committed artifacts; the companion script's stdout is the figure, and
it is reproduced verbatim in §8. The one thing that *would* deserve a figure — `max|PF|` vs `ell` overlaid on
the `F_pred ∝ ell` line and the `sqrt(kBT·k(ell))` thermal band — needs the missing comparable-state
measurements from §9.4 first, and drawing it from the confounded 3-point series would be exactly the kind of
visual that manufactures confidence.
