# H.7 SF — Route-B NMII force-scaling: the SF generation-limit is REAL (HALT→PI)

**Date:** 2026-06-09
**Session:** (i) "CLOSE THE FLOOR" (PLATFORM_PI_QUEUE SESSION SPLIT), worktree
`ffn_cellsim-nmii`, branch `h7/sf-nmii-forcescale`.
**Mandate:** find the native NMII content/density of a single ventral stress fiber +
per-minifilament stall force (deep-research, cross-verified), DERIVE the Route-B
mesoscale force-scale factor, apply to `sf_myosin_`, test the Kumar 10-30 nN single-SF
tension band — or HALT→PI if no usable native-NMII-per-SF datum exists.

> ⚠️ **Citation-integrity discipline** (3 confirmed hallucinated SourceEvidence rows on
> record). Every value below is tied to a real, retrievable paper with a measurement
> method and a SOLID / ORDER / SPARSE confidence label. No value is back-solved to make a
> gate pass — the HARD RULE for this session ("force-scaling factor는 밴드 통과를 위해
> 고르면 안 됨 — 밀도 datum에서 유도만").

---

## TL;DR (≤8 lines)

1. **OUTCOME = HALT→PI (REFUTE).** The Kumar 10-30 nN single-SF tension is **NOT closable**
   from a measured native-NMII density datum via Route B. The one multiplier Route B
   needs that I could NOT find is the **parallel minifilament count per SF cross-section**.
2. **The molecular content IS solid** (Anchor 2a): native minifilament ≈ 30 heads/side,
   ~300 nm, 28-30 molecules (Billington 2013 / Hu 2017 / Melli 2018 / Niederman 1975).
3. **The per-minifilament FORCE is ORDER/extrapolated** (Anchor 2b): only primary estimate
   **fs ≈ 17 pN** (Stachowiak 2009, muscle-myosin extrapolation, self-flagged as not a
   direct NMII measurement). → a **DERIVABLE per-minifilament molecular correction ≈ 3.4×**
   over the runtime's 5 pN minifilament (10 heads × 0.5 pN).
4. **The parallel cross-section count is MISSING** (Anchor 1): the only explicit literature
   estimate (~50/cross-section) was **REFUTED 0-3**; back-solving 10-30 nN ÷ 17 pN ⇒
   ~590-1760, **geometrically impossible** (cross-section 50-250 nm r, ~10-30 actin → O(5-15)).
5. **REFRAME (HIGH conf):** Kassianidou/Schwarz/**Kumar** 2017 (PNAS) — single-fiber
   **ACTIVE** myosin force ≈ **5-6 nN**; the 10-30 nN is mostly **network/prestress**.
   Even the refuted-generous 50 × 17 pN ≈ **0.85 nN** is ~7× under the ~6 nN active target
   and ~12-35× under Kumar.
6. **This is the SF instance of the cortical γ-floor** — same ½·n·f·ℓ generation budget,
   same MISSING motor-density datum, same insight that the literature band includes
   passive/prestress the active motors alone do not supply.
7. **PI decision needed:** accept generation-bound (→ session ii) vs. re-scope the gate to
   the ~6 nN single-fiber ACTIVE target vs. open the load-dependent-duty-ratio lever.

---

## 1. What the runtime models vs. what the literature says

| quantity | runtime (`phase1_h3.yaml` + `myosin.py`) | native (deep-research) | confidence |
|---|---|---|---|
| heads per side | **10** (brief literal) | **~30** (≈30 molecules, ~60 heads) | SOLID |
| per-head stall | **0.5 pN** | ~1.7 pN (muscle extrap.); NMII span 0.7-10 pN | ORDER |
| per-minifilament force | φ·10·0.5 ≈ **5 pN** | **fs ≈ 17 pN** (Stachowiak 2009) | ORDER |
| minifilament length | 700 nm rigid rod (mesoscale) | ~300 nm | SOLID |
| parallel / cross-section | 1 effective (series tension) | **NO DATUM** (~50 refuted 0-3) | MISSING |
| single-fiber ACTIVE force | — | **~5-6 nN** (Kassianidou 2017) | HIGH |
| single-fiber NETWORK total | (Kumar band target) | ~25 nN center; 10-30 nN composite | HIGH |

## 2. The decomposition (why the gate can't be closed honestly)

A 1-D contractile cable's tension at a cross-section = (parallel motor units crossing that
plane) × (force per unit). Minifilaments in **series** along the fiber share the SAME
tension (they do not add) — so the tension ceiling is set by the **parallel cross-section
count**, exactly the datum that does not exist.

The Route-B factor therefore factors as
`factor = factor_mini (molecular) × N_parallel (cross-section)`:

* **`factor_mini` — DERIVABLE, lit-anchored (≈ 3.4×).** Bring the model minifilament
  (φ·10·0.5 pN = 5 pN) up to the literature per-minifilament estimate (17 pN, Stachowiak
  2009). Grid-invariant, not band-chased. Encoded in `h7_basal_sf_force_budget.py`.
* **`N_parallel` — MISSING.** The ~50/cross-section estimate was REFUTED 0-3 (its
  ~100-actin ÷ 2 derivation did not survive). The geometric ceiling is O(5-15) (50-250 nm
  radius, ~10-30 actin filaments across, 2 actin/minifilament). The back-solve to Kumar
  (~590-1760) exceeds geometry by ~40-350×. **No grid-invariant count datum exists.**

**Bounding calculation (force-budget output, `production/h7_basal_sf_force_budget.json`):**

| parallel count | × 17 pN/minifilament | vs ~6 nN ACTIVE | vs Kumar 10 nN |
|---|---|---|---|
| geometric 5-15 | 0.085-0.255 nN | ~24-70× under | ~40-120× under |
| **refuted-generous 50** | **0.85 nN** | **~7× under** | **~12-35× under** |
| back-solve 590-1760 | 10-30 nN (by construction) | — | **geometrically impossible** |

Even the **most generous** (refuted) parallel count, with the molecular correction
applied, lands ~7× under the single-fiber **active** target — so the gap is **not** an
artifact of the runtime's small per-minifilament force; it is irreducible without the
missing count, and the count the band requires cannot physically fit the cross-section.

## 3. The reframe — active generation vs. network prestress (HIGH confidence)

Kassianidou, Brand, Schwarz & **Kumar** 2017 (PNAS 114:2622, U2OS, the Kumar lab's own
active-Kelvin-Voigt + active-cable-network analysis):

* single-fiber aggregate **motor stall force** `Fs = k·Lo ≈ 6 nN` (k≈3 nN/µm, Lo≈2 µm);
* a **connecting** SF adds only **~5 nN of ACTIVE myosin force**; the length-defined SF
  reaches **~25 nN at center** — the remainder is **network connectivity / prestress**.

Kumar 2006 itself reports **stress (Pa)**, not a clean per-SF nN tension; the "10-30 nN"
is a composite community figure (Deguchi 2006 microcantilever; per-adhesion). So the
Kumar band is **not a pure single-fiber active-generation target** — it bundles the
~5-6 nN active component with ~20 nN of network/prestress. In our platform the network/
prestress component is exactly what the **connected basal mesh + FA anchoring + whole-cell
turgor prestress** supply — NOT `sf_myosin_` generation alone.

## 4. Unification with the cortical γ-floor

This is the **same generation-limit, surfaced from the SF side**:

* Same ½·n·f·ℓ active-gel budget; same outcome that the active motors fall ~1-3 orders
  under the literature band.
* Same **missing motor-density datum**: cortically, no MCF7 minifilament areal density
  (Nie 2015 HeLa 0.6/µm² is the only proxy; `H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07`);
  for the SF, no parallel-cross-section count (refuted). The SF session was meant to
  *fill* the cortical density gap from the SF side — **it cannot**, because the SF count
  datum is itself missing.
* Same structural insight: the literature target (γ band; Kumar band) **includes passive/
  prestress contributions** the active motors alone do not generate. The active component
  is consistent with the model's order; the band is not a pure-active target.

See `gamma-floor-layered-resolution` memory + `H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09`
(main repo). Session (ii) "ACCEPT THE FLOOR" should fold this in as one bounded story.

## 5. What was done (auditable)

* **Deep-research** (104-agent harness, 5 angles, 22 sources, 25 claims 3-vote verified,
  18 confirmed / 7 killed). The ~50/cross-section parallel estimate and 3 weaker SF-tension
  claims were among the 7 KILLED — the refutation is the central result, not a side note.
* **`h7_basal_sf_force_budget.py` extended** to encode (a) the DERIVABLE molecular
  correction (≈3.4×, Anchor 2a/2b), (b) the MISSING parallel count (None-gated; refuted /
  geometric / back-solved bounds, never picked), (c) the Kassianidou active/network split.
  Verdict **REFUTE**, `halt_to_pi: true`. No `myosin.py` `mesoscale_force_scaling` edit —
  a full SF factor cannot be derived honestly (the parallel multiplier is missing), so
  wiring one would be a back-solve (forbidden).

## 6. SourceEvidence candidates (PI-gated — NOT auto-registered)

For a future PI-gated SE registration (citation-integrity verified here; method + ID):

| key | claim | method / cell | conf |
|---|---|---|---|
| Billington2013_JBC | NMII minifilament 29-30 molecules / 58-60 heads | EM, recombinant human | SOLID |
| Hu2017_NatCellBiol (PMID 28114270) | ~30 molecules, head doublet 300±20 nm | 3D-SIM + PREM, U2OS | SOLID |
| Melli2018_eLife (32871) | 30 molecules (NM2A/B), 16 (NM2C), ~300 nm | negative-stain EM, Sf9 | SOLID |
| Niederman1975_JCB (67:72) | 28 molecules / 56 heads | EM, platelet | SOLID |
| Stachowiak2009_BiophysJ (PMC2711311) | per-minifilament stall fs ≈ 17 pN (muscle extrap.) | model | ORDER |
| Kassianidou2017_PNAS (114:2622) | single-fiber active ~5-6 nN; network total ~25 nN | active-KV / cable-net model, U2OS | HIGH |
| Kumar2006_BiophysJ (PMID 16500961) | single-living-SF mechanics (reports stress Pa) | laser nanoscissor + TFM | SOLID(method) |

⚠️ The "~50 minifilaments/cross-section" value (Stachowiak-derived) is **REFUTED** —
do NOT register it as a count. The ~17 pN is an extrapolation, **not** a direct NMII
measurement — register with the ORDER label and the "muscle-myosin extrapolation" note.

## 7. HALT → PI — decision request

The mandate's "no usable native-NMII-per-SF datum exists → HALT to PI" branch is the
honest outcome, with a richer finding than "not found":

**(Q-SF-1) Gate scope.** Is the SF gate target the **single-fiber ACTIVE** myosin force
(~5-6 nN, Kassianidou) or the **network-total** 10-30 nN (Kumar/Deguchi)? The active
component is the only thing `sf_myosin_` generates; the rest is mesh + FA + turgor prestress
(emergent from the full build, not a myosin dial). **Recommend:** scope the `sf_myosin_`
active gate to ~5-6 nN and let the network/prestress total emerge from the integrated cell.

**(Q-SF-2) The residual ~7-70× even at the active target.** Closing even ~6 nN at 17 pN/
minifilament needs ~350 parallel (vs geometric ~5-15). The only unclosed physics that
could bridge this is the **load-dependent duty ratio** (engaged-head fraction rising from
~0.05 unloaded toward ~0.5 under resistive load — Kovacs 2007), which would raise the
per-minifilament force well above 17 pN. **Recommend:** treat this as the next mechanistic
lever (a per-head load-dependent engagement term), NOT a density factor — but only on PI
direction, since it touches the Hill/Bell-Evans head kinetics in `myosin.py` (frozen path).

**(Q-SF-3) Accept generation-bound now?** If PI prefers, the honest fallback is session
(ii): document SF + cortical γ as ONE bounded generation-limit and report with the bounds.
This result already supplies that conclusion.

**Default if no PI response:** keep `ventral_stress_fibers` LIVE-but-passive-active-bounded
(the build gates all PASS; only the Kumar active band is unmet), do NOT wire a back-solved
factor, and hand the unified generation-limit story to session (ii).
