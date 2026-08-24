# Cortex filament count (70,686) — independent cross-check + PI issues

**Date:** 2026-07-21 · **Trigger:** PI question — "how was the cortex filament count derived, is it
MCF7-based, and can it be found another way?" · **Status:** PI-issue (documentation only; no code/count
change this session, per PI 2026-07-21).

## 1. What 70,686 actually is

It is **not a count** — it is a product of two factors:

```
70,686 = ρ_areal × 4πR²  =  100 µm⁻²  ×  4π(7.5 µm)²  =  100 × 706.86
```
Source lines: [aleph/laws/cortex_assembly.py:48](../../ff/cortex_assembly.py#L48),
[aleph/components/incumbent/assemble.py:144](../../ac/cell/assemble.py#L144).

- **ρ_areal = 100 µm⁻²** — a **GENERIC mammalian cortex** areal density (KB-3.1 / KB-3.18, citing
  Salbreux2012_TCB, ChughPaluch2018_JCS, Murrell2015_NRMCB, Charras2008_BJ). **No MCF7 measurement.**
- **R = 7.5 µm** — MCF7 suspended radius (Wagner2011, see §3). MCF7-specific.

So "70,686 cortical filaments" is an **MCF7-sized geometric estimate using a generic cortex density**, not a
measured MCF7 cortical filament count. No paper reports this number for MCF7.

## 2. How the count is ACTUALLY determined — literature + Cytosim (2026-07-21 correction)

**First, an important correction to an earlier draft of this doc.** A prior version anchored the cross-check on
a "100 nm actin mesh" and concluded 70,686 was ~9–15× too dense. **That was wrong** — 100 nm is the *perturbed*
cortex mesh, not the native one (see Bovellan below). Corrected below.

### 2a. Nobody counts cortical filaments directly

The real cortex is too dense for microscopy to resolve individual filament ends. Fritzsche/Baranov super-res
(Jurkat, [PMC5390943](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5390943/)) measured cortex thickness
230 (+105/−125) nm and membrane–cortex distance (bimodal 50 / 120 nm) but state the **high actin density
precludes filament-end identification** → no filament length / precise density. **There is no measured MCF7
cortical filament census.** Every "count" (ours included) is a density/concentration/mesh estimate.

### 2b. Cytosim (the reference engine) SETS the number explicitly

Belmonte, Leptin & Nédélec 2017 ([PMC5615920](https://pmc.ncbi.nlm.nih.gov/articles/PMC5615920/)) simulate the
cortex as a 2D patch with the filament count **chosen directly**: ~1,500 filaments over a circular area of
radius 15 µm (≈707 µm²) ⇒ **~2 filaments/µm²** — a coarse representative-filament density picked to reproduce
network mechanics, NOT derived from concentration/geometry. So "100 µm⁻² × 4πR²" is *a* convention, not the
field standard; the areal density is a modeling choice, not a measured value.

### 2c. Measured native mesh (Bovellan et al. 2014, Curr Biol, [PMC4110400](https://pmc.ncbi.nlm.nih.gov/articles/PMC4110400/))

- **Control (native) cortex: homogeneous mesh, gap ~30 nm.**
- mDia1 depletion → high-density patches + 100–200 nm gaps; both nucleators perturbed → mesh **>100 nm**.
- ⇒ The KB "actin mesh ~100 nm" (KB-3.1/3.18) is the **perturbed / loose** value, not native. Native ≈ 30 nm.

### 2d. Count vs the CORRECT native mesh (reproducer: `scripts/cortex_count_crosscheck.py`)

Cortex shell A = 4πR² = 706.9 µm², V = A·h = 141.4 µm³ (h = 200 nm). Actin linear density = 741 monomers/µm.
Filament count via the Schmidt/Käs ξ–[c] relation (ξ[µm]=0.3/√c[mg/mL]) at the model L = 3 µm:

| mesh ξ | regime | implied [F-actin] | **N (L=3 µm)** |
|---|---|---|---|
| 100 nm | *perturbed* cortex (the value the earlier draft wrongly used) | ~214 µM | ~8,200 |
| 50 nm | intermediate | ~860 µM | ~32,800 |
| **30 nm** | **native cortex (Bovellan)** | ~2.4 mM | **~91,200** |

**70,686 sits in the native 30–50 nm regime — it is order-correct, if anything slightly conservative vs a
30 nm mesh.** The earlier "9–15× too many" claim is **withdrawn**; it came from anchoring on the perturbed mesh.

> **All "[F-actin]" values here are LOCAL, polymerized F-actin inside the 200 nm cortex shell** — computed
> purely from filament contour length (monomers-in-filaments ÷ 741 mon·µm⁻¹). The free G-actin monomer pool
> never enters the count. So the ~1.8–2.4 mM figures are NOT comparable to the cell-averaged total-actin
> (G+F) ~100–500 µM; a thin dense shell is *expected* to be mM-scale F-actin, and this is consistent with the
> 30 nm native mesh. (An earlier "concentration looks high" caveat was a total-vs-F-actin category error,
> retracted.) The count is a pure F-actin interpretation; nothing depends on the G-actin pool.

### 2e. The residual issue is filament LENGTH, not the count magnitude

Native cortical filaments are **short and high-turnover (~hundreds of nm)**, not the 3 µm the model uses. At a
fixed native mesh, N ∝ 1/L, so a faithful short-filament census would be **~10× larger** (912k at L=0.3 µm,
30 nm mesh). The model's 70,686 × 3 µm are therefore **coarse representative filaments** (Cytosim-style), which
is defensible — but it should be *labelled* as such, and per-filament kinetics/turnover tuned accordingly. This
is the real fidelity gap, more than the head-count.

## 3. Wagner2011 citation conflict — RESOLVED (real, was un-registered)

- **Wagner BA et al. 2011**, *Free Radic Biol Med* 51(3):700, **PMC3147247**. Coulter Counter, MCF7 suspended
  **dia 14.8 µm / vol 1.70 pL — MEASURED, open-access** → R ≈ 7.4 µm ≈ 7.5 µm.
- Used consistently: dcm/geometry.py, dcm/nondim.py, tests, docs/LAYER2_ANCHORS_2026-06-02.md,
  MCF7_CONVERSION_SPEC_2026-06-02.md. **Not fabricated.**
- The conflict: absent from the KB registry (`source_evidence`/`paper_refs`/`knowledge_claim` all empty on
  Wagner), so [aleph/components/incumbent/assemble.py](../../ac/cell/assemble.py) over-warned *"do NOT assert Wagner 2011"*. It
  was already flagged as an un-registered candidate in
  [SE_REGISTRATION_CANDIDATES_2026-06-30.md:131](../../references/SE_REGISTRATION_CANDIDATES_2026-06-30.md#L131).
- **PI 2026-07-21 decision:** register Wagner2011 (SourceEvidence + suspended-radius KnowledgeClaim) in Notion,
  refresh, and reconcile the assemble.py comment (done). R = 7.5 µm is thereby anchored, not soft.

## 4. Open PI issues (this doc = the surface, no count change made)

1. **The count MAGNITUDE is defensible** — 70,686 sits in the native 30–50 nm mesh regime (§2d). No downward
   revision is warranted; the earlier "reduce to ~5–8k" suggestion is retracted (it was anchored on the
   perturbed 100 nm mesh).
2. **Filament LENGTH is the real gap**: native cortical filaments are ~hundreds of nm, not 3 µm. The 70,686 are
   coarse *representative* filaments (Cytosim-style, ~2/µm² there vs 100/µm² here). Decide: (a) label + tune
   per-filament turnover/kinetics to match short-filament reality, or (b) move to a short-filament census
   (~10× more filaments) — a large native-run cost. PI call.
3. **KB mesh label fix**: KB-3.1/3.18 "actin mesh ~100 nm" is the *perturbed* value; native ≈ 30 nm
   (Bovellan2014). Register the native mesh and re-label the 100 nm as perturbed. (Register Bovellan2014 +
   Belmonte2017 as SourceEvidence.)
4. **Generic vs MCF7 density**: the 100 µm⁻² areal density has no MCF7 anchor (its magnitude happens to match a
   native mesh, so it is not disqualifying). If an MCF7-specific cortical F-actin density exists, prefer it.

## 5. DECISION REQUIRED (PI) — cortical filament length / count representation

The model uses **N = 70,686 filaments each L = 3 µm**. Native cortex = short (~hundreds of nm), high-turnover
filaments in a ~30 nm mesh. At a fixed native mesh, N ∝ 1/L, so the two are coupled. Three ways forward:

| Option | What | N (native 30 nm mesh) | GPU cost | Fidelity |
|---|---|---|---|---|
| **A — coarse representative** | keep 70,686 × 3 µm, LABEL as representative, tune per-filament turnover so aggregate mechanics match | ~70,686 (unchanged) | low | ⚠️ lumps length distribution → **conflicts with the CLAUDE.md hard rule** (fine-grained over lumped); length-dependent buckling/severing/turnover misrepresented |
| **B — short-filament census** | realistic L≈200–400 nm short filaments at 30 nm mesh | **~300k–900k cortical** | high (~5–10× population; likely > A5000 16 GB with full compartment stack → larger/multi-GPU) | ✅ true fine-grained; correct length-dependent mechanics + turnover |
| **C — length distribution** | sample cortical L from a cortex-specific short-mean exponential (KB-3.18 gives 1–10 µm *linear*; cortex needs a shorter mean) | intermediate | medium | ✅ distribution-faithful, cheaper than uniform-short |

**Recommendation: B (or C as a staged step toward B).** The CLAUDE.md architecture principle is explicit —
"pick the full-fidelity, fine-grained, mechanistic option over lumped even at higher cost … larger/multi-GPU is
the response to memory pressure; never lower biological density to fit memory." Option A is the lumped path and
violates that rule; it is acceptable only as an *interim, explicitly-labelled* placeholder, not production.
This is a native-run-cost decision → **PI sign-off before any population change.**
