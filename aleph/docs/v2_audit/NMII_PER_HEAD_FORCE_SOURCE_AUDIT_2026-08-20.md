# The 0.5–2.0 pN bracket is not two measurements. Neither end is one.

**Status:** SOURCE AUDIT. Decides nothing, selects no value. Written 2026-08-20 by session `d380041d`
at PI request — *"한 번만 다시 자료 조사해보자"* — before banding `F_stall_head` as a test axis.

Literature retrieved from **PubMed** during this audit; every article is cited with its DOI below.

---

## 1. What the repository currently claims

`aleph/components/motor/params_i0b3.yaml`, `F_stall_head` — `value: null`, `evidence_status: "GAP — PI"`:

| | value | recorded source |
|---|---:|---|
| `claim_a` | **0.5 pN** | *"AFINES (Freedman 2017 PLoS Comput Biol); Kovacs 2003; Tam 2021; ff/hand_kmc NMIIA preset"* |
| `claim_b` | **2.0 pN** | *"Billington structural interpretation; ff/myosin_linear.F_HEAD_PN=2.0 (KB-3.18 ~2 pN)"* |

The file already flags the conflict: *"0.5 pN attribution not yet an audited SourceEvidence row; 2.0 pN /
Billington also not closed for head-level use."* This audit closes that flag in the negative.

---

## 2. What the sources actually are

### 2.1 Kovács/Wang 2003 measured no force

**Wang F, Kovács M, Hu A, Limouze J, Harvey EV, Sellers JR (2003)**, *"Kinetic mechanism of non-muscle
myosin IIB: functional adaptations for tension generation and maintenance"*, J Biol Chem 278(30):27439-48
— [10.1074/jbc.M302510200](https://doi.org/10.1074/jbc.M302510200).

It is a **transient-kinetics** study of an NMIIB subfragment-1 construct: ATPase cycle rates, ADP affinity,
ADP release rate. The reported finding is that *"non-muscle myosin IIB subfragment-1 spends a significantly
higher proportion of its kinetic cycle strongly attached to actin than do the muscle myosins."*

⚠ **There is no force measurement in it.** A solution-kinetics paper cannot be the provenance of a per-head
isometric stall force in pN. The citation is real; it does not support the quantity it is attached to.

### 2.2 AFINES is a simulation parameterisation, not an assay

`claim_a`'s primary is AFINES (Freedman et al. 2017), and the note carried beside it in our own ledger says
so: *"small per-head stall; 'small relative to 3-4 pN cardiac' (AFINES_ALGORITHM_NOTES L516)"*. That is a
modelling choice with a comparative justification, not a measured number. Using it as a literature anchor
makes another model's parameter into our evidence.

### 2.3 `claim_b`'s 2.0 pN traces to a kernel this project already retired

Our own ledger records it: *"higher per-head stall; **used by the RETIRED aggregate linear kernel**"*, and
`ff/myosin_linear.F_HEAD_PN = 2.0`. A constant inside a retired lumped kernel is not a head-level
measurement, and "Billington structural interpretation" is an inference from filament architecture rather
than a force assay.

### 2.4 The field said the measurement had not been made

**Melli L, Billington N, Sun SA, Bird JE, Nagy A, Friedman TB, Takagi Y, Sellers JR (2018)**, *"Bipolar
filaments of human nonmuscle myosin 2-A and 2-B have distinct motile and mechanical properties"*, eLife 7
— [10.7554/eLife.32871](https://doi.org/10.7554/eLife.32871). Full text retrieved; **it contains the string
"pN" zero times**, and states directly:

> *"The effect of force on the kinetics of NM2 paralogs **has not been measured using optical trapping**."*

That is the Sellers lab, in 2018, on their own subject.

### 2.5 The one optical-trap study on NMIIB reports displacements, not a stall force

**Norstrom MF, Smithback PA, Rock RS (2010)**, *"Unconventional processive mechanics of non-muscle myosin
IIB"*, J Biol Chem 285(34):26326-34 — [10.1074/jbc.M110.123851](https://doi.org/10.1074/jbc.M110.123851).
A three-bead optical trap on NMIIB. Its reported quantities are a **forward displacement of 5.4 nm** and
**backward steps of −5.9 nm**, with forward steps and detachment *"weakly force-dependent at all forces"*
and the conclusion that *"NMIIB can readily move in both directions at stall."*

⚠ **Scope note, and it cuts against a single stall value:** a motor that steps backwards as readily as
forwards at stall does not have "the" isometric per-head force in the sense our `f_head` uses it. Full text
is not available through PubMed Central for this article, so **no pN value from it is quoted here.**

---

## 3. What this audit concludes

**The bracket is not a disagreement between two measurements. It is the gap between a simulation parameter
and a retired kernel's constant.** One more literature pass does not produce a better number, because the
number does not exist in the form the engine asks for.

⚠ **Two corroborations found along the way**, recorded because they are stronger than what they replace:

* `N_side` — Melli 2018's own text has *"30 motors per filament end"*, and Tripathi/Bond/Sellers/Billington/
  Takagi (2021) J Vis Exp — [10.3791/62180](https://doi.org/10.3791/62180) — describes NM2-B as *"a bipolar
  filament containing approximately 30 molecules."* Both support `claim_b` (28–30) over `claim_a` (10,
  AFINES). **Not a decision** — recorded so the PI decides on evidence rather than on symmetry.
* **Duty ratio is a measured, load-dependent quantity.** Kovács M, Thirumurugan K, Knight PJ, Sellers JR
  (2007) PNAS 104(24):9994-9 — [10.1073/pnas.0701181104](https://doi.org/10.1073/pnas.0701181104) — reports
  that ADP release is *"increased 4-fold by the assisting load on one head and decreased 5-fold (for 2A) or
  12-fold (for 2B) by the resisting load on the other"*, and that *"resistive load increases duty ratio to
  favor tension maintenance by two-headed attachment."* Melli 2018 adds that NM2-B has *"a four-fold higher
  duty ratio than NM2-A."*

  ⚠ **This bears directly on the resting bound fraction.** The engine's own position (`hand.py:187`, `:230`)
  is that the duty ratio EMERGES from `k_on`/`k_off`, never imposed — and the measured literature is about
  exactly those load-dependent rates, not about a fraction. So the resting fraction is the quantity for
  which sourced kinetics DO exist, one layer down. Recorded; not acted on, because the PI selected the
  direct-band option for now.

---

## 4. What follows for the ask

The PI's instruction — *"우리는 어떤 것에 대한 기준을 잡지 않고 생리 범위만을 일단 테스트용으로 잡는다"* —
is not a shortcut past a sourced value here. **It is the only honest option available**, because this audit
finds no measured per-head isometric stall force for non-muscle myosin II to be had.

What the `source` field should therefore say is what is true: a test axis over 0.5–2.0 pN, with the
provenance of each end named as what it is, and the run record carrying that sentence rather than a
citation that does not support the number.

## 5. What this does NOT do

It selects no value, narrows no bracket, and does not retire `claim_a` or `claim_b` — both stay recorded,
which is what makes the conflict auditable. It does not read Norstrom 2010's full text (unavailable) and
quotes no pN from it. It does not touch `N_side`, whose corroboration above is evidence for a PI decision,
not the decision.
