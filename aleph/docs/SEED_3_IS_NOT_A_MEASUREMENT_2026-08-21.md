# seed 3 is not a measurement, and the 1.30 ratio was a property of the cut

**2026-08-21, session `d380041d`. A correction to my own commit `f2a…` of forty minutes earlier.**

When seed 3 landed I reported two things that are wrong:

> *"Two of three seeds have windows above the bar that pass every gate."*
> *"The first genuine seed-to-seed tau comparison exists … cut 22,000: seed 1 tau 115.8, seed 3 tau
> 89.2, ratio 1.30."*

Session 21 refuted both. **Reproduced here independently before accepting**, at cuts the earlier table
did not carry.

## τ is not stable across cuts for seed 3

```
   cut    s1 tau   s3 tau    ratio
21,000      97.7     87.9     1.11
22,000     115.8     89.2     1.30     <-- the one I quoted
23,000     103.2     64.7     1.59
24,000     108.1     51.3     2.11
25,000     112.3     51.5     2.18
26,000     111.8     49.0     2.28

spread over cuts    seed 1  16.7%      seed 3  61.2%
```

**Seed 1 converges. Seed 3 falls monotonically by a factor of 1.8 over 5,000 samples.** A τ that
depends this strongly on where the window starts is not that window's τ; it is a description of a
correlation structure still changing along the series.

## ⚠ The ratio moves 2.1× with the cut, so it is not a property of the seeds

`1.11 … 2.28`. **A property of two seeds cannot depend on where the observer chooses to cut.** What I
reported as the first seed-to-seed τ comparison was a reading of one arbitrary column of that table.

## ⚠ And seed 3's n_eff passing is NOT independent evidence — this is the error under the error

```
n_eff = N / (2 tau)

s3 cut 22,000   tau 89.2  ->  n_eff 44.8
s3 cut 26,000   tau 49.0  ->  n_eff 40.8      (N fell from 8,000 to 4,000 and n_eff barely moved)
```

**A falling τ mechanically raises n_eff.** I treated `n_eff ≥ 25` as an independent gate confirming
that seed 3 had settled, when n_eff is algebraically a function of the very quantity in question. Seed
3's high n_eff is a **consequence of its unconverged τ**, not evidence against it.

⚠ **This is the sharpest form of tonight's family and I committed it after cataloguing nineteen of
them.** I spent the night writing that a check must not be entangled with its own subject — the
tolerance that is exact about the wrong thing, the assert that ratifies the value its docstring
forbids, the query that returns zero for every input. **`n_eff` is entangled with `tau` by definition,
and I used one to certify the other.**

## What actually stands

**One seed of three produced a τ. There is no τ scatter.**

⚠ **And the end-state column may not be used either.** Seeds 1/2/3 end at 540.33 / 544.41 / 540.04 —
tantalisingly close — but those are γ magnitudes, the set is `REFUSED`, and computing a spread from
them **is** the seed-scatter measurement, which may only be done on a set that passed. `(c) 3` and
`(c) 17` stand.

**What the three runs establish is that settling behaviour differs between seeds at fixed
configuration.** That is not an input to D-2 (the scatter of a settled observable); it is an input to
**D-3, initial-condition independence** — *"without this, stationary is compatible with stuck where it
started."* Which of the two it belongs to is the PI's.

## The practice note, which is session 21's and worth keeping

They observed that sending the numbers **without interpretation, at three cuts**, is what made the
refutation possible: had I sent only cut 22,000 — the one that supported my reading — **1.30 would have
entered the record as a seed property.** The habit that caught this was not skepticism; it was handing
over more columns than my conclusion needed.
