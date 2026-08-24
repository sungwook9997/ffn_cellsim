# Run 1 — the criterion could not have shown descent. Recorded as VOID, not as a negative result.

**Slurm job 69, RTX 4090-1, 2026-08-20 14:54–15:06 KST.** Full native, 70,686 filaments / 494,802 actin
nodes / 511,114 total, 4,000 inner iterations, 10 declared points, all 10 attempted, all rc=0.
Code: `89a3f38a` + `c7651281` + `c47803c6` (`COMMIT_STAMP` on the run host; the host is not a git
checkout, so this is `source: declared`).

## The declared criterion was unfalsifiable, and that is this run's finding

Declared before the run in `ac_l0_balance_point_exists.py`:

    descent_ratio = residual_end / residual_start  < 1

`residual_end` is read **after the outer step's accept/reject**. Every point in this run was rejected
(`outer_accepted: false`, `outer_rolled_back: true`), and the rollback is bit-exact — so `residual_end`
was restored to `residual_start` in all ten arms, to every digit:

    residual_start = residual_end = 47.405172   in 10/10 points, ratio exactly 1.000000

**A criterion that reads a post-rollback field cannot show descent on a rejected step.** It is the mirror
of the defect this session spent the day documenting in `balance_ok` — a gate that cannot fail — and it
was written by the same session, one day later, into its own experiment.

⚠ **The summary's verdict string `NO_DESCENT_ANYWHERE` is therefore VOID and may not be quoted.** It is
what the criterion said; the criterion said nothing.

## The instrument was fine — that is why this is a criterion failure and not a run failure

| | |
|---|---|
| `inner_iters` | **4,000 of 4,000** — the solve ran its full budget |
| `n_actin` / `n_total` | **494,802 / 511,114** — full native, confirmed from the run's own ledger |
| `myosin_bound_after_step` | **88 / 442 / 1,768** at fraction 0.01 / 0.05 / 0.20 — the setpoint bound heads and scaled with the axis |
| positive control | job 68 aborted 5 min in when the instrument check read the wrong field; fixed in `c47803c6`. **The check did its job.** |

## Post-hoc observation, NOT the declared test, and not quotable as one

`residual_candidate` — what the solve reached before accept/reject — was in every artifact all along:

    control_OFF        cand/start = 1.00139957     0 bound heads
    all nine ON points  cand/start = 1.0012372x    88 / 442 / 1,768 bound heads

Two things are visible and **both need a re-declared run before either is a result**:

1. Nothing descends. Every arm rises; ON rises slightly less than OFF (0.124% vs 0.140%).
2. **A 20× change in bound heads (88 → 1,768) moves the ratio by nothing** — the nine ON points agree to
   the 7th decimal. The tension source is on and the residual does not respond to how much of it there is.

## The probable reason, and what run 2 must therefore vary

`residual_start` is **47.405172 at both 3,000 and 70,686 filaments** — identical. The residual maximum is
therefore not on the cortex; it is on the membrane carrying the turgor. And this run did **not** pass
`--preload-erm-balance`, so the ERM tethers were force-free at build.
`RESTING_BASELINE_DIAGNOSIS_2026-07-23c.md` §2 already says what that means: *a force-free ERM transmits
nothing*. Cortical tension with no transmission path cannot move a membrane residual, however much of it
there is — which is exactly the insensitivity observed.

So run 1 tested **the setpoint alone**, and the honest statement of its outcome is that the setpoint alone
does not couple to the residual. Run 2 must vary the setpoint **and** the ERM preload together.
