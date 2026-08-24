r"""Cortical tension γ on the RESTING path, reported against SEED SCATTER.

**What this is for.** `STATE.md` (e) 1 — replace stage 1's acceptance criterion *convergence* with
*stationarity* — has been HELD since 2026-07-28. `PI_DECISION_e1_STATIONARITY_2026-08-16.md` §2.2 says
why, and it is not a wording problem: *"A stationarity gate cannot be written against an observable the
run does not produce."* γ is the observable the resting state most needs and the resting driver does not
emit it; γ exists only in the GATE-B slice drivers. §5 item 4 lists emitting it as the one prerequisite
that is engineering rather than decision. This module is that prerequisite and nothing more.

**⚠ NO MAGNITUDE IS CLAIMED HERE.** `STATE.md` (c) 3 retires *"γ = 3.72 pN/µm and every SF/cortex motor
tension"* and (c) 17 retires *"every cortex-motor γ at dt=0.01 — magnitudes AND ratios"*. This module
adds no γ to either list and removes nothing from them. It builds the PATH that computes γ where the
resting run can reach it; whether any number it produces is quotable is a separate PI decision that this
file must not pre-empt. Every record it emits carries :data:`NO_MAGNITUDE_CLAIMED` verbatim.

It also writes **no stationarity verdict about the cell**. :func:`observe_gamma_seed_scatter` reports
whether a γ *mean* is statistically supportable across replicates, which is the input to (e) 1's D-2, not
D-2's answer. Choosing the observable set (D-1), the bands, and the window is the PI's.

The one thing this module refuses to let the caller do
------------------------------------------------------
`STATE.md` (f) carries a standing measurement rule, and §3 D-2 of the decision document calls violating it
*"the single most likely way to write the amendment and still be measuring nothing"*:

    compare γ against the SEED SCATTER, not the within-run ``sem``.

Three replicates measured S = 0.0668 pN/µm = 1.95% of γ — **67× the within-run sem**. A stationarity test
built on within-run variance is wrong by that factor, in the direction that certifies noise as a result.
So :func:`observe_gamma_seed_scatter` takes a mapping **keyed by seed** and there is no single-run entry
point that returns an error bar. The within-run ``sem`` is still reported — labelled, beside
``scatter_inflation``, so the size of the discrepancy is visible in every record rather than only here.
That is the same device `aleph/observe/stationarity.py` uses for ``sem_naive``/``sem_inflation``, one
level up: naive → correlated → across-seed, each ratio recorded.

What is ported and what is new
-------------------------------
The γ arithmetic is **ported, not re-derived**. The method of planes is
:func:`aleph.laws.gamma_estimator.method_of_planes_gamma` — imported, not copied; it is the same
instrument the GATE-B drivers use. The two family extractions read the same physics as
``components/incumbent/assembled_cortical_stress.py`` (which `world/` may not import — see
`tests/architecture/test_layer_directions.py`; the arena may bind `laws/` and may not reach into the
trees it replaces). ``actin_axial_tension`` here is the vectorised form of that module's per-node Python
loop and is pinned equal to it, term for term, in ``tests/world/test_observe_gamma.py``. The loop is
O(n_actin) in Python — at 494,802 actin nodes that is a per-sample cost paid once per recorded step, and
a stationarity window needs hundreds of samples.

⚠ **Two incumbent callers disagree on the cut centre**, and it changes γ.
``measure_assembled_cortical_stress`` passes the FULL node array (so the centre is the mean over membrane
and nucleus too); ``ac_gate_b_cortex_motor_native._measure_gamma`` passes ``pos[:n_actin]``. This module
takes the actin-only centre — γ is a property of the cortex — and exposes ``centre`` for a caller who
must reproduce the other. Surfaced rather than silently chosen.

⚠ **``force_mask`` is required and has no default.** It must be the omit mask the inner solve ran with.
Measured 2026-07-29: with the mask applied to the solve and not to the measurement, the two arms' dynamics
were bit-identical while γ differed by a constant 0.41526 pN/µm — the estimator was reading a different
force field from the one being integrated. A default here would make that silent again.

⚠ **The step this γ is measured on is not JUDGED today.** `aleph/world/step.py` stands on
``UndefinedAcceptance``: the step runs, the verdict is ``ACCEPTANCE_UNDEFINED``, and ``is_accepted`` is
falsy for it on purpose. So "measured between accepted steps" is, right now, "measured between steps
whose acceptance no criterion has been applied to". Every record here carries
:data:`STEP_ACCEPTANCE_NOTE` saying so, because a γ series silently presented as sampled from accepted
physics is the promotion the charter's evidence-class rule exists to stop. That module's own self-check
already writes a criterion of the form *"accepted iff ``gamma_pn_per_um`` is in the observables"* — this
module is what makes such a criterion writable at all, and it is still the PI's to declare.

Sanity Gate (charter, before first execution)
---------------------------------------------
* **Dimensions** — ``f_ext`` [pN], ``pos`` [µm], ``k_xl`` [pN/µm], ``R_um`` [µm] ⇒ γ [pN/µm]. 1 pN/µm =
  1e-6 N/m; the literature bands in `laws/gamma_estimator.py` are hundreds of pN/µm, not fractions.
* **Boundary cases** — no filament of ≥ 2 nodes, no crosslink, zero-length segment, a constant series, a
  single replicate: each returns a stated refusal or an empty family, never a silent zero passed off as a
  measurement. Covered in :func:`_demo`.
* **Conservation** — the actin extraction is a running sum of the external force from a FREE filament end,
  so a filament whose net external force is ≈ 0 closes; ``net_force_vector_pn`` is carried as that
  diagnostic and is NOT a tension.
* **Sign sense** — ``T = −(Σ_{j≤k} f_ext[j])·û_k``: positive is tension, negative compression, matching
  the ported module and the estimator's convention.
* **Precision** — float64 throughout; the cumulative sum is over ~5e5 terms of like magnitude, which is
  where the loop and the vectorised form can differ in the last bits. The pinning test asserts a relative
  tolerance, not bit equality, and says so.
* **Measurement protocol** — host post-processing of a synchronised device read, BETWEEN accepted steps
  (I0-A). Nothing here runs in the hot loop and nothing here mutates state.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.gamma_estimator import fibonacci_plane_normals, method_of_planes_gamma
from aleph.observe.stationarity import (
    DEFAULT_DRIFT_SIGMA,
    DEFAULT_MIN_TAU_WINDOWS,
    DEFAULT_SOKAL_C,
    DEFAULT_TRANSIENT_SIGMA,
    assess_stationarity,
)

__all__ = [
    "NO_MAGNITUDE_CLAIMED",
    "STEP_ACCEPTANCE_NOTE",
    "DEFAULT_N_PLANES",
    "DEFAULT_MIN_REPLICATES",
    "DEFAULT_RESOLVE_K",
    "RestingGamma",
    "SeedScatterReport",
    "actin_axial_tension",
    "plane_force_sums",
    "segment_offsets",
    "gamma_planes_device",
    "crosslink_tension",
    "gamma_from_resting_readback",
    "observe_gamma_seed_scatter",
]

#: Stamped into every record this module emits. `STATE.md` (c) 3 and (c) 17 stand; emitting γ on the
#: resting path does not lift either, and a record that did not say so would be quoted as if it had.
NO_MAGNITUDE_CLAIMED: str = (
    "no magnitude claimed — this record establishes that gamma is COMPUTED on the resting path and how "
    "it scatters across seeds. STATE.md (c) 3 and (c) 17 retire the gamma magnitudes and ratios measured "
    "so far and are NOT lifted by this measurement. Quoting any number here as a physical cortical "
    "tension requires a separate PI decision."
)

#: Stamped beside :data:`NO_MAGNITUDE_CLAIMED`. Deliberately a STATEMENT plus a pointer rather than a
#: value read out of `aleph/world/step.py` at call time: that module is another lane's and is moving,
#: and a stale copy of somebody else's verdict is worse than a sentence telling the reader where the
#: verdict actually lives.
STEP_ACCEPTANCE_NOTE: str = (
    "the step this was sampled on is NOT judged: aleph/world/step.py stands on UndefinedAcceptance, so "
    "the verdict is ACCEPTANCE_UNDEFINED and is_accepted is False for it. Read the run's own verdicts "
    "there; do not read this series as sampled from accepted physical steps."
)

#: Diametral cut orientations averaged by the method of planes. 64 is what the GATE-B cortex-motor gate
#: uses; kept identical so a resting γ and a GATE-B γ differ by the cell state and not by the instrument.
DEFAULT_N_PLANES: int = 64

#: Elements per chunk in the plane projection — large enough for BLAS to matter, small enough to bound
#: scratch. A pure performance knob: `plane_force_sums` returns the same numbers for any chunk size, and
#: `tests/world/test_observe_gamma.py` asserts exactly that.
#:
#: ⚠ This said "one chunk's (elements x planes) intermediate is ~64 MB". That is one array, and about
#: SIX are live at once (`side_a`, `side_b`, the raw and abs-ed `u @ normals.T`, `crossing`, and the
#: product being summed). Measured with `tracemalloc` at 400,000 elements: **371 MB peak transient**,
#: 5.8x the figure that stood here. Corrected rather than quietly dropped.
#:
#: The same measurement explains why an external minor-fault counter cannot count steps, which cost the
#: Lead an hour and nearly a 1.3 h run: at native this path allocates and frees ~3.9 GB per family per
#: call — roughly 7x the 576 MB readback — so the allocator traffic is dominated by THESE temporaries,
#: and their fault count depends on chunking and allocator reuse rather than on step boundaries. I told
#: the Lead the faults were the readback. They are mostly this, and this is my code.
_PLANE_CHUNK: int = 131_072

#: Replicates below which no scatter is reported. Two seeds give a std with one degree of freedom, which
#: is a number rather than a measurement; the standing rule in `STATE.md` (f) was measured on three.
DEFAULT_MIN_REPLICATES: int = 3

#: Coverage factor turning a scatter into a resolution floor. Declared here, before any run, and copied
#: into every report — it is a contract, not a knob to be picked once the spread is known.
DEFAULT_RESOLVE_K: float = 2.0

_ZERO_LENGTH_UM: float = 1e-12


@dataclass(frozen=True, slots=True)
class RestingGamma:
    """γ for ONE resting configuration, with the provenance needed to re-judge it.

    Attributes:
        gamma_total_pn_per_um: method-of-planes γ over every load-bearing family [pN/µm].
        gamma_network_pn_per_um: γ over the passive network families (actin constraint + crosslink).
            ``None`` on an empty group, for the same reason as ``gamma_source_pn_per_um``.
        gamma_source_pn_per_um: γ over the active families (crossbridge). ⚠ **``None`` when the group
            has NO ELEMENTS, which is not the same as measuring zero.** The PHASE 4 cell the 2026-08-21
            τ runs were measured on stands 11 populations and never calls `build_nmii` — there is no
            motor in it — so a reported ``0.0`` there would say "the motors contributed nothing" when the
            truth is "there are no motors", and a reader would have to know the population list by heart
            to tell. ``None`` cannot be summed, averaged or plotted as a zero; it breaks loudly instead.
            Read it beside ``n_source_elements``.
        gamma_by_family_pn_per_um: per-family γ, so a total can be attributed rather than trusted.
        net_force_vector_pn: Σ of the per-node force vectors [pN]. A Newton-closure DIAGNOSTIC, ≈ 0 for
            internal forces. **Not a tension** — it is on this record so the closure can be checked, and
            reading it as γ is the legacy error `components/incumbent/cortical_tension.py` warns about.
        n_planes / R_um / centre_um: the instrument, recorded because γ depends on all three.
        force_mask: the omit mask the inner solve ran with, as supplied by the caller.
        n_actin / n_crosslinks / n_source_elements: the population the number was measured on.
        magnitude_claim: :data:`NO_MAGNITUDE_CLAIMED`.
    """

    gamma_total_pn_per_um: float
    gamma_network_pn_per_um: float | None
    gamma_source_pn_per_um: float | None
    gamma_by_family_pn_per_um: dict[str, float]
    net_force_vector_pn: tuple[float, float, float]
    n_planes: int
    R_um: float
    centre_um: tuple[float, float, float]
    force_mask: tuple[str, ...]
    n_actin: int
    n_crosslinks: int
    n_source_elements: int
    magnitude_claim: str = NO_MAGNITUDE_CLAIMED

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form, for a run record."""
        return {
            "gamma_total_pn_per_um": self.gamma_total_pn_per_um,
            "gamma_network_pn_per_um": self.gamma_network_pn_per_um,
            "gamma_source_pn_per_um": self.gamma_source_pn_per_um,
            "gamma_by_family_pn_per_um": dict(self.gamma_by_family_pn_per_um),
            "net_force_vector_pn": list(self.net_force_vector_pn),
            "n_planes": self.n_planes,
            "R_um": self.R_um,
            "centre_um": list(self.centre_um),
            "force_mask": list(self.force_mask),
            "n_actin": self.n_actin,
            "n_crosslinks": self.n_crosslinks,
            "n_source_elements": self.n_source_elements,
            "magnitude_claim": self.magnitude_claim,
            "step_acceptance": STEP_ACCEPTANCE_NOTE,
        }


@dataclass(frozen=True, slots=True)
class SeedScatterReport:
    """γ across REPLICATE SEEDS — the error bar `STATE.md` (f) says is the only defensible one.

    Attributes:
        observable: which series was judged, for the record.
        status: ``"MEASURED"`` when a cross-seed mean and scatter were computed, ``"REFUSED"`` otherwise.
            **Neither value is a claim that the cell is at rest, and neither makes a γ quotable.**
        n_replicates / seeds: how many independent seeds, and which.
        mean_pn_per_um: mean over replicates of each replicate's steady-state mean. ``None`` if refused.
        seed_scatter_pn_per_um: **S** — the sample standard deviation (``ddof=1``) of those replicate
            means. This is the scale the mean must be judged against.
        sem_across_seeds_pn_per_um: ``S / sqrt(n_replicates)``.
        mean_within_run_sem_pn_per_um: the average of the per-replicate correlation-aware ``sem``.
            Reported so it can be COMPARED, never so it can be quoted as the error bar.
        scatter_inflation: ``S / mean_within_run_sem``. The measured version of the standing 67×. One
            would mean the two agree; anything far above it means a within-run error bar is a fiction.
        seed_scatter_fraction: ``S / |mean|`` — the standing rule's 1.95%, re-measured per run.
        resolvable_fraction_one_seed: ``resolve_k * S / |mean|``. Below this, one seed cannot tell two γ
            apart — the standing "~5%" limit, recomputed rather than inherited.
        resolvable_fraction_n_seeds: ``resolve_k * sem_across_seeds / |mean|``, the same for this ensemble.
        window_s: **W** — the SHORTEST analysed window across replicates. The window a joint statement can
            be made over is the worst one, not the mean one.
        equilibration_s: the LATEST transient end across replicates: no replicate may be judged before it.
        max_tau_int_s / min_n_effective: the worst correlation time and the worst effective sample count.
        all_sokal_windows_closed: ``False`` if any replicate's Sokal window never closed, in which case its
            ``tau_int`` is a LOWER bound and its ``n_eff`` an upper one.
        required_window_samples: how long the analysed window must be, IN SAMPLES, for the declared
            length bar to be cleared — ``min_tau_windows * tau_samples``, worst replicate. On a
            ``TOO_SHORT`` refusal this is the answer to "how much longer", in the only unit that is
            sourced when the run's time axis is not. ⚠ It is computed from the CURRENT ``tau`` estimate,
            so if ``min_tau_reliability`` is low it is itself a lower bound: a longer run can measure a
            longer ``tau`` and move this target further out. That is not a defect in the number, it is
            what an unconverged correlation time means.
        max_tau_relative_error: the worst Madras-Sokal relative uncertainty on ``tau_int`` itself,
            ``sqrt(2*(2M+1)/N)`` with ``M`` the Sokal window lag. It is what turns *"tau grew between
            runs"* into *"tau grew by more than its own error bar"*, which is the only form of that
            observation that decides anything. Validated against the realised scatter of ``tau`` over 200
            AR(1) realisations: predicted 0.217 / 0.349 against realised 0.171 / 0.257 at rho = 0.9 / 0.98
            — conservative by ~30%, i.e. erring toward refusing to call a change significant.
        min_tau_reliability: the worst ``T / tau_int`` across replicates — how many correlation times the
            series spans. Below roughly 50 the ``tau_int`` estimate is materially biased LOW, so
            ``n_effective`` is an OVER-estimate and any error bar built on it is a lower bound. Surfaced
            at the top level because it decides whether ``min_n_effective`` may be read at face value.
        windows_open_on_transient: seeds whose RETAINED window still begins on a rise. ⚠⚠ **SINCE PI
            RULING 18 (b), 2026-08-22, EMPTY IS THE NORMAL CASE AND CARRIES ALMOST NO INFORMATION.**
            `equilibration_start` now SKIPS contaminated candidates instead of labelling the winner
            afterwards, so this field is non-empty only when NO candidate was clean and the search fell
            back — `fell_back_to_contaminated` per replicate. When that happens the flag means what it
            always meant, and it is the stronger signal it ever was. Otherwise read
            `windows_skipped_on_transient` beside it, which is where the warning moved.
            ⚠ **EMPTY is not
            evidence of settledness on its own** — the detector loses power as the window shrinks, so an
            unflagged short window means "could not see a rise", not "there was none". Read
            `min_n_effective` beside it: below the length bar a window is not a measurement, and two such
            windows agreeing that they see no transient is agreement about nothing. Measured 2026-08-21
            on contract III: at a 4,000-sample cut BOTH seeds read unflagged (1.75 and 1.43 sigma) with
            n_eff 17.9 and 14.5 against a bar of 25, while at the one cut where seed 1 IS a measurement
            (n_eff 54.4) seed 2 fails this test, the drift test and the length gate together. ⚠ **When this is
            non-empty every tau-derived field above describes the residual transient, not a correlation
            time, and must not be read** — `max_tau_samples`, `min_tau_reliability`,
            `max_tau_relative_error`, `min_n_effective` and `required_window_samples` alike. Measured
            2026-08-21 on contract III seed 1: the report carried `max_tau_samples` 1,280.2 and
            `min_tau_reliability` 12.01 from a window the equilibration search itself flagged at 9.28
            sigma, while tau on the clean part of the same series is ~102 — a factor 12.6. Those numbers
            were read at face value because nothing at this level said they had been invalidated, which
            is the defect this field closes.
        windows_skipped_on_transient: ``(seed, n_candidate_windows_rejected)`` for every replicate whose
            search had to walk around a rise. ⚠ **This exists because the fix for one blind spot opened
            another.** PI ruling 18 (b) filters contaminated candidates out of the selection, which is
            right — before it, widening the search bound alone let a seed take the highest reliability of
            three on a window that opened on a transient. But filtering also means the contamination
            stops being announced: `windows_open_on_transient` goes empty and the report reads clean.
            ⚠ **Non-empty is NOT the signal; the COUNT is.** A settled series also skips a few, measured
            at (12, 9, 3) against 28 each for a contaminated series of the same length. The reason is
            structural rather than noise-dependent, and it is the peer session's: **an early cut opens
            near the start of ANY series, and the start is the easiest place in a series to find a
            rise** — so the first handful of candidates are the ones most likely to trip the detector
            whatever the shape underneath. A floor of a few skips is therefore the expected reading, not
            a warning. A LARGE count says the series
            had a long dirty region, and it is the one field that distinguishes "settled early" from
            "the search skipped 36 candidates to find somewhere it could stand". Measured 2026-08-22 on contract III: seed 3 skipped 36 and took the HIGHEST
            reliability of the three replicates, at a cut where tau is still falling (89.2 -> 51.7).
            ⚠ **The filter does not catch that seed, and neither does this field on its own** — both see
            *where the window opened*, not *whether tau has settled*. What this field buys is that the
            "2 of 3 seeds gave a tau, 105.4 against 51.7" reading is visibly not a seed scatter.
            ⚠ **Read `window_samples` beside it against `stationarity.MIN_RETAINED_SAMPLES`.** The filter
            pushes the window later, so on a series it has to walk a long way through, the retained
            length can come down near that floor — and there the floor stops being a safety rail and
            starts deciding the answer. Measured 2026-08-22 on the opposed-drift fixture: at floors 512 /
            256 / 128 the chosen start for one seed moved 5,480 -> 5,434 -> 5,548 and for another
            5,206 -> 5,720 -> 5,694, i.e. **it moves and it is not even monotone**. No field was added for
            this, deliberately: `window_samples` already carries the length, and on the native gamma
            series the optimum retains 6,000-11,000 samples, nowhere near the floor. If a native run ever
            reports a window within a few hundred of it, this note is the reason to distrust the cut.
        max_tau_samples: the worst ``tau_int`` as a COUNT of samples, and ``equilibration_samples`` the
            transient's end likewise. ⚠ These exist because the ``*_s`` fields beside them are named for
            seconds and hold ``samples * dt_s``: a caller with no physical time axis passes ``dt_s = 1.0``
            — which is the correct usage — and then ``max_tau_int_s: 188`` is a correct value of the
            wrong quantity, indistinguishable from 188 seconds to anyone reading the record later. Every
            entry in the intended read order now has a count twin, so the reading never depends on
            knowing what ``dt_s`` was.
        equilibration_samples: see ``max_tau_samples``.
        window_samples: **W in samples** — the shortest analysed window as a COUNT. Always sourced. The
            seconds above are ``samples * dt_s`` and are only as sourced as the caller's ``dt_s``, which
            for a run whose mobility is unsourced is not sourced at all.
        drift_per_s_by_seed: each replicate's **SIGNED** least-squares slope, keyed by seed. ⚠ A RATE:
            it scales as ``1 / dt_s``, the OPPOSITE direction to ``window_s``, ``max_tau_int_s`` and
            every other ``_s`` field on this report. At ``dt_s = 1.0`` it is per SAMPLE, and a reader who
            has learned that the ``_s`` fields hold counts will correct this one the wrong way — the same
            suffix, the reciprocal trap. Measured: ``dt_s`` 1.0 → 0.05 moves the others ×0.05 and this
            ×20. Only the MAGNITUDE is affected; the sign (:attr:`drift_directions_agree`) and the
            significance (:attr:`drift_sigma_by_seed`) are ``dt_s``-invariant, and those are the two the
            verdict rests on.
            `aleph/observe/stationarity.py` keeps this signed on purpose — *"a quantity falling toward
            its steady state and one rising toward it are different physics, and ``abs`` hides which
            happened from the person who has to decide what to do next"* — and a summary that dropped it
            applied that ``abs`` one level up.
        drift_sigma_by_seed: each replicate's ``drift_sigma_measured`` — the correlation-aware
            significance the ``DRIFTING`` verdict is actually thresholded on, keyed by seed. The slope
            beside it carries direction and physical size; this carries whether it is distinguishable
            from zero, which is the number that decides. Per-replicate only until now, and on 2026-08-21
            that meant a reader comparing the transient detector against a plain linear fit had no way to
            see that the module's own drift test was already answering the same question — 8.75 and 6.52
            on seed 2's longer windows, 0.11-0.97 on seed 1's.
        drift_directions_agree: whether every replicate's slope shares one sign. ``None`` when no slope
            was measurable. ⚠ **This is an input to D-3, not an answer to it.** Seeds that drift the SAME
            way are one system still equilibrating; seeds that drift OPPOSITE ways are seeds going to
            different places, which is `PI_DECISION_e1_STATIONARITY_2026-08-16.md` §3 D-3 —
            initial-condition independence, *"without this, stationary is compatible with stuck where it
            started"*. Both read as ``REFUSED, DRIFTING`` and nothing else here told them apart. Whether
            a disagreement means the resting state is conditional is the PI's call and is NOT decided by
            this field; different seeds are also not yet different CONSTRUCTIONS, which is what D-3 asks
            for.
        degenerate_seeds: seeds whose analysed window is constant to float64 round-off. See
            :func:`observe_gamma_seed_scatter` for why these refuse rather than pass.
        per_replicate: each seed's full :class:`~aleph.observe.stationarity.StationarityReport` dict.
        refusals: why, in words. Empty iff ``status == "MEASURED"``.
        contract: the thresholds, declared before the run and copied here.
        magnitude_claim: :data:`NO_MAGNITUDE_CLAIMED`.
    """

    observable: str
    status: str
    n_replicates: int
    seeds: tuple[int, ...]
    mean_pn_per_um: float | None
    seed_scatter_pn_per_um: float | None
    sem_across_seeds_pn_per_um: float | None
    mean_within_run_sem_pn_per_um: float | None
    scatter_inflation: float | None
    seed_scatter_fraction: float | None
    resolvable_fraction_one_seed: float | None
    resolvable_fraction_n_seeds: float | None
    window_s: float | None
    equilibration_s: float | None
    max_tau_int_s: float | None
    min_n_effective: float | None
    all_sokal_windows_closed: bool | None
    window_samples: int | None = None
    windows_open_on_transient: tuple[int, ...] = ()
    windows_skipped_on_transient: tuple[tuple[int, int], ...] = ()
    max_tau_samples: float | None = None
    equilibration_samples: int | None = None
    required_window_samples: int | None = None
    min_tau_reliability: float | None = None
    max_tau_relative_error: float | None = None
    drift_per_s_by_seed: dict[int, float | None] = field(default_factory=dict)
    drift_sigma_by_seed: dict[int, float | None] = field(default_factory=dict)
    drift_directions_agree: bool | None = None
    degenerate_seeds: tuple[int, ...] = ()
    per_replicate: tuple[dict[str, Any], ...] = ()
    refusals: tuple[str, ...] = ()
    contract: dict[str, Any] = field(default_factory=dict)
    magnitude_claim: str = NO_MAGNITUDE_CLAIMED

    @property
    def refused(self) -> bool:
        """``True`` when no cross-seed mean may be formed."""
        return self.status != "MEASURED"

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form, for a run record."""
        return {
            # ⚠ ORDER IS LOAD-BEARING. Everything that says "do not read the rest" comes first, because
            # a reader that truncates long fields for legibility swallowed exactly this on 2026-08-21:
            # `opens_on_transient: True` and `transient_ratio: 9.28` reached the terminal as the string
            # "per_replicate (4,712 chars)". Pinning individual keys to the front is what had already
            # failed there — two were pinned, the third was not thought of — so the invalidation flags
            # are top-level keys AND first.
            "status": self.status,
            "refusals": list(self.refusals),
            "windows_open_on_transient": list(self.windows_open_on_transient),
            "windows_skipped_on_transient": [list(t) for t in self.windows_skipped_on_transient],
            "degenerate_seeds": list(self.degenerate_seeds),
            "magnitude_claim": self.magnitude_claim,
            "step_acceptance": STEP_ACCEPTANCE_NOTE,
            "observable": self.observable,
            "n_replicates": self.n_replicates,
            "seeds": list(self.seeds),
            "mean_pn_per_um": self.mean_pn_per_um,
            "seed_scatter_pn_per_um": self.seed_scatter_pn_per_um,
            "sem_across_seeds_pn_per_um": self.sem_across_seeds_pn_per_um,
            "mean_within_run_sem_pn_per_um": self.mean_within_run_sem_pn_per_um,
            "scatter_inflation": self.scatter_inflation,
            "seed_scatter_fraction": self.seed_scatter_fraction,
            "resolvable_fraction_one_seed": self.resolvable_fraction_one_seed,
            "resolvable_fraction_n_seeds": self.resolvable_fraction_n_seeds,
            "window_s": self.window_s,
            "equilibration_s": self.equilibration_s,
            "max_tau_int_s": self.max_tau_int_s,
            "min_n_effective": self.min_n_effective,
            "all_sokal_windows_closed": self.all_sokal_windows_closed,
            "window_samples": self.window_samples,
            "max_tau_samples": self.max_tau_samples,
            "equilibration_samples": self.equilibration_samples,
            "required_window_samples": self.required_window_samples,
            "min_tau_reliability": self.min_tau_reliability,
            "max_tau_relative_error": self.max_tau_relative_error,
            "drift_per_s_by_seed": {int(k): v for k, v in self.drift_per_s_by_seed.items()},
            "drift_sigma_by_seed": {int(k): v for k, v in self.drift_sigma_by_seed.items()},
            "drift_directions_agree": self.drift_directions_agree,
            "per_replicate": [dict(r) for r in self.per_replicate],
            "contract": dict(self.contract),
        }


# ---------------------------------------------------------------------------------------------
# Load-bearing element families — ported from `components/incumbent/assembled_cortical_stress.py`
# ---------------------------------------------------------------------------------------------


def actin_axial_tension(
    f_ext: npt.ArrayLike, pos_actin: npt.ArrayLike, fiber_offsets: npt.ArrayLike,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    r"""Per-segment actin axial tension from the external-force running sum along each filament.

    Actin is an inextensible bead-rod (NF2007 reshape projection), so its axial tension is NOT a spring
    ``k·Δx`` — it is the constraint force holding each segment at rest length. At the overdamped balance
    the constraint force on each node balances the external force accumulated BEFORE the reshape, and
    force balance from a free filament end gives, for segment ``k``::

        T_k = − ( Σ_{j ≤ k} f_ext[node_j] ) · û_k

    so the actin tension carries everything the motors inject and the crosslinks transmit, with no
    double count and no dependence on the motor's internal anchor representation.

    Vectorised form of the ported per-node Python loop; ``tests/world/test_observe_gamma.py`` pins the two
    equal. The cumulative sum is taken over the whole array and the per-filament offset subtracted, which
    is exact for float64 up to the reassociation the pinning test tolerances allow.

    Args:
        f_ext: ``(n_actin, 3)`` external force per actin node [pN] — the accumulator output, pre-reshape.
        pos_actin: ``(n_actin, 3)`` actin node positions [µm].
        fiber_offsets: ``(n_fibers + 1,)`` node-range offsets; fiber ``i`` owns ``[off[i], off[i+1])``.

    Returns:
        ``(seg_a, seg_b, tension)`` — the two node indices of each segment and its SIGNED axial tension
        [pN] (+ tension, − compression). Filaments with fewer than 2 nodes contribute nothing.

    Raises:
        ValueError: if the shapes disagree or the offsets are not non-decreasing and in range. These are
            caller bugs — an offset array that is silently clipped measures a different cell.
    """
    f = np.asarray(f_ext, dtype=np.float64)
    pos = np.asarray(pos_actin, dtype=np.float64)
    off = np.asarray(fiber_offsets, dtype=np.int64).ravel()
    if f.ndim != 2 or f.shape[1] != 3 or pos.shape != f.shape:
        raise ValueError(f"f_ext {f.shape} and pos_actin {pos.shape} must both be (n_actin, 3)")
    if off.size < 2:
        raise ValueError("fiber_offsets needs at least 2 entries to describe one filament")
    if off[0] < 0 or off[-1] > f.shape[0] or np.any(np.diff(off) < 0):
        raise ValueError("fiber_offsets must be non-decreasing and within [0, n_actin]")

    starts, ends = off[:-1], off[1:]
    keep = (ends - starts) >= 2
    starts, ends = starts[keep], ends[keep]
    n_seg = ends - starts - 1
    total = int(n_seg.sum())
    empty_i, empty_f = np.zeros(0, np.int64), np.zeros(0, np.float64)
    if total == 0:
        return empty_i, empty_i, empty_f

    # Ragged arange: index within its own filament, for every segment of every filament.
    within = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(n_seg) - n_seg, n_seg)
    seg_a = np.repeat(starts, n_seg) + within
    seg_b = seg_a + 1

    csum = np.cumsum(f, axis=0)
    base = np.zeros((starts.size, 3), dtype=np.float64)
    nonzero_start = starts > 0
    base[nonzero_start] = csum[starts[nonzero_start] - 1]
    running = csum[seg_a] - np.repeat(base, n_seg, axis=0)

    d = pos[seg_b] - pos[seg_a]
    length = np.linalg.norm(d, axis=1)
    ok = length > _ZERO_LENGTH_UM
    tension = np.zeros(total, dtype=np.float64)
    tension[ok] = -np.einsum("ij,ij->i", running[ok], d[ok] / length[ok, None])
    return seg_a, seg_b, tension


def crosslink_tension(
    pos: npt.ArrayLike, xl: npt.ArrayLike, kxl: npt.ArrayLike, r0xl: npt.ArrayLike,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Hookean crosslink tensions ``T = k (|r| − r0)`` [pN] → ``(nodeA, nodeB, tension)``.

    Args:
        pos: ``(n_total, 3)`` node positions [µm].
        xl: ``(n_xl, 2)`` crosslink node-index pairs.
        kxl: ``(n_xl,)`` stiffnesses [pN/µm].
        r0xl: ``(n_xl,)`` rest lengths [µm].

    Returns:
        ``(node_a, node_b, tension)``, empty arrays when there is no crosslink.
    """
    p = np.asarray(pos, dtype=np.float64)
    pairs = np.asarray(xl, dtype=np.int64).reshape(-1, 2)
    if pairs.shape[0] == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float64)
    a, b = pairs[:, 0], pairs[:, 1]
    length = np.linalg.norm(p[b] - p[a], axis=1)
    tension = np.asarray(kxl, dtype=np.float64).ravel() * (
        length - np.asarray(r0xl, dtype=np.float64).ravel()
    )
    return a, b, tension


# ---------------------------------------------------------------------------------------------
# The method of planes, one pass per family
# ---------------------------------------------------------------------------------------------


def plane_force_sums(
    rA: npt.ArrayLike,
    rB: npt.ArrayLike,
    tension: npt.ArrayLike,
    *,
    R_um: float,
    centre: npt.ArrayLike | None,
    n_planes: int = DEFAULT_N_PLANES,
    chunk: int = _PLANE_CHUNK,
) -> npt.NDArray[np.float64]:
    r"""One family's per-plane cut tension [pN/µm], BEFORE the absolute value and the orientation mean.

    :func:`aleph.laws.gamma_estimator.method_of_planes_gamma` computes, for each of ``n_planes``
    Fibonacci normals ``n̂``::

        gamma_n = ( Σ_{elements crossing n̂} T · |û·n̂| ) / (2πR)     and returns  mean_n |gamma_n|

    This returns the vector ``gamma_n`` instead of collapsing it, and that is the whole point: the plane
    sum is **additive over element families**, while ``|·|`` and ``mean`` are not. So the per-family
    vectors can be added to get any GROUP's γ — source, network, total — for the cost of ONE pass over
    each family's elements, instead of one pass per group per family.

    **Why this is re-expressed rather than called.** The module docstring says the instrument is imported
    so a resting γ and a GATE-B γ differ by the cell state and not by the estimator, and that still holds:
    the normals come from :func:`~aleph.laws.gamma_estimator.fibonacci_plane_normals`, the formula is the
    one above, and ``test_the_plane_sums_reproduce_the_shared_estimator`` asserts agreement with
    ``method_of_planes_gamma`` itself. The imported function is the reference this is pinned to, which is
    a stronger guarantee than the earlier one — it is now checked rather than argued.

    **Why it is worth doing at all.** Measured on the dev Mac at 2 M elements: the shared estimator's
    per-plane Python loop costs 0.786 s per call. At the 4,361,496-node native cell the Lead measured
    8.34 s/step with γ on against 0.16 s/step with it off — a 50× that a stationarity window thousands of
    steps long pays every step, three times over for three seeds. Moving the cut onto the device is still
    the real fix and is not this module's to make.

    ⚠ **The realized factor depends on the FAMILY COUNT, and is not the 4.0× the bench prints.** With
    ``F`` families the old path made ``F + 2`` passes (one per family, one for the network group, one for
    the total) and this one makes ``F``, so the γ portion speeds up by ``(F + 2) / F``: 3× at one family,
    **2× at the two the resting path has** (actin + crosslink), 1.67× at three. The 4.0× benchmark is
    four estimator calls against one pass, i.e. the ``F = 1`` case with the source group also populated.
    Native confirms the right number: 8.34 → 4.06 s/step is 8.18 → 3.90 s of γ = **2.10× against a
    predicted 2.00×**. Quote ``(F+2)/F``, not 4×.

    Args:
        rA: ``(M, 3)`` element endpoint-A positions [µm].
        rB: ``(M, 3)`` element endpoint-B positions [µm].
        tension: ``(M,)`` signed scalar element tensions [pN].
        R_um: cell radius [µm] for the cut circumference ``2πR``.
        centre: ``(3,)`` cut centre [µm], or ``None`` to derive it ON DEVICE from the actin segment
            endpoints — the same definition as :func:`gamma_from_resting_readback`'s default, computed
            by the same fixed-order reduction, so the two agree by construction rather than by the
            caller passing the same thing to both. ⚠ This function used to REQUIRE it, and the
            2026-08-21 A/B is what that cost: the host defaulted to its centroid, the driver passed the
            device side the origin, and the comparison measured two centres rather than two kernels —
            7.5e-06 at step 0 growing to 15%. A default that only works if every caller remembers is
            the same defect one level out.
        n_planes: diametral cut orientations.
        chunk: elements per projection block. Performance only.

    Returns:
        ``(n_planes,)`` signed γ contribution of this family, in pN/µm. All zeros for an empty family,
        which is the additive identity and therefore the right answer rather than a special case.
    """
    normals = fibonacci_plane_normals(int(n_planes))
    out = np.zeros(int(n_planes), dtype=np.float64)
    rA_arr = np.asarray(rA, dtype=np.float64)
    if rA_arr.size == 0:
        return out
    rB_arr = np.asarray(rB, dtype=np.float64)
    c = (
        0.5 * (rA_arr.mean(axis=0) + rB_arr.mean(axis=0))
        if centre is None
        else np.asarray(centre, dtype=np.float64).ravel()
    )
    a_all = rA_arr - c
    b_all = rB_arr - c
    t_all = np.asarray(tension, dtype=np.float64).ravel()
    step = max(1, int(chunk))
    for lo in range(0, a_all.shape[0], step):
        a, b = a_all[lo : lo + step], b_all[lo : lo + step]
        d = b - a
        length = np.linalg.norm(d, axis=1)
        u = np.zeros_like(d)
        ok = length > 0.0
        u[ok] = d[ok] / length[ok, None]
        side_a = a @ normals.T                       # (m, n_planes)
        side_b = b @ normals.T
        crossing = (side_a * side_b) < 0.0
        weight = np.abs(u @ normals.T)
        out += (t_all[lo : lo + step, None] * weight * crossing).sum(axis=0)
    return out / (2.0 * np.pi * float(R_um))


# ---------------------------------------------------------------------------------------------
# γ for one resting configuration
# ---------------------------------------------------------------------------------------------


def gamma_from_resting_readback(
    *,
    R_um: float,
    pos: npt.ArrayLike,
    f_ext: npt.ArrayLike,
    n_actin: int,
    fiber_offsets: npt.ArrayLike,
    force_mask: Sequence[str],
    crosslinks: npt.ArrayLike | None = None,
    k_xl: npt.ArrayLike | None = None,
    r0_xl: npt.ArrayLike | None = None,
    extra_families: Mapping[str, tuple[npt.ArrayLike, npt.ArrayLike, npt.ArrayLike]] | None = None,
    source_families: Sequence[str] = ("crossbridge",),
    n_planes: int = DEFAULT_N_PLANES,
    centre: npt.ArrayLike | None = None,
) -> RestingGamma:
    """Method-of-planes γ from ONE synchronised host readback of a resting configuration.

    This takes host arrays, not a device cell, on purpose. `world/` may not import `engine/` or
    `components/` (`tests/architecture/test_layer_directions.py`), the resting driver lives in the
    latter, and the only thing crossing that line here is float64 NumPy the driver already has in hand:
    the accumulator output ``f``, ``pos_d``, ``foff_d``, ``xl_d``. Keeping it that way also means the
    whole path type-checks and self-tests on a machine with no CUDA, while the number itself is still
    measured on the GPU host by whoever calls it.

    Args:
        R_um: cortex shell radius [µm] for the cut circumference ``2πR``.
        pos: ``(n_total, 3)`` node positions [µm].
        n_actin: EXCLUSIVE UPPER BOUND of the actin nodes in ``pos`` — the slice handed to the
            extraction is ``[0, n_actin)``, and ``fiber_offsets`` index into it. ⚠ It does NOT assert
            that actin STARTS at 0, and in the arena it does not: a cortex claimed at ``[lo, lo+N)`` is
            addressed by absolute offsets, so a caller correctly passes ``lo + N`` here. This used to
            double as the population whose centroid became the default cut centre, which silently made
            that centroid the mean over membrane and envelope nodes too. It no longer does — see
            ``centre``.
        f_ext: ``(n_total, 3)`` accumulated external force per node [pN], PRE-reshape.
        n_actin: number of actin nodes.
        fiber_offsets: ``(n_fibers + 1,)`` actin filament node-range offsets.
        force_mask: the omit mask the inner solve ran with. **Required, no default** — see the module
            docstring for the constant 0.41526 pN/µm this prevents. Pass ``()`` to state "nothing omitted".
        crosslinks / k_xl / r0_xl: the crosslink family, or ``None`` for a cell with none.
        extra_families: further ``name -> (rA, rB, tension)`` families, for a caller that has them —
            the crossbridge family when a resting myosin setpoint makes one exist. Its extraction lives
            with the motor, not here.
        source_families: which family names roll into ``gamma_source``; everything else is network.
        n_planes: diametral cut orientations averaged.
        centre: cut centre [µm]. Defaults to the centroid of the ELEMENTS BEING CUT — the mean of the
            actin segment endpoints, falling back to all families when there is no actin. Derived from
            the elements rather than from an index range, because an index range encodes an assumption
            about where the population sits and the arena does not honour it. See the module docstring
            on why the two incumbent callers disagree and why this one chooses.

    Returns:
        A :class:`RestingGamma`, carrying :data:`NO_MAGNITUDE_CLAIMED`.

    Raises:
        ValueError: on a non-positive ``R_um`` or ``n_planes``, or on ``n_actin`` out of range.
    """
    R = float(R_um)
    if not np.isfinite(R) or R <= 0.0:
        raise ValueError(f"R_um must be finite and positive, got {R_um!r}")
    if int(n_planes) < 1:
        raise ValueError(f"n_planes must be >= 1, got {n_planes!r}")

    p = np.asarray(pos, dtype=np.float64)
    f = np.asarray(f_ext, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3 or f.shape != p.shape:
        raise ValueError(f"pos {p.shape} and f_ext {f.shape} must both be (n_total, 3)")
    n_a = int(n_actin)
    if not 0 <= n_a <= p.shape[0]:
        raise ValueError(f"n_actin={n_a} outside [0, {p.shape[0]}]")

    families: dict[str, tuple[npt.NDArray, npt.NDArray, npt.NDArray]] = {}
    sa, sb, t_actin = actin_axial_tension(f[:n_a], p[:n_a], fiber_offsets)
    if sa.size:
        families["actin"] = (p[sa], p[sb], t_actin)
    n_xl = 0
    if crosslinks is not None:
        xa, xb, t_xl = crosslink_tension(p, crosslinks, k_xl, r0_xl)
        n_xl = int(xa.size)
        if n_xl:
            families["crosslink"] = (p[xa], p[xb], t_xl)
    for name, fam in (extra_families or {}).items():
        rA = np.asarray(fam[0], dtype=np.float64)
        if rA.size:
            families[name] = (rA, np.asarray(fam[1], np.float64), np.asarray(fam[2], np.float64))

    if centre is None:
        # ⚠ Derived from the ELEMENTS, never from an index range. This was `p[:n_a].mean(axis=0)`, which
        # is the actin centroid only if actin starts at index 0. In the arena it does not: a cortex
        # claimed at [lo, lo+N) is addressed by ABSOLUTE offsets, so a correct caller passes lo+N as
        # `n_actin` and that slice then averages every membrane and envelope node below `lo` as well.
        # Measured 2026-08-21 in a driver A/B: γ differed from a device path passing the origin by
        # 7.5e-06 at step 0 — against 3.865e-12 for the same code on identical inputs — and grew to 15%
        # as the cell moved, because a centre offset does not perturb the sum, it FLIPS WHICH ELEMENTS
        # CROSS a plane. Both arms were wrong; neither centre was the cortex centroid.
        endpoints = [fam for fam in (families.get("actin"),) if fam is not None] or list(
            families.values()
        )
        c = (
            0.5
            * (
                np.vstack([fam[0] for fam in endpoints]).mean(axis=0)
                + np.vstack([fam[1] for fam in endpoints]).mean(axis=0)
            )
            if endpoints
            else np.zeros(3)
        )
    else:
        c = np.asarray(centre, dtype=np.float64).ravel()
    if c.shape != (3,):
        raise ValueError(f"centre must be 3 components, got {c.shape}")

    src = tuple(source_families)
    net = tuple(name for name in families if name not in src)

    # ONE pass per family. Every γ below is a sum of these vectors — see :func:`plane_force_sums`.
    planes = {
        name: plane_force_sums(fam[0], fam[1], fam[2], R_um=R, centre=c, n_planes=int(n_planes))
        for name, fam in families.items()
    }

    def _gamma(names: Sequence[str]) -> float:
        present = [planes[name] for name in names if name in planes]
        if not present:
            return 0.0
        return float(np.mean(np.abs(np.sum(present, axis=0))))

    def _group(names: Sequence[str]) -> float | None:
        """``None`` for a group with no elements — absent is not zero. See the dataclass docstring."""
        return _gamma(names) if any(name in planes for name in names) else None

    by_family = {name: _gamma([name]) for name in families}
    # Over the nodes the actin segments actually touch, for the same reason as the centre: `f[:n_a]`
    # summed every foreign node below the cortex claim into a diagnostic labelled "actin".
    if sa.size:
        touched = np.unique(np.concatenate([sa, sb]))
        net_vec = f[touched].sum(axis=0)
    else:
        net_vec = np.zeros(3)
    return RestingGamma(
        gamma_total_pn_per_um=_gamma(tuple(families)),
        gamma_network_pn_per_um=_group(net),
        gamma_source_pn_per_um=_group(src),
        gamma_by_family_pn_per_um=by_family,
        net_force_vector_pn=(float(net_vec[0]), float(net_vec[1]), float(net_vec[2])),
        n_planes=int(n_planes),
        R_um=R,
        centre_um=(float(c[0]), float(c[1]), float(c[2])),
        force_mask=tuple(str(name) for name in force_mask),
        n_actin=n_a,
        n_crosslinks=n_xl,
        n_source_elements=sum(
            int(families[name][2].size) for name in src if name in families
        ),
    )


# ---------------------------------------------------------------------------------------------
# γ across replicate seeds — (e) 1 D-2
# ---------------------------------------------------------------------------------------------


def observe_gamma_seed_scatter(
    replicates: Mapping[int, npt.ArrayLike],
    dt_s: float,
    *,
    observable: str = "gamma_total_pn_per_um",
    min_replicates: int = DEFAULT_MIN_REPLICATES,
    resolve_k: float = DEFAULT_RESOLVE_K,
    min_tau_windows: float = DEFAULT_MIN_TAU_WINDOWS,
    drift_sigma: float = DEFAULT_DRIFT_SIGMA,
    transient_sigma: float = DEFAULT_TRANSIENT_SIGMA,
    min_span_s: float | None = None,
    sokal_c: float = DEFAULT_SOKAL_C,
) -> SeedScatterReport:
    """Judge γ across independent SEEDS, and report the error bar the seeds support.

    Each replicate is first judged on its own by
    :func:`aleph.observe.stationarity.assess_stationarity` — Sokal automatic windowing for ``tau_int``,
    an equilibration search for where the transient ends, a correlation-aware drift significance. That
    step is what measures the window **W**. Only then are the replicate means combined, and the spread
    that is reported is the spread ACROSS them.

    A replicate that is not ``STATIONARY`` refuses the whole report rather than being dropped: dropping
    the replicates that did not settle and averaging the rest is selection on the outcome.

    ⚠ **Investigating a refusal: do NOT call** :func:`~aleph.observe.stationarity.assess_stationarity`
    **on a slice.** It runs its OWN equilibration search on whatever it is handed, so slicing a series
    and re-calling it applies the search TWICE and the second pass hides what the first did. Measured
    2026-08-21: probing cuts that way returned ``STATIONARY, tau ~ 102`` at 8,000 and 12,000 samples,
    which reads as "the retained window was fine" — while the single-pass search on the full series was
    in fact returning a cut 6,000 samples too early, and the apparent counter-evidence was the
    instrument being applied to its own output. Use
    :func:`~aleph.observe.stationarity.integrated_autocorrelation_time` and
    :func:`~aleph.observe.stationarity.window_opens_on_a_transient` directly on the slice: those compute,
    they do not search.

    ⚠ **A CONSTANT series is refused too, and that is the more dangerous case.** Measured 2026-08-21 on
    the resting path: with ``relax = 0`` the engine emits a γ that is constant to float64 round-off, and
    a long enough run of it is judged ``STATIONARY`` with ``mean = 0.18860334`` and ``sem = 0.0`` — a
    perfect steady state, from a series that says nothing. At 20 samples it refuses only because
    ``n_eff`` is under the length bar; at 300 it passes. "The solver stopped moving" being certified as
    "this is the resting cell" is objection (b) of `PI_DECISION_e1_STATIONARITY_2026-08-16.md` §1
    arriving through the statistics instead of through the solver, and more samples make it worse rather
    than better. So degeneracy is judged INDEPENDENTLY of the verdict and of the length.

    The test for it is float64 round-off — ``std <= eps * |mean|`` — and that is deliberately a machine
    property rather than a threshold. A tuning constant here would be a number chosen to decide which
    runs pass, which the charter forbids; ``eps`` is not chosen and is grid-invariant.

    **Why there is NO separate refusal for an unclosed Sokal window**, though it looks like the same
    shape of hole. When the window never closes, ``tau_int`` is a LOWER bound and ``n_eff`` an UPPER one,
    i.e. the length gate would be tested against a number biased in the flattering direction. But the two
    cannot co-occur: an unclosed window means the loop reached ``m = n-1`` without ``m >= c*tau_samples``,
    so ``tau_samples > (n-1)/c`` and therefore::

        n_eff = n / (2*tau_samples) < c*n / (2*(n-1))  <=  (2/3)*c     for n >= 4

    At ``sokal_c = 5`` that is ``n_eff < 3.34`` against a required 25 — excluded with 7.5x to spare, and
    searching AR(1) up to rho = 0.99999 and random walks produced no unclosed window at all. So
    ``all_sokal_windows_closed`` is a diagnostic to read on a REFUSED report; in a ``MEASURED`` one it
    cannot be ``False``. That guarantee is arithmetic between two independently declared thresholds, so
    :func:`observe_gamma_seed_scatter` REFUSES A CONTRACT that breaks it rather than letting a lowered
    length bar reopen the hole in silence.

    Args:
        replicates: ``seed -> gamma series``, one entry per independent run. A mapping keyed by seed
            rather than a list, so two runs of the same seed cannot be passed off as two replicates.
        dt_s: SAMPLING interval in the SIMULATION's time [s] — the recording stride times the integrator
            step, not the integrator step alone. Handing the step understates ``tau_int`` by exactly that
            stride, which flatters ``n_eff`` and the error bar.
            ⚠ **Three different "seconds per step" are in reach of a caller and only one is this one.**
            Wall-clock seconds per step (0.2173 on the 2026-08-21 device run) is throughput and has
            nothing to do with a correlation time; the integrator's physical ``dt`` is this parameter;
            and a step COUNT is neither but is the only one always sourced. Passing wall clock produces a
            complete, plausible report whose every ``*_s`` field is a throughput measurement wearing a
            physics label — nothing in here can detect it, because both are positive finite floats.
            When the run's own time axis is unsourced — as it is wherever mobility is — read
            ``window_samples``, ``min_n_effective`` and ``required_window_samples``, which are counts.
        observable: name recorded in the report.
        min_replicates: below this, no scatter is reported. Contract, recorded.
        resolve_k: coverage factor turning the scatter into a resolution floor. Contract, recorded.
        min_tau_windows, drift_sigma, transient_sigma, min_span_s, sokal_c: forwarded verbatim to
            :func:`~aleph.observe.stationarity.assess_stationarity` and recorded in each replicate's own
            contract block. Declare them before the run.

    Returns:
        A :class:`SeedScatterReport`. ``status`` is ``"MEASURED"`` or ``"REFUSED"``; **neither is a
        stationarity verdict about the cell and neither makes a γ quotable.**

    Raises:
        ValueError: if ``min_replicates`` < 2, ``resolve_k`` is not positive and finite, or a series is
            unusable (fewer than 4 samples, non-finite). Those are caller bugs, not findings.
    """
    if int(min_replicates) < 2:
        raise ValueError("min_replicates must be at least 2 — a scatter needs two numbers to exist")
    if not np.isfinite(resolve_k) or resolve_k <= 0.0:
        raise ValueError(f"resolve_k must be finite and positive, got {resolve_k!r}")
    # An unclosed Sokal window bounds n_eff below (2/3)*sokal_c (module docstring). The length bar is
    # n_eff >= min_tau_windows/2, so a bar under (4/3)*sokal_c would let a series whose tau_int is only a
    # LOWER bound clear a gate built on the UPPER bound of its own n_eff. Not a threshold of ours — a
    # relation between two the caller declared, and it raises rather than refusing because a contract
    # that cannot mean what it says is a caller bug, not a finding about the data.
    if float(min_tau_windows) < (4.0 / 3.0) * float(sokal_c):
        raise ValueError(
            f"min_tau_windows={float(min_tau_windows):g} is below (4/3)*sokal_c="
            f"{(4.0 / 3.0) * float(sokal_c):g}. Under that bar an unclosed Sokal window can pass, and "
            "there tau_int is a LOWER bound and n_effective an UPPER one — the gate would be scored on "
            "the flattering side of its own estimate. Raise the bar or lower sokal_c."
        )

    seeds = tuple(int(s) for s in replicates)
    contract: dict[str, Any] = {
        "min_replicates": int(min_replicates),
        "resolve_k": float(resolve_k),
        "min_tau_windows": float(min_tau_windows),
        "drift_sigma": float(drift_sigma),
        "transient_sigma": float(transient_sigma),
        "min_span_s": None if min_span_s is None else float(min_span_s),
        "sokal_c": float(sokal_c),
        "sample_dt_s": float(dt_s),
        "error_bar_basis": (
            "SEED SCATTER (std over replicate means, ddof=1), never the within-run sem — STATE.md (f) "
            "and PI_DECISION_e1_STATIONARITY_2026-08-16.md D-2"
        ),
        "degeneracy_basis": (
            "a replicate whose analysed window is constant to float64 round-off (std <= eps*|mean|) is "
            "REFUSED regardless of verdict or length. Not a tuning constant: eps is a machine property"
        ),
        "time_axis": (
            "seconds here are samples * sample_dt_s and are only as sourced as the caller's dt_s, which "
            "must be the SIMULATION's sampling interval — not wall-clock seconds per step, which would "
            "make every *_s field a throughput measurement wearing a physics label. A caller with no "
            "physical time axis should pass dt_s = 1.0, and then EVERY *_s field holds a COUNT despite "
            "its name: read window_samples / max_tau_samples / equilibration_samples / "
            "required_window_samples / min_n_effective, which are counts under any dt_s. "
            "⚠ drift_per_s_by_seed is a RATE and scales as 1/dt_s, the OPPOSITE way to every other _s "
            "field here: at dt_s = 1.0 it is per SAMPLE, and applying the counts-despite-their-name "
            "correction to it corrects in the wrong direction. Its sign and its significance "
            "(drift_sigma_by_seed) are dt_s-invariant; only its magnitude is not"
        ),
        "declared": "before the run; recorded so the report is judged against its own criterion",
    }

    reports = {
        seed: assess_stationarity(
            replicates[seed], dt_s, observable=f"{observable}[seed={seed}]",
            min_tau_windows=min_tau_windows, drift_sigma=drift_sigma,
            transient_sigma=transient_sigma, min_span_s=min_span_s, c=sokal_c,
        )
        for seed in seeds
    }
    per_replicate = tuple(reports[seed].as_dict() | {"seed": seed} for seed in seeds)

    refusals: list[str] = []
    if len(seeds) < int(min_replicates):
        refusals.append(
            f"{len(seeds)} replicate(s) against a declared minimum of {int(min_replicates)}; a seed "
            "scatter cannot be measured from fewer, and a within-run sem is not a substitute"
        )
    eps = float(np.finfo(np.float64).eps)
    degenerate: list[int] = []
    for seed in seeds:
        rep = reports[seed]
        # Judged BEFORE the verdict and independently of it: a constant series PASSES the verdict once it
        # is long enough, which is the failure this catches. The scale is the raw series mean — used as a
        # magnitude to compare round-off against, not quoted, and available even when the verdict
        # withholds a mean.
        scale = float(abs(np.mean(np.asarray(replicates[seed], dtype=np.float64))))
        if (rep.tau is not None and rep.tau.zero_variance) or rep.std <= eps * scale:
            degenerate.append(seed)
            refusals.append(
                f"seed {seed}: DEGENERATE — the analysed window is constant to float64 round-off "
                f"(std {rep.std:.3g} <= eps*|mean| {eps * scale:.3g}). A constant series has no "
                f"correlation time and nothing to be stationary about, so MORE SAMPLES WILL NOT HELP: "
                f"they make it pass. Verdict on its own was {rep.verdict.value}"
            )
        elif rep.refused:
            flag = (
                " ⚠ its retained window OPENS ON A TRANSIENT, so every tau-derived field in this report "
                "describes that transient and not a correlation time — do not read them"
                if rep.equilibration is not None and rep.equilibration.opens_on_transient
                else ""
            )
            refusals.append(
                f"seed {seed}: {rep.verdict.value} — "
                f"{rep.notes.get('reason', 'see its own report')}{flag}"
            )

    # Window bookkeeping is reported whether or not the mean is, because "how long did we look" is the
    # thing (e) 1 has to fix and it is answerable even when the answer is "not long enough".
    spans = [rep.analysed_time_s for rep in reports.values()]
    starts = [
        rep.equilibration.start_time_s for rep in reports.values() if rep.equilibration is not None
    ]
    start_idx = [
        rep.equilibration.start_index for rep in reports.values() if rep.equilibration is not None
    ]
    opens = tuple(
        seed for seed in seeds
        if reports[seed].equilibration is not None and reports[seed].equilibration.opens_on_transient
    )
    # ⚠ `opens` became a near-empty field on 2026-08-22 and this is where the signal it used to carry
    # now lives. PI ruling 18 (b) makes `equilibration_start` SKIP candidates whose window opens on a
    # rise, so the winner is almost never contaminated and `opens` reads clean on a series that is not.
    # The contamination did not go away; the search walked around it. This counts what it walked around.
    skipped = tuple(
        (seed, reports[seed].equilibration.n_contaminated_skipped) for seed in seeds
        if reports[seed].equilibration is not None
        and reports[seed].equilibration.n_contaminated_skipped
    )
    taus = [rep.tau.tau_int_s for rep in reports.values() if rep.tau is not None]
    n_effs = [rep.tau.n_effective for rep in reports.values() if rep.tau is not None]
    closed = [rep.tau.window_closed for rep in reports.values() if rep.tau is not None]
    reliabilities = [rep.tau.reliability for rep in reports.values() if rep.tau is not None]
    slopes = {seed: reports[seed].drift_per_s for seed in seeds}
    drift_sigmas = {seed: reports[seed].drift_sigma_measured for seed in seeds}
    measured_slopes = [v for v in slopes.values() if v is not None]
    tau_samples = [rep.tau.tau_samples for rep in reports.values() if rep.tau is not None]
    # Madras & Sokal 1988: sigma(tau_int)/tau_int ~ sqrt(2*(2M+1)/N). M is the window the Sokal sum was
    # cut at, which `TauEstimate` already carries; nothing else in the tree computes this.
    tau_errors = [
        float(np.sqrt(2.0 * (2.0 * rep.tau.window_lag + 1.0) / rep.n_samples))
        for rep in reports.values()
        if rep.tau is not None and rep.n_samples > 0
    ]

    common = dict(
        observable=observable,
        n_replicates=len(seeds),
        seeds=seeds,
        window_s=min(spans) if spans else None,
        equilibration_s=max(starts) if starts else None,
        max_tau_int_s=max(taus) if taus else None,
        min_n_effective=min(n_effs) if n_effs else None,
        all_sokal_windows_closed=all(closed) if closed else None,
        window_samples=min(rep.n_samples for rep in reports.values()) if reports else None,
        windows_open_on_transient=opens,
        windows_skipped_on_transient=skipped,
        max_tau_samples=max(tau_samples) if tau_samples else None,
        equilibration_samples=max(start_idx) if start_idx else None,
        required_window_samples=(
            int(np.ceil(float(min_tau_windows) * max(tau_samples))) if tau_samples else None
        ),
        min_tau_reliability=min(reliabilities) if reliabilities else None,
        max_tau_relative_error=max(tau_errors) if tau_errors else None,
        drift_per_s_by_seed=slopes,
        drift_sigma_by_seed=drift_sigmas,
        drift_directions_agree=(
            len({float(np.sign(v)) for v in measured_slopes}) == 1 if measured_slopes else None
        ),
        degenerate_seeds=tuple(degenerate),
        per_replicate=per_replicate,
        contract=contract,
    )

    if refusals:
        return SeedScatterReport(
            status="REFUSED", mean_pn_per_um=None, seed_scatter_pn_per_um=None,
            sem_across_seeds_pn_per_um=None, mean_within_run_sem_pn_per_um=None,
            scatter_inflation=None, seed_scatter_fraction=None,
            resolvable_fraction_one_seed=None, resolvable_fraction_n_seeds=None,
            refusals=tuple(refusals), **common,
        )

    means = np.array([reports[seed].mean for seed in seeds], dtype=np.float64)
    within = np.array([reports[seed].sem for seed in seeds], dtype=np.float64)
    mean = float(means.mean())
    scatter = float(means.std(ddof=1))
    sem_across = scatter / np.sqrt(len(seeds))
    mean_within = float(within.mean())

    # A zero denominator is a real case (a perfectly reproducible series, or a γ that is zero), and it is
    # reported as unknown rather than as an infinity that would read like a huge measured inflation.
    inflation = scatter / mean_within if mean_within > 0.0 else None
    frac = scatter / abs(mean) if mean != 0.0 else None
    return SeedScatterReport(
        status="MEASURED",
        mean_pn_per_um=mean,
        seed_scatter_pn_per_um=scatter,
        sem_across_seeds_pn_per_um=float(sem_across),
        mean_within_run_sem_pn_per_um=mean_within,
        scatter_inflation=inflation,
        seed_scatter_fraction=frac,
        resolvable_fraction_one_seed=None if frac is None else float(resolve_k) * frac,
        resolvable_fraction_n_seeds=(
            None if mean == 0.0 else float(resolve_k) * float(sem_across) / abs(mean)
        ),
        refusals=(),
        **common,
    )


# ---------------------------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------------------------


def _ar1(rho: float, n: int, seed: int, mean: float, sigma: float = 1.0) -> np.ndarray:
    """A correlated series with a known mean — the case a within-run sem gets wrong."""
    rng = np.random.default_rng(seed)
    out = np.zeros(n)
    for i in range(1, n):
        out[i] = rho * out[i - 1] + rng.standard_normal()
    return mean + sigma * out


def _demo() -> None:
    """Self-check: γ is computed, the window is measured, and the error bar comes from the seeds."""
    # -- γ on a configuration whose answer is known by hand ------------------------------------
    # One straight filament of 3 nodes along +x, its two free ends pushed TOWARD each other by 1 pN.
    # The running sum from the first end gives −1 pN in both segments: uniform compression. Reversing
    # the end forces pulls the rod apart and gives uniform +1 pN tension. Sign sense, both ways.
    pos = np.array([[-1.0, 0, 0], [0.0, 0, 0], [1.0, 0, 0]])
    f = np.array([[1.0, 0, 0], [0.0, 0, 0], [-1.0, 0, 0]])
    sa, sb, T = actin_axial_tension(f, pos, [0, 3])
    assert sa.tolist() == [0, 1] and sb.tolist() == [1, 2], (sa, sb)
    assert np.allclose(T, [-1.0, -1.0]), T
    assert np.allclose(actin_axial_tension(-f, pos, [0, 3])[2], [1.0, 1.0])

    # A filament of one node is not a filament; a zero-length segment is not a division by zero.
    assert actin_axial_tension(f, pos, [0, 1, 3])[0].tolist() == [1]
    degenerate = np.array([[0.0, 0, 0], [0.0, 0, 0]])
    assert actin_axial_tension(np.ones((2, 3)), degenerate, [0, 2])[2].tolist() == [0.0]

    # Against the ported per-node loop, on a ragged multi-filament case.
    rng = np.random.default_rng(7)
    off = np.array([0, 5, 6, 12, 20])
    px, fx = rng.normal(size=(20, 3)), rng.normal(size=(20, 3))
    loop_a, loop_t = [], []
    for i in range(off.size - 1):
        s, e = int(off[i]), int(off[i + 1])
        if e - s < 2:
            continue
        cum = np.zeros(3)
        for k in range(e - s - 1):
            cum = cum + fx[s + k]
            u = px[s + k + 1] - px[s + k]
            loop_a.append(s + k)
            loop_t.append(-float(np.dot(cum, u / np.linalg.norm(u))))
    vec_a, _, vec_t = actin_axial_tension(fx, px, off)
    assert vec_a.tolist() == loop_a, (vec_a, loop_a)
    assert np.allclose(vec_t, loop_t, rtol=1e-12, atol=1e-12), np.max(np.abs(vec_t - np.array(loop_t)))

    # -- the whole readback path, with the crosslink family and the required mask ----------------
    theta = np.linspace(0, 2 * np.pi, 33)[:-1]
    ring = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1) * 5.0
    inward = -ring / np.linalg.norm(ring, axis=1, keepdims=True)
    rep = gamma_from_resting_readback(
        R_um=5.0, pos=ring, f_ext=inward, n_actin=ring.shape[0],
        fiber_offsets=[0, ring.shape[0]], force_mask=(),
        crosslinks=[[0, 16]], k_xl=[100.0], r0_xl=[9.0],
    )
    assert rep.gamma_total_pn_per_um > 0.0
    assert rep.gamma_source_pn_per_um is None, "no crossbridge family: ABSENT, not zero"
    assert rep.gamma_network_pn_per_um is not None
    assert set(rep.gamma_by_family_pn_per_um) == {"actin", "crosslink"}
    assert rep.n_crosslinks == 1 and rep.force_mask == () and rep.magnitude_claim is NO_MAGNITUDE_CLAIMED
    assert "gamma_total_pn_per_um" in rep.as_dict()
    assert "ACCEPTANCE_UNDEFINED" in rep.as_dict()["step_acceptance"]
    try:
        gamma_from_resting_readback(R_um=0.0, pos=ring, f_ext=inward, n_actin=32,
                                    fiber_offsets=[0, 32], force_mask=())
    except ValueError as exc:
        assert "finite and positive" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a non-positive radius must refuse")

    # -- seed scatter: the 67x failure mode, reproduced and refused ------------------------------
    # Three strongly correlated replicates that each settle TIGHTLY about their OWN level. That is the
    # measured situation: within-run fluctuation small, seed-to-seed offset large. A within-run sem here
    # is a confident error bar around the wrong number, which is exactly the 67x STATE.md (f) records.
    dt = 0.01
    series = {s: _ar1(0.97, 6000, s, mean=100.0 + 0.5 * i, sigma=0.02)
              for i, s in enumerate((11, 12, 13))}
    got = observe_gamma_seed_scatter(series, dt)
    assert got.status == "MEASURED", got.refusals
    assert got.window_s is not None and got.window_s > 0.0, "the window W must be measured"
    assert got.max_tau_int_s is not None and got.min_n_effective is not None
    assert got.scatter_inflation is not None and got.scatter_inflation > 10.0, got.scatter_inflation
    assert got.sem_across_seeds_pn_per_um > got.mean_within_run_sem_pn_per_um
    assert got.resolvable_fraction_one_seed > got.resolvable_fraction_n_seeds > 0.0

    # Two seeds is not three, and a drifting replicate is not dropped — it refuses the report.
    two = observe_gamma_seed_scatter({11: series[11], 12: series[12]}, dt)
    assert two.status == "REFUSED" and two.mean_pn_per_um is None
    assert two.window_s is not None, "a refusal still reports how long we looked"
    # A constant series is judged STATIONARY with sem 0 once it is long enough. It must not reach a mean.
    flat = dict(series)
    flat[15] = np.full(6000, 0.18860334)
    flat_report = observe_gamma_seed_scatter(flat, dt)
    assert flat_report.status == "REFUSED" and flat_report.degenerate_seeds == (15,)
    assert "MORE SAMPLES WILL NOT HELP" in " ".join(flat_report.refusals)
    assert flat_report.window_samples is not None and flat_report.window_samples > 0

    drifting = dict(series)
    drifting[14] = series[13] + np.linspace(0.0, 50.0, series[13].size)
    assert observe_gamma_seed_scatter(drifting, dt).status == "REFUSED"

    # Seeds drifting OPPOSITE ways read identically to seeds drifting together unless the sign survives.
    opposed = observe_gamma_seed_scatter(
        {s: _ar1(0.97, 6000, s, mean=100.0, sigma=0.02) + np.linspace(0.0, sign * 3.0, 6000)
         for s, sign in ((21, +1.0), (22, -1.0), (23, +1.0))}, dt)
    assert opposed.drift_directions_agree is False, opposed.drift_per_s_by_seed
    assert set(opposed.drift_per_s_by_seed) == {21, 22, 23}
    together = observe_gamma_seed_scatter(
        {s: _ar1(0.97, 6000, s, mean=100.0, sigma=0.02) + np.linspace(0.0, 3.0, 6000)
         for s in (24, 25, 26)}, dt)
    assert together.drift_directions_agree is True
    assert together.status == opposed.status == "REFUSED", "both refuse; only the SIGNS tell them apart"

    if wp.get_cuda_device_count() > 0:
        _device_demo()
    else:
        print(
            f"observe_gamma: DEVICE GATES NOT RUN — warp {wp.config.version} sees no CUDA here. "
            "This is NOT 'not applicable': gamma_planes_device is UNVERIFIED until this module is run "
            "in the PRODUCTION runtime on the GPU host as: srun --jobid=<N> --overlap env "
            "CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<n> python -m aleph.world.observe_gamma "
            "(all three parts required — the device pin alone does not make the process legitimate)."
        )

    print(
        "observe_gamma self-check OK — "
        f"gamma_total={rep.gamma_total_pn_per_um:.4g} pN/um (NO MAGNITUDE CLAIMED), "
        f"W={got.window_s:.4g} s, tau_max={got.max_tau_int_s:.4g} s, "
        f"S/sem_within={got.scatter_inflation:.4g}x"
    )


def _device_demo() -> None:
    """The device gates, as plain asserts — runnable with NO test runner.

    ⚠ **This exists because `skipped` is not `verified`.** The equivalence and reproducibility gates for
    :func:`gamma_planes_device` were written as CUDA-only pytest tests, which skip on the dev Mac
    (correctly) and cannot be invoked in the production runtime at all: `ffn_sim` — the environment every
    physics number this project has was measured in — has no pytest, and the GPU host has never had a
    pytest process run on it. So the count of CUDA-gated tests ever executed in the runtime that produces
    the numbers was ZERO, while the suite printed them as `skipped`, which reads as *not applicable here*
    when it means *not verified anywhere* (Lead, 2026-08-21,
    `NO_CUDA_TEST_HAS_EVER_RUN_ON_THE_RUN_HOST_2026-08-21.md`, PI queue item 15).

    The fix needs no dependency. This module already carries a framework-free self-check that the repo
    invokes as `python -m aleph.world.observe_gamma`, and that runs under any interpreter — including the
    production one, with pytest absent. ⚠ The invocation on the GPU host is::

        srun --jobid=<N> --overlap env CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<n> \
            python -m aleph.world.observe_gamma

    All three parts are load-bearing and the shorter forms are each half of it. The `CUDA_VISIBLE_DEVICES`
    pin confines the process to the granted card but does NOT make it legitimate: `gpu-bypass-watch` has
    run as root since June, reads each GPU process's environ for `SLURM_JOB_ID`, validates it against
    `scontrol`, and kills anything without one (INTERVAL 5, GRACE 10, HITS 2) — so a bare process is
    admitted for about twenty seconds, and a probe that finishes faster reads as proof there is no
    enforcement. `srun --jobid --overlap` alone supplies the job id but leaves `CUDA_VISIBLE_DEVICES`
    unset and all three cards visible: attaching to an allocation does not inherit the wrapper's
    confinement. (Lead, 2026-08-21; PI queue item 16 asks for this line to be documented, not built.) So the gate moves here and the pytest test calls it, rather than
    a second copy of the same arithmetic living where it cannot run.

    It also PRINTS the interpreter's Warp version and device. Two Warp versions on one host with no
    record naming either is how an equivalence test between a device kernel and a host estimator gets
    scored against the wrong compiler — and for that test the compiler is the variable under test.

    Raises:
        RuntimeError: if called without CUDA. Refusing is the point; a CPU fallback here would be the
            path the charter says must not exist.
    """
    if wp.get_cuda_device_count() < 1:
        raise RuntimeError("_device_demo needs a CUDA device; there is no CPU path and must not be one")
    dev = "cuda:0"
    rng = np.random.default_rng(2026)
    # ⚠ The fixture must have MANY MORE elements than `_PLANE_BLOCKS`, and that is not a size
    # preference. At 320 segments against 1024 blocks — the first version of this check, and the one
    # that passed on 2026-08-21 — 704 blocks receive nothing and the rest receive exactly one element,
    # so `_plane_partial_sums_kernel`'s `while e < n_elements: ... e += n_blocks` body runs at most ONCE
    # per thread and the striding accumulation, which is the part that can actually be wrong, is never
    # executed. A gate can be green and still have missed its own subject. Asserted below rather than
    # left to whoever next edits these numbers.
    n_fibers, per_fiber, n_xl = 4_000, 9, 12_000
    off = (np.arange(n_fibers + 1, dtype=np.int64) * per_fiber)
    n_actin = int(off[-1])
    # ⚠ CHAINS, not random points grouped by index. The earlier fixture drew independent shell points
    # and called every 9 of them a filament; `actin_axial_tension` walks that happily and the
    # equivalence arithmetic is unaffected, but anything CENTRE-related is then unrepresentative — the
    # segment-endpoint mean and the node mean differ by ~sigma/sqrt(n_fibers) for random grouping and
    # coincide to 1e-16 for uniformly spaced chains, which is what a cortex actually is.
    seed_pt = rng.normal(size=(n_fibers, 3))
    seed_pt /= np.linalg.norm(seed_pt, axis=1, keepdims=True)
    tangent = rng.normal(size=(n_fibers, 3))
    tangent -= (tangent * seed_pt).sum(1, keepdims=True) * seed_pt
    tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
    arc = np.arange(per_fiber)[None, :, None]
    pos = (5.0 * seed_pt[:, None, :] + 0.05 * arc * tangent[:, None, :]).reshape(n_actin, 3)
    f_ext = rng.normal(size=(n_actin, 3)) * 0.5
    bonds = rng.integers(0, n_actin, size=(n_xl, 2)).astype(np.int32)
    k_xl, r0_xl = np.full(n_xl, 8.2e5), np.full(n_xl, 0.05)
    centre = pos.mean(axis=0)

    sa, sb, t_actin = actin_axial_tension(f_ext, pos, off)
    xa, xb, t_xl = crosslink_tension(pos, bonds, k_xl, r0_xl)
    host = {
        "actin": plane_force_sums(pos[sa], pos[sb], t_actin, R_um=5.0, centre=centre, n_planes=64),
        "crosslink": plane_force_sums(pos[xa], pos[xb], t_xl, R_um=5.0, centre=centre, n_planes=64),
    }

    seg_off = segment_offsets(off)
    assert int(seg_off[-1]) == sa.size, (seg_off[-1], sa.size)
    per_block = sa.size // _PLANE_BLOCKS
    assert per_block >= 8, (
        f"only {per_block} element(s) per reduction block: the striding accumulation in "
        f"`_plane_partial_sums_kernel` is barely exercised. Raise n_fibers, do not lower _PLANE_BLOCKS."
    )
    static = dict(
        centre=centre,
        fiber_offsets_d=wp.array(off.astype(np.int32), dtype=wp.int32, device=dev),
        seg_offsets_d=wp.array(seg_off.astype(np.int32), dtype=wp.int32, device=dev),
        n_segments=int(seg_off[-1]), n_actin=n_actin, R_um=5.0,
        crosslinks_d=wp.array(bonds, dtype=wp.int32, device=dev),
        k_xl_d=wp.array(k_xl, dtype=wp.float64, device=dev),
        r0_xl_d=wp.array(r0_xl, dtype=wp.float64, device=dev),
        n_planes=64, device=dev,
    )

    runs = [
        gamma_planes_device(
            pos_d=wp.array(pos, dtype=wp.vec3d, device=dev),
            f_ext_d=wp.array(f_ext, dtype=wp.vec3d, device=dev), **static)
        for _ in range(5)
    ]

    # (1) EQUIVALENCE with the host estimator, which is itself pinned to `laws.method_of_planes_gamma`.
    #     Relative, not bitwise: the device reduces in (plane, block) order and the host in element
    #     order — the same difference the host's own chunk-invariance check already allows.
    worst = 0.0
    for name in ("actin", "crosslink"):
        rel = np.max(np.abs(runs[0][name] - host[name]) / (np.abs(host[name]) + 1e-300))
        worst = max(worst, float(rel))
        assert np.allclose(runs[0][name], host[name], rtol=1e-11, atol=1e-12), (name, rel)

    # (2) The chain closes on the shared estimator itself, not only on this module's host copy.
    total_dev = float(np.mean(np.abs(runs[0]["actin"] + runs[0]["crosslink"])))
    total_ref = float(method_of_planes_gamma(
        np.vstack([pos[sa], pos[xa]]), np.vstack([pos[sb], pos[xb]]),
        np.concatenate([t_actin, t_xl]), 5.0, n_planes=64, centre=centre))
    assert abs(total_dev - total_ref) <= 1e-11 * abs(total_ref), (total_dev, total_ref)

    # (3) BITWISE reproducibility. Not a nicety: `observe_gamma_seed_scatter` refuses a series whose
    #     analysed window is constant to within `eps*|mean|`, and an order-varying reduction would
    #     manufacture noise above that floor and switch the refusal off.
    for other in runs[1:]:
        for name in ("actin", "crosslink"):
            assert np.array_equal(runs[0][name], other[name]), name

    # (4) The two DEFAULTS must agree, by construction rather than by the caller passing the same thing
    #     to both. This is the check the 2026-08-21 A/B did not have.
    #     ⚠ It must compare DEFAULT against DEFAULT. Its first version compared the device default
    #     against `runs[0]`, which was computed with an explicit and different centre — i.e. it made the
    #     very mistake it exists to catch, and reported 193% on the host for it.
    host_default = plane_force_sums(pos[sa], pos[sb], t_actin, R_um=5.0, centre=None, n_planes=64)
    derived = gamma_planes_device(
        pos_d=wp.array(pos, dtype=wp.vec3d, device=dev),
        f_ext_d=wp.array(f_ext, dtype=wp.vec3d, device=dev),
        **{**static, "centre": None})
    centre_rel = float(
        np.max(np.abs(derived["actin"] - host_default) / (np.abs(host_default) + 1e-300))
    )
    assert centre_rel <= 1e-11, centre_rel

    print(
        f"observe_gamma DEVICE gates PASS — warp {wp.config.version} on {wp.get_device(dev)}; "
        f"{sa.size} actin segments + {n_xl} crosslinks over {_PLANE_BLOCKS} blocks "
        f"({per_block} elements/block); equivalence worst rel {worst:.3e} vs host; "
        f"chain to method_of_planes_gamma closed; 5 identical launches bitwise equal; "
        f"device-derived centre matches the host default to {centre_rel:.3e}"
    )



# ---------------------------------------------------------------------------------------------
# The same measurement on the device
# ---------------------------------------------------------------------------------------------
#
# γ dominates a native step, and the device→host copy is not why. Measured at native by the Lead
# (2026-08-21, RTX 4090-1, production runtime, inside Slurm job 73; 68,814 filaments × 61 nodes =
# 4,197,654 nodes / 4,128,840 segments):
#
#     gamma_planes_device        10.77 ms/call   (n=8)
#     host readback             178.66 ms/call   (n=8)   — 16.6x the whole device computation
#     host estimator          5,964.91 ms/call   (n=3, median of 6156 / 5965 / 5815)
#     host path total         6,143.57 ms/step                          ratio 570x
#
# So "avoid the roundtrip" was never the answer: a free transfer still leaves ~6 s of host NumPy per
# step. The work is the estimator, and the estimator lives here, so the kernel does too.
#
# ⚠ Two earlier figures quoted in this header are RETRACTED, both the Lead's and both corrected by the
# same measurement. The readback is **576.0 MB, not 209.4** — the arrays return at full arena CAPACITY
# (12,000,000 nodes × 3 × float64) while ~4.59 M are live, so 62% of every transfer is unallocated
# arena; the original figure was computed from the live count rather than measured. And the copy's share
# of the step follows from the corrected size.
#
# ⚠ An open discrepancy, left open. The run difference says host γ costs **3.95 s/step** (4.1212 with
# `--emit-gamma` against 0.1677 without, 5 runs); the isolated timing above says **6.14 s/step** on the
# same arrays, same offsets, same `force_mask=()`. Nobody knows why yet and no reason has been invented.
# Both are ≫ 10.77 ms, so nothing here depends on the resolution — but it is not written down as settled.
#
# ⚠ **A device-γ step cost is still ARITHMETIC.** "≈0.18 s/step" is 0.1677 measured plus 0.0108
# measured, added by hand. The honest number comes from one driver run with this path wired in, and it
# has not happened. Not quoted as a measurement here or anywhere.

# **MEASURED 2026-08-21**, in the PRODUCTION runtime (`ffn_sim`, warp 1.14.0) on the granted RTX 4090,
# pinned with `srun --jobid=<N> --overlap env CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<n>`.
# Fixture: 32,000 actin segments (CHAINS) + 12,000 crosslinks over 1,024 blocks = 31 elements per block.
#
#     equivalence with the host estimator      worst rel  2.442e-13
#     chain closed on method_of_planes_gamma   yes
#     five identical launches                  bitwise equal
#     device-derived centre vs host default    2.392e-13
#
# The last line is the check the 2026-08-21 A/B did not have, and it caught a real fault on its FIRST
# run — 1.935 against a 1e-11 bar, from three causes none of which was the reduction.
#
# ⚠ **Not comparable to the 6.626e-13 that stood here before.** That was the same code on a fixture whose
# "filaments" were independent shell points bundled in nines; this one builds chains, which is what a
# cortex is. Different geometry, not a smaller error on the same problem.
#
# ⚠ **What that PASS covers and does not.** A code-equivalence result at the self-check's own scale. It
# is not a physics result and not a magnitude: no γ value is claimed, `STATE.md` (c) 3 and (c) 17 stand,
# and the fixture is not the native population — see NATIVE EQUIVALENCE below for that.
#
# ⚠ Native equivalence has since been measured for the `actin` family — see below. The `crosslink` family
# has NOT been checked at native, and deliberately is not queued: it reaches the same
# `_plane_partial_sums_kernel` through the same index arrays, so there is nothing family-specific in the
# summation to test, and `_crosslink_tension_kernel` is three lines already covered at self-check scale.
# The driver's resting γ path does not pass crosslinks today either.
#
# An earlier run of the same code passed at 4.521e-14 on a fixture of 320 segments over 1,024 blocks —
# one element per block, so the striding accumulation never executed and the number bounded the indexing
# rather than the loop. Kept as the reason the fixture is what it is, not as a result.
#
# **NATIVE EQUIVALENCE, measured 2026-08-21** on the same cortex — 68,814 filaments / 4,128,840 segments,
# 4,032 elements per block, `actin` family: worst relative error **3.865e-12** against the check's 1e-11,
# WITHIN. `_PLANE_BLOCKS` stays 1024 and the tolerance was never touched, which is what declaring the
# contingency in advance was for. The arena's force field is zero before any step, so both paths would
# have summed to zero and the ratio would have been 0/0 — the run filled the cortex range with a
# deterministic `default_rng(20260821)` normal field, scale 0.5. Nothing about that field is physical and
# no γ is claimed; the host median |plane sum| of 3.740259 is recorded so the ratio is visibly not
# flattered by a near-zero denominator.
#
# ⚠ **My extrapolation that produced that contingency was built on a regime that did not apply, and
# bracketing the answer does not make it right.** The crossing guard admits only ~0.3% of elements to any
# one plane (measured on shell geometry, host-side), so the terms actually ACCUMULATED per thread are:
#
#     fixture            elements/block    accumulated/block    worst rel error
#     320 segments              1               ~0.001            4.521e-14
#     32,000 segments          31               ~0.10             6.626e-13
#     4,128,840 segments    4,032               ~12.9             3.865e-12
#
# The first two threads sum ZERO OR ONE term. There is no accumulation there to have a scaling law about
# — those points measure indexing and single-term rounding — so fitting sqrt-versus-linear across all
# three was never valid, and the constant crossing fraction cannot explain the apparent change in law
# either (it scales both intervals identically). **Only the native point is in the accumulation regime,
# so there is ONE point and no established law.** If the population or `_PLANE_BLOCKS` changes
# materially, MEASURE — do not extrapolate from this, including with the numbers above.
#
# ⚠ **And half of the declared contingency is now doubtful.** Raising `_PLANE_BLOCKS` shortens the
# per-thread stride but LENGTHENS the second stage, which sums `n_blocks` partials sequentially — at
# native that reduce is already 1024 deep against a stride of ~12.9 accumulated terms, so the sequential
# depth is dominated by the reduce and raising blocks moves the larger term the wrong way. The half that
# stands unconditionally is the other one: **never loosen the tolerance.** If a future measurement does
# cross it, the fix to investigate first is a pairwise/tree second stage, not more blocks — and that is a
# measurement to make, not a conclusion reached here.

#: Element-strided blocks per plane. Sets the reduction's summation order, so it is part of what makes a
#: result reproducible: the same `n_blocks` gives the same bits, a different one gives a different
#: rounding. Not a tuning knob — change it and re-pin, do not sweep it.
_PLANE_BLOCKS: int = 1024


@wp.kernel
def _actin_axial_tension_kernel(
    pos: wp.array(dtype=wp.vec3d),                  # (n_actin,) node positions [µm]
    f_ext: wp.array(dtype=wp.vec3d),                # (n_actin,) external force per node [pN], pre-reshape
    fiber_offsets: wp.array(dtype=wp.int32),        # (n_fibers+1,) node ranges
    seg_offsets: wp.array(dtype=wp.int32),          # (n_fibers+1,) output ranges (host, topology-static)
    seg_a: wp.array(dtype=wp.int32),                # out
    seg_b: wp.array(dtype=wp.int32),                # out
    tension: wp.array(dtype=wp.float64),            # out (n_segments,) signed axial tension [pN]
):
    """One thread per FILAMENT, walking its own nodes — the running sum is serial by nature and each
    filament's is independent, so the parallelism is over fibers and no scan primitive is needed."""
    i = wp.tid()
    s = fiber_offsets[i]
    e = fiber_offsets[i + 1]
    if e - s < 2:
        return
    o = seg_offsets[i]
    cum = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    for k in range(e - s - 1):
        cum = cum + f_ext[s + k]
        d = pos[s + k + 1] - pos[s + k]
        length = wp.length(d)
        t = wp.float64(0.0)
        if length > wp.float64(0.0):
            t = -wp.dot(cum, d / length)
        seg_a[o + k] = s + k
        seg_b[o + k] = s + k + 1
        tension[o + k] = t


@wp.kernel
def _crosslink_tension_kernel(
    pos: wp.array(dtype=wp.vec3d),                  # (n_total,) node positions [µm]
    bonds: wp.array(dtype=wp.int32, ndim=2),        # (n_xl, 2) node-index pairs
    k_xl: wp.array(dtype=wp.float64),               # (n_xl,) [pN/µm]
    r0_xl: wp.array(dtype=wp.float64),              # (n_xl,) [µm]
    idx_a: wp.array(dtype=wp.int32),                # out (n_xl,) — endpoint indices for the plane pass
    idx_b: wp.array(dtype=wp.int32),                # out (n_xl,)
    tension: wp.array(dtype=wp.float64),            # out (n_xl,) T = k(|r| - r0) [pN]
):
    """Hookean crosslink tension, one thread per bond."""
    b = wp.tid()
    i = bonds[b, 0]
    j = bonds[b, 1]
    idx_a[b] = i
    idx_b[b] = j
    tension[b] = k_xl[b] * (wp.length(pos[j] - pos[i]) - r0_xl[b])


@wp.kernel
def _plane_partial_sums_kernel(
    pos: wp.array(dtype=wp.vec3d),                  # (n_total,) node positions [µm]
    idx_a: wp.array(dtype=wp.int32),                # (M,) element endpoint-A node index
    idx_b: wp.array(dtype=wp.int32),                # (M,) element endpoint-B node index
    tension: wp.array(dtype=wp.float64),            # (M,) signed element tension [pN]
    normals: wp.array(dtype=wp.vec3d),              # (n_planes,) Fibonacci unit normals
    centre: wp.vec3d,                               # cut centre [µm]
    n_elements: wp.int32,
    n_blocks: wp.int32,
    partial: wp.array(dtype=wp.float64, ndim=2),    # out (n_planes, n_blocks) — NOT reduced here
):
    """Per-(plane, block) partial cut force. One thread strides the element list; nothing is atomic.

    The formula is :func:`aleph.laws.gamma_estimator.method_of_planes_gamma`'s inner sum, term for term:
    an element crosses plane ``n̂`` when its endpoints sit on opposite sides of it, and contributes
    ``T * |û·n̂|``.
    """
    p, k = wp.tid()
    n_hat = normals[p]
    acc = wp.float64(0.0)
    e = k
    while e < n_elements:
        ra = pos[idx_a[e]] - centre
        rb = pos[idx_b[e]] - centre
        if wp.dot(ra, n_hat) * wp.dot(rb, n_hat) < wp.float64(0.0):
            d = rb - ra
            length = wp.length(d)
            if length > wp.float64(0.0):
                acc += tension[e] * wp.abs(wp.dot(d / length, n_hat))
        e += n_blocks
    partial[p, k] = acc


@wp.kernel
def _plane_reduce_kernel(
    partial: wp.array(dtype=wp.float64, ndim=2),    # (n_planes, n_blocks)
    n_blocks: wp.int32,
    inv_circumference: wp.float64,                  # 1 / (2*pi*R) [1/µm]
    out: wp.array(dtype=wp.float64),                # out (n_planes,) [pN/µm]
):
    """Sum each plane's blocks in FIXED index order. One thread per plane; deterministic by construction."""
    p = wp.tid()
    total = wp.float64(0.0)
    for k in range(n_blocks):
        total = total + partial[p, k]
    out[p] = total * inv_circumference


@wp.kernel
def _endpoint_sum_partial_kernel(
    pos: wp.array(dtype=wp.vec3d),                  # (n_total,) node positions [µm]
    idx_a: wp.array(dtype=wp.int32),                # (M,) element endpoint-A node index
    idx_b: wp.array(dtype=wp.int32),                # (M,) element endpoint-B node index
    n_elements: wp.int32,
    n_blocks: wp.int32,
    partial: wp.array(dtype=wp.float64, ndim=2),    # out (3, n_blocks) — NOT reduced here
):
    """Per-(component, block) partial sum of both endpoints, for the cut centre.

    Shaped so :data:`_plane_reduce_kernel` finishes it unchanged: pass ``inv = 1/(2*M)`` and the three
    outputs are ``sum(a_c + b_c) / (2M)``, which is exactly the host default's
    ``0.5*(mean(rA) + mean(rB))``. Same striding, same fixed-order second stage, same reason — an
    order-varying centre would jitter which elements cross a plane, and that jitter is what the
    degeneracy refusal in `observe_gamma_seed_scatter` must not see manufactured.
    """
    c, k = wp.tid()
    acc = wp.float64(0.0)
    e = k
    while e < n_elements:
        acc += pos[idx_a[e]][c] + pos[idx_b[e]][c]
        e += n_blocks
    partial[c, k] = acc


def segment_offsets(fiber_offsets: npt.ArrayLike) -> npt.NDArray[np.int64]:
    """Output ranges for :data:`_actin_axial_tension_kernel`, one entry per fiber plus a total.

    Topology-static: compute once at build and keep the device copy, never per step. Fibers with fewer
    than two nodes own an EMPTY range rather than being dropped, so fiber ``i``'s segments are always at
    ``[out[i], out[i+1])`` and the kernel's thread index is the fiber index with no compaction map.
    """
    off = np.asarray(fiber_offsets, dtype=np.int64).ravel()
    per_fiber = np.maximum(np.diff(off) - 1, 0)
    return np.concatenate([[0], np.cumsum(per_fiber)]).astype(np.int64)


def gamma_planes_device(
    *,
    pos_d: wp.array,
    f_ext_d: wp.array,
    fiber_offsets_d: wp.array,
    seg_offsets_d: wp.array,
    n_segments: int,
    n_actin: int,
    R_um: float,
    centre: npt.ArrayLike,
    crosslinks_d: wp.array | None = None,
    k_xl_d: wp.array | None = None,
    r0_xl_d: wp.array | None = None,
    n_planes: int = DEFAULT_N_PLANES,
    n_blocks: int = _PLANE_BLOCKS,
    device: str | None = None,
) -> dict[str, npt.NDArray[np.float64]]:
    """Per-family plane sums computed ON THE DEVICE, returning ``n_planes`` float64 per family.

    Same quantity as :func:`plane_force_sums`, same units, same sign convention; the equivalence is
    asserted by ``test_device_matches_host``, which is CUDA-only. What crosses the bus is ``n_planes``
    doubles per family — 512 bytes at the default — instead of the 209.4 MB/step the host path reads
    back, and none of the O(n_elements) arithmetic happens on the host at all.

    Verified equivalent to :func:`plane_force_sums` and bitwise reproducible on the production runtime —
    see the section header for what that PASS covers and what it does not. It is a code-equivalence
    result at the self-check's scale; it says nothing about any γ magnitude and nothing about native.

    ⚠ Scratch is allocated per call. At a stationarity window of thousands of steps a caller should hoist
    it; that is a caller-side change and it does not alter a single returned digit, so it is deliberately
    not hidden inside this function.

    Args:
        pos_d: ``(n_total,)`` ``wp.vec3d`` node positions [µm], device-resident.
        f_ext_d: ``(n_total,)`` ``wp.vec3d`` accumulated external force [pN], PRE-reshape.
        fiber_offsets_d: ``(n_fibers+1,)`` ``wp.int32`` actin node ranges.
        seg_offsets_d: ``(n_fibers+1,)`` ``wp.int32`` from :func:`segment_offsets`, uploaded once.
        n_segments: ``seg_offsets[-1]``.
        n_actin: number of actin nodes (the kernel launch dimension is fibers, this is for the record).
        R_um: cortex shell radius [µm].
        centre: ``(3,)`` cut centre [µm]. Required — this function does not choose it; see the module
            docstring on why the two incumbent callers disagree.
        crosslinks_d / k_xl_d / r0_xl_d: the crosslink family, or ``None``.
        n_planes: diametral cut orientations.
        n_blocks: reduction blocks per plane. Fixes the summation order — see the section header.
        device: Warp device. Must resolve to CUDA; a non-CUDA device RAISES rather than proceeding.

    Returns:
        ``{family: (n_planes,) plane sums [pN/µm]}``, ready to add across families exactly as
        :func:`plane_force_sums`' output is.

    Raises:
        RuntimeError: if the resolved device is not CUDA. The charter has no CPU path and a fallback
            here would be one.
        ValueError: if ``centre`` is ``None`` and there are no actin segments to derive it from. The
            fallback the host default uses there needs the crosslink endpoints, which is a different
            reduction; refusing is better than silently centring on something else.
    """
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(
            f"gamma_planes_device resolved a non-CUDA device ({dev}). There is no CPU simulation path "
            "and this must not silently become one — run inside a Slurm allocation on the GPU host."
        )
    n_fibers = int(fiber_offsets_d.shape[0]) - 1
    normals_d = wp.array(fibonacci_plane_normals(int(n_planes)), dtype=wp.vec3d, device=dev)
    partial_d = wp.zeros((int(n_planes), int(n_blocks)), dtype=wp.float64, device=dev)
    out_d = wp.zeros(int(n_planes), dtype=wp.float64, device=dev)
    inv_circ = wp.float64(1.0 / (2.0 * np.pi * float(R_um)))

    def _planes(idx_a_d, idx_b_d, tension_d, count: int) -> npt.NDArray[np.float64]:
        if count == 0:
            return np.zeros(int(n_planes), dtype=np.float64)
        partial_d.zero_()
        wp.launch(_plane_partial_sums_kernel, dim=(int(n_planes), int(n_blocks)),
                  inputs=[pos_d, idx_a_d, idx_b_d, tension_d, normals_d, centre_w,
                          wp.int32(count), wp.int32(n_blocks), partial_d], device=dev)
        wp.launch(_plane_reduce_kernel, dim=int(n_planes),
                  inputs=[partial_d, wp.int32(n_blocks), inv_circ, out_d], device=dev)
        return out_d.numpy().copy()

    families: dict[str, npt.NDArray[np.float64]] = {}

    if centre is None and not int(n_segments):
        raise ValueError(
            "centre=None needs actin segments to derive the cut centre from, and there are none. Pass "
            "an explicit centre, or use gamma_from_resting_readback, whose fallback covers this case."
        )

    seg_a_d = wp.zeros(max(1, int(n_segments)), dtype=wp.int32, device=dev)
    seg_b_d = wp.zeros(max(1, int(n_segments)), dtype=wp.int32, device=dev)
    seg_t_d = wp.zeros(max(1, int(n_segments)), dtype=wp.float64, device=dev)
    wp.launch(_actin_axial_tension_kernel, dim=n_fibers,
              inputs=[pos_d, f_ext_d, fiber_offsets_d, seg_offsets_d, seg_a_d, seg_b_d, seg_t_d],
              device=dev)

    if centre is None:
        centre_partial_d = wp.zeros((3, int(n_blocks)), dtype=wp.float64, device=dev)
        centre_out_d = wp.zeros(3, dtype=wp.float64, device=dev)
        wp.launch(_endpoint_sum_partial_kernel, dim=(3, int(n_blocks)),
                  inputs=[pos_d, seg_a_d, seg_b_d, wp.int32(int(n_segments)), wp.int32(n_blocks),
                          centre_partial_d], device=dev)
        wp.launch(_plane_reduce_kernel, dim=3,
                  inputs=[centre_partial_d, wp.int32(n_blocks),
                          wp.float64(1.0 / (2.0 * int(n_segments))), centre_out_d], device=dev)
        centre_v = centre_out_d.numpy()
    else:
        centre_v = np.asarray(centre, dtype=np.float64).ravel()
    if centre_v.shape != (3,):
        raise ValueError(f"centre must be 3 components, got {centre_v.shape}")
    centre_w = wp.vec3d(centre_v[0], centre_v[1], centre_v[2])

    if int(n_segments):
        families["actin"] = _planes(seg_a_d, seg_b_d, seg_t_d, int(n_segments))

    if crosslinks_d is not None and int(crosslinks_d.shape[0]):
        n_xl = int(crosslinks_d.shape[0])
        xa_d = wp.zeros(n_xl, dtype=wp.int32, device=dev)
        xb_d = wp.zeros(n_xl, dtype=wp.int32, device=dev)
        xt_d = wp.zeros(n_xl, dtype=wp.float64, device=dev)
        wp.launch(_crosslink_tension_kernel, dim=n_xl,
                  inputs=[pos_d, crosslinks_d, k_xl_d, r0_xl_d, xa_d, xb_d, xt_d], device=dev)
        families["crosslink"] = _planes(xa_d, xb_d, xt_d, n_xl)

    return families


if __name__ == "__main__":
    _demo()
