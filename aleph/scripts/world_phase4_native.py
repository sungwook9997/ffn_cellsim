"""PHASE 4 — stand the cell, bind the laws that exist, take steps, and judge nothing.

**What this driver claims: the step machinery runs on a native arena with real forces from bound
laws. What it does NOT claim: that any step was accepted, or any magnitude.**

Both of those are deliberate and neither is a hedge.

⚠ **No step here is accepted, because this engine has no acceptance criterion.** Both predicates it
has shipped are one-sided and both were measured: ``ledger.py:239-259`` accepts inside the Higham
summation bound ``n_terms*eps64`` — the rounding error of the summation, which ``sf_motor_slice.py:65``
says tests *adjoint closure* — so **it cannot fail for physics**; and ``descent_ratio`` read a
post-rollback residual pinned to its own start, so **it cannot succeed** (job 69, VOID). The amendment
meant to replace them, ``STATE.md`` (e) 1, is still undecided — ⚠ **but not for the reason this line
gave until 2026-08-22.** Cortical tension IS emitted on the resting path now (``--emit-gamma``, below);
what is open is WHICH statistic, over WHAT window, against WHAT variance, and all three are the PI's.
So the run uses :class:`~aleph.world.step.UndefinedAcceptance` and every step records
``ACCEPTANCE_UNDEFINED``.

⚠ **No magnitude is claimed.** The forces are real — ``laws/`` kernels over arena-addressed ranges,
the binding measured force-identical to private arrays at 1.108e-16 relative (``78942fc4``) — but the
COEFFICIENTS they are launched with are declared axes and test points, not sourced values for this
cell, and the residual they produce is therefore not a physical residual. It is reported so the
machinery can be seen working, and the record says which it is.

**What the run is FOR.** Three questions that do not need a verdict to answer:

1. Does the arena take a step at native population without the partition invariant breaking?
2. What does one step cost, in wall time, at 4.6 M nodes?
3. Does the residual move at all when a force is applied — i.e. is anything actually wired?

Sanity Gate:
    * dimensions — positions µm, forces pN, ``dt_phys`` s; no conversion happens in this driver.
    * boundary cases — a run of zero steps is refused; a population that fails
      :func:`~aleph.world.geometry.assert_inside_membrane` aborts before any step.
    * conservation — ``assert_partitioned`` after every step, via :class:`~aleph.world.step.WorldStep`.
    * CFL/precision — float64 throughout; the inner relaxation's step size is DECLARED, never derived
      from a stability estimate this driver would then be free to adjust.
    * sign sense — a positive surface tension must SHRINK a closed surface; asserted against the
      built membrane before stepping, so a sign flip aborts rather than producing a plausible run.
    * measurement protocol — one ``StepVerdict`` per step, appended in order, each carrying the
      predicate that judged it.

Runtime: NVIDIA Warp on CUDA. Refuses a non-CUDA device.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--dt-phys-s", type=float, default=0.05,
                    help="outer physical timestep [s]. This driver's clock, NOT a physiological "
                         "quantity; it is recorded as such.")
    ap.add_argument("--cortex-seg-um", type=float, default=0.05)
    ap.add_argument("--sigma-pn-per-um", type=float, default=1.0,
                    help="membrane surface tension [pN/um]. ⚠ A TEST POINT, not a sourced value for "
                         "this cell. It sets the force scale, so no residual from this run is a "
                         "physical residual and the record says so.")
    ap.add_argument("--kappa-tilde-pn-um", type=float, default=20.0,
                    help="calibrated dihedral coupling [pN*um]. ⚠ Also a test point.")
    ap.add_argument("--heartbeat-s", type=float, default=30.0, metavar="SECONDS",
                    help="wall-clock cadence for the progress file. TIME-based, not step-based: the "
                         "step cost is the unknown a heartbeat exists to measure, so a step interval "
                         "would write at an unknown rate. Costs a few hundred bytes and no device "
                         "work; it is not the snapshot and does not share its budget.")
    ap.add_argument("--snapshot-every", type=int, default=0, metavar="N",
                    help="write a position frame every N steps into an .alephcell sequence beside "
                         "--out, openable by aleph.viz.cell_app. TOPOLOGY IS WRITTEN ONCE — a frame "
                         "of 4.59 M positions is 55 MB against 39 MB of topology, so re-writing the "
                         "topology per frame would triple the file AND imply the connectivity is "
                         "time-dependent, which it is not. 0 disables.")
    ap.add_argument("--core-only", action="store_true",
                    help="stand only cortex, membrane and nuclear envelope. Default is the FULL cell, "
                         "the same twelve populations the PHASE 1 driver measures — stepping a "
                         "different cell from the one that was measured is how two numbers stop being "
                         "comparable. The record says which arm ran.")
    ap.add_argument("--envelope-kappa-tilde-pn-um", type=float, default=None,
                    help="bind Helfrich bending to the NUCLEAR ENVELOPE [pN·um]. "
                         "laws.nucleus_envelope.KAPPA_NE_PN_UM = 0.0828 (~20 kBT, lamina+bilayer). "
                         "⚠ BENDING ONLY, and the omission is the point: the envelope's AREAL "
                         "TENSION is bilinear in the areal strain (chromatin below a ~10%% knee, "
                         "lamin-A/C above), so a CONSTANT sigma would be an imposed prestress on a "
                         "shell built at rest — which is what the GATE-B contract forbids when it "
                         "requires tension to EMERGE. Binding it needs the bilinear law and a "
                         "per-step strain, and that is scoped separately.")
    ap.add_argument("--bind", action="append", default=[], metavar="POP=k_axial[,kappa[,rest_um]]",
                    help="bind a strand population's axial and bending laws, e.g. "
                         "--bind microtubule=1e5,20.0. Repeatable. BOTH coefficients are required for "
                         "each population and neither has a default: k_axial is EA/L_seg and kappa is "
                         "kB*T*L_p, and a persistence length is per filament TYPE — actin ~17 um "
                         "against a microtubule's ~5000 um, which are not interchangeable. A "
                         "population named here that PHASE 1 did not build is refused rather than "
                         "skipped.")
    ap.add_argument("--cortex-k-axial", type=float, default=None,
                    help="cortex axial stiffness [pN/um]. NO DEFAULT. For actin this is EA/L_seg and "
                         "both factors are declared axes. Omitted, the cortex is not bound and the "
                         "record says the run carried membrane forces only.")
    ap.add_argument("--cortex-kappa", type=float, default=None,
                    help="cortex flexural rigidity kappa [pN*um^2] = kB*T*L_p. NO DEFAULT — the "
                         "persistence length is per filament type (actin ~17 um, MT ~5000 um) and "
                         "they are not interchangeable. Enters as kappa/seg^3, so it is CUBED in the "
                         "step.")
    ap.add_argument("--seed", type=int, default=0,
                    help="cortex build seed. Recorded, because gamma is judged against SEED SCATTER "
                         "and a replicate without its seed is not a replicate.")
    ap.add_argument("--gamma-device", action="store_true",
                    help="compute gamma with the device kernel instead of the host estimator. Same "
                         "quantity — verified to 3.87e-12 worst relative error at NATIVE scale — and "
                         "570x faster, because the host path reads back 576 MB and then spends 5.96 s "
                         "of NumPy per step. Off by default: a run in flight was sized on the host "
                         "path and its contract names that instrument.")
    ap.add_argument("--emit-gamma", action="store_true",
                    help="compute cortical tension on this path and hand it to the predicate through "
                         "StepContext.observables. ⚠ Its ABSENCE is what made STATE.md (e) 1 "
                         "unwritable rather than merely undecided, and this flag is what ended that "
                         "— so passing it is what makes a run usable as (e) 1 input at all. NO gamma "
                         "VALUE is quotable: (c) 3 and "
                         "(c) 17 stand and the record carries NO_MAGNITUDE_CLAIMED. Requires the "
                         "cortex to be bound, since gamma is a property of the cortex.")
    ap.add_argument("--mobility", type=float, default=None,
                    help="per-node mobility 1/gamma [um/(pN*s)] for an overdamped LANGEVIN step. "
                         "⚠ UNSOURCED — the arena's mobility is zero because a mobility belongs to a "
                         "drag law nobody has written, and the one this repo has is on "
                         "FORBIDDEN_DRAG_SYMBOLS. Given, the run is DRIVEN AND DAMPED and can produce "
                         "a fluctuating series; the TIMESCALE of that series is unsourced and the "
                         "record says so. Mutually exclusive with --relax-um-per-pn, which is a "
                         "deterministic relaxation and not a thermostat.")
    ap.add_argument("--temperature-k", type=float, default=310.0,
                    help="[K]. 310 is physiological and is the one sourced quantity in the step.")
    ap.add_argument("--relax-um-per-pn", type=float, default=0.0,
                    help="explicit relaxation gain [um/pN] for the inner solve. DEFAULT 0.0 = the "
                         "geometry does not move, and that is the honest default: a mobility is a "
                         "drag law's parameter and the builders leave it at zero deliberately "
                         "(membrane.py:49-50). A nonzero value here is a DECLARED numerical "
                         "relaxation, not a drag law, and the record labels it that way.")
    args = ap.parse_args()

    if args.steps < 1:
        raise SystemExit("refused: a run of zero steps measures nothing")

    import numpy as np
    import warp as wp
    import warp.utils

    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.build import build_all, builders_not_standing
    from aleph.world.geometry import assert_inside_membrane, footprint
    from aleph.world.laws_bind import (
        accumulate_area_tension,
        accumulate_helfrich,
        upload_surface_topology,
    )
    # ⚠ `segment_offsets` and `gamma_planes_device` were used below without being imported, so
    # --gamma-device died with a NameError on its first real invocation while --emit-gamma was
    # fine. Same shape as observe_gamma.py's __main__ ordering bug this morning: a branch that
    # nothing had executed yet. The A/B run that found it is the only reason it did not ship.
    from aleph.world.observe_gamma import (gamma_from_resting_readback, gamma_planes_device,
                                          segment_offsets)
    from aleph.world.populations import build_remaining_populations
    from aleph.world.step import UndefinedAcceptance, WorldStep
    from aleph.world.thermostat import langevin_step, stability_bound_s
    from aleph.world.strand_bind import (
        accumulate_axial,
        accumulate_bending,
        upload_strand_topology,
    )

    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:
        raise RuntimeError(f"PHASE 4 is native-only; {args.device!r} resolved to {dev}, not CUDA")

    t_build = time.perf_counter()
    cap = {Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000, Kind.ANGLE3: 12_000_000,
           Kind.ANGLE4: 2_000_000, Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
           Kind.BOND: 1_000_000, Kind.GRID_CELL: 4_000_000}
    arena = WorldArena(capacity=cap, device=args.device)
    # ⚠ The seed reaches build_all only if it forwards one. If it does not, this run is NOT an
    # independent replicate however many times it is repeated, and the record says which it is —
    # three identical runs offered as a seed scatter would be the single most misleading artifact
    # this driver could produce.
    import inspect as _inspect
    seed_reaches_build = "seed" in _inspect.signature(build_all).parameters
    cell = (build_all(args.device, arena=arena, seed=args.seed) if seed_reaches_build
            else build_all(args.device, arena=arena))
    fp = footprint(7.5)
    wp.synchronize_device(dev)
    # ⚠ The full cell, from the SAME builder the PHASE 1 driver and the renderer use. Standing only
    # cortex/membrane/envelope here would step a different cell from the one PHASE 1 measures, and
    # --bind could only ever name three populations. `--core-only` keeps the old three-population
    # arm available for a cheap timing comparison, and says which it ran in the record.
    built: dict[str, object] = ({} if args.core_only
                                else build_remaining_populations(arena, fp, seg_um=args.cortex_seg_um))
    n_live_nodes = int(arena.census()["live"]["node"])   # every claim, not just the last
    build_s = time.perf_counter() - t_build

    pops = {"cortex": cell.cortex, "membrane": cell.membrane,
            "nuclear_envelope": cell.envelope, **built}
    placement = assert_inside_membrane(arena.node_arrays["position"].numpy(), pops, fp.r_cell_um)
    arena.assert_partitioned()

    # Bind the surface laws that exist. The binding itself is not re-derived here: `78942fc4` measured
    # production kernels over arena ranges at 1.108e-16 relative against private arrays.
    topo = upload_surface_topology(arena, cell.membrane)

    # ── the nuclear envelope, if the caller states its bending rigidity ───────────────────────────
    # ⚠ 81,920 faces stood in every run to date and carried NO law: `upload_surface_topology` was
    # called for the membrane alone. Measured 2026-08-24 on the tau record — two of eleven populations
    # carried force, and this is one of the nine that did not.
    env_topo = None
    if args.envelope_kappa_tilde_pn_um is not None:
        env_topo = upload_surface_topology(arena, cell.envelope)

    # ── the cortex, if the caller states its coefficients ─────────────────────────────────────────
    # ⚠ Both have NO default. k_xl is the constant this engine has already been bitten by — production
    # ran 8.2e5 pN/um against a field band of 0.1-100, because a Bell-Evans BINDING-BARRIER CURVATURE
    # was read as a structural spring constant (CROSSLINK_STIFFNESS_AXIS_2026-08-16.md). A default here
    # would put it back under a new name, in the binding that reaches every filament in the cell.
    ctx = None
    if (args.cortex_k_axial is None) != (args.cortex_kappa is None):
        raise SystemExit(
            "refused: --cortex-k-axial and --cortex-kappa go together. Binding an axial law without "
            "bending gives a cortex of floppy strings and binding bending without axial gives one that "
            "cannot hold its length; either alone is a different material, not a partial cortex.")
    if args.cortex_k_axial is not None:
        ctx = upload_strand_topology(arena, "cortex", cell.cortex)
    force_d = arena.node_arrays["force"]
    pos_d = arena.node_arrays["position"]

    # ── any other strand population the caller states coefficients for ────────────────────────────
    # Generic on purpose: strand_bind takes uniform-chain dicts and device-array objects alike, so
    # there is nothing population-specific left to write. What IS population-specific is the physics,
    # and that arrives from the caller or the bind refuses.
    extra_binds: list = []
    for spec in args.bind:
        if "=" not in spec:
            raise SystemExit(f"refused: --bind {spec!r} is not POP=k_axial[,kappa[,rest_um]]")
        pop, coeffs = spec.split("=", 1)
        parts = [float(x) for x in coeffs.split(",")]
        if not 1 <= len(parts) <= 3:
            raise SystemExit(f"refused: --bind {spec!r} takes k_axial[,kappa[,rest_um]]")
        # ⚠ kappa is OPTIONAL, and only for a population with no bending triple. `chromatin` is the
        # case: `build/chromatin` claims NO ANGLE3 range because its source models chromatin as
        # extensible springs of ZERO bending modulus, and `laws` therefore has no bending term to
        # give it. Until 2026-08-24 this parser required two numbers, so binding chromatin at all
        # meant passing a kappa that nothing would read — a magnitude in the record that governs
        # nothing, which is the shape this repository keeps finding at the bottom of a wrong result.
        # A kappa supplied for a population that HAS triples is still required; the check is below,
        # after the topology is known, so the refusal can name what it saw.
        k_ax = parts[0]
        kap = parts[1] if len(parts) >= 2 else None
        rest = parts[2] if len(parts) == 3 else None
        if pop not in built:
            raise SystemExit(
                f"refused: --bind names {pop!r}, which PHASE 1 did not build. Built: "
                f"{sorted(built)}. Skipping it silently would let a run report a binding it does "
                "not have.")
        t = upload_strand_topology(arena, pop, built[pop])
        if t.n_angles and kap is None:
            raise SystemExit(
                f"refused: --bind {pop}={coeffs} gave no kappa, but {pop} carries {t.n_angles:,} "
                f"bending triples. Binding the axial law alone would leave a population of floppy "
                f"strings whose bending force is silently zero — say the number or say why not.")
        if kap is not None and not t.n_angles:
            raise SystemExit(
                f"refused: --bind {pop}={coeffs} supplied kappa={kap}, but {pop} has NO bending "
                f"triple to launch it over. The value would be recorded and read by nothing.")
        extra_binds.append((pop, t, k_ax, kap, rest))

    def _zero_forces() -> None:
        force_d.zero_()

    def _accumulate() -> None:
        accumulate_area_tension(arena, topo, args.sigma_pn_per_um)
        accumulate_helfrich(arena, topo, args.kappa_tilde_pn_um)
        if env_topo is not None:
            # Bending only — see the flag's help for why no areal tension is applied here.
            accumulate_helfrich(arena, env_topo, args.envelope_kappa_tilde_pn_um)
        for _pop, _t, _k, _kap, _rest in extra_binds:
            # rest_um is not passed for the same reason as the cortex: where a builder supplied the
            # BUILT chord, overriding it pre-strains every segment. Where it did not, strand_bind
            # refuses rather than guessing, and the caller has to say.
            accumulate_axial(arena, _t, k_axial_pn_per_um=_k, rest_um=_rest)
            if _t.n_angles:
                accumulate_bending(arena, _t, kappa_pn_um2=_kap, seg_um=args.cortex_seg_um)
        if ctx is not None:
            # rest_um is NOT passed: cortex.py supplies the BUILT chord, which is shorter than the
            # requested arc step by the curvature of the shell. Overriding it would pre-strain every
            # segment, and strand_bind refuses the override for exactly that reason.
            accumulate_axial(arena, ctx, k_axial_pn_per_um=args.cortex_k_axial)
            accumulate_bending(arena, ctx, kappa_pn_um2=args.cortex_kappa,
                               seg_um=args.cortex_seg_um)

    # ── Sanity Gate, sign sense: a positive surface tension must pull a closed surface INWARD ──────
    # Run before any step. A sign flip here produces a run that looks fine and means the opposite, so
    # it aborts rather than being recorded.
    _zero_forces()
    _accumulate()
    wp.synchronize_device(dev)
    lo, hi = cell.membrane.nodes.lo, cell.membrane.nodes.hi
    p0 = pos_d.numpy()[lo:hi]
    f0 = force_d.numpy()[lo:hi]
    radial = np.einsum("ij,ij->i", f0, p0 / np.linalg.norm(p0, axis=1, keepdims=True))
    mean_radial_pn = float(radial.mean())
    if not mean_radial_pn < 0.0:
        raise RuntimeError(
            f"SIGN GATE FAILED: mean radial force on the membrane is {mean_radial_pn:+.6g} pN under a "
            f"POSITIVE surface tension of {args.sigma_pn_per_um} pN/um. Tension must pull a closed "
            "surface inward. A run past this point would look plausible and mean the opposite."
        )

    if args.mobility is not None and args.relax_um_per_pn:
        raise SystemExit(
            "refused: --mobility and --relax-um-per-pn together. One is a thermostat whose noise is "
            "fixed by fluctuation-dissipation and the other is a deterministic relaxation; running "
            "both makes a series that is neither, and its stationarity would mean nothing.")

    if args.emit_gamma and ctx is None:
        raise SystemExit(
            "refused: --emit-gamma without a bound cortex. Cortical tension is a property of the "
            "cortex, and computing it from a membrane-only force field would return a number with "
            "nothing cortical in it.")

    gamma_trace: list[float] = []
    thermo: list = []

    # ⚠ Hoisted, not built per step. `gamma_planes_device`'s docstring says its scratch is per-call and
    # that a caller running a long window should hoist it; these index arrays are the caller's half of
    # that, and rebuilding them each step would put back an allocation the device path exists to remove.
    _GD: dict = {}
    if args.gamma_device:
        _n_str, _n_per = int(cell.cortex.n_strands), int(cell.cortex.nodes_per_strand)
        _lo = int(cell.cortex.nodes.lo)
        _off = _lo + np.arange(_n_str + 1, dtype=np.int64) * _n_per
        _seg = segment_offsets(_off)
        _GD = {"fiber_offsets_d": wp.array(_off.astype(np.int32), dtype=wp.int32, device=dev),
               "seg_offsets_d": wp.array(_seg.astype(np.int32), dtype=wp.int32, device=dev),
               "n_segments": int(_seg[-1]), "n_actin": _lo + _n_str * _n_per,
               "R_um": float(cell.cortex.radius_um)}
        print(f"[phase4] gamma on DEVICE: {_GD['n_segments']:,} segments, "
              f"{_GD['n_segments'] / 1024:.0f} elements/block", flush=True)

    def _cortex_centroid(n_a: int) -> np.ndarray:
        """The cut centre, reduced on the device and derived from THIS population, not an index range.

        `wp.utils.array_sum` over 4.2 M vec3d costs 0.474 ms, against 101 MB of transfer per step to
        compute three numbers. Returned as host float64 because both estimators take it that way.
        """
        return np.asarray(wp.utils.array_sum(pos_d[:n_a]), np.float64) / float(n_a)

    def _gamma() -> dict[str, float]:
        """Cortical tension from the step's own force field, for StepContext.observables.

        ⚠ This is the seat (e) 1 needs and it is why the observables hook exists. Emitting gamma here
        does not decide the amendment, it makes the amendment DECIDABLE — and as of 2026-08-21 it HAS,
        so the sentence that used to stand here ("gamma is not emitted on the resting path") described
        the state before this function ran rather than after. The statistic, its window and the
        variance it is judged against remain the PI's.

        ⚠ And NO gamma value from this run is quotable. STATE.md (c) 3 and (c) 17 stand; the
        coefficients are test points; and session C's module reports against SEED SCATTER, not a
        single run's error bar, because (e) 1's D-2 calls within-run variance "the single most likely
        way to write the amendment and still be measuring nothing". One run cannot supply that.
        """
        if args.gamma_device:
            # ⚠ SAME QUANTITY, computed where the data already is. `gamma_from_resting_readback`
            # forms gamma_total as mean(|sum over families of plane sums|); this forms it from the
            # device plane sums the same way. The two agree to 3.87e-12 worst relative error at
            # NATIVE scale (4,128,840 segments, 4,032 elements/block) —
            # `world_gamma_native_equiv.py`, and to 6.63e-13 in the module's own gate.
            #
            # ⚠ WHY IT IS HERE AT ALL: measured 2026-08-21, the host path costs 6,143.57 ms/step
            # (178.66 readback + 5,964.91 estimator) against 10.77 ms for this one — 570x — and the
            # readback alone moves 576.0 MB, of which ~62% is unallocated arena capacity. A 9,000-step
            # run at 6.12 s/step does not fit a twelve-hour grant; at the device rate the same grant
            # buys thirty times the window.
            # ⚠ THE CENTRE IS NOT THE ORIGIN AND PASSING THE ORIGIN CHANGED THE ANSWER BY 15%.
            # `gamma_from_resting_readback`'s centre defaults to the actin centroid — measured here at
            # |c| = 2.05e-02 um, twenty nanometres off origin — and a centre offset does not perturb
            # the sum, it flips which elements CROSS a plane. Discrete, not round-off, which is why an
            # A/B on 12 steps showed 7.5e-06 at step 0 growing to 1.5e-01, against an identical-input
            # kernel equivalence of 3.9e-12.
            #
            # Reduced on the device — `wp.utils.array_sum` over 4.2 M vec3d is 0.474 ms — rather than
            # read back, because a 101 MB transfer per step to compute three numbers would put back
            # most of what the device path exists to remove.
            n_a = _GD["n_actin"]
            centre = _cortex_centroid(n_a)
            res = gamma_planes_device(
                pos_d=pos_d, f_ext_d=force_d,
                fiber_offsets_d=_GD["fiber_offsets_d"], seg_offsets_d=_GD["seg_offsets_d"],
                n_segments=_GD["n_segments"], n_actin=n_a,
                R_um=_GD["R_um"], centre=centre, device=dev)
            total = float(np.mean(np.abs(sum(np.asarray(v, np.float64) for v in res.values()))))
            gamma_trace.append(total)
            return {"gamma_pn_per_um": total}

        n_str, n_per = int(cell.cortex.n_strands), int(cell.cortex.nodes_per_strand)
        lo = int(cell.cortex.nodes.lo)
        offsets = lo + np.arange(n_str + 1, dtype=np.int64) * n_per
        g = gamma_from_resting_readback(
            R_um=float(cell.cortex.radius_um),
            pos=pos_d.numpy(), f_ext=force_d.numpy(),
            n_actin=lo + n_str * n_per, fiber_offsets=offsets,
            # ⚠ EXPLICIT, and neither path takes a default. The estimator's default names the centre
            # population by an INDEX RANGE (`p[:n_actin]`), which is the cortex only because
            # build_all happens to claim the cortex at 0; a different claim order silently averages
            # someone else's nodes into an "actin centroid" — measured elsewhere at 68% wrong. Two
            # instruments each taking their own default is a comparison of defaults, not of kernels,
            # and that is what made the first A/B read as a 15% kernel disagreement.
            centre=_cortex_centroid(lo + n_str * n_per),
            # ⚠ Explicit and empty, never defaulted. 2026-07-29 measured that masking the SOLVE and
            # not the MEASUREMENT left the dynamics bit-identical while gamma differed by a constant
            # 0.41526 pN/um — the estimator was reading a different force field from the one being
            # integrated. A default would make that silent again.
            force_mask=())
        gamma_trace.append(float(g.gamma_total_pn_per_um))
        return {"gamma_pn_per_um": float(g.gamma_total_pn_per_um)}

    def _residual() -> tuple[float, float]:
        """Assembled |sum F| over every live node [pN], and the tolerance it is judged against.

        ⚠ The tolerance returned is the Higham summation bound, and it is returned so a caller CAN
        compare — not so this driver calls the comparison convergence. Nothing here reads it.
        """
        f = force_d.numpy()[: n_live_nodes]
        resultant = float(np.linalg.norm(f.sum(axis=0)))
        n_terms = int(f.size)
        return resultant, n_terms * float(np.finfo(np.float64).eps) * max(
            float(np.abs(f).sum()), 1.0)

    def _solve() -> None:
        """The inner mechanical solve: accumulate, then an explicitly DECLARED relaxation.

        With ``--relax-um-per-pn 0`` nothing moves, which is the honest default — the builders leave
        ``mobility`` at zero because a mobility is a drag law's parameter and a builder filling it in
        would be inventing physics. A nonzero gain is a numerical relaxation and is labelled one.
        """
        _zero_forces()
        _accumulate()
        if args.relax_um_per_pn:
            wp.launch(_relax_kernel, dim=n_live_nodes,
                      inputs=[pos_d, force_d, wp.float64(float(args.relax_um_per_pn))], device=dev)
        elif args.mobility is not None:
            # Driven AND damped, which is the only arm that can produce a series a stationarity
            # criterion has anything to say about — see GAMMA_SEED_SCATTER_REFUSED_TWICE, where the
            # frozen arm is degenerate and the relaxed arm drifts at 4e+05 sigma.
            thermo.append(langevin_step(
                arena, 0, n_live_nodes, mobility_um_per_pn_s=args.mobility,
                temperature_k=args.temperature_k, dt_s=args.dt_phys_s,
                base_seed=args.seed, step_index=len(thermo)))
        wp.synchronize_device(dev)

    # ── optional position sequence, for the native viewer ─────────────────────────────────────────
    snap_frames: list[tuple[int, dict]] = []

    def _snapshot(step_i: int) -> None:
        """Capture positions for one frame, and WRITE A HEARTBEAT so the run is observable.

        ⚠ **The frames still accumulate in memory and the record is still written at the end**, which
        means a run killed at the wall clock loses them. That was true before this change and it is
        true after it; what changes is that the heartbeat lands on disk immediately, so *how far a
        running job has got* is answerable.

        **Why that mattered enough to fix mid-flight.** On 2026-08-21 a 9,000-step run was launched
        against a grant sized from a 4.06 s/step measurement, and an hour later there was **no way to
        tell whether it was on that rate** — no output, no file, nothing but a process. A second
        timing of the same quantity said 6.14 s/step, under which the run would be killed at ~6,672
        steps having written nothing at all. The contract for that run stated that snapshots were
        written periodically. **They were not, and I had written that sentence without checking it.**
        """
        host = pos_d.numpy()
        snap_frames.append((step_i, {n: host[lo:hi].astype(np.float32) for n, (lo, hi) in (
            (nm, (o["claims"]["node"] if isinstance(o, dict) else (o.nodes.lo, o.nodes.hi)))
            for nm, o in pops.items())}))
        _heartbeat(step_i)

    def _heartbeat(step_i: int) -> None:
        """Write the progress file. Cheap, atomic, and DECOUPLED from snapshots.

        ⚠ These were one flag until 2026-08-21, and the coupling made the instrument useless for the
        thing it was added for. A snapshot costs 55 MB of host RAM, so its interval is set by memory —
        250 steps for a 5,250-step run. **That put the first rate reading twenty-seven minutes into a
        run whose sizing was the open question.** A heartbeat costs a few hundred bytes; there is no
        reason for it to inherit the snapshot's budget.

        Time-based rather than step-based on purpose: the step cost is what is unknown, so a
        step-based interval writes at an unknown rate, which is the same defect one level down.
        """
        beat = args.out.with_suffix(".progress.json")
        done = time.perf_counter() - t_steps if step_i else 0.0
        tmp = beat.with_suffix(".tmp")
        tmp.parent.mkdir(parents=True, exist_ok=True)
        tmp.write_text(json.dumps({
            "steps_done": step_i, "steps_requested": args.steps,
            "elapsed_s": round(done, 1),
            "s_per_step": round(done / step_i, 4) if step_i else None,
            "eta_s": round(done / step_i * (args.steps - step_i), 1) if step_i else None,
            "frames_held_in_memory": len(snap_frames),
            "⚠": ("frames are held in RAM and written with the record at the end; a run killed at "
                  "the wall clock loses them. This file is a heartbeat, not a checkpoint."),
        }, indent=1, ensure_ascii=False))
        tmp.replace(beat)                    # atomic: a reader never sees a half-written heartbeat

    ws = WorldStep(arena=arena, predicate=UndefinedAcceptance())
    t_steps = time.perf_counter()
    _last_beat = [t_steps]
    if args.snapshot_every:
        _snapshot(0)
    for _step in range(args.steps):
        ws.step(dt_phys_s=args.dt_phys_s, solve=_solve, assemble_residual=_residual,
                observables=(_gamma if args.emit_gamma
                             else (lambda: {"mean_radial_pn": mean_radial_pn})))
        if args.snapshot_every and (_step + 1) % args.snapshot_every == 0:
            _snapshot(_step + 1)
        elif time.perf_counter() - _last_beat[0] >= args.heartbeat_s:
            _last_beat[0] = time.perf_counter()
            _heartbeat(_step + 1)
    wp.synchronize_device(dev)
    steps_s = time.perf_counter() - t_steps

    from aleph.scripts.run_provenance import stamp

    record = {
        "schema": "ffn-world-phase4-native@1",
        "kind": "diagnostic",
        # ⚠ Stamp the code that actually RAN, not the code someone believes ran. The GPU tree is
        # hand-deployed, so a commit alone can be hours stale — 2026-08-21 measured exactly that,
        # and found it only by hashing the remote tree afterwards. The closure digest is over the
        # imported bytes, so it cannot go stale the way a commit can.
        "provenance": stamp(__file__),
        "quantitative_claim_status": "BLOCKED",
        "device": str(dev),
        "build_s": round(build_s, 3),
        "steps_s": round(steps_s, 3),
        "s_per_step": round(steps_s / args.steps, 4),
        "n_nodes": n_live_nodes,
        "cell": ("core only: cortex, membrane, nuclear_envelope" if args.core_only else
                 f"FULL: cortex, membrane, nuclear_envelope + {len(built)} more, from "
                 "world/populations.py — the same builder PHASE 1 and the renderer use"),
        "populations": sorted(pops),
        # ⚠ What is NOT in this cell, derived from the builder set rather than remembered. The tau
        # runs this driver produced stand eleven populations and `build/nmii.py` is not among the
        # builders that made them -- THERE IS NO MOTOR IN THIS CELL -- and `build/cytosol.py` never
        # ran, which is why `census.live.grid_cell` is 0. Neither fact was in the record; a reader had
        # to know the population list by heart. `populations` alone cannot carry an absence.
        "builders_not_standing": list(builders_not_standing(pops)),
        "builders_not_standing_note": (
            "Builder modules under world/build/ that contributed NO population to this cell. A name "
            "here means the module did not run, never that a population is empty -- absent and "
            "measured-zero are different findings and may not be summed, averaged or plotted alike."),
        "census": arena.census(),
        "placement_check": placement,
        "placement_envelope": fp.as_record(),
        "sign_gate": {
            "mean_radial_force_pn": mean_radial_pn,
            "asserted": "a POSITIVE surface tension pulls a closed surface INWARD",
            "verdict": "PASS",
        },
        "envelope_binding": (
            None if env_topo is None else {
                "population": "nuclear_envelope",
                "n_faces": int(env_topo.n_faces) if hasattr(env_topo, "n_faces") else None,
                "kappa_tilde_pn_um": float(args.envelope_kappa_tilde_pn_um),
                "law": "Helfrich bending only",
                "areal_tension": "NOT BOUND — bilinear in areal strain (chromatin below a ~10% knee, "
                                 "lamin-A/C above). A constant sigma on a shell built at rest is an "
                                 "imposed prestress, and tension must EMERGE.",
                "coefficient_basis": "laws.nucleus_envelope.KAPPA_NE_PN_UM = 0.0828 pN·um, ~20 kBT "
                                     "(lamina + bilayer). ⚠ not measured on MCF7.",
            }),
        "extra_bindings": [
            {**t.record(), "k_axial_pn_per_um": k, "kappa_pn_um2": kap,
             "rest_um": rest,
             "coefficient_basis": "TEST POINTS stated by the caller; no default exists for any of them"}
            for pop, t, k, kap, rest in extra_binds],
        "cortex_binding": (None if ctx is None else {
            **ctx.record(),
            "k_axial_pn_per_um": args.cortex_k_axial,
            "kappa_pn_um2": args.cortex_kappa,
            "coefficient_basis": "TEST POINTS stated by the caller. NOT sourced values for this cell "
                                 "— see CROSSLINK_STIFFNESS_AXIS_2026-08-16.md, where the production "
                                 "k_xl of 8.2e5 pN/um is a Bell-Evans binding-barrier curvature read "
                                 "as a structural spring constant against a field band of 0.1-100.",
            "rest_length_basis": "the BUILT chord from cortex.py, not the requested arc step",
        }),
        "coefficients": {
            "sigma_pn_per_um": args.sigma_pn_per_um,
            "kappa_tilde_pn_um": args.kappa_tilde_pn_um,
            "basis": "TEST POINTS, not sourced values for this cell. They set the force scale, so no "
                     "residual from this run is a physical residual.",
            "relax_um_per_pn": args.relax_um_per_pn,
            "relax_basis": "a DECLARED numerical relaxation, not a drag law. The builders leave the "
                           "arena's per-node mobility at zero deliberately (membrane.py:49-50).",
            "dt_phys_s": args.dt_phys_s,
            "dt_phys_basis": "this driver's clock, not a physiological quantity",
        },
        "thermostat": (None if not thermo else {
            **thermo[-1].record(),
            "n_steps": len(thermo),
            "stability_bound_s": (stability_bound_s(args.mobility, args.cortex_k_axial)
                                  if args.cortex_k_axial else None),
            "stability_note": "REPORTED, never enforced. dt < 1/(mobility * k_max) is what the "
                              "stiffest bond implies; a module that clamped the caller's dt would be "
                              "choosing a timestep, and a timestep chosen so a run behaves is what "
                              "the charter forbids.",
        }),
        "gamma": (None if not args.emit_gamma else {
            "trace_pn_per_um": gamma_trace,
            "NO_MAGNITUDE_CLAIMED": True,
            "why_emitted": "A stationarity gate cannot be written against an observable the run does "
                           "not produce, and cortical tension is that observable. Emitting it does "
                           "not decide STATE.md (e) 1 — it makes the amendment decidable, which it "
                           "was not before 2026-08-21. What is still open is the PI's: WHICH "
                           "statistic, over WHAT window, judged against WHAT variance.",
            "why_not_quotable": "(c) 3 and (c) 17 retire every gamma this engine has produced; the "
                                "coefficients here are declared test points; and one run cannot "
                                "supply the SEED SCATTER (e) 1's D-2 requires — within-run variance "
                                "is called 'the single most likely way to write the amendment and "
                                "still be measuring nothing'.",
            "seed": args.seed,
            "is_an_independent_replicate": seed_reaches_build,
            "replicate_note": ("build_all forwards the seed, so runs at different --seed are "
                               "independent replicates" if seed_reaches_build else
                               "⚠ build_all does NOT forward a seed, so repeating this run at a "
                               "different --seed produces the SAME cell. These are not replicates "
                               "and must not be combined as a seed scatter."),
            "next": "three seeds through observe_gamma.observe_gamma_seed_scatter, then the PI fixes "
                    "the statistic, the window and the variance it is judged against",
        }),
        "step": ws.as_record(),
        "not_a_claim": (
            "NO step in this run was accepted, and none was rejected — this engine has no acceptance "
            "criterion. Both it has shipped are one-sided and both were measured: the balance gate is "
            "the Higham summation bound and cannot fail for physics; descent_ratio read a "
            "post-rollback residual and cannot succeed. STATE.md (e) 1 is UNDECIDED — no longer "
            "unwritable, since cortical tension is emitted on this path as of 2026-08-21, but the "
            "statistic, the window and the variance it is judged against are all still the PI's. "
            "Every step here therefore records ACCEPTANCE_UNDEFINED, and is_accepted is False for it."
        ),
    }
    if snap_frames:
        # The container is world_export_cell's, with a `frame` field on each positions block. A
        # single-frame file is byte-identical in shape to a static export, so an older viewer still
        # opens one and a newer one plays a sequence — the format did not fork.
        from aleph.scripts.world_export_cell import write_sequence

        seq = args.out.with_suffix(".alephcell")
        n_bytes = write_sequence(seq, pops, snap_frames, header_extra={
            "source": "world_phase4_native.py",
            "dt_phys_s": args.dt_phys_s,
            "snapshot_every": args.snapshot_every,
            "time_axis": ("⚠ frame index * dt_phys. With an unsourced mobility the dt is a step "
                          "index wearing a unit, so read FRAMES, not seconds."),
            "physics_claim": "NONE. Positions moved under a declared relaxation or thermostat; no "
                             "step in this sequence was accepted.",
        })
        record["sequence"] = {"path": str(seq), "n_frames": len(snap_frames),
                              "bytes": n_bytes, "steps": [i for i, _ in snap_frames]}
        print(f"[phase4] wrote {len(snap_frames)} frames, {n_bytes / 1e6:.0f} MB -> {seq}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1))
    print(f"[phase4] {args.steps} steps, {steps_s / args.steps:.4f} s/step, "
          f"{n_live_nodes:,} nodes, sign gate PASS, "
          f"{ws.as_record()['n_accepted']} accepted (of {args.steps}), wrote {args.out}")
    return 0


try:                                          # pragma: no cover - import-time on the GPU host only
    import warp as wp

    @wp.kernel
    def _relax_kernel(pos: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
                      gain: wp.float64):
        """One explicit overdamped move, ``x += gain * F``.

        Not a drag law and not claiming to be one: ``gain`` is a declared numerical relaxation and the
        record labels it that way. A real mobility is ``1/gamma`` from a drag law nobody has bound
        here yet, which is why the arena's per-node ``mobility`` is still zero.
        """
        i = wp.tid()
        pos[i] = pos[i] + gain * force[i]
except Exception:                             # pragma: no cover - dev machine without warp
    _relax_kernel = None  # type: ignore[assignment]


if __name__ == "__main__":
    raise SystemExit(main())
