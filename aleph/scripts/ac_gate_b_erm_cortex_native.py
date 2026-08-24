#!/usr/bin/env python
r"""GATE-B NATIVE driver — emergent membrane--cortex ERM tether tension on the FULL native cell.

The SECOND dynamic-runtime participant (after the NMII cortex-motor slice): the ERM (ezrin/radixin/moesin)
membrane--cortex clutch.  This demonstrates the GATE-B MECHANISM on the real composed cell (CUDA / gbook
A5000):

    resting turgor pushes the membrane outward  →  the ERM tethers (membrane↔cortex) stretch and develop
    tension that rides on the actomyosin cortex  →  the most-loaded tethers SLIP-detach (Bell 1978; ERM is a
    pure slip bond) while shed pairs rebind within capture  →  a dynamic steady-state bound fraction + an
    EMERGENT membrane--cortex tether tension.

This extends the actomyosin cortex OUTWARD toward the membrane.  The tether tension EMERGES from the bound
population (no imposed prestress), and under load the Bell-slip shedding IS the mechanistic bleb-onset — it
replaces the hard ``f_rupt`` membrane-energy threshold with the force-dependent single-molecule detachment.
The slice (:class:`aleph.engine.erm_cortex_slice.ERMCortexSlice`) drives
:class:`aleph.engine.transaction.CellTransaction`: each accepted step snapshots the ERM SoA + surface
positions → (near-no-op propose: pairs are fixed) → inner mechanical relax (the membrane rides on the cortex
through the tethers under turgor) → commit the accepted-gated ERM Bell KMC on the converged geometry → advance
the device clock.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (this needs a CUDA GPU; it will NOT run on the dev Mac):

    ssh gbook
    cd ~/ffn_ac_native
    # DELIVERABLE — force-balanced resting baseline (preload ON, default) + turgor ramp so tension EMERGES:
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_erm_cortex_native.py --steps 40 --dt 0.01 --outer 40 \
      --membrane-subdiv 8 --turgor-ramp 1.0

    # OLD un-preloaded control (non-converged, geometry rolled back, tension≡0) kept for the record:
    #   … ac_gate_b_erm_cortex_native.py --steps 40 --dt 0.01 --outer 40 --no-preload

    # full RESTING-production ERM density (membrane subdiv 8 = 655,362 tethers) + a k_off0 sweep:
    for KO in 0.1 0.3 1.0; do PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_erm_cortex_native.py --steps 40 --dt 0.01 --outer 40 \
      --membrane-subdiv 8 --k-off0 $KO ; done

Long run: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env python).
────────────────────────────────────────────────────────────────────────────────────────────────────────

WHAT TO EXPECT (with the PROVISIONAL Bell kinetics below — a MECHANISM demonstration, NOT quantitative
production; every ERM on/off magnitude is a PI decision, absent from the Contract-Graph 2026-07-21 audit):
  * with the RESTING PRELOAD ON (default), each ERM rest length is pre-stretched so its pre-tension holds the
    membrane node's turgor: the resting baseline is FORCE-BALANCED (the inner solve converges, accepted geometry
    PERSISTS) and every bound tether already carries its physiological resting tension ~turgor/node (≈ 0.6 pN at
    subdiv 6).  This is the baseline the ramp perturbs FROM (PI physiological-baseline rule);  --no-preload skips
    it and reproduces the OLD non-converged behaviour (geometry rolled back → length≡rest → tension≡0);
  * under the --turgor-ramp the membrane is pushed further outward, the tethers stretch (length > rest), and the
    per-tether tension ``k_erm·(length − rest)`` RISES from the resting baseline (k_erm = 4600 pN/µm is the
    SOURCED Braunger single-bond stiffness);
  * `bound%` falls from 100% and settles to a dynamic steady state as the Bell SLIP off-rate ``k_off0·
    exp(F/F0)`` sheds the most-loaded tethers, balanced by capture-gated rebinding — an EMERGENT duty, not an
    imposed fraction.  Where local tension is highest the shedding is fastest: that is the bleb-onset signature;
  * the emergent membrane--cortex tether tension (mean per-tether [pN], total [pN], areal traction [pN/µm²])
    rises with the load and then rides on the bound population;
  * max|PF| (the projected residual) stays finite and the membrane tracks the cortex each step.

FORCE-ACCEPT: the demo force-accepts every step (`accepted_d = ones`) so the KMC accumulates and the mechanism
is visible; a resting-baseline PRODUCTION run would instead gate on the force-balance predicate.

NO DOUBLE COUNT: the built `MembraneCompartment` owns an identical ERM force+KMC path.  This driver EXTRACTS
the tether SoA into the slice's connector, then sets `cell.membrane.n_erm = 0` so the membrane compartment
stops touching ERM, and wraps `cell.membrane` in a thin adapter so the driver's own inner relax adds the
Helfrich+area force (real membrane) + the slice's ERM tether force (one SoA, single owner) — exactly how the
NMII GATE-B driver routes its motor through `_accumulate_all` via `_SliceMyosinAdapter`.

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac (no CUDA) — the Lead runs it on gbook.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import R_CELL_UM, CellConfig, build_cell
from aleph.components.incumbent.driver import (
    _accumulate_all,
    _preload_cortex_pretension,
    _preload_erm_resting_balance,
    make_inner_solve,
)
from aleph.engine.erm_cortex_slice import (
    ERMSliceParams,
    build_erm_cortex_connector,
    build_erm_cortex_slice,
    build_native_surface_owner,
)
from aleph.components.motor.bell_kinetics_analytic import bell_f0_from_x_beta

# ── PROVISIONAL ERM Bell kinetics — NOT SOURCED, PI/KB GAP ─────────────────────────────────────────────
# The Contract-Graph has NO MCF7 ERM on/off datum (k_on / k_off0 / F0 / capture are all absent — 2026-07-21
# audit).  These proxy values only give the physiological SHAPE (a bound tether sheds faster the harder it is
# pulled — Bell 1978 slip) so the mechanism path runs; they are MECHANISM-DEMO PROXIES, never sourced
# magnitudes, and must be sourced by PI before any quantitative claim.  ``k_erm`` is DIFFERENT: it is the
# SOURCED Braunger 2014 single-bond stiffness (4.6 pN/nm = 4600 pN/µm, PI-ratified 2026-07-21), read from the
# built membrane, not invented here.  (k in 1/s, x/capture in µm, F0 in pN.)
ERM_KON_PROVISIONAL = 1.0          # ERM rebind rate [1/s]  — ezrin residence ~seconds (slower than myosin)
ERM_KOFF0_PROVISIONAL = 0.3        # ERM zero-force off-rate [1/s]
ERM_XBETA_PROVISIONAL_UM = 1.0e-3  # ERM Bell bond length x_beta [µm] → F0 = kBT/x_beta ≈ 4.28 pN
ERM_CAPTURE_PROVISIONAL_UM = 0.05  # ERM rebind capture radius [µm]  (~ the tether length scale)

_FLUID_CODE = 1                    # FieldGrid mask: OUTSIDE=0 / FLUID=1 / NUCLEUS=2 (ac.fluid.field_grid)


@wp.kernel
def _ramp_turgor_fluid_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    d_pa: wp.float64,
) -> None:
    """Raise the interior osmotic-turgor setpoint by ``d_pa`` [Pa == pN/µm²] on every FLUID cell (device-only).

    This is the LOAD RAMP.  ``p`` is the same conservative field the membrane traction reads at each live
    triangle centroid (``MembranePressureTraction``, mask-aware), so bumping the FLUID interior raises the
    outward pressure the inner solve must balance — the membrane genuinely bulges against the ERM tethers,
    which then develop ``k_erm·(length − rest)`` tension.  Nothing is injected into the tether; the tension is
    the real reaction of the membrane moving.  OUTSIDE/NUCLEUS cells are untouched (no spurious body gradient).
    """
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID_CODE:
        p[i, j, k] = p[i, j, k] + d_pa


def _banner(params: ERMSliceParams, x_beta_um: float, n_erm: int) -> None:
    """Print the loud PROVISIONAL-PARAMS banner listing every value + its PI-GAP status."""
    gap = "PROVISIONAL PROXY — NOT SOURCED, PI/KB GAP (mechanism-shape only; source before any claim)"
    rows = [
        ("k_erm  [pN/µm]", params.k_erm, "SOURCED — Braunger 2014 single-bond 4.6 pN/nm (PI-ratified 2026-07-21)"),
        ("k_on   [1/s]  ", params.k_on, gap),
        ("k_off0 [1/s]  ", params.k_off0, gap),
        ("F0     [pN]   ", params.bell_force_pn, f"{gap}  (= kBT/x_beta, x_beta={x_beta_um*1e3:.2f} nm — GAP)"),
        ("capture[µm]   ", params.capture_radius_um, gap),
    ]
    print("=" * 104)
    print("  GATE-B PROVISIONAL PARAMETERS  —  ERM (ezrin/radixin/moesin) membrane–cortex clutch, pure Bell SLIP")
    print("  MECHANISM DEMONSTRATION, NOT quantitative production — every Bell on/off magnitude is a PI decision.")
    print(f"  {n_erm:,} ERM tethers (one per membrane node), all BOUND at t0 (force-free formation).")
    print("-" * 104)
    for label, value, status in rows:
        print(f"    {label} = {float(value):>10.5g}   {status}")
    print("=" * 104, flush=True)


def _tether_stats(cell, slc) -> dict:
    """Emergent membrane--cortex tether tension telemetry (out-of-hot-loop host readback; I0-A honoured).

    Per bound tether the tensile load is ``k_erm·max(length − rest, 0)`` [pN] (the same unilateral law the
    device force kernel scatters).  Reports the bound fraction plus the emergent tension as mean-per-tether
    [pN], total [pN], and an areal traction [pN/µm²] over the membrane area — the tension EMERGES from the
    bound population, no imposed prestress.
    """
    pos = cell.pos_d.numpy()
    mi = slc.connector.membrane_idx_d.numpy()
    ci = slc.connector.cortex_idx_d.numpy()
    bound = slc.connector.bound_d.numpy().astype(bool)
    rest = slc.connector.rest_d.numpy()
    k_erm = float(slc.connector.k_erm)
    length = np.linalg.norm(pos[mi] - pos[ci], axis=1)
    tension = np.where(bound, np.maximum(k_erm * (length - rest), 0.0), 0.0)  # per-tether [pN]
    n_erm = int(mi.shape[0])
    n_bound = int(bound.sum())
    area = 4.0 * np.pi * R_CELL_UM * R_CELL_UM
    return {
        "n_erm": n_erm, "n_bound": n_bound,
        "bound_pct": 100.0 * n_bound / max(n_erm, 1),
        "mean_pn": float(tension[bound].mean()) if n_bound else 0.0,
        "total_pn": float(tension.sum()),
        "areal_pn_per_um2": float(tension.sum() / area),
        "detach_events": int(slc.connector.detach_events_d.numpy()[0]),
        "attach_events": int(slc.connector.attach_events_d.numpy()[0]),
    }


class _SliceMembraneAdapter:
    """Replaces ``cell.membrane`` so the driver's inner relax adds the SLICE's ERM force (one SoA, no double).

    ``_accumulate_all`` calls ``cell.membrane.accumulate(pos, f)`` each inner iteration; this adapter runs the
    REAL membrane compartment's Helfrich + area force (its ERM is disabled via ``n_erm = 0``) and then the
    slice's two-array ERM tether over the SAME global ``pos``/``f`` (GLOBAL membrane/cortex indices — the split
    is logical).  Any other MembraneCompartment attribute the driver/solver reads is delegated to the original.
    """

    def __init__(self, slc, original) -> None:
        self._original = original           # real MembraneCompartment, its own ERM disabled (n_erm = 0)
        self._connector = slc.connector

    def __getattr__(self, name: str):
        # only reached when normal lookup fails (accumulate/_original/_connector are found first).
        return getattr(object.__getattribute__(self, "_original"), name)

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        self._original.accumulate(pos, f)                       # Helfrich + area (n_erm == 0 → no ERM here)
        self._connector.accumulate_pair(pos, f, pos, f)         # slice ERM tether: membrane −f / cortex +f


def _params(args) -> ERMSliceParams:
    """Assemble the ERM slice params: sourced k_erm (from CLI/Braunger) + the provisional Bell proxies."""
    f0 = bell_f0_from_x_beta(args.x_beta)  # F0 = kBT/x_beta [pN]; x_beta is the ERM GAP proxy
    return ERMSliceParams(
        k_erm=args.k_erm, k_on=args.k_on, k_off0=args.k_off0,
        bell_force_pn=f0, capture_radius_um=args.capture,
        source=("k_erm=Braunger2014(4.6pN/nm,PI-ratified); "
                "k_on/k_off0/F0/capture=PROVISIONAL PI-GAP proxy (no MCF7 ERM on/off datum in Contract-Graph)"),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="GATE-B: emergent membrane–cortex ERM tether tension (native).")
    ap.add_argument("--steps", type=int, default=40, help="number of ACCEPTED physical steps to run")
    ap.add_argument("--dt", type=float, default=0.01, help="physical timestep dt_phys [s] (attach/detach tick)")
    ap.add_argument("--outer", type=int, default=40, help="inner mechanical relax iterations per step (n_inner)")
    ap.add_argument("--filaments", type=int, default=70686, help="cortex F-actin count (native = 70686)")
    ap.add_argument("--membrane-subdiv", type=int, default=6, dest="membrane_subdiv",
                    help="membrane icosphere level = ERM tether count (6=40,962 dynamic; 8=655,362 resting prod)")
    ap.add_argument("--k-erm", type=float, default=4.6e3, dest="k_erm",
                    help="ERM linker stiffness [pN/µm] (SOURCED Braunger 2014 = 4600; overrides the membrane's)")
    ap.add_argument("--k-on", type=float, default=ERM_KON_PROVISIONAL, dest="k_on",
                    help="ERM rebind rate [1/s] (PROVISIONAL PI-GAP proxy)")
    ap.add_argument("--k-off0", type=float, default=ERM_KOFF0_PROVISIONAL, dest="k_off0",
                    help="ERM zero-force off-rate [1/s] (PROVISIONAL PI-GAP proxy)")
    ap.add_argument("--x-beta", type=float, default=ERM_XBETA_PROVISIONAL_UM, dest="x_beta",
                    help="ERM Bell bond length x_beta [µm] → F0 = kBT/x_beta (PROVISIONAL PI-GAP proxy)")
    ap.add_argument("--capture", type=float, default=ERM_CAPTURE_PROVISIONAL_UM, dest="capture",
                    help="ERM rebind capture radius [µm] (PROVISIONAL PI-GAP proxy)")
    ap.add_argument("--turgor-ramp", type=float, default=0.0, dest="turgor_ramp",
                    help="LOAD RAMP: raise the interior turgor Π by this many Pa PER ACCEPTED STEP "
                         "(0.0 = constant-load baseline, bit-identical to before). >0 bulges the membrane "
                         "against the ERM tethers so k_erm·(length−rest) tension EMERGES and the most-loaded "
                         "tethers shed by the force-driven Bell slip = bleb-onset. e.g. 1.0 ramps Π 40→80 Pa "
                         "(~2× resting) over 40 steps.")
    ap.add_argument("--preload", action=argparse.BooleanOptionalAction, default=True, dest="preload",
                    help="establish the FORCE-BALANCED resting baseline before the loop via the SAME preloads "
                         "the production driver uses (run_from_resting, driver.py:1011-1014): "
                         "_preload_erm_resting_balance (pre-stretch each ERM rest so its pre-tension cancels the "
                         "membrane node's turgor → the inner solve CONVERGES → accepted geometry PERSISTS instead "
                         "of being rolled back) + _preload_cortex_pretension(--cortex-prestrain). DEFAULT ON — this "
                         "is the correct physiological baseline (PI physiological-baseline rule). --no-preload "
                         "reproduces the OLD un-preloaded behaviour (relaxed cortex + positive turgor = no "
                         "equilibrium → non-converged → geometry rolled back → length≡rest → tension≡0) for the record.")
    ap.add_argument("--cortex-prestrain", type=float, default=0.0, dest="cortex_prestrain",
                    help="cortex crosslink pre-strain fed to _preload_cortex_pretension [dimensionless]. DEFAULT "
                         "0.0 = the value the GATE-A native clean-turgor baseline uses (ac_resting_go_probe.py "
                         "--native: preload=True, prestrain=0.0): with myosin OFF the structural-fix cortex "
                         "SELF-BALANCES the ERM reaction (clean cortex max|F|=0.036 pN, 0 nodes over the 0.21 gate) "
                         "so NO Laplace hoop pre-stress is applied. A non-zero value is the resting-MYOSIN-substitute "
                         "static prestress — a GATE-A P0.3 PI modeling decision, NOT used for this myosin-off baseline.")
    ap.add_argument("--seed", type=int, default=0, help="base RNG seed")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"GATE-B ERM native driver needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    print(f"[gate-b-erm] device={dev}  steps={args.steps}  dt={args.dt}s  inner={args.outer}  "
          f"membrane_subdiv={args.membrane_subdiv}  (ERM = pure Bell SLIP; shedding-under-load = bleb-onset)",
          flush=True)

    # 1. Build the FULL native cell — membrane + cortex, radial ERM pairing (tethers transmit the turgor load).
    t0 = time.time()
    cfg = CellConfig(
        n_filaments=args.filaments, with_myosin=False, overlap_free_cortex=True,
        with_membrane=True, membrane_subdivisions=args.membrane_subdiv, nucleus_subdivisions=3,
        with_pressure=True, erm_radial_pairing=True, seed=args.seed,
    )
    cell = build_cell(cfg)
    membrane = cell.membrane
    if membrane is None or int(membrane.n_erm) == 0:
        raise RuntimeError("no ERM tethers: build with with_membrane=True and erm_radial_pairing=True")
    n_erm = int(membrane.n_erm)
    print(f"[gate-b-erm] built cell in {time.time() - t0:.1f}s: n_actin={cell.n_actin} n_total={cell.n_total} "
          f"n_erm={n_erm} (membrane subdiv {args.membrane_subdiv}) k_erm={float(membrane.k_erm):.1f} pN/µm",
          flush=True)

    # 1b. RESTING BASELINE — the PHYSIOLOGICAL-BASELINE fix (PI rule).  build_cell alone leaves a RELAXED cortex
    #     under POSITIVE turgor: every ERM rest length == its force-free formation length, so the membrane carries
    #     the full ~turgor/node UNBALANCED — there is no force-balanced equilibrium, the inner solve NEVER converges
    #     (converged=0, residual ~pN), and the driver ROLLS BACK every accepted step's geometry
    #     (conditional_rollback_vec3_kernel).  Membrane then can't move → length ≡ rest → tether tension ≡ 0.
    #     These are the SAME two preloads the production driver runs (run_from_resting, driver.py:1011-1014).  Must
    #     run BEFORE the connector extracts erm_rest_d and BEFORE membrane.n_erm := 0 (the preload needs n_erm > 0
    #     and the live membrane ERM in _accumulate_all).  Order matches production: cortex pretension, then ERM.
    if args.preload:
        pre_cortex = _preload_cortex_pretension(cell, args.cortex_prestrain)   # cortex hoop pre-stress (0.0 = off)
        pre_erm = _preload_erm_resting_balance(cell)                           # ERM pre-tension holds the turgor
        print(f"[gate-b-erm] RESTING PRELOAD ON (physiological baseline): cortex_prestrain={args.cortex_prestrain:g} "
              f"({pre_cortex}); ERM balance n_erm={pre_erm['n_erm']} residual "
              f"{pre_erm['res_before']:.4g}→{pre_erm['res_after']:.4g} pN/node — the ERM tethers now pre-stretch to "
              f"hold turgor (each carries ~turgor/node ≈ its resting tension), the inner solve CONVERGES, and "
              f"accepted geometry PERSISTS (no rollback).  With --turgor-ramp the tension rises FROM this baseline.",
              flush=True)
    else:
        print("[gate-b-erm] RESTING PRELOAD OFF (--no-preload): relaxed cortex + positive turgor = NO force-balanced "
              "equilibrium → inner solve does NOT converge → geometry rolled back each step → length≡rest → "
              "tether tension≡0.  This reproduces the OLD un-preloaded behaviour for the record; not a valid baseline.",
              flush=True)

    # 2. Extract the ERM SoA → the connector (aliasing the global pos_d/f_d with GLOBAL membrane/cortex indices).
    params = _params(args)
    _banner(params, args.x_beta, n_erm)
    connector = build_erm_cortex_connector(
        membrane_idx_d=membrane.erm_m_d, cortex_idx_d=membrane.erm_c_d,
        bound_d=membrane.erm_bound_d, rest_d=membrane.erm_rest_d,
        membrane_pos_d=cell.pos_d, cortex_pos_d=cell.pos_d, params=params, device=dev,
        rng_epoch_d=membrane.erm_rng_epoch_d,
        detach_events_d=membrane.erm_detach_events_d, attach_events_d=membrane.erm_attach_events_d,
    )
    membrane_owner = build_native_surface_owner(
        name="membrane", position_d=cell.pos_d, force_d=cell.f_d, device=dev)
    cortex_owner = build_native_surface_owner(
        name="cortex", position_d=cell.pos_d, force_d=cell.f_d, device=dev)
    slc = build_erm_cortex_slice(
        membrane_owner=membrane_owner, cortex_owner=cortex_owner, connector=connector,
        base_seed=args.seed, device=dev)

    # 3. Move ERM ownership to the slice (NO double count) + route its force through the driver's inner relax.
    membrane.n_erm = 0                                    # the membrane compartment stops touching ERM
    cell.membrane = _SliceMembraneAdapter(slc, original=membrane)
    inner_solve = make_inner_solve(cell, n_inner=args.outer, reshape_every=20, inner_solver="explicit")

    def membrane_relax() -> None:
        # bending + inextensibility + crosslink + steric + turgor + nucleus + (via the adapter) Helfrich/area
        # of the membrane + the slice's ERM tether — the membrane rides on the cortex through the tethers.
        inner_solve(float(args.dt))

    accepted_ones = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=dev)

    # LOAD RAMP setpoint: bumping the interior FLUID turgor each accepted step is the physical perturbation that
    # makes the membrane bulge against the tethers (no tether-force injection).  Needs the Biot pressure field.
    ramp_pa = float(args.turgor_ramp)
    if ramp_pa > 0.0 and (cell.grid is None or cell.membrane_pressure is None):
        raise RuntimeError("--turgor-ramp needs the Biot pressure field + membrane traction (with_pressure=True)")
    if ramp_pa > 0.0:
        print(f"[gate-b-erm] LOAD RAMP ON: +{ramp_pa:g} Pa/step on the FLUID interior "
              f"(Π {40.0:g}→{40.0 + ramp_pa * args.steps:g} Pa over {args.steps} steps) — tether tension will "
              f"EMERGE and the most-loaded tethers shed by force-driven Bell slip (bleb-onset).", flush=True)
    elif args.preload:
        print("[gate-b-erm] constant-load baseline (--turgor-ramp 0) with preload ON: Π held at resting 40 Pa; "
              "the pre-stretched tethers carry their EMERGENT resting tension (~turgor/node) — the force-balanced "
              "control the ramp perturbs from.", flush=True)
    else:
        print("[gate-b-erm] constant-load baseline (--turgor-ramp 0, --no-preload): Π held at resting 40 Pa, "
              "tethers sit at rest → tension ~0 AND geometry rolled back (non-converged control).", flush=True)

    # 4. Accepted-step event loop — the tether tension EMERGES and the Bell-slip shedding begins under load.
    print(f"\n{'step':>4} {'bound':>8} {'bound%':>7} {'mean_pN':>9} {'total_pN':>11} "
          f"{'areal_pN/µm2':>12} {'detach':>7} {'attach':>7} {'maxPF':>10}", flush=True)
    for step in range(args.steps):
        if ramp_pa > 0.0:
            # Raise the turgor setpoint the inner solve balances BEFORE this step's relax, so the membrane moves
            # outward against the tethers this step; the tether tension is then the k_erm·(length−rest) reaction.
            wp.launch(_ramp_turgor_fluid_kernel, dim=cell.grid.shape,
                      inputs=[cell.grid.p, cell.grid.mask, wp.float64(ramp_pa)], device=dev)
        slc.step(membrane_relax, dt_phys=float(args.dt), accepted_d=accepted_ones)  # force-accept: watch shedding

        # out-of-hot-loop telemetry (device→host readback BETWEEN steps only; I0-A honoured inside the loop)
        st = _tether_stats(cell, slc)
        max_pf = float(inner_solve.convergence_force_d.numpy()[0])
        print(f"{step:>4d} {st['n_bound']:>8d} {st['bound_pct']:>6.1f}% {st['mean_pn']:>9.4g} "
              f"{st['total_pn']:>11.4g} {st['areal_pn_per_um2']:>12.4g} {st['detach_events']:>7d} "
              f"{st['attach_events']:>7d} {max_pf:>10.4g}", flush=True)

    final = _tether_stats(cell, slc)
    print(f"\n[gate-b-erm] DONE: {args.steps} accepted steps. Bound {final['n_bound']}/{n_erm} tethers "
          f"({final['bound_pct']:.1f}%); {final['detach_events']} Bell detachments, {final['attach_events']} "
          f"rebinds. Emergent tether tension total={final['total_pn']:.4g} pN "
          f"(areal {final['areal_pn_per_um2']:.4g} pN/µm²) — EMERGED from the bound population, no imposed "
          f"prestress. Bell-slip shedding under load = bleb-onset. (Bell magnitudes PROVISIONAL — a PI decision.)",
          flush=True)


if __name__ == "__main__":
    main()
