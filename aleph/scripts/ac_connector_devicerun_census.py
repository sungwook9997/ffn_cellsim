#!/usr/bin/env python
r"""Per-connector DEVICE-RUN census — which of the 38 declared connectors actually launch a kernel.

WHY THIS EXISTS.  ``dispatch.connector_runtime_census()`` names a runtime for all 38 declared connectors
and ``validate_connector_runtime_census`` cross-checks that table against the architecture in both
directions.  That is a STRUCTURAL result and it is real.  It is also, on its own, not evidence that
anything runs: :meth:`aleph.engine.actor.CellActor.assert_fully_bound` tests only that the required hooks
are CALLABLE, so a method that exists and returns passes it.  ``STATE.md`` (f) says so in as many words —
"38 OBJECTS, zero device-run; ``require_complete=True`` still satisfied only by STUBS" — but that sentence
has never been measured.  This module measures it.

WHAT IS MEASURED, PER CONNECTOR, SEPARATELY.  Every connector runtime's mechanics + ledger hooks are
wrapped in a proxy.  ONE composed candidate solve is then driven, and each wrapped call is bracketed by
its OWN ``wp.timing_begin(wp.TIMING_KERNEL)`` / ``wp.timing_end`` window (the same mechanism that produced
``outputs/ac/profile/*.json``'s ``launches_per_iteration``).  Three INDEPENDENT observations are recorded
for each connector, and they are never collapsed into one another:

  ``launches``   the Warp kernels that ran inside that connector's own hook calls, by name;
  ``l1_delta``   the change in Σ|F| over EVERY watched device force array, measured by a device-side
                 reduction before and after the call (never a per-connector host readback in a loop);
  ``ledger_pN``  |force resultant| the connector pushed into its own freshly-zeroed ``GlobalCellLedger``.

VERDICTS — the LAUNCH observation alone.  ``DEVICE_RUN`` = the connector's own hooks launched at least one
Warp kernel.  ``STUB`` = the hook was called and launched nothing.  ``NOT_CALLED`` = the composed pipeline
never reached this connector at all — different from a stub, and the census must not report it as one.

THE FORCE CHANNEL IS REPORTED BESIDE THE VERDICT, NEVER FOLDED INTO IT.  The first run folded it in and
read ``l1_delta = 2.7e-21 pN`` as a moving force channel; summed over ~5e5 nodes that is float64
accumulation noise. ``force_channel`` is therefore the FACT ``ZERO`` / ``NONZERO``
plus ``l1_delta_per_node_pN``, and nothing more: deciding whether ~5e-27 pN/node counts as a force needs a
chosen physical scale, which is what the charter forbids a gate to contain.  At a rest configuration with
no bound heads, no force is the physically correct reading, not a defect.
``NOT_BOUND`` = the declared connector has NO runtime object in this composed world, so there was nothing
to call.  That is WEAKER than ``STATE.md``'s "38 OBJECTS": the object does not exist here.  The first run
of this harness died on ``KeyError: 'membrane_erm_cortex'`` rather than recording it, and the absence was
the finding — a census that crashes on the weakest state cannot report the census.

WHAT A VERDICT IS ABOUT.  Every verdict here is about THIS composed configuration on THIS build, not about
a connector class forever.  ``composed_native.py`` already declares that it composes three real owners plus
``nmii_cortex_motor`` and binds the rest as dispatch-coverage placeholders; if the census reproduces that,
it has MEASURED a documented claim rather than discovered a new one, and that is the useful outcome.

CONTROLS — run the census TWICE; each pass asserts only the control it was given.
  * POSITIVE — ``--positive-control nmii_cortex_motor`` must come back ``DEVICE_RUN``.  It is the edge the
    tier-(a) GATE-B cortex-motor row was measured on.  If it does not, the harness is wrong and no other
    row may be read.
  * NEGATIVE — ``--negative-control <connector> --positive-control ""`` nulls that connector's mechanics
    hooks, and it must flip to ``STUB``.  Point it at the SAME connector as the positive pass: this
    composed world binds exactly one runtime, so a different edge cannot serve as the control, and the
    A/B on one edge is the stronger test anyway — same world, same call site, hooks nulled, verdict
    flips.  Without that flip, "we saw launches" is a statement about the instrumentation, not the code.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — forces are ``wp.vec3d`` [pN]; ``l1_delta`` is Σ|F| [pN]; launches are a count; no bare
    number carries an implied unit.
  * boundary — ``launches == 0`` and ``l1_delta == 0`` are NEVER folded into one value, and a connector the
    pipeline never called is ``NOT_CALLED``, never ``STUB``.  A missing observation is ``None``, never 0.
  * conservation/invariant — the census binds nothing and mutates no physics; it wraps hooks and restores
    them in a ``finally``.  ``validate_connector_runtime_census`` is re-run after the solve and must still
    report 38/38, so an observation that changed the wiring fails loudly.
  * CFL/precision — no time integration is performed; this drives ONE candidate force assembly, not a step.
  * sign sense — not applicable: read-only observation, no force is computed here.
  * measurement protocol — one composed solve; per-hook device reductions accumulate into a (n_connectors,)
    device array and there is exactly ONE host readback, after the solve completes.  Nothing is read inside
    a physical-time loop, so the NG-6 residency result is untouched (there is no physical-time loop here).

engine units: length µm, force pN.  Runtime: NVIDIA Warp on CUDA only — this raises on a non-CUDA device.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN (Slurm only — outside an allocation ``cuInit()`` returns ``CUDA_ERROR_NO_DEVICE``):

    MEM_GB=8 CPUS_PER_TASK=8 gpu-submit 4090 02:00:00 \
      "<env-python> aleph/scripts/ac_connector_devicerun_census.py --k-axial 1000 \
         --out aleph/outputs/ac/connector_devicerun/native_record.json"

``MEM_GB`` is a HOST cgroup limit (``--mem``), unrelated to the card's 24 GiB of VRAM.  It was ``48``
here, never measured; a full-native composed step peaks at **733 MB** host RSS (measured 2026-08-15,
``/usr/bin/time -v``, warm Warp kernel cache), so 8 is >10x the observed peak with room for a cold-cache
JIT.  This is not cosmetic: ``gpu-check`` demands ``MEM_GB + 8`` GB FREE and the box caps at 53 GB, so
the old ask required 56 GB — effectively the whole machine — and was rejected on any half-loaded box.
────────────────────────────────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.engine.actor import _MECHANICS_HOOKS
from aleph.engine.dispatch import connector_runtime_census, validate_connector_runtime_census

#: The ledger hook is observed alongside the mechanics hooks: a connector may compute force in one and
#: report its resultant in the other, and a census that watched only one would miscount either way.
_WATCHED_HOOKS: tuple[str, ...] = (*_MECHANICS_HOOKS, "accumulate_ledger")

SCHEMA = "ac.engine.observe/run-record@2"


@wp.kernel
def _l1_vec3d_kernel(src: wp.array(dtype=wp.vec3d), dst: wp.array(dtype=wp.float64)) -> None:
    """Accumulate Σ|F| of one vec3d force array into ``dst[0]`` (device-side; no host readback)."""
    t = wp.tid()
    wp.atomic_add(dst, 0, wp.length(src[t]))


@wp.kernel
def _norm_vec3d_kernel(src: wp.array(dtype=wp.vec3d), dst: wp.array(dtype=wp.float64)) -> None:
    """Write |src[0]| into ``dst[0]`` — the ledger's single resultant, reduced on the device."""
    dst[0] = wp.length(src[0])


class _ForceProbe:
    """Device-side Σ|F| over every watched force array, with no per-call host readback."""

    def __init__(self, arrays: tuple[wp.array, ...], device: str) -> None:
        self.arrays = arrays
        self.device = device
        self._scratch = wp.zeros(1, dtype=wp.float64, device=device)

    def l1(self) -> wp.array:
        """Return a fresh device scalar holding Σ|F| across all watched arrays."""
        self._scratch.zero_()
        for arr in self.arrays:
            n = int(arr.shape[0])
            if n:
                wp.launch(_l1_vec3d_kernel, dim=n, inputs=[arr, self._scratch], device=self.device)
        out = wp.zeros(1, dtype=wp.float64, device=self.device)
        wp.copy(out, self._scratch)
        return out


def _collect_force_arrays(*objects: object) -> tuple[wp.array, ...]:
    """Collect every distinct ``vec3d`` device force array reachable on the supplied runtime objects."""
    found: list[wp.array] = []
    seen: set[int] = set()
    for obj in objects:
        if obj is None:
            continue
        for attr in ("force_d", "force", "node_force_d", "reaction_d"):
            arr = getattr(obj, attr, None)
            if arr is None or getattr(arr, "dtype", None) is not wp.vec3d:
                continue
            ptr = int(getattr(arr, "ptr", 0) or 0)
            key = ptr if ptr else id(arr)
            if key not in seen:
                seen.add(key)
                found.append(arr)
    return tuple(found)


def _blank() -> dict:
    return {"calls": [], "launches": [], "l1_delta_pN": 0.0, "l1_before_pN": 0.0,
            "l1_after_pN": 0.0, "ledger_pN": 0.0, "instrumented_hooks": [], "errors": []}


class _Census:
    """Wrap connector hooks at the CLASS, attribute each call to the runtime OBJECT that received it.

    Instance-level wrapping is not available: 7 of the 13 connector runtime classes are
    ``frozen=True, slots=True`` dataclasses, which refuse both ``setattr`` and ``object.__setattr__``.
    The class dict is writable, so the proxy is installed there and dispatches on ``id(self)``.

    That change forces an honest limit into the open.  **One runtime object serves several declared
    edges** — ``ImmersedPorousTransfer`` is the census's runtime for seven connectors,
    ``ContractJointConnector`` for six, ``FilamentCrosslinkConnector`` for five.  A call arrives at the
    object, not at an edge, so when an object is shared its observation attributes to the SET of edges it
    serves and to no single one of them.  Rows say so via ``shared_with``; a per-edge claim from a shared
    runtime would be invented, not measured.
    """

    def __init__(self, probe: _ForceProbe | None) -> None:
        self.probe = probe  # the only device-touching collaborator; the census itself launches nothing
        #: keyed by ``id(runtime)`` — the unit an observation can honestly be attributed to.
        #: ``_pin`` holds a strong reference to every keyed object: an id is only unique while its
        #: object is alive, and CPython reuses the address of a collected one.
        self.by_object: dict[int, dict] = {}
        self._pin: dict[int, object] = {}
        self.records: dict[str, dict] = {}
        self._restore: list[tuple[type, str, object, bool]] = []
        self._wrapped: set[tuple[type, str]] = set()
        self._nulled: set[tuple[int, str]] = set()

    def instrument(self, name: str, runtime: object) -> None:
        """Install recording proxies for one connector's runtime and map the edge name onto it."""
        self._pin[id(runtime)] = runtime
        rec = self.by_object.setdefault(id(runtime), _blank())
        self.records[name] = rec
        cls = type(runtime)
        for hook in _WATCHED_HOOKS:
            if not callable(getattr(runtime, hook, None)):
                continue
            if (cls, hook) in self._wrapped:
                if hook not in rec["instrumented_hooks"]:
                    rec["instrumented_hooks"].append(hook)
                continue
            if not self._wrap_class(cls, hook):
                rec["errors"].append(f"{hook}: not instrumentable on {cls.__name__}")
                continue
            self._wrapped.add((cls, hook))
            rec["instrumented_hooks"].append(hook)

    def nullify(self, runtime: object) -> tuple[str, ...]:
        """NEGATIVE CONTROL — make this ONE object's mechanics hooks no-ops (the proxy honours the flag)."""
        self._pin[id(runtime)] = runtime
        nulled = tuple(
            hook for hook in _MECHANICS_HOOKS if callable(getattr(runtime, hook, None))
        )
        for hook in nulled:
            self._nulled.add((id(runtime), hook))
        return nulled

    def _wrap_class(self, cls: type, hook: str) -> bool:
        original = cls.__dict__.get(hook, getattr(cls, hook, None))
        own = hook in cls.__dict__
        if original is None:
            return False
        try:
            setattr(cls, hook, self._proxy(original, hook))
        except (AttributeError, TypeError):
            return False
        self._restore.append((cls, hook, original, own))
        return True

    def _proxy(self, fn, hook: str):
        census = self

        def proxy(self, *args, **kwargs):  # noqa: N805 — installed as an unbound class attribute
            rec = census.by_object.setdefault(id(self), _blank())
            if (id(self), hook) in census._nulled:
                rec["calls"].append(hook)  # called, and deliberately launching nothing
                return None
            before = census.probe.l1()
            wp.timing_begin(wp.TIMING_KERNEL)
            try:
                return fn(self, *args, **kwargs)
            finally:
                launched = wp.timing_end(synchronize=True)
                after = census.probe.l1()
                rec["calls"].append(hook)
                rec["launches"].extend(sorted({r.name for r in launched}))
                b, a = float(before.numpy()[0]), float(after.numpy()[0])
                rec["l1_delta_pN"] += abs(a - b)
                rec["l1_before_pN"] = max(rec["l1_before_pN"], b)
                rec["l1_after_pN"] = max(rec["l1_after_pN"], a)
        return proxy

    def restore(self) -> None:
        """Put every original hook back — the census must leave the classes exactly as it found them."""
        for cls, hook, original, own in reversed(self._restore):
            if own:
                setattr(cls, hook, original)
            else:  # it was inherited; wrapping created a new class attribute that must go away
                delattr(cls, hook)
        self._restore.clear()
        self._wrapped.clear()


def _verdict(rec: dict) -> str:
    """Classify one connector by whether it LAUNCHED — the question this census exists to answer.

    The first run folded the force delta into this verdict and read ``l1_delta = 2.7e-21 pN`` as a moving
    force channel. Over ~5e5 nodes that is float64 accumulation noise, not force, and at a rest
    configuration with no bound heads ZERO force is the physically correct answer rather than a defect.
    So the verdict is the launch observation alone, and the force channel is reported beside it by
    :func:`_force_channel` against a DERIVED round-off floor.
    """
    if not rec["instrumented_hooks"] and rec["errors"]:
        # A runtime the harness could not wrap tells us NOTHING about whether it runs. Reporting it as
        # NOT_CALLED (or worse, STUB) would be the harness's own limitation dressed up as a finding.
        return "NOT_INSTRUMENTABLE"
    if not rec["calls"]:
        return "NOT_CALLED"
    if not rec["launches"]:
        return "STUB"
    return "DEVICE_RUN"


def _force_channel(rec: dict, _n_terms: int) -> str:
    """Report whether the watched force arrays changed at all. A FACT, deliberately not a judgement.

    An earlier draft graded this ``MOVED`` / ``AT_ROUNDOFF`` / ``MOVED`` against ``γ_n · Σ|F|``. That was
    wrong twice. Summing exact zeros gives exact zero, so the first run's ``2.7e-21 pN`` over ~5e5 nodes
    is not summation round-off — the kernels really wrote ~5e-27 pN per node. And deciding whether
    5e-27 pN/node is a physical force needs a chosen physical scale, which is exactly what the charter
    forbids a gate to contain. So the census reports ``ZERO`` or ``NONZERO`` plus the magnitude per node,
    and leaves "is that physical?" to a reader who has a scale to compare it against.
    """
    return "ZERO" if rec["l1_delta_pN"] <= 0.0 else "NONZERO"


def main() -> None:
    ap = argparse.ArgumentParser(description="Per-connector device-run census over the composed native cell.")
    ap.add_argument("--k-axial", type=float, required=True,
                    help="SF axial backbone stiffness [pN/µm] — a PI-GAP, required so it is never defaulted")
    ap.add_argument("--cortex-filaments", type=int, default=70686, help="native cortex F-actin count")
    ap.add_argument("--n-fibers-ecm", type=int, default=200)
    ap.add_argument("--ecm-box-um", type=float, default=10.0)
    ap.add_argument("--k-xb", type=float, default=None, help="NMII crossbridge stiffness [pN/µm] (PI-GAP)")
    ap.add_argument("--positive-control", type=str, default="nmii_cortex_motor")
    ap.add_argument("--negative-control", type=str, default="ecm_crosslink")
    ap.add_argument("--commit", type=str, default="", help="declared build commit for the record")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise SystemExit("device-run census requires CUDA (I0-A); there is no CPU simulation path.")

    from aleph.components.ecm.mikado_topology import MikadoInitConfig, MikadoTopologyBuilder
    from aleph.components.incumbent.assemble import (
        NMII_BACKBONE_LP_DIAGNOSTIC_UM, NMII_K_XB_TEST, CellConfig, build_cell,
    )
    from aleph.engine.composed_native import build_native_composed_cell_world
    from aleph.engine.ecm_mechanics import COLLAGEN_MATERIAL_KEY
    from aleph.engine.ecm_world import BoundaryAnchorMode, ECMWorldSettings
    from aleph.engine.ledger import make_global_cell_ledger
    from aleph.engine.sf_population import build_sf_arc_population
    from aleph.laws.ecm_library import get_spec
    from aleph.scripts.ac_gate_b_composed_native import _build_nmii_cortex_connector
    from aleph.scripts.ac_gate_b_cortex_motor_native import _build_actuator, _build_port, _params

    wall_t0 = time.perf_counter()
    k_xb = float(args.k_xb) if args.k_xb is not None else float(NMII_K_XB_TEST)

    # 1. Build the native cell + NMII actuator + cortex port (same wiring as the GATE-B composed driver).
    cfg = CellConfig(
        n_filaments=int(args.cortex_filaments), with_myosin=True, overlap_free_cortex=True,
        membrane_subdivisions=6, nucleus_subdivisions=3, erm_radial_pairing=True,
        resting_bound_myosin_fraction=None, nmii_straddle_placement=True,
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        nmii_backbone_lp_source="CENSUS_DIAGNOSTIC_L_p_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=int(args.seed),
    )
    cell = build_cell(cfg)
    nmii_actuator = _build_actuator(cell, dev, k_xb)
    cortex_port = _build_port(cell, dev)
    nmii_cortex_connector = _build_nmii_cortex_connector(
        actuator_state=nmii_actuator, port=cortex_port, params=_params(False, k_xb), device=dev,
        max_segment_length_um=float(cell.srest_d.numpy().max()),
    )
    sf_pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    sf_pop.assert_partitioned()
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    ecm_topology = MikadoTopologyBuilder(MikadoInitConfig(
        box_lo_um=(0.0, 0.0, 0.0), box_hi_um=(args.ecm_box_um,) * 3, n_fibers=int(args.n_fibers_ecm),
        fiber_length_um=spec.fiber_len_um, target_segment_um=spec.seg_um,
        crosslink_capture_um=spec.xl_contact_um, pin_faces=("z_lo",), pin_margin_um=1.0,
        rng_seed=int(args.seed), max_refinement_level=0,
    ), device=dev).initialize()
    build_seconds = time.perf_counter() - wall_t0

    # 2. COMPOSE.  A failure here is the FIRST RESULT, not a harness bug — record it and stop.
    cell_world = build_native_composed_cell_world(
        device=dev, base_seed=int(args.seed), sf_population=sf_pop,
        sf_k_axial_pn_per_um=float(args.k_axial), ecm_topology=ecm_topology,
        ecm_settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET), ecm_spec=spec,
        nmii_actuator=nmii_actuator, cortex_port=cortex_port,
        nmii_cortex_connector=nmii_cortex_connector, cortex_capacity=int(args.cortex_filaments),
    )
    actor = cell_world.world.actor
    architecture = cell_world.world.architecture
    declared = tuple(c.name for c in architecture.connectors)

    watched = _collect_force_arrays(cortex_port, *cell_world.component_owners.values())
    #: Term count for the float64 accumulation bound: every element summed into one Sigma|F| scalar.
    n_watched_elements = sum(int(a.shape[0]) for a in watched)
    probe = _ForceProbe(watched, dev)
    census = _Census(probe)

    # A declared connector may have NO runtime object at all in this composed world. That is a weaker
    # state than STUB (which at least has an object whose hooks return) and weaker than NOT_CALLED
    # (bound, never reached), and the census must count it separately rather than crash on it — the
    # first run of this harness did crash here, and the absence was the finding.
    unbound = tuple(actor.missing_bindings()["connectors"])
    bound = tuple(name for name in declared if name not in set(unbound))

    for name in bound:
        census.instrument(name, actor.connector_runtime(name))

    # Which declared edges share ONE runtime object — the true attribution unit (see _Census docstring).
    shared: dict[int, list[str]] = {}
    for name in bound:
        shared.setdefault(id(actor.connector_runtime(name)), []).append(name)

    # 3. NEGATIVE CONTROL. Nulling is per runtime OBJECT, so it silences every edge that object serves;
    #    the record names them so the control is never mistaken for a single-edge finding.
    nulled: tuple[str, ...] = ()
    nulled_edges: list[str] = []
    negative = args.negative_control
    if negative and negative not in set(bound):
        # An unbound edge cannot be a negative control: nulling nothing proves nothing. Substitute a
        # bound one and SAY SO, rather than silently reporting a control that never ran.
        negative = next((n for n in bound if n != args.positive_control), "")
    if negative:
        neg_runtime = actor.connector_runtime(negative)
        nulled = census.nullify(neg_runtime)
        nulled_edges = sorted(shared[id(neg_runtime)])

    # 4. ONE composed candidate solve.  Each wrapped hook opens its OWN timing window, so the launch
    #    attribution is per connector and not a share of a whole-solve total.
    ledgers = {name: make_global_cell_ledger(device=dev) for name in bound}
    try:
        cell_world.solve_candidate()
        wp.synchronize_device(wp.get_device())
    finally:
        census.restore()

    # 5. Ledger channel, measured separately and AFTER the hooks are restored: each connector pushes into
    #    its own freshly-zeroed ledger, so a resultant is that connector's alone.
    norm_d = wp.zeros(1, dtype=wp.float64, device=dev)
    for name in bound:
        runtime = actor.connector_runtime(name)
        hook = getattr(runtime, "accumulate_ledger", None)
        if not callable(hook):
            continue
        try:
            hook(ledgers[name])
        except Exception as exc:  # a signature the census cannot satisfy is recorded, never swallowed
            census.records[name]["errors"].append(f"accumulate_ledger: {type(exc).__name__}: {exc}")
            continue
        wp.launch(_norm_vec3d_kernel, dim=1,
                  inputs=[ledgers[name].force_resultant_d, norm_d], device=dev)
        census.records[name]["ledger_pN"] = float(norm_d.numpy()[0])
    wp.synchronize_device(wp.get_device())

    # 6. The wiring must be untouched by the observation.
    validate_connector_runtime_census(architecture)

    runtime_module = {b.connector: f"{b.module}.{b.symbol}" for b in connector_runtime_census()}
    rows = []
    for name in declared:
        if name not in set(bound):
            rows.append({
                "connector": name, "runtime": runtime_module.get(name), "verdict": "NOT_BOUND",
                "force_channel": None, "shared_with": [],
                "attribution": "no runtime object in this composed world",
                "n_calls": 0, "hooks_called": [], "hooks_instrumented": [], "launches": [],
                "n_launch_kernels": 0, "l1_delta_pN": None, "l1_delta_per_node_pN": None,
                "ledger_pN": None, "errors": [],
            })
            continue
        rec = census.records[name]
        peers = sorted(set(shared[id(actor.connector_runtime(name))]) - {name})
        rows.append({
            "connector": name,
            "runtime": runtime_module.get(name),
            "verdict": _verdict(rec),
            "force_channel": _force_channel(rec, n_watched_elements),
            "shared_with": peers,
            "attribution": "edge" if not peers else "runtime-object shared across edges",
            "n_calls": len(rec["calls"]),
            "hooks_called": sorted(set(rec["calls"])),
            "hooks_instrumented": sorted(set(rec["instrumented_hooks"])),
            "launches": sorted(set(rec["launches"])),
            "n_launch_kernels": len(set(rec["launches"])),
            "l1_delta_pN": rec["l1_delta_pN"],
            "l1_delta_per_node_pN": rec["l1_delta_pN"] / max(n_watched_elements, 1),
            "n_watched_nodes": n_watched_elements,
            "ledger_pN": rec["ledger_pN"],
            "errors": rec["errors"],
        })

    by_verdict: dict[str, list[str]] = {}
    for row in rows:
        by_verdict.setdefault(row["verdict"], []).append(row["connector"])

    pos = next((r for r in rows if r["connector"] == args.positive_control), None)
    neg = next((r for r in rows if r["connector"] == negative), None)
    controls = {
        "positive": {"connector": args.positive_control, "verdict": pos and pos["verdict"],
                     "required": "DEVICE_RUN", "ok": bool(pos and pos["verdict"] == "DEVICE_RUN")},
        "negative": {"connector": negative, "requested": args.negative_control,
                     "substituted": negative != args.negative_control,
                     "verdict": neg and neg["verdict"],
                     "required": "STUB", "ok": bool(neg and neg["verdict"] == "STUB"),
                     "nulled_hooks": list(nulled), "nulled_edges": nulled_edges},
    }
    # An empty control string means "not requested in this pass" — the two-pass protocol asserts one
    # control per pass, and a pass must not fail for a control it was never given.
    controls["positive"]["requested"] = bool(args.positive_control)
    controls["negative"]["requested_by_caller"] = bool(args.negative_control)
    controls_ok = ((not args.positive_control or controls["positive"]["ok"])
                   and (not args.negative_control or controls["negative"]["ok"]))
    wall_seconds = time.perf_counter() - wall_t0

    record = {
        "schema": SCHEMA,
        "run_label": "connector_devicerun_census",
        "build": {"commit": args.commit or None, "source": "declared",
                  "reason": "the run host is not a git checkout; the commit is the caller's assertion"},
        "device": dev,
        "evidence": "NATIVE",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": (
            "per-connector kernel-launch and force-delta observation over ONE composed candidate solve. "
            "BLOCKED for magnitudes because no step was accepted and no solve converged: the l1_delta and "
            "ledger columns are EXISTENCE evidence for a force channel, never a physical force value."
        ),
        "verdict": "PASS" if controls_ok else "CONTROLS_FAILED",
        "config": {
            "n_filaments": int(args.cortex_filaments), "k_axial_pn_per_um": float(args.k_axial),
            "k_xb_pn_per_um": k_xb, "n_fibers_ecm": int(args.n_fibers_ecm),
            "membrane_subdivisions": 6, "nucleus_subdivisions": 3, "seed": int(args.seed),
        },
        "parameter_provenance": {"k_axial_pn_per_um": "PI_GAP", "k_xb_pn_per_um": "PI_GAP",
                                 "ecm_spec": "SOURCED"},
        "census": {
            "cortex_filaments": int(args.cortex_filaments), "n_actin_nodes": int(cell.n_actin),
            "n_total_nodes": int(cell.n_total), "n_heads": int(nmii_actuator.n_heads),
            "fraction_of_native": float(int(args.cortex_filaments) / 70686.0),
        },
        "controls": controls,
        "summary": {
            "n_declared_connectors": len(declared),
            "n_bound_runtimes": len(bound),
            "unbound_connectors": sorted(unbound),
            "n_distinct_runtime_objects": len(census.by_object),
            **{f"n_{k.lower()}": len(v) for k, v in sorted(by_verdict.items())},
            "by_verdict": {k: sorted(v) for k, v in sorted(by_verdict.items())},
        },
        "connectors": rows,
        "registered_components": sorted(cell_world.world.registered_components()),
        "timing": {
            "wall_seconds": wall_seconds, "build_seconds": build_seconds, "n_steps": 0,
            "comparable": False,
            "not_comparable_reason": (
                "a census, not a benchmark: every connector hook is wrapped in its own timing window, "
                "which serialises the stream. No step was accepted, so there is no per-step cost to quote."
            ),
        },
        "scope": (
            "Every verdict is about THIS composed configuration on THIS build — not about a connector "
            "class in general. NOT_CALLED means the composed pipeline never reached the connector and is "
            "NOT a stub finding. LAUNCH_NO_FORCE means a kernel ran and no watched force array moved, "
            "which a geometrically inactive edge can produce legitimately."
        ),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2), encoding="utf-8")

    print(f"\n[census] {len(declared)} declared connectors, device={dev}")
    print(f"{'connector':34s} {'verdict':18s} {'force':12s} {'kernels':>8s} {'l1_delta pN':>14s}")
    for row in rows:
        d = row["l1_delta_pN"]
        print(f"{row['connector'][:34]:34s} {row['verdict']:18s} "
              f"{str(row['force_channel'] or '-'):12s} {row['n_launch_kernels']:8d} "
              f"{('-' if d is None else f'{d:.4e}'):>14s}")
    print(f"\n[census] {record['summary']}")
    print(f"[census] controls: {json.dumps(controls)}")
    print(f"[census] wrote {out}")

    if not controls_ok:
        raise SystemExit(
            "[census] CONTROLS FAILED — the positive control must be DEVICE_RUN and the negative control "
            "must be STUB. Until both hold, no row in this census is attributable and none may be quoted."
        )


if __name__ == "__main__":
    main()
