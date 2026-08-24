r"""The barbed-end channel: what makes a lamellipodium and a filopodium PUSH, bound to arena ranges.

WHAT WAS MISSING.  ``world/build/lamellipodium.py`` and ``world/build/filopodium.py`` lay down two
populations of straight chains with uniform polarity, and that is all they do — a comb and a set of
fingers, geometrically correct and completely inert.  Nothing in the arena makes a barbed end grow, so
neither population can protrude, and "ACTIVE_LOAD_PATH" is at present a role name with no active term
behind it.  This module supplies the term.

**IT WRITES NO LAW, AND THAT IS THE DESIGN.**  The Brownian ratchet (KB-3.6, Mogilner-Oster) is already
in this repository twice — :func:`aleph.laws.polymerization_warp.resolve_polymerization` resolves
``v0`` from Pollard's barbed-end kinetics and :func:`aleph.laws.motility_warp.barbed_end_growth_kernel`
applies ``v(f) = v0 * exp(-f*delta/kBT)`` to the tip segment of a uniform-stride filament population.
That kernel's addressing (``fiber_off``, ``seg_off``, tip = last node) is EXACTLY the arena's layout,
so the binding is two index arrays and a launch.  ``world/`` does not own laws; re-deriving the ratchet
here would put a second copy of a force law under a directory that has no authority to hold one.

WHY ``fiber_off``/``seg_off`` ARE ARITHMETIC.  A built population is uniform-stride: strand ``s`` owns
nodes ``[node_lo + s*n_per, node_lo + (s+1)*n_per)`` and segments ``[seg_lo + s*(n_per-1), ...)``.  So
the two offset arrays are ``arange``s, uploaded ONCE at bind time and never again.  This matters: the
alternative (walk the strand list on the host each step) is precisely the host mirror the world port
exists to delete.

THE SIGN, WHICH IS THE WHOLE CLAIM OF THIS SLICE.  Three things have to point the same way or the
"protrusion" runs backwards, and each is checked rather than assumed:

  1. *tip is the LAST node.*  Both builders write bead ``i`` at ``root + i*seg*d``, so the last bead is
     the outermost and :attr:`aleph.world.strand.Strand.polarity` is ``+1``.  Checked on the host from
     the census, per kind, in :func:`plan_barbed_ratchet`.
  2. *the tip tangent points OUT.*  For the lamellipodium ``d = cos(mode)*axis + -sin(mode)*lateral``
     with ``cos(mode) > 0``, so the mean tip tangent is the protrusion axis; for the filopodium ``d``
     is the outward finger direction and the root sits at ``root_R_um * d``.  Both are host-checkable
     facts about the census and both are asserted.
  3. *load OPPOSES growth.*  The kernel takes ``f_load = max(-force[tip] . that, 0)``, so a force
     pushing the tip inward is positive load and SLOWS the ratchet; a tip under tension grows at
     ``v0`` and no faster.  That is the law's own sign and it is not touched here.

**NO VELOCITY AND NO FORCE IS CLAIMED BY THIS MODULE.**  What is claimed is that the channel exists,
that it is addressed to the right nodes, and that its sign is the one above.  ``v0`` follows from a
declared G-actin concentration through the law's own resolver; whether that concentration is this
cell's is a separate question this module does not answer and must not appear to.

WHAT IS DELIBERATELY ABSENT — and it is half of KB-3.8.  The filopodial length balance is
``dL/dt = v_p - v_retro - k_cap*L``.  This module supplies ``v_p`` and NOTHING ELSE: there is no
retrograde-flow term and no capping term in :mod:`aleph.laws` that is addressed to a free filopodial
bundle (``motility_warp.fiber_treadmill_kernel`` and ``pointed_end_depoly_kernel`` are both gated on a
cortex polarity cap and are not that object).  Writing them would be writing a NEW LAW, which belongs
in ``laws/`` and is a Lead decision, not a ``world/`` file.  So the binding here is named for what it
is — a barbed-end ratchet — and never for the balance it is one term of.  A filopodium driven by this
module alone only ever gets longer, and that is a stated gap, not a result.

⚠ **THE LAMELLIPODIAL COUNT IS 2.5-10x UNDER** and this module does not fix it.  ``KB-3.7`` is
``verified`` with *"~100 fil per um leading edge"*; at a 5-20 um edge that is 500-2,000 filaments
against the build's 200.  The load a tip feels is divided across the barbed ends sharing it, so the
count is not cosmetic — it sets the per-filament load.  The leading-edge width is a PI decision
(``PER_CELL_STRUCTURE_COUNTS_2026-08-20.md`` §5 item 1) and until it lands, a run through this module
is a run at the placeholder count.  :func:`plan_barbed_ratchet` reports the count it was handed
alongside the KB band so the artifact carries the discrepancy rather than hiding it.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [um]; ``force`` [pN]; ``seg_rest`` and ``seg_max_um`` [um]; ``dt`` [s];
    ``v0`` [um/s]; ``delta`` [um]; ``kBT`` [pN.um]; ``G_actin_uM`` [uM].  ``grown`` is [um] summed.
  * boundary — a population of fewer than 3 beads per strand is refused (no bending triple, so the tip
    tangent has no interior neighbour to be measured against); a non-uniform census is refused rather
    than strided over; a ``seg_max_um`` at or below the realised segment length is refused, because a
    tip with zero room reads downstream as "stalled by load" when it is stalled by bookkeeping.
  * conservation/invariant — every index the kernel can reach is asserted on the host to lie inside
    this population's own NODE and SEGMENT claims before any launch.  Growth is rest-length only: no
    node is moved here, so this module conserves the centre of mass exactly and the tip advance is the
    reshape's to make.  ``grown[0]`` accumulates the total drawn length for the G-actin budget.
  * CFL/precision — float64 throughout.  ``dt`` is bounded by the tip cap rather than by stability:
    the update is kinematic and monotone, and ``wp.min(r + v*dt, seg_max)`` cannot overshoot.
  * sign sense — the three-part argument above; asserted host-side in :func:`_demo` against the law's
    own reference :func:`aleph.laws.polymerization_warp.ratchet_velocity_np`.
  * measurement protocol — two int32 offset arrays and one float64 rest-length array are uploaded ONCE
    at bind time.  :meth:`BarbedRatchet.step` launches and returns; nothing is read back.

engine units: length um, force pN, time s.  Runtime: NVIDIA Warp on CUDA; :func:`plan_barbed_ratchet`
touches no device and :func:`_demo` needs no card.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.laws.motility_warp import barbed_end_growth_kernel
from aleph.laws.polymerization_warp import ResolvedPolymerization, resolve_polymerization
from aleph.world.arena import Kind, WorldArena

__all__ = ["BarbedRatchet", "plan_barbed_ratchet", "bind_barbed_ratchet", "KB37_FIL_PER_UM_EDGE"]

#: ``KB-3.7``, status ``verified``: *"fil density ~100 per um leading edge"*.  Carried here ONLY so the
#: plan can report how far the census's count sits from it.  It is never used to change a count — the
#: leading-edge width is an open PI decision and a count derived from it would be a build detail
#: answering a decision.
KB37_FIL_PER_UM_EDGE = 100.0

#: Per kind, the census key whose POSITIVITY witnesses that the tip tangent points outward, and why.
#: A kind absent from this table is refused: "the tip is the last node" is a property of the builder
#: that wrote the population, and assuming it for an unknown builder is how a protrusion runs backwards.
_OUTWARD_WITNESS: dict[str, tuple[str, str]] = {
    "lamellipodium": (
        "cos_mode",
        "bead i sits at seed + i*seg*(cos_mode*axis + -sin_mode*lateral); cos_mode > 0 makes the mean "
        "tip tangent the protrusion axis, so the patch grows FORWARD rather than across itself",
    ),
    "filopodium": (
        "root_R_um",
        "the finger roots at root_R_um*d on the cell surface and bead i sits at root + i*seg*d, so the "
        "last bead is the outermost and growth is away from the cell",
    ),
}


@dataclass(frozen=True, slots=True)
class BarbedRatchet:
    """One population's bound barbed-end channel: the law's parameters plus its arena addressing.

    Attributes:
        population: the population these ranges belong to.
        plan: the host-side census this was bound from, including the count discrepancy report.
        spec: the resolved ratchet, from :mod:`aleph.laws.polymerization_warp`.
        seg_max_um: the tip-segment rest-length cap [um]; growth past it needs bead INSERTION, which is
            a dynamic topology change and is not this module's.
        fiber_off_d / seg_off_d: ``(n_strands + 1,)`` and ``(n_strands,)`` GLOBAL offsets, int32.
        seg_rest_d: globally indexed segment rest lengths [um].  Mutated in place by :meth:`step`.
        grown_d: ``(1,)`` accumulated drawn length [um], for the G-actin budget (KB-3.21).
    """

    population: str
    plan: dict[str, object]
    spec: ResolvedPolymerization
    seg_max_um: float
    fiber_off_d: wp.array
    seg_off_d: wp.array
    seg_rest_d: wp.array
    grown_d: wp.array

    @property
    def n_tips(self) -> int:
        """Number of barbed ends this channel drives — one per strand."""
        return int(self.plan["n_strands"])

    def step(self, arena: WorldArena, dt: float) -> None:
        """Advance every barbed end by ``v(f)*dt`` of REST length.  Launches; reads nothing back.

        Args:
            arena: the world holding ``position`` and ``force``.  ``force`` must already hold this
                step's assembled load — the ratchet reads the load the network actually applied, which
                is the whole point of calling it load-dependent.
            dt: physical timestep [s].  Must be positive; a zero step is refused rather than launched,
                since "grew by nothing" and "was never launched" are different facts downstream.
        """
        if not (dt > 0.0):
            raise ValueError(
                f"dt={dt} must be positive. A zero-length step launches a kernel that changes nothing, "
                "which is indistinguishable downstream from a channel that was never bound."
            )
        wp.launch(
            barbed_end_growth_kernel,
            dim=self.n_tips,
            inputs=[
                arena.node_arrays["position"], self.fiber_off_d, self.seg_off_d, self.seg_rest_d,
                self.spec.v0_um_s * float(dt), self.spec.delta_um, self.spec.kBT_pN_um,
                arena.node_arrays["force"], float(self.seg_max_um), self.grown_d,
            ],
            device=arena.device,
        )


def plan_barbed_ratchet(
    census: dict[str, object], *, G_actin_uM: float | None = None, seg_max_um: float | None = None,
) -> dict[str, object]:
    """Resolve the channel from a built population's census.  No device is touched.

    Args:
        census: what ``build_lamellipodium`` / ``build_filopodia`` returned.  Read, never modified.
        G_actin_uM: the free G-actin concentration [uM] the barbed-end kinetics run at.  **REQUIRED**
            — it sets ``v0`` through ``(k_on*G - k_off)*delta_full``, so defaulting it would choose the
            protrusion speed silently.  The law bands it to [1, 100] and refuses below the critical
            concentration.
        seg_max_um: cap on the tip segment's rest length [um].  **REQUIRED** — sustained growth past it
            needs bead insertion, and a cap picked here would be an unsourced contour limit.

    Returns:
        The plan: addressing, the resolved ratchet's numbers, the sign witness, and — for a
        lamellipodium — how far the built count sits from ``KB-3.7``'s verified density.

    Raises:
        ValueError: on an unknown kind, a missing declaration, a census that is not uniform-stride, a
            strand of fewer than 3 beads, a non-positive sign witness, or a cap at or below the
            realised segment length.
    """
    kind = str(census.get("kind", ""))
    if kind not in _OUTWARD_WITNESS:
        raise ValueError(
            f"kind={kind!r} has no outward-tip witness in this module. That the barbed end is the LAST "
            "node is a property of the builder that wrote the population, not a convention; binding a "
            "ratchet to an unknown layout is how a protrusion runs inward while reporting a speed. "
            f"Known: {sorted(_OUTWARD_WITNESS)}."
        )
    if G_actin_uM is None or seg_max_um is None:
        missing = "G_actin_uM" if G_actin_uM is None else "seg_max_um"
        raise ValueError(
            f"{missing} has no default and must be declared. A default here would be an unsourced "
            "constant governing the protrusion under a new name."
        )

    claims = census.get("claims")
    if not isinstance(claims, dict) or "node" not in claims or "segment" not in claims:
        raise ValueError(
            "census carries no node/segment claims, so it is a PLAN and not a BUILD. A channel bound to "
            "a plan would address node 0 of the arena, which belongs to whichever population claimed "
            "first."
        )
    n_strand = int(census["n_strands"])
    n_per = int(census["nodes_per_strand"])
    if n_per < 3:
        raise ValueError(
            f"nodes_per_strand={n_per}: the ratchet measures the tip tangent as pos[tip]-pos[tip-1] and "
            "reports the load along it, so a 2-bead strand has a tangent but no interior neighbour to "
            "carry the reaction. That is a different object, not a coarser one."
        )
    node_lo, node_hi = (int(v) for v in claims["node"])
    seg_lo, seg_hi = (int(v) for v in claims["segment"])
    if node_hi - node_lo != n_strand * n_per:
        raise ValueError(
            f"the node claim spans {node_hi - node_lo} but {n_strand} strands x {n_per} beads is "
            f"{n_strand * n_per}. This population is not uniform-stride, and striding over it anyway "
            "would address one strand's tip inside another strand's body."
        )
    if seg_hi - seg_lo != n_strand * (n_per - 1):
        raise ValueError(
            f"the segment claim spans {seg_hi - seg_lo} but {n_strand} strands x {n_per - 1} segments "
            f"is {n_strand * (n_per - 1)}; the rest-length addressing would be off by a strand."
        )

    witness_key, witness_why = _OUTWARD_WITNESS[kind]
    witness = float(census[witness_key])
    if not (witness > 0.0):
        raise ValueError(
            f"{kind}: the outward-tip witness {witness_key}={witness} is not positive, so the last bead "
            f"is not the outermost. {witness_why}. Growing the tip segment of this population would "
            "protrude it INWARD while the artifact called it protrusion."
        )

    seg_realised = float(census["seg_um_realised"])
    cap = float(seg_max_um)
    if cap <= seg_realised:
        raise ValueError(
            f"seg_max_um={cap} is at or below the realised segment length {seg_realised} um, so every "
            "tip starts with zero room and stops on the first step. Downstream that reads as 'stalled "
            "by load', which is the one thing this channel exists to measure."
        )

    spec = resolve_polymerization(G_actin_uM=float(G_actin_uM))

    plan: dict[str, object] = {
        "channel": "barbed_end_ratchet",
        "law": "aleph.laws.motility_warp.barbed_end_growth_kernel (KB-3.6 Mogilner-Oster)",
        "kind": kind,
        "population": str(census.get("population", kind)),
        "n_strands": n_strand,
        "nodes_per_strand": n_per,
        "node_claim": (node_lo, node_hi),
        "segment_claim": (seg_lo, seg_hi),
        "tip_node_first": node_lo + n_per - 1,
        "tip_node_last": node_lo + (n_strand - 1) * n_per + n_per - 1,
        "tip_segment_first": seg_lo + n_per - 2,
        "seg_um_realised": seg_realised,
        "seg_max_um": cap,
        "sign_witness": {"key": witness_key, "value": witness, "why": witness_why},
        "ratchet": {
            "G_actin_uM": spec.G_actin_uM, "v0_um_s": spec.v0_um_s, "delta_um": spec.delta_um,
            "kBT_pN_um": spec.kBT_pN_um, "efold_force_pN": spec.efold_force_pN,
        },
        # The two terms of KB-3.8 this channel is NOT. Recorded in the plan so the artifact carries the
        # gap rather than a reader inferring a length balance from a growth-only run.
        "kb38_terms_absent": ["v_retro (retrograde flow)", "k_cap*L (capping)"],
        "kb38_note": (
            "dL/dt = v_p - v_retro - k_cap*L; only v_p is bound. Neither shrink term exists in laws/ "
            "addressed to a free bundle, and writing one is a new law and a Lead decision."
        ),
    }
    if kind == "lamellipodium":
        # KB-3.7 is verified and the build's count is explicit, so the two can be compared but only one
        # of them is a measurement. Report the ratio; change nothing.
        implied_edge_um = n_strand / KB37_FIL_PER_UM_EDGE
        plan["kb37_count_check"] = {
            "built_n_filaments": n_strand,
            "kb37_fil_per_um_edge": KB37_FIL_PER_UM_EDGE,
            "implied_leading_edge_um": implied_edge_um,
            "kb37_edge_band_um": (5.0, 20.0),
            "under_by": (5.0 / implied_edge_um, 20.0 / implied_edge_um),
            "note": (
                "the leading-edge width is an OPEN PI decision (PER_CELL_STRUCTURE_COUNTS_2026-08-20 "
                "§5 item 1); the load a tip feels is shared across the barbed ends at the edge, so an "
                "under-count raises the per-filament load and the ratchet reports it as stall"
            ),
        }
    return plan


def bind_barbed_ratchet(
    arena: WorldArena, census: dict[str, object], *,
    G_actin_uM: float | None = None, seg_max_um: float | None = None,
    seg_rest_d: wp.array | None = None,
) -> BarbedRatchet:
    """Upload the channel's addressing once and return the handle that steps it.

    Args:
        arena: the world the population was built into.  Must hold a CUDA allocation.
        census / G_actin_uM / seg_max_um: see :func:`plan_barbed_ratchet`.
        seg_rest_d: globally indexed segment rest lengths [um].  **Pass the arena-wide array once
            PHASE 3 bonds own one** — rest length is a bond property and this module is only its first
            consumer.  ``None`` allocates one sized to the arena's SEGMENT capacity, with this
            population's slice at its realised segment length and everything else zero; that is a
            stand-in, and it is the one piece of state this module owns that it would rather not.

    Returns:
        The bound channel.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_barbed_ratchet`, or a ``seg_rest_d`` too short to
            address this population's segment claim.
    """
    plan = plan_barbed_ratchet(census, G_actin_uM=G_actin_uM, seg_max_um=seg_max_um)
    n_strand = int(plan["n_strands"])
    n_per = int(plan["nodes_per_strand"])
    node_lo, _ = plan["node_claim"]
    seg_lo, seg_hi = plan["segment_claim"]

    if not arena.node_arrays:
        raise RuntimeError(
            f"{plan['population']}: this arena holds no device allocation. Warp CUDA is the only "
            f"runtime and there is no CPU path. The channel it would have bound: {n_strand} barbed "
            f"ends at v0={plan['ratchet']['v0_um_s']} um/s."
        )

    # Uniform stride, so both offset tables are aranges. Uploaded once; never rebuilt.
    fiber_off = node_lo + np.arange(n_strand + 1, dtype=np.int64) * n_per
    seg_off = seg_lo + np.arange(n_strand, dtype=np.int64) * (n_per - 1)
    if int(fiber_off[-1]) > np.iinfo(np.int32).max:
        raise ValueError(
            f"node index {int(fiber_off[-1])} exceeds int32, which is the width the laws/ kernels take. "
            "A multi-cell arena reaches this; narrowing silently is how it would wrap."
        )

    if seg_rest_d is None:
        capacity = int(arena.capacity.get(Kind.SEGMENT, 0))
        host = np.zeros(max(capacity, seg_hi), dtype=np.float64)
        host[seg_lo:seg_hi] = float(plan["seg_um_realised"])
        seg_rest_d = wp.array(host, dtype=wp.float64, device=arena.device)
    elif int(seg_rest_d.shape[0]) < seg_hi:
        raise ValueError(
            f"seg_rest_d holds {int(seg_rest_d.shape[0])} segments but this population's claim ends at "
            f"{seg_hi}. Segment rest lengths are addressed GLOBALLY here, the same as nodes."
        )

    return BarbedRatchet(
        population=str(plan["population"]),
        plan=plan,
        spec=resolve_polymerization(G_actin_uM=float(G_actin_uM)),
        seg_max_um=float(plan["seg_max_um"]),
        fiber_off_d=wp.array(fiber_off.astype(np.int32), dtype=wp.int32, device=arena.device),
        seg_off_d=wp.array(seg_off.astype(np.int32), dtype=wp.int32, device=arena.device),
        seg_rest_d=seg_rest_d,
        grown_d=wp.zeros(1, dtype=wp.float64, device=arena.device),
    )


def _demo() -> None:
    """Self-check: the addressing, every refusal, and the SIGN — none of which needs a device."""
    from aleph.laws.polymerization_warp import ratchet_velocity_np

    # A lamellipodium census in the shape build_lamellipodium returns, at the PLACEHOLDER count.
    lam = {
        "kind": "lamellipodium", "population": "lamellipodium", "n_strands": 200,
        "nodes_per_strand": 21, "seg_um_realised": 0.05, "cos_mode": 0.8191520442889918,
        "claims": {"node": (1000, 1000 + 200 * 21), "segment": (500, 500 + 200 * 20)},
    }
    plan = plan_barbed_ratchet(lam, G_actin_uM=20.0, seg_max_um=0.2)

    # Addressing: strand 0's tip is the LAST bead of the first stride, strand 199's the last of the last.
    assert plan["tip_node_first"] == 1000 + 20
    assert plan["tip_node_last"] == 1000 + 199 * 21 + 20
    assert plan["tip_segment_first"] == 500 + 19          # segment 19 joins beads 19 and 20
    assert plan["sign_witness"]["key"] == "cos_mode"

    # The count discrepancy is REPORTED, not repaired: 200 filaments implies a 2 um edge, and KB-3.7's
    # band starts at 5 um. Nothing here changes the count.
    chk = plan["kb37_count_check"]
    assert abs(chk["implied_leading_edge_um"] - 2.0) < 1e-12
    assert chk["under_by"][0] == 2.5 and chk["under_by"][1] == 10.0

    # SIGN 1 - load OPPOSES growth, and tension does not speed it past v0. The law's own reference.
    spec = resolve_polymerization(G_actin_uM=20.0)
    v = ratchet_velocity_np([0.0, 1.0, 5.0, 20.0], spec)
    assert v[0] == spec.v0_um_s                                   # no load -> the full rate
    assert np.all(np.diff(v) < 0.0)                               # heavier load -> slower, monotone
    assert np.all(v > 0.0)                                        # a ratchet stalls asymptotically
    assert ratchet_velocity_np(-100.0, spec) == spec.v0_um_s      # tension is clamped, not a boost

    # SIGN 2 - the tip tangent points along the protrusion axis, which is what cos_mode > 0 witnesses.
    # Rebuild one lamellipodial filament's beads from the builder's own formula and check it directly.
    axis, lat = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    d = lam["cos_mode"] * axis - np.sqrt(1.0 - lam["cos_mode"] ** 2) * lat
    beads = np.array([i * 0.05 * d for i in range(21)])
    assert float(np.dot(beads[-1] - beads[-2], axis)) > 0.0       # the tip grows FORWARD
    assert float(np.dot(beads[-1] - beads[0], axis)) > 0.0        # and the whole filament leans forward

    # A filopodium binds through the same door, on its own witness.
    filo = {
        "kind": "filopodium", "population": "filopodium", "n_strands": 50 * 20,
        "nodes_per_strand": 41, "seg_um_realised": 0.05, "root_R_um": 7.5,
        "claims": {"node": (0, 1000 * 41), "segment": (0, 1000 * 40)},
    }
    fplan = plan_barbed_ratchet(filo, G_actin_uM=20.0, seg_max_um=0.2)
    assert fplan["sign_witness"]["key"] == "root_R_um"
    assert "kb37_count_check" not in fplan                        # that band is the lamellipodium's

    # And the two shrink terms of KB-3.8 are declared ABSENT rather than silently unimplemented.
    assert len(fplan["kb38_terms_absent"]) == 2
    assert "v_retro" in fplan["kb38_terms_absent"][0]

    # ── every refusal ───────────────────────────────────────────────────────────────────────────
    def refuses(fragment: str, census: dict[str, object], **kwargs: object) -> None:
        try:
            plan_barbed_ratchet(census, **kwargs)                      # type: ignore[arg-type]
        except ValueError as exc:
            assert fragment in str(exc), f"expected {fragment!r} in: {exc}"
        else:  # pragma: no cover
            raise AssertionError(f"must refuse: {fragment}")

    ok = {"G_actin_uM": 20.0, "seg_max_um": 0.2}
    refuses("no outward-tip witness", lam | {"kind": "cortex"}, **ok)
    refuses("no default", lam, seg_max_um=0.2)
    refuses("no default", lam, G_actin_uM=20.0)
    refuses("a PLAN and not a BUILD", {k: v for k, v in lam.items() if k != "claims"}, **ok)
    refuses("different object", lam | {"nodes_per_strand": 2}, **ok)
    refuses("not uniform-stride", lam | {"n_strands": 199}, **ok)
    refuses("off by a strand", lam | {"claims": {"node": lam["claims"]["node"],   # type: ignore[index]
                                                 "segment": (500, 500 + 199 * 20)}}, **ok)
    refuses("not positive", lam | {"cos_mode": 0.0}, **ok)
    refuses("zero room", lam, G_actin_uM=20.0, seg_max_um=0.05)
    # The law's own band still binds through this module - it is not re-implemented, so it cannot drift.
    try:
        plan_barbed_ratchet(lam, G_actin_uM=1000.0, seg_max_um=0.2)
    except ValueError as exc:
        assert "physiological band" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an out-of-band G-actin must refuse")

    # A device-less arena refuses the BIND while still reporting what it would have bound.
    arena = WorldArena(capacity={Kind.NODE: 10, Kind.SEGMENT: 10})
    try:
        bind_barbed_ratchet(arena, lam, **ok)                          # type: ignore[arg-type]
    except RuntimeError as exc:
        assert "only" in str(exc) and "200 barbed ends" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less arena must refuse the bind")

    print("protrusion: addressing, sign and refusals OK (no device, no claim about any speed)")


if __name__ == "__main__":
    _demo()
