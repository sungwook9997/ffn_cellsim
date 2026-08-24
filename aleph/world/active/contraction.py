r"""The contractile channel: what makes a stress fibre PULL, bound to arena ranges — and why it is empty.

WHAT WAS MISSING, AND WHAT STILL IS.  ``world/build/stress_fiber.py`` lays down ventral bundles of
straight chains, and a stress fibre with no motor in it is a passive cable: it can resist being pulled
but it cannot pull.  ``KB-DRAFT-7-08`` says the object is a PRESTRESSED actomyosin cable, and the
prestress has one source — bound NMII.  This module supplies the channel that source acts through.

**IT DOES NOT SUPPLY THE STATIONS, AND THAT IS THE POINT.**  A2's builder returns
``n_stations_placed: 0`` and records ``sarcomere_um`` without using it, for a stated MODEL reason
(``STATE.md`` (e) item 5): the sarcomeric banding along a ventral fibre is CURVED, and a rigid bipolar
minifilament cannot straddle two anti-parallel filaments across a curved band.  That is a physics
question, not a missing loop.  Placing stations here from ``sarcomere_um`` would answer an open
decision with a build detail — the exact thing the builder's own docstring refuses to do.  So the
stations arrive as an ARGUMENT, this module validates them, and with none supplied it refuses the bind
and says why.  **A zero-station launch is not offered**: a kernel that runs over nothing reads
downstream as "the motors produced no contraction", which is the opposite of "there are no motors".

**IT WRITES NO LAW.**  :func:`aleph.laws.myosin_linear.minifilament_kernel` already is the Stam-Hocky
bipolar minifilament with the sourced non-muscle force-velocity, and its addressing — a ``(M, 2)``
table of node indices — is exactly what a station list is.  Binding is a validation pass and a launch.

⚠ **THE FORCE-VELOCITY LAW IS LINEAR, NOT HILL, AND THAT IS RATIFIED HERE.**  A Hill hyperbola for
non-muscle myosin IIA does not exist in the literature — ``a/F0 ~ 0.25`` is muscle (Hill 1938).  The
sourced law is affine, ``v(F) = v0*(1 - F/F_stall)``, PI-ratified 2026-07-07 off the lit-ingest and
recorded in :mod:`aleph.laws.myosin_linear`.  Any request phrased as "Hill w~2" for this population is
superseded by that ruling, and this module inherits the ratified law rather than restating it.

THE SIGN, WHICH IS THE WHOLE CLAIM OF THIS SLICE.  The kernel inverts the force-velocity law: at the
current inter-anchor sliding velocity it exerts ``F = F_stall*(1 - v_slide/v0)``, applies ``+f`` to
``i`` along ``(pos[j] - pos[i])`` and ``-f`` to ``j``.  So each anchor is pulled TOWARD the other:
the station SHORTENS, never lengthens, and past stall the magnitude is clamped at zero rather than
going negative — a motor does not actively push its own anchors apart.  The pair is equal-and-opposite
by construction, so a station applies no net force to the cell.  All three are asserted host-side in
:func:`_demo` against a numpy twin of the kernel.

WHY ``v_slide`` STARTS AT ZERO, AND WHY THAT IS THE PHYSIOLOGICAL BASELINE RATHER THAN A CONVENIENT
NULL.  ``v_slide = 0`` is the ISOMETRIC condition, and the law returns the full ensemble stall there.
That is the resting prestress of a ventral fibre — the state the cell is actually in before anything
is done to it — not a zeroed-out motor.  A quasi-static solve relaxes toward ``v -> 0``, so the
equilibrium this channel drives toward IS the stall prestress.  Starting at ``v0`` instead would start
the fibre force-free, which is the unphysical baseline the charter forbids.

**NO FORCE MAGNITUDE IS CLAIMED BY THIS MODULE.**  ``F_stall`` follows from ``KB-3.18``'s composition
through the law's own resolver; whether that composition is this cell's, and where the stations sit,
are the two open questions this module is deliberately unable to answer on its own.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [um]; ``force`` [pN]; ``v_slide`` and ``v0`` [um/s]; ``f_stall`` [pN].
    Station indices carry no unit.
  * boundary — an empty station table is REFUSED, not launched at dimension zero; a station whose two
    nodes lie on the SAME strand is refused, because a bipolar minifilament bridges two anti-parallel
    filaments and a link inside one filament is a spring wearing a motor's name; a station bridging
    two CO-ORIENTED strands is refused for the same reason one layer up, and a census that cannot say
    which way a strand points is refused rather than trusted (see :func:`_anti_parallel_witness` — no
    built population can say it today, and that is the finding, not this module's gap); a self-link is
    refused; a duplicate station is refused, since it doubles a force by bookkeeping.
  * conservation/invariant — every index is asserted on the host to lie inside this population's own
    NODE claim before any launch.  The kernel's two ``atomic_add``s are equal and opposite, so the
    channel contributes exactly zero net force and zero net torque about the station's own axis; the
    force is internal, which is what makes it prestress rather than propulsion.
  * CFL/precision — nothing is integrated here; the channel accumulates into ``force`` and the outer
    solver owns the step.  float64 throughout; the int32 narrowing of node indices is range-checked.
  * sign sense — the three-part argument above; asserted in :func:`_demo` against a numpy twin.
  * measurement protocol — the station table and its sliding-velocity array are uploaded ONCE at bind
    time.  :meth:`ContractileStations.accumulate` launches and returns; nothing is read back.

engine units: length um, force pN, time s.  Runtime: NVIDIA Warp on CUDA; :func:`plan_contraction`
touches no device and :func:`_demo` needs no card.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.myosin_linear import MyosinMinifilament, minifilament_kernel, resolve_myosin
from aleph.world.arena import WorldArena

__all__ = ["ContractileStations", "plan_contraction", "bind_contraction", "STATION_GAP"]

#: The standing MODEL reason there are no stations to bind, quoted where a caller will meet it.  It is
#: an open PI item, not a bug, and the correct response to it is a decision rather than a placement.
STATION_GAP = (
    "STATE.md (e) item 5 — SF geometry REFUSES a motor station (N_stations: 0) for a MODEL reason: the "
    "sarcomeric banding along a ventral fibre is curved and a rigid bipolar minifilament cannot "
    "straddle two anti-parallel filaments across it. Do not force it. Placing stations from the "
    "recorded sarcomere_um would answer an open decision with a build detail."
)


def _anti_parallel_witness(census: dict[str, object]) -> tuple[np.ndarray | None, dict[str, object]]:
    """Resolve the per-strand polarity a station table must be checked against, or say what is missing.

    A bipolar minifilament SLIDES two anti-parallel filaments past each other.  Between two filaments
    that point the same way it translates the pair instead — it is a spring, and the force it produces
    would be recorded as tension.  So "different filaments" is only half the invariant; the other half
    needs a polarity PER STRAND.

    ⚠ **No population in this arena carries one today, and that is the finding rather than a gap in
    this module.**  ``build/cortex.py``'s ``StrandPopulation`` records ONE ``polarity`` for the whole
    population, so no two of its filaments are anti-parallel — the pair is not unverifiable, it is
    unrepresentable.  ``build/stress_fiber.py`` writes no polarity at all, deliberately (*"writing a
    polarity now would fix it before its law exists"*), while its own docstring quotes *"mixed polarity
    does NOT contract"*.  Filopodium and lamellipodium are uniform by construction and correctly so.

    Making polarity per-strand is not a free edit: it selects ``Strand.tip_node``, which every
    barbed-end law uses as an ADDRESS, so an array would propagate into the growth kernels.  It is an
    A2/G build decision and a PI item, not something a consumer may assume its way past.
    """
    per_strand = census.get("strand_polarity")
    if per_strand is not None:
        pol = np.asarray(per_strand, dtype=np.int64).reshape(-1)
        n_strand = int(census.get("n_strands", pol.size))
        if pol.size != n_strand:
            raise ValueError(
                f"strand_polarity holds {pol.size} entries but the census declares {n_strand} strands. "
                "A polarity table shorter than its population silently orients the tail by index wrap."
            )
        if not bool(np.isin(pol, (-1, 1)).all()):
            raise ValueError(
                "strand_polarity must be +1 or -1 per strand; a third value is not a third orientation."
            )
        return pol, {"kind": "per_strand", "n_strands": n_strand,
                     "n_plus": int((pol > 0).sum()), "n_minus": int((pol < 0).sum())}
    if "polarity" in census:
        return None, {
            "kind": "population_wide",
            "value": int(census["polarity"]),
            "why": (
                f"census kind={census.get('kind')!r} records ONE polarity "
                f"({int(census['polarity'])}) for the whole population, so NO two of its filaments are "
                "anti-parallel — the pair this channel needs is not unverified, it is unrepresentable"
            ),
        }
    return None, {
        "kind": "absent",
        "why": (
            f"census kind={census.get('kind')!r} carries no polarity at all, per-strand or otherwise, "
            "so nothing in it says which way any filament points"
        ),
    }


@dataclass(frozen=True, slots=True)
class ContractileStations:
    """One population's bound contractile channel: the motor's parameters plus its arena addressing.

    Attributes:
        population: the population these stations were validated against.
        plan: the host-side plan this was bound from.
        motor: the resolved Stam-Hocky minifilament, from :mod:`aleph.laws.myosin_linear`.
        links_d: ``(M, 2)`` GLOBAL node indices, int32 — the two anchors each station bridges.
        v_slide_d: ``(M,)`` inter-anchor sliding velocity [um/s].  Zero is ISOMETRIC and returns the
            full stall; the outer loop updates it if it tracks sliding.
    """

    population: str
    plan: dict[str, object]
    motor: MyosinMinifilament
    links_d: wp.array
    v_slide_d: wp.array

    @property
    def n_stations(self) -> int:
        """Number of bipolar minifilaments this channel drives."""
        return int(self.plan["n_stations"])

    def accumulate(self, arena: WorldArena) -> None:
        """Add every station's contractile pair into ``arena``'s force array.  Launches; reads nothing.

        Args:
            arena: the world holding ``position`` and ``force``.  Each station adds ``+f`` to one
                anchor and ``-f`` to the other, so the contribution to the assembled total is zero and
                what is left behind is prestress.
        """
        wp.launch(
            minifilament_kernel,
            dim=self.n_stations,
            inputs=[
                arena.node_arrays["position"], self.links_d, self.v_slide_d,
                self.motor.f_stall_pn, self.motor.v0_um_s,
            ],
            outputs=[arena.node_arrays["force"]],
            device=arena.device,
        )


def plan_contraction(
    census: dict[str, object], stations: npt.ArrayLike | None = None, *,
    motor: MyosinMinifilament | None = None,
) -> dict[str, object]:
    """Validate a station table against a built population's census.  No device is touched.

    Args:
        census: what ``build_stress_fibers`` returned.  Read, never modified.
        stations: ``(M, 2)`` GLOBAL node indices, the two anchors of each bipolar minifilament.
            ``None`` or empty is ACCEPTED here and reported as ``bindable: False`` with
            :data:`STATION_GAP` — so an assembly driver can record the gap in its own artifact — but
            :func:`bind_contraction` will refuse it.
        motor: the resolved minifilament.  Defaults to :func:`aleph.laws.myosin_linear.resolve_myosin`,
            which is a defaulted LAW rather than a defaulted constant: its composition comes from
            ``KB-3.18`` and lives in ``laws/``, and overriding it here is what a sweep does.

    Returns:
        The plan, bindable or not, with the motor's resolved numbers and the gap recorded either way.

    Raises:
        ValueError: on a census that is not a BUILD, a station table of the wrong shape, an index
            outside the population's node claim, a self-link, a duplicate, or a station whose two
            nodes lie on the same strand.
    """
    claims = census.get("claims")
    if not isinstance(claims, dict) or "node" not in claims:
        raise ValueError(
            "census carries no node claim, so it is a PLAN and not a BUILD. Stations validated against "
            "a plan would be validated against nobody's addresses."
        )
    node_lo, node_hi = (int(v) for v in claims["node"])
    if "nodes_per_strand" not in census:
        raise ValueError(
            f"census kind={census.get('kind')!r} carries no nodes_per_strand, so it is not a "
            "uniform-stride strand population and this module cannot tell which filament a node "
            "belongs to. That is not a missing convenience: the ONE invariant a station table has to "
            "satisfy is that it straddles two ANTI-PARALLEL filaments, and without a stride the check "
            "silently passes for a pair inside a single filament. Refusing beats binding a channel "
            "whose only defence has been skipped."
        )
    n_per = int(census["nodes_per_strand"])
    motor = motor or resolve_myosin()
    witness, witness_report = _anti_parallel_witness(census)

    plan: dict[str, object] = {
        "channel": "bipolar_minifilament_contraction",
        "law": "aleph.laws.myosin_linear.minifilament_kernel (KB-3.18, LINEAR FV, PI-ratified 2026-07-07)",
        "kind": str(census.get("kind", "")),
        "population": str(census.get("population", census.get("kind", ""))),
        "node_claim": (node_lo, node_hi),
        "nodes_per_strand": n_per,
        "motor": {
            "n_heads": motor.n_heads, "n_side": motor.n_side, "f_head_pn": motor.f_head_pn,
            "duty": motor.duty, "v0_um_s": motor.v0_um_s,
            "f_stall_pn": motor.f_stall_pn, "f_resting_pn": motor.f_resting_pn,
        },
        "force_velocity": "linear v(F) = v0*(1 - F/F_stall) — NOT Hill; Hill a/F0~0.25 is muscle-only",
        "v_slide_initial": 0.0,
        "v_slide_note": "zero is ISOMETRIC and returns the full stall — the resting prestress, not a null",
        "n_stations_placed_by_builder": int(census.get("n_stations_placed", 0)),
        "sarcomere_um_declared": census.get("sarcomere_um_declared"),
        "anti_parallel_witness": witness_report,
    }

    links = np.zeros((0, 2), dtype=np.int64) if stations is None else np.asarray(stations, dtype=np.int64)
    if links.size == 0:
        links = links.reshape(0, 2)
    if links.ndim != 2 or links.shape[1] != 2:
        raise ValueError(
            f"stations must be (M, 2) node-index pairs; got shape {links.shape}. A minifilament bridges "
            "exactly two anchors — one index is a site and three is a network."
        )
    if links.shape[0] == 0:
        return plan | {"n_stations": 0, "bindable": False, "gap": STATION_GAP}

    inside = (links >= node_lo) & (links < node_hi)
    if not bool(inside.all()):
        bad = links[~inside.all(axis=1)][:3].tolist()
        raise ValueError(
            f"stations reach outside this population's node claim [{node_lo}, {node_hi}); first "
            f"offenders {bad}. Node indices are GLOBAL here, so an out-of-range station silently "
            "contracts whichever population claimed that address."
        )
    if bool((links[:, 0] == links[:, 1]).any()):
        raise ValueError(
            "a station links a node to itself. The kernel returns early on a zero-length axis, so this "
            "does not crash — it produces a motor that is present in the count and absent in the force."
        )
    strand = (links - node_lo) // n_per
    same = strand[:, 0] == strand[:, 1]
    if bool(same.any()):
        bad = links[same][:3].tolist()
        raise ValueError(
            f"{int(same.sum())} station(s) link two nodes of the SAME filament; first offenders {bad}. "
            "A bipolar minifilament bridges two ANTI-PARALLEL filaments and slides them past each "
            "other; inside one filament the same pull is a spring wearing a motor's name, and it would "
            "report as contraction while only shortening a segment."
        )
    if witness is None:
        raise ValueError(
            f"{plan['population']}: {plan['anti_parallel_witness']['why']}\n\nDifferent filaments is "
            "only HALF the invariant. The other half is that they run OPPOSITE ways, and a station "
            "between two co-oriented filaments is the same spring wearing a motor's name that the "
            "same-strand check just refused — one layer up, where the check cannot see it. Refusing "
            "here rather than binding: a co-oriented pair does not fail, it produces a force, and that "
            "force gets recorded as contractile tension."
        )
    opposed = witness[strand[:, 0]] != witness[strand[:, 1]]
    if not bool(opposed.all()):
        bad = links[~opposed][:3].tolist()
        raise ValueError(
            f"{int((~opposed).sum())} station(s) bridge two CO-ORIENTED filaments; first offenders "
            f"{bad}. Their barbed ends point the same way, so the minifilament translates the pair "
            "instead of sliding them past each other — it is a spring, and the artifact would call the "
            "force it produces cortical/fibre tension."
        )
    key = np.sort(links, axis=1)
    if len(np.unique(key, axis=0)) != len(key):
        raise ValueError(
            "the station table holds duplicate anchor pairs. Two stations on one pair is not two "
            "minifilaments, it is one force counted twice — and the count is what the artifact reports."
        )
    if int(links.max()) > np.iinfo(np.int32).max:
        raise ValueError(
            f"node index {int(links.max())} exceeds int32, which is the width the laws/ kernel takes. A "
            "multi-cell arena reaches this; narrowing silently is how it would wrap."
        )

    return plan | {
        "n_stations": int(links.shape[0]),
        "bindable": True,
        "strands_bridged": int(len(np.unique(strand))),
        "prestress_pn_isometric": links.shape[0] * motor.f_stall_pn,
        "gap": None,
    }


def bind_contraction(
    arena: WorldArena, census: dict[str, object], stations: npt.ArrayLike | None = None, *,
    motor: MyosinMinifilament | None = None,
) -> ContractileStations:
    """Upload a validated station table once and return the handle that accumulates it.

    Args:
        arena: the world the population was built into.  Must hold a CUDA allocation.
        census / stations / motor: see :func:`plan_contraction`.

    Returns:
        The bound channel.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_contraction`, and on an EMPTY station table — which
            is where a stress fibre built by A2 lands today, on purpose.  See :data:`STATION_GAP`.
    """
    plan = plan_contraction(census, stations, motor=motor)
    if not plan["bindable"]:
        raise ValueError(
            f"{plan['population']}: no motor stations were supplied, so there is nothing to contract. "
            "This channel refuses to launch over an empty table rather than reporting zero contraction "
            "— 'the motors produced nothing' and 'there are no motors' are different facts, and only "
            f"the second one is true here.\n\n{STATION_GAP}\n\nSupply stations explicitly once the "
            "placement decision lands, or leave this population passive and say so in the artifact."
        )
    if not arena.node_arrays:
        raise RuntimeError(
            f"{plan['population']}: this arena holds no device allocation. Warp CUDA is the only "
            f"runtime and there is no CPU path. The channel it would have bound: {plan['n_stations']} "
            f"minifilaments at F_stall={plan['motor']['f_stall_pn']} pN."
        )

    links = np.asarray(stations, dtype=np.int64).reshape(-1, 2).astype(np.int32)
    return ContractileStations(
        population=str(plan["population"]),
        plan=plan,
        motor=motor or resolve_myosin(),
        links_d=wp.array(links, dtype=wp.int32, ndim=2, device=arena.device),
        v_slide_d=wp.zeros(links.shape[0], dtype=wp.float64, device=arena.device),
    )


def _twin(pos: np.ndarray, links: np.ndarray, v_slide: np.ndarray,
          m: MyosinMinifilament) -> np.ndarray:
    """Host twin of ``minifilament_kernel`` — the arbiter for the sign check, and nothing else.

    It is not a second implementation of the law: it exists so :func:`_demo` can assert the DIRECTION
    of the pair on a machine with no card, and it is checked against the kernel's own source above.
    """
    force = np.zeros_like(pos)
    for t, (i, j) in enumerate(links):
        d = pos[j] - pos[i]
        L = float(np.linalg.norm(d))
        if L < 1e-12:
            continue
        mag = max(m.f_stall_pn * (1.0 - v_slide[t] / m.v0_um_s), 0.0)
        f = (mag / L) * d
        force[i] += f
        force[j] -= f
    return force


def _demo() -> None:
    """Self-check: the SIGN, the empty-table refusal, and every validation — none needs a device."""
    # A stress-fibre census in the shape build_stress_fibers returns, INCLUDING its n_stations_placed: 0.
    sf = {
        "kind": "stress_fiber", "population": "sf_arc", "n_strands": 40, "nodes_per_strand": 61,
        "sarcomere_um_declared": 1.0, "n_stations_placed": 0,
        "claims": {"node": (0, 40 * 61), "segment": (0, 40 * 60)},
    }

    # ── the state the repository is actually in: geometry stands, the channel has nothing to drive ──
    empty = plan_contraction(sf)
    assert empty["bindable"] is False and empty["n_stations"] == 0
    assert "MODEL reason" in str(empty["gap"])
    assert empty["sarcomere_um_declared"] == 1.0     # recorded by the builder, used by nobody
    arena = WorldArena(capacity={})
    try:
        bind_contraction(arena, sf)
    except ValueError as exc:
        assert "different facts" in str(exc) and "N_stations: 0" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an empty station table must refuse the bind, not launch over nothing")

    # ── SIGN — asserted on a station that straddles two filaments, which is the only legal kind ──
    m = resolve_myosin()
    # Filament 0's bead 30 and filament 1's bead 30: global 30 and 61+30, 0.05 um apart in y.
    pos = np.zeros((2 * 61, 3), dtype=np.float64)
    pos[30] = [1.0, 0.00, 0.0]
    pos[61 + 30] = [1.0, 0.05, 0.0]
    links = np.array([[30, 61 + 30]], dtype=np.int64)
    # ⚠ NO built population carries this today (see _anti_parallel_witness). It is supplied here so the
    # accept path and the sign can be pinned at all; that it has to be synthesised IS the finding.
    mixed = sf | {"strand_polarity": [(+1 if k % 2 == 0 else -1) for k in range(40)]}

    iso = _twin(pos, links, np.zeros(1), m)
    sep = pos[links[0, 1]] - pos[links[0, 0]]
    assert float(np.dot(iso[links[0, 0]], sep)) > 0.0     # anchor i is pulled TOWARD j
    assert float(np.dot(iso[links[0, 1]], sep)) < 0.0     # anchor j is pulled TOWARD i  => SHORTENS
    assert np.allclose(iso.sum(axis=0), 0.0)              # equal and opposite => internal, no propulsion
    assert abs(float(np.linalg.norm(iso[links[0, 0]])) - m.f_stall_pn) < 1e-12   # isometric = full stall

    # LINEAR, not Hill: sample the FV law at v/v0 = 0, 1/4, 1/2, 3/4 and require exact affinity. A Hill
    # hyperbola would show curvature here; the second difference of an affine sequence is exactly zero.
    fracs = np.array([0.0, 0.25, 0.5, 0.75])
    mags = np.array([float(np.linalg.norm(_twin(pos, links, np.array([f * m.v0_um_s]), m)[30]))
                     for f in fracs])
    assert np.allclose(np.diff(mags, n=2), 0.0, atol=1e-12)
    assert mags[0] > mags[-1] > 0.0                        # faster sliding => less force

    # Past stall the magnitude CLAMPS at zero: a motor does not actively push its anchors apart.
    over = _twin(pos, links, np.array([2.0 * m.v0_um_s]), m)
    assert np.allclose(over, 0.0)

    # ── every validation ────────────────────────────────────────────────────────────────────────
    def refuses(fragment: str, stations: object) -> None:
        try:
            plan_contraction(mixed, stations)                         # type: ignore[arg-type]
        except ValueError as exc:
            assert fragment in str(exc), f"expected {fragment!r} in: {exc}"
        else:  # pragma: no cover
            raise AssertionError(f"must refuse: {fragment}")

    good = plan_contraction(mixed, links)
    assert good["bindable"] is True and good["n_stations"] == 1 and good["strands_bridged"] == 2
    assert good["prestress_pn_isometric"] == m.f_stall_pn
    assert good["anti_parallel_witness"]["kind"] == "per_strand"

    refuses("exactly two anchors", np.zeros((2, 3), dtype=np.int64))
    refuses("outside this population's node claim", np.array([[30, 40 * 61]]))
    refuses("links a node to itself", np.array([[30, 30]]))
    refuses("SAME filament", np.array([[30, 31]]))         # both on filament 0 — a spring, not a motor
    refuses("counted twice", np.array([[30, 91], [91, 30]]))

    # A census that is a plan and not a build has no addresses to validate against.
    try:
        plan_contraction({k: v for k, v in mixed.items() if k != "claims"}, links)
    except ValueError as exc:
        assert "a PLAN and not a BUILD" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a plan census must refuse")

    # ── the OTHER populations this channel will meet, now that G's minifilaments landed on cortex ──
    # A cortical shell IS uniform-stride (build/cortex.py StrandPopulation.record), so it binds through
    # the same door on its own claim — the same-strand check stays meaningful because "which filament"
    # is still answerable.
    cortex = {
        "kind": "strand_population", "population": "cortex", "n_strands": 60000,
        "nodes_per_strand": 70, "polarity": +1, "claims": {"node": (0, 60000 * 70)},
    }
    # The stride IS there, so the same-strand half of the invariant still works on cortex.
    try:
        plan_contraction(cortex, np.array([[5, 6]]))      # inside one cortical filament
    except ValueError as exc:
        assert "SAME filament" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("the straddle check must survive on cortex too")
    # But the OTHER half does not: StrandPopulation records one polarity for 60,000 filaments, so a
    # station between two DIFFERENT cortical filaments is still a co-oriented pair. This is the case
    # that read `bindable: True` until 2026-08-20 and would have recorded a spring as cortical tension.
    assert plan_contraction(cortex)["anti_parallel_witness"]["kind"] == "population_wide"
    try:
        plan_contraction(cortex, np.array([[5, 70 + 5]]))
    except ValueError as exc:
        assert "unrepresentable" in str(exc) and "only HALF the invariant" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a population-wide polarity cannot witness an anti-parallel pair")

    # sf_arc writes no polarity at all, deliberately — so it is refused on this layer too, and the
    # artifact can tell the two reasons apart: (e) 5 is about STATIONS, this is about ORIENTATION.
    assert plan_contraction(sf)["anti_parallel_witness"]["kind"] == "absent"
    try:
        plan_contraction(sf, links)
    except ValueError as exc:
        assert "no polarity at all" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a census with no polarity cannot witness an anti-parallel pair")

    # A polarity table that does not cover its population wraps by index; refused rather than trusted.
    for bad, frag in ((sf | {"strand_polarity": [1, -1]}, "shorter than its population"),
                      (sf | {"strand_polarity": [0] * 40}, "third orientation")):
        try:
            plan_contraction(bad, links)
        except ValueError as exc:
            assert frag in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"must refuse: {frag}")

    # A closed surface does NOT, and the refusal names why rather than raising KeyError: without a
    # stride the anti-parallel straddle check silently passes, and that check is the only defence here.
    try:
        plan_contraction({"kind": "closed_surface", "population": "membrane",
                          "claims": {"node": (0, 163842)}}, np.array([[5, 900]]))
    except ValueError as exc:
        assert "only defence has been skipped" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a surface census must refuse, and say why")

    # And a device-less arena refuses a BINDABLE table while reporting what it would have bound.
    try:
        bind_contraction(WorldArena(capacity={}), mixed, links)
    except RuntimeError as exc:
        assert "1 minifilaments" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less arena must refuse the bind")

    print("contraction: sign, linearity and refusals OK — and sf_arc stays unbindable, by decision")


if __name__ == "__main__":
    _demo()
