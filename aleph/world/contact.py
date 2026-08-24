"""Non-penetration constraints — the destination for the CONTACT edges that are not bonds.

PI decision 1, 2026-08-21 (``docs/v2_audit/PI_DECISIONS_2026-08-21.md``): the cardless CONTACT edges
are **not bonds — they are constraints**, and ``world/`` gets a place for them. This module is that
place. It holds the DATA TYPE and the CONTRACT and nothing else: no kernel, no neighbour query, no
force, no penetration depth, no magnitude. What it can say is *"these edges arrive here, and this type
can express them"*; what it cannot say is what any of them is worth.

WHY A SECOND PRIMITIVE AT ALL, RATHER THAN A FLAG ON ``BondFamily``.  Three independent reasons, and
each one alone is sufficient:

* **The missing chemistry card is the identity, not the gap.** ``BondFamily`` takes ``chemistry_card``
  as a required field because that is what a family IS. A steric exclusion has no molecule — nothing
  binds and nothing unbinds — so the field cannot be filled in, ever, and a type whose required field
  is permanently unanswerable is the wrong type.
* **``BondCount`` demands a population, and a contact set is not one.** It asks *how many, against
  what, for which cell, on whose authority* and refuses without all four. A contact count is a function
  of the CONFIGURATION: it changes every step, it is an OUTPUT of a broad phase, and answering "how
  many" for it would produce one number per membrane quadrature point — the icosphere vertex count,
  which is the exact defect ``bond.py`` is shaped around, re-entered through a different door.
* **The law is unilateral and the sign is structural.** ``bond.py`` stores one ``rest_um`` and one
  stiffness per bond and has no side; that is a bilateral spring by construction. A bilateral spring at
  the contact distance lets the membrane PULL the cortex back, and ``engine/connector_joints.py:27-30``
  names that: *"a contact that pulls is an adhesion nobody declared. The sign is therefore structural,
  gated on the gap, and cannot be turned off by a stiffness value."*

THE FORMULATION, AND WHY THE STIFFNESS IS OPTIONAL.  aLENS (``docs/v2_audit/
CROSSLINK_STIFFNESS_AXIS_2026-08-16.md`` §4) puts this physics at **unilateral complementarity**::

    0 <= Phi_u   ⊥   gamma_u >= 0

— the gap may not be negative, the contact multiplier may not be negative, and their product is zero:
no interpenetration, no adhesion, and no force at a distance. **There is no stiffness in that
statement.** ``gamma`` is a Lagrange multiplier the solve produces, not ``k`` times an overlap.

The incumbent runtime is instead the REGULARISED form: a one-sided Hookean penalty below the contact
distance, whose stiffness ``contracts.py:950-953`` describes as *"a numerical exclusion stiffness
derived from the existing WCA / derived-mobility pattern"*. That is an honest description of a
numerical parameter and a dishonest one of a physiological magnitude, so this type keeps them apart:
:attr:`ContactFamily.regularisation_stiffness_pn_per_um` is ``None`` for the complementarity form and
carries the penalty constant otherwise, where ``None`` is the physics and a number is an approximation
of it. A result standing on a number there is standing on a regularisation and must say so.

**The one physical magnitude a contact has is its CONTACT DISTANCE** ``d0`` — the separation at which
the force vanishes — and :class:`ContactDistance` demands the same scope and provenance ``BondCount``
does, for the same reason. It is not currently available for any of the edges below; see them.

THE IDENTITY QUESTION, WHICH IS THE REAL ONE.  ``families.DUPLICATE_CRITERION`` says a bond family's
identity is its CHEMISTRY CARD, and that a shared runtime class or a shared force law is not a
duplicate. Cardless edges cannot be judged by it — the criterion returns nothing to compare, which is
itself evidence they are not bonds. The analogue for a constraint is :data:`CONTACT_IDENTITY`:

    **A constraint's identity is its GAP FUNCTION** — which primitive kinds are measured, and by what
    metric. The surface pair is SCOPE, not identity. The force law is downstream and decides nothing.

Read that against the bond criterion and the two are the same rule pointed at different objects. For a
bond, the physics lives in the chemistry and the geometry is incidental, so identity is the card. For a
constraint there is no chemistry, and the physics IS the geometry — *what is being kept out of what,
measured how* — so identity is Phi. In both cases what is rejected is the same thing: the declared
component pair, which is bookkeeping in one and scope in the other, and identity in neither.

⚠ AND THE ANTI-SYMMETRY IS DELIBERATE.  ``bond.py`` refuses to store the component pair, because
declaring it is what let ONE declaration stand in for a POPULATION. A :class:`ContactScope` declares
its pair openly, and that is safe here for a reason specific to constraints: there is no population to
hide. The active set is regenerated per accepted step and its size is reported by the step, so the
declaration cannot absorb a count — **this type has no field a count fits in**, and that is enforced by
its shape rather than by a comment.

WHAT THIS MEANS FOR THE THREE EDGES PI DECISION 1 NAMES.  ``connector_joints.py:7-13`` lists FIVE
declared CONTACT edges and splits them 3/2 on ``kinetics``, and the split does not fall where the
decision's wording does:

===============================  ========  =========================================
edge                             kinetics  chemistry card
===============================  ========  =========================================
``membrane_cortex_contact``      False     none
``membrane_ecm_contact``         False     none
``nucleus_cortex_contact``       False     none
``lamellipodium_membrane_contact``  True   ``actin_membrane_brownian_ratchet_contact``
``filopodium_membrane_tip``      True      ``actin_membrane_brownian_ratchet_contact``
===============================  ========  =========================================

The decision names ``lamellipodium_membrane_contact`` as one of *"the three cardless CONTACT edges"*;
by the frozen source the third cardless edge is ``membrane_ecm_contact``. **The ruling still holds for
all four**, and it holds for a reason better than a correction: the unilateral gap law is SHARED by all
five, and ``DUPLICATE_CRITERION`` says verbatim that a shared force law is not a duplicate. So the
ratchet edges are BOTH — a non-penetration constraint (a barbed end may push a membrane and may never
tow it, under every reading) AND, if the ratchet is the *tethered* one, a bond population on top of it.
The constraint half is unconditional and lands here. The bond half stays open, where its own module
already left it. That is not the decision being wrong; it is the decision being about a shape, and one
edge having two.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``d0`` [µm]; ``Phi`` [µm]; ``gamma`` [pN]; a regularisation stiffness [pN/µm]; the
    complementarity product [pN·µm]. :func:`check_complementarity` carries the two tolerances with
    their units in the argument names and DERIVES the product tolerance from them rather than taking a
    third — a chosen product tolerance would be a threshold picked after seeing a residual.
  * boundary — an empty constraint set is legal and is not a failure: a configuration in which nothing
    touches satisfies the complementarity statement vacuously, and that is the physical answer. A
    family with zero scopes is rejected, because that is a declaration of nothing.
  * conservation/invariant — no arena claim is taken and none may be. A constraint is not an
    allocation: ``arena.py:76-78`` closed the ``Kind`` set on exactly this argument, that a ``SITE`` is
    *"an address resolved from live geometry, not an allocation"*. Asserted in :func:`_demo`.
  * CFL/precision — no integration here. Recorded for whoever writes the solve: the complementarity
    form is what removes the contact stiffness from the CFL bound entirely, which is the whole reason
    aLENS reports timesteps *"two or more orders of magnitude larger"*; the regularised form does not,
    and its stiffness enters the bound through the row sum of its endpoints exactly as a bond's does.
  * gate soundness — **A TOLERANCE MUST HAVE THE SAME EXTENSIVITY AS THE QUANTITY IT BOUNDS.**
    :func:`check_complementarity` is written to satisfy that rule and this is the reason, recorded
    here rather than left in the test that proves it. Every quantity in the row is a MAX-NORM and the
    derived product tolerance is ``gap_tol·max|gamma| + force_tol·max|Phi|`` — magnitudes only, with
    **no term in the constraint count**. So all of them are INTENSIVE in ``C``, and a step that
    activates a million contacts is judged by the same numbers as one that activates ten.

    ⚠ That matters here more than anywhere else in the engine, because a contact set's SIZE IS A
    FUNCTION OF THE CONFIGURATION BEING JUDGED. A tolerance carrying a count term would be enlarged
    by the very events the gate exists to rule on, and the gate would loosen exactly when it is
    being tested hardest. The engine has a measured instance of that failure: ``ledger.py``'s
    balance gate bounds a max-norm residual with a Higham SUM over force terms, and on 2026-08-21
    its tolerance grew **626x** across configurations while the residual it bounds did not move at
    all. Same engine, two gates, opposite extensivity — and the intensive one is the correct design.
    Copy this one, not that one; and note the property is INVISIBLE IN A DIFF, so it is held by
    ``test_contact.py``'s replication check rather than by anyone remembering it.
  * sign sense — the sign is the contract. ``Phi >= 0`` separation is admissible, ``Phi < 0``
    interpenetration is not; ``gamma >= 0`` pushes apart, ``gamma < 0`` is an undeclared adhesion.
    Both are checked, separately, so a violation says which one happened.
  * measurement protocol — host-side declaration only; no device is touched and nothing is uploaded.
    :func:`check_complementarity` reads arrays a caller supplies and never fetches them itself.

engine units: length µm, force pN, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.world.bond import SourceClass

__all__ = [
    "CONTACT_IDENTITY",
    "GAP_LAW_POINT_SURFACE",
    "ContactDistance",
    "ContactScope",
    "ContactFamily",
    "check_complementarity",
]

#: The duplicate criterion for edges that have no chemistry card, and the answer to the question
#: ``families.DUPLICATE_CRITERION`` cannot be asked about them. See the module docstring for the
#: argument; the short form is that a bond's physics is its chemistry and a constraint's physics is its
#: geometry, so each is identified by the one it actually has.
CONTACT_IDENTITY = (
    "gap function — which primitive kinds are measured and by what metric. The surface pair is SCOPE, "
    "not identity; the force law is downstream and decides nothing."
)

#: The only gap law ``world/`` has a use for today: the signed Euclidean separation between a point and
#: the surface it may not pass through, minus the contact distance. Every one of the five declared
#: CONTACT edges resolves to this one law (``engine/connector_joints.py:7-13`` lists them as *"the same
#: unilateral gap law"*), which under :data:`CONTACT_IDENTITY` makes them ONE family with several
#: scopes rather than several families — the merge ``bond.py`` performs with a shared card, performed
#: here with a shared Phi.
GAP_LAW_POINT_SURFACE = "unilateral_point_surface_gap"


@dataclass(frozen=True, slots=True)
class ContactDistance:
    """The separation at which a contact force vanishes — the ONE physical magnitude a contact has.

    Every field is required, for the reason :class:`~aleph.world.bond.BondCount` gives: a magnitude
    labelled physiological is physiological FOR SOME PARTICULAR cell, and without saying which, a
    sourced value cannot be told apart from a proxy.

    ⚠ This is NOT a rest length and must not be read as one. ``engine/connector_joints.py:56-58``
    records that a contact joint packs the contact distance into the SoA slot the bilateral law reads
    as ``rest_um``, which is a storage decision and not a claim that the two are the same quantity: a
    rest length is where a bilateral force is zero and either sign is admissible around it, while a
    contact distance is where a unilateral force STARTS and only one side of it exists at all.

    Attributes:
        value_um: the contact distance [µm]. Must be finite and strictly positive — a contact distance
            of zero says two surfaces may touch at a point with no force, which is not a regularised
            contact but the absence of one.
        scope: the cell line, state, assay and temperature the value is true FOR.
        source_class: see :class:`~aleph.world.bond.SourceClass`, imported rather than re-spelled.
        provenance: the citation or the derivation. Required and non-empty.
    """

    value_um: float
    scope: str
    source_class: SourceClass
    provenance: str

    def __post_init__(self) -> None:
        if not np.isfinite(self.value_um) or self.value_um <= 0.0:
            raise ValueError(
                f"contact distance must be finite and strictly positive; got {self.value_um!r}. A zero "
                "contact distance is the absence of a contact, not a contact with no gap."
            )
        if not str(self.scope).strip():
            raise ValueError(
                "scope is required: the membrane/cortex apposition gap and the perinuclear gap are "
                "different lengths of different cells, and a value with no scope cannot say which."
            )
        if not str(self.provenance).strip():
            raise ValueError(
                "provenance is required and must be non-empty. A magnitude with no stated origin is "
                "how a fixture value acquires a physiological label."
            )


@dataclass(frozen=True, slots=True)
class ContactScope:
    """One surface pair a gap law is enforced over, and what is still unanswered about it.

    A scope is a DECLARATION of where the broad phase looks. It is not a population and holds no pair
    list, no count and no capacity: those are properties of the configuration at a step, not of the
    declaration, and this class deliberately has nowhere to put them.

    **THE TWO SIDES ARE NOT INTERCHANGEABLE, AND THEY ARE NAMED.** ``GAP_LAW_POINT_SURFACE`` says so
    in its own name: one side supplies the SURFACE that may not be passed through, the other supplies
    the POINTS that may not pass through it. The first version of this class stored an unordered-
    looking ``populations`` pair and recorded neither role, and the three declared scopes promptly
    disagreed about the order — ``(membrane, cortex)`` and ``(nuclear_envelope, cortex)`` were written
    surface-first while ``(lamellipodium, membrane)`` was written point-first, with nothing noticing.
    Naming the roles is what makes an illegal endpoint sayable at all.

    ⚠ **AND THE SURFACE SIDE MUST OWN FACES — this class cannot check it, so it states it.** A
    non-penetration constraint is about a REGION a point may not enter, and ``surface.py`` records
    that a FACE is *"the only primitive whose measure is an AREA"*. A population claiming ``NODE``
    only bounds no region and cannot be a surface. The live instance is the nuclear lamina:
    ``build/lamina.py`` claims ``Kind.NODE`` and nothing else — 2,043 points whose
    ``n_filaments_implied`` is ``0.5 · n_nodes · LAMIN_CONNECTIVITY``, an implied count rather than a
    built topology — so it cannot be a surface here, while ``build/envelope.py`` claims ``NODE`` +
    ``FACE`` + ``ANGLE4`` and can. Checking it needs a population's primitive census, which belongs to
    the builders and not to a declaration type; this is why the requirement is written and not coded.

    ⚠ **A point set is no better as the POINT side, and that half is the stronger one.** A multiplier
    enters at a node. A node with no incident SEGMENT and no incident FACE has nothing to transmit it
    to, so the force is delivered to a free particle and leaves the structure: the constraint would
    push individual lamin nodes and the nucleus would feel nothing. An unconnected population is a
    load SINK, not a load PATH, and is therefore outside the contact system in BOTH roles — for two
    different reasons. Recorded here because PI decision 9's LINC/SUN coupling assumes a nuclear side
    that can carry load. See PI queue item 17.

    Attributes:
        connector: the declared connector name, exactly as the device-run census spells it.
        point_population: the population whose points are kept OUT of the other. Any population with
            addressable nodes may serve, subject to the load-path caveat above.
        surface_population: the population that supplies the surface. Must own FACEs.
        d0: the contact distance, or ``None`` when it is not yet available — which is the state all
            four contact edges are in today. ``None`` is a census entry, not a placeholder value.
        blocked_by: what is unanswered. Empty is only legal once ``d0`` exists.
    """

    connector: str
    point_population: str
    surface_population: str
    d0: ContactDistance | None = None
    blocked_by: tuple[str, ...] = ()

    @property
    def populations(self) -> tuple[str, str]:
        """``(point, surface)`` — ordered by ROLE, never by declaration convenience."""
        return (self.point_population, self.surface_population)

    def __post_init__(self) -> None:
        if not str(self.connector).strip():
            raise ValueError("a contact scope needs the connector name it answers for")
        if not str(self.point_population).strip() or not str(self.surface_population).strip():
            raise ValueError(
                f"{self.connector}: both sides must be named. A gap is measured FROM points TO a "
                "surface, and a side left blank is a role nobody assigned."
            )
        if self.point_population == self.surface_population:
            raise ValueError(
                f"{self.connector}: {self.point_population!r} is named as both the point side and the "
                "surface side. A population cannot be kept out of itself — that is self-avoidance, "
                "which is a property of one object's own law, not a contact between two."
            )
        if self.d0 is None and not self.blocked_by:
            raise ValueError(
                f"{self.connector}: no contact distance and nothing named as missing. A scope with no "
                "d0 is blocked, and a blocked scope that will not say on what is a refusal with no "
                "content — the queue would carry a count instead of a question."
            )
        if self.d0 is not None and self.blocked_by:
            raise ValueError(
                f"{self.connector}: a d0 is supplied AND blockers are named. Say which: an answered "
                "scope with open questions is how a partial answer gets read as a whole one."
            )

    @property
    def enforceable(self) -> bool:
        """True when the scope has its contact distance. Says nothing about whether it is RIGHT."""
        return self.d0 is not None


@dataclass(frozen=True, slots=True)
class ContactFamily:
    """One gap law, enforced over one or more surface pairs. The constraint counterpart of a family.

    Attributes:
        gap_law: the identity — see :data:`CONTACT_IDENTITY`. Two declared edges computing the same Phi
            over different surfaces are ONE family with two scopes.
        scopes: the surface pairs this law is enforced over. At least one.
        regularisation_stiffness_pn_per_um: ``None`` for the complementarity form
            ``0 <= Phi ⊥ gamma >= 0``, in which no stiffness exists because ``gamma`` is a multiplier.
            A number is the REGULARISED (one-sided penalty) approximation of it, and is a numerical
            parameter rather than a physiological one — a result standing on it carries that.
    """

    gap_law: str
    scopes: tuple[ContactScope, ...]
    regularisation_stiffness_pn_per_um: float | None = None

    def __post_init__(self) -> None:
        if not str(self.gap_law).strip():
            raise ValueError(
                "a contact family is identified by its gap law and cannot be declared without one; "
                f"the identity criterion is: {CONTACT_IDENTITY}"
            )
        if not self.scopes:
            raise ValueError(f"{self.gap_law}: a family enforced over no surface pair declares nothing")
        seen: set[frozenset[str]] = set()
        for s in self.scopes:
            key = frozenset(s.populations)
            if key in seen:
                raise ValueError(
                    f"{self.gap_law}: two scopes over the same population pair {sorted(key)}. Under "
                    "the identity criterion those are the same constraint declared twice, which is "
                    "the duplication this shape exists to make visible."
                )
            seen.add(key)
        k = self.regularisation_stiffness_pn_per_um
        if k is not None and (not np.isfinite(k) or k <= 0.0):
            raise ValueError(
                f"{self.gap_law}: a regularisation stiffness must be finite and positive; got {k!r}. "
                "Use None for the complementarity form, where no stiffness exists at all — zero is "
                "neither, and would silently disable the constraint."
            )

    @property
    def is_complementarity(self) -> bool:
        """True when the family is stated as ``0 <= Phi ⊥ gamma >= 0`` with no penalty stiffness."""
        return self.regularisation_stiffness_pn_per_um is None

    @property
    def enforceable(self) -> bool:
        """True when every scope has its contact distance."""
        return all(s.enforceable for s in self.scopes)

    def blocked(self) -> dict[str, tuple[str, ...]]:
        """The open questions per connector, for Lead's queue. Empty when the family is enforceable."""
        return {s.connector: s.blocked_by for s in self.scopes if s.blocked_by}

    def provenance_row(self) -> dict[str, object]:
        """The artifact row — the law, its scopes, and whether a regularisation is in play.

        Together on purpose, and the last field is why: a contact result read without knowing whether a
        penalty stiffness was in the loop cannot be told apart from one that was not.
        """
        return {
            "gap_law": self.gap_law,
            "identity_criterion": CONTACT_IDENTITY,
            "formulation": "complementarity" if self.is_complementarity else "regularised_penalty",
            "regularisation_stiffness_pn_per_um": self.regularisation_stiffness_pn_per_um,
            "scopes": [
                {
                    "connector": s.connector,
                    "populations": list(s.populations),
                    "d0_um": None if s.d0 is None else s.d0.value_um,
                    "d0_scope": None if s.d0 is None else s.d0.scope,
                    "d0_source_class": None if s.d0 is None else str(s.d0.source_class),
                    "d0_provenance": None if s.d0 is None else s.d0.provenance,
                    "blocked_by": list(s.blocked_by),
                }
                for s in self.scopes
            ],
        }


def check_complementarity(
    phi_um: npt.ArrayLike,
    gamma_pn: npt.ArrayLike,
    *,
    gap_tol_um: float,
    force_tol_pn: float,
) -> dict[str, float]:
    """Score ``0 <= Phi ⊥ gamma >= 0`` on a supplied configuration, and raise if it does not hold.

    This is the CONTRACT, not a solve: it says what "the constraint is satisfied" means, so a later
    step machinery has something written before the run to be scored against. It computes no force and
    proposes no correction.

    The product tolerance is DERIVED from the two supplied tolerances rather than taken as a third
    argument. ``|Phi·gamma|`` can be non-zero at first order for two admissible reasons — a gap within
    ``gap_tol_um`` carrying a real multiplier, or a multiplier within ``force_tol_pn`` at a real gap —
    so the admissible product is ``gap_tol_um·max|gamma| + force_tol_pn·max|Phi|`` and nothing else. A
    tolerance chosen after seeing a residual would be a threshold edited to fit the data.

    Args:
        phi_um: ``(C,)`` gaps [µm], one per active constraint. Positive is separation.
        gamma_pn: ``(C,)`` contact multipliers [pN], one per active constraint. Positive pushes apart.
        gap_tol_um: how far below zero a gap may sit and still be called non-penetrating [µm].
        force_tol_pn: how far below zero a multiplier may sit and still be called non-adhesive [pN].

    Returns:
        The artifact row a gate records: the three residuals — ``penetration_um`` (the deepest
        negative gap, as a positive depth), ``adhesion_pn`` (the most negative multiplier, as a
        positive magnitude) and ``product_pn_um`` (the largest ``|Phi·gamma|``) — **each beside the
        threshold that judged it** (``gap_tol_um``, ``force_tol_pn``, ``product_tol_pn_um``), plus
        ``n_constraints``. The thresholds travel with the values because two of them are ARGUMENTS
        rather than properties of this file: a provenance stamp fixes *which code ran* and cannot
        recover *what it was given*, so a row without them turns a PASS into an unfalsifiable claim.

    Raises:
        ValueError: on mismatched shapes, a non-finite entry, or a negative tolerance.
        AssertionError: if any of the three residuals exceeds its tolerance, naming which one — the
            three fail for different physical reasons and a single verdict would hide which happened.
    """
    phi = np.asarray(phi_um, np.float64).reshape(-1)
    gamma = np.asarray(gamma_pn, np.float64).reshape(-1)
    if phi.shape != gamma.shape:
        raise ValueError(f"Phi and gamma must pair up one-to-one; got {phi.shape} and {gamma.shape}")
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(gamma))):
        raise ValueError("Phi and gamma must be finite; a non-finite entry is a solve that diverged")
    if gap_tol_um < 0.0 or force_tol_pn < 0.0:
        raise ValueError(f"tolerances must be nonnegative; got {gap_tol_um!r}, {force_tol_pn!r}")

    # An empty active set satisfies the statement vacuously — nothing touches, and that is an answer.
    penetration = float(max(0.0, -phi.min())) if phi.size else 0.0
    adhesion = float(max(0.0, -gamma.min())) if gamma.size else 0.0
    product = float(np.max(np.abs(phi * gamma))) if phi.size else 0.0
    product_tol = (
        gap_tol_um * float(np.max(np.abs(gamma))) + force_tol_pn * float(np.max(np.abs(phi)))
        if phi.size
        else 0.0
    )

    # EVERY threshold a verdict rests on, beside the value it judges. The two sign tolerances are
    # ARGUMENTS, not properties of this file, so no amount of provenance stamping recovers them: a
    # row carrying `penetration_um` without `gap_tol_um` cannot distinguish a PASS at 1e-9 um from a
    # PASS at 1e-3, and a PASS that carries no trace of what it was measured against is worse than a
    # wrong number, because the number can at least be checked. Caught by Lead reading this function
    # rather than its description, 2026-08-21.
    row = {
        "n_constraints": float(phi.size),
        "penetration_um": penetration,
        "gap_tol_um": float(gap_tol_um),
        "adhesion_pn": adhesion,
        "force_tol_pn": float(force_tol_pn),
        "product_pn_um": product,
        "product_tol_pn_um": product_tol,
    }
    if penetration > gap_tol_um:
        raise AssertionError(
            f"non-penetration violated: deepest gap {-phi.min():.6g} um exceeds the tolerance "
            f"{gap_tol_um:.6g} um. Two surfaces are inside each other."
        )
    if adhesion > force_tol_pn:
        raise AssertionError(
            f"unilaterality violated: multiplier {gamma.min():.6g} pN is below -{force_tol_pn:.6g} pN. "
            "A contact that pulls is an adhesion nobody declared."
        )
    if product > product_tol:
        raise AssertionError(
            f"complementarity violated: max |Phi*gamma| = {product:.6g} pN*um exceeds the derived "
            f"tolerance {product_tol:.6g} pN*um. A separated pair is carrying force at a distance."
        )
    return row


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Self-check: a count has nowhere to live, the identity merges by gap law, and the sign has teeth."""
    import dataclasses

    # 1. THE COUNT IS UNREPRESENTABLE, and that is the design. If a field for one ever appears, this
    #    fails — which is the only durable way to keep the ERM defect out of a type that has no builder.
    fields = {f.name for f in dataclasses.fields(ContactScope)} | {
        f.name for f in dataclasses.fields(ContactFamily)
    }
    assert not {"count", "n_contacts", "pairs", "capacity", "density"} & fields, sorted(fields)

    blocked = ContactScope(
        connector="membrane_cortex_contact",
        point_population="cortex",
        surface_population="membrane",
        blocked_by=("contact distance d0 [um] is registered nowhere for any cell line",),
    )
    assert not blocked.enforceable
    assert blocked.populations == ("cortex", "membrane"), "populations reads (point, surface)"

    # 1b. THE ROLES ARE THE POINT. A population may not be kept out of itself, and a side may not be
    #     left unnamed — both were expressible before the roles existed.
    for kwargs, expect in (
        ({"surface_population": "cortex"}, "cannot be kept out of itself"),
        ({"surface_population": "  "}, "both sides must be named"),
    ):
        try:
            dataclasses.replace(blocked, **kwargs)
        except ValueError as exc:
            assert expect in str(exc), (kwargs, exc)
        else:  # pragma: no cover
            raise AssertionError(f"ContactScope({kwargs}) must refuse")

    # A scope must say what is missing, and may not say both at once.
    for kwargs, expect in (
        ({"blocked_by": ()}, "refusal with no content"),
        ({"d0": ContactDistance(0.03, "fixture", SourceClass.CONVENIENCE, "test_connector_joints.py")},
         "Say which"),
    ):
        try:
            dataclasses.replace(blocked, **kwargs)
        except ValueError as exc:
            assert expect in str(exc), (kwargs, exc)
        else:  # pragma: no cover
            raise AssertionError(f"ContactScope({kwargs}) must refuse")

    # A contact distance refuses to be a bare number, exactly as BondCount does.
    for kwargs, expect in (
        ({"value_um": 0.0}, "absence of a contact"),
        ({"scope": " "}, "scope is required"),
        ({"provenance": ""}, "provenance is required"),
    ):
        base = dict(value_um=0.03, scope="MCF7 interphase adherent 37C",
                    source_class=SourceClass.PI_GAP, provenance="no Contract-Graph datum")
        try:
            ContactDistance(**{**base, **kwargs})
        except ValueError as exc:
            assert expect in str(exc), (kwargs, exc)
        else:  # pragma: no cover
            raise AssertionError(f"ContactDistance({kwargs}) must refuse")

    # 2. THE IDENTITY MERGES BY GAP LAW. Four declared edges, one Phi, one family, four scopes.
    #    Read as (point, surface): the surface side is the one that owns FACEs in every case.
    edges = (
        ("membrane_cortex_contact", "cortex", "membrane"),
        ("membrane_ecm_contact", "ecm", "membrane"),
        ("nucleus_cortex_contact", "cortex", "nuclear_envelope"),
        ("lamellipodium_membrane_contact", "lamellipodium", "membrane"),
    )
    fam = ContactFamily(
        gap_law=GAP_LAW_POINT_SURFACE,
        scopes=tuple(
            ContactScope(connector=c, point_population=pt, surface_population=sf,
                         blocked_by=("d0 unsourced",))
            for c, pt, sf in edges
        ),
    )
    assert fam.is_complementarity, "no stiffness declared means the multiplier form, not a free penalty"
    assert not fam.enforceable and len(fam.blocked()) == 4
    row = fam.provenance_row()
    assert row["formulation"] == "complementarity" and row["regularisation_stiffness_pn_per_um"] is None

    # The same pair declared twice is the duplication this shape makes visible.
    try:
        ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=fam.scopes + (fam.scopes[0],))
    except ValueError as exc:
        assert "declared twice" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a repeated population pair must refuse")

    # A regularisation stiffness of zero is neither form, and would disable the constraint silently.
    try:
        ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=fam.scopes[:1],
                      regularisation_stiffness_pn_per_um=0.0)
    except ValueError as exc:
        assert "silently disable" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a zero regularisation stiffness must refuse")

    # 3. THE CONTRACT HAS TEETH. Complementarity holds exactly when each pair is separated OR loaded.
    ok = check_complementarity([0.5, 0.0, 0.0], [0.0, 12.0, 0.0], gap_tol_um=1e-9, force_tol_pn=1e-9)
    assert ok["penetration_um"] == 0.0 and ok["product_pn_um"] == 0.0 and ok["n_constraints"] == 3.0
    # An empty active set is a physical answer, not a failure: nothing is touching.
    assert check_complementarity([], [], gap_tol_um=0.0, force_tol_pn=0.0)["n_constraints"] == 0.0

    for phi, gamma, expect in (
        ([-0.02], [5.0], "non-penetration violated"),   # inside each other
        ([0.10], [-5.0], "unilaterality violated"),     # a contact that pulls
        ([0.10], [5.0], "complementarity violated"),    # force at a distance
    ):
        try:
            check_complementarity(phi, gamma, gap_tol_um=1e-9, force_tol_pn=1e-9)
        except AssertionError as exc:
            assert expect in str(exc), (phi, gamma, exc)
        else:  # pragma: no cover
            raise AssertionError(f"({phi}, {gamma}) must be rejected as {expect}")

    # 4. NO ARENA CLAIM, EVER. A constraint is an address resolved from live geometry (arena.py:76-78),
    #    so declaring one may not consume capacity a real population will need.
    from aleph.world.arena import Kind, WorldArena

    arena = WorldArena(capacity={Kind.NODE: 8, Kind.SEGMENT: 8, Kind.ANGLE3: 8, Kind.BOND: 8})
    ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=fam.scopes)
    assert arena.n_live(Kind.BOND) == 0 and arena.n_live(Kind.NODE) == 0
    arena.assert_partitioned()

    print(f"contact self-check OK — {fam.gap_law}: {len(fam.scopes)} scopes, "
          f"{len(fam.blocked())} blocked, formulation={row['formulation']}")


if __name__ == "__main__":
    _demo()
