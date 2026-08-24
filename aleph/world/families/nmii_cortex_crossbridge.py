r"""``nmii_cortex_motor`` — the head->actin crossbridge, this engine's FIRST kinetic bond family.

PI decision 2 of 2026-08-21 chose this edge to open the kinetics machinery, **overruling this repo's
own recommendation of ERM**, and recorded why: a processive motor needs a rate law, a duty ratio and
catch-slip detachment all written for the first time, so ``bond.py``'s deferred kinetic interface gets
designed against the most demanding consumer rather than the easiest. That is the standard this module
should be read against — not "does it build", but "what did it MEASURE about the interface".

⚠ FILE NAME vs CONNECTOR NAME.  The declared census name is ``nmii_cortex_motor``
(``engine/contracts.py:883``); the file is named for the physics because three other declared edges
carry the same one. :data:`CONNECTOR` is the census name, so Lead's queue matches.

WHAT THIS MODULE MEASURED ABOUT ``bond.py``'S DEFERRED INTERFACE — the actual product
-----------------------------------------------------------------------------------
``bond.py`` defers "no kinetic law, no attach/detach, no free list, no snapshot twins ... They arrive
with the first family that has kinetics, which is also when the right shape for them becomes
measurable rather than guessable." This is that family, and three of the four answers are now
measured rather than guessed:

1. **A free list is unnecessary here, and would be wrong.** The bond population IS the head
   population: a head has exactly one crossbridge slot and that slot is its identity. Detaching frees
   no storage and rebinding allocates none, so a free list would be bookkeeping over a static array.
   ⚠ This answer is about THIS family. ERM sheds tethers and does not obviously share it.
2. **Snapshot twins are unnecessary if the kernel cannot write committed state.**
   ``laws/crossbridge_kmc.propose_crossbridge_transitions_kernel`` reads ``bound`` and writes a
   separate ``proposal``, so a rejected step is handled by NOT copying rather than by rolling back.
   That is the cheapest possible form of "commits ONLY on an accepted physical step"
   (``CLAUDE.md`` §Architectural principle) and it needs no second copy of the state.
3. ⛔ **``BondFamily`` COULD NOT EXPRESS THIS FAMILY — raised, and now ANSWERED in ``bond.py``.**
   ``node_i``/``node_j`` are built once and never change, which is correct for a crosslink or a
   tether. For a crossbridge the actin partner is **STATE, NOT TOPOLOGY**: a head detaches, the
   filament slides, and the head rebinds to a DIFFERENT site. The incumbent has this defect already
   and names it — ``components/motor/INTEGRATION.md:142``: the crossbridge kernels "pulled the head
   toward a FIXED anchor and never READ the abscissa", so the runtime realised a passive spring
   network while the active contraction was never generated.
   This module did not work around it: it raised the question to Lead, who owns ``bond.py``, and
   ``f99400a6`` added ``partner_is_dynamic`` + ``live_partner_owner`` — paired, so a family cannot
   declare one without the other — and ``build_bond_family`` forwards them. **This family DECLARES
   them**, so ``node_j`` can never appear in an artifact as the partner in force at that step. A
   re-targeting search written here instead would have put a kinetic degree of freedom inside a
   connector module where nothing could see it.

WHAT IS ANSWERED — the count, and it is answerable for a structural reason
--------------------------------------------------------------------------
:class:`~aleph.world.bond.BondCount` asks *how many, against what, for which cell, on whose
authority*. For a KINETIC family that first question splits in two, and only one half is a population:

* **CAPACITY** — how many crossbridge SITES exist. One per head, ``2H`` per bipolar minifilament.
  Structural, static, and answerable. This is what :data:`SPEC` declares and what the family holds.
* **OCCUPANCY** — how many are bound at an instant. ``params_i0b3.yaml`` is explicit that the engaged
  fraction "is EMERGENT from k_on/k_off ... NOT imposed as this duty", so declaring it as a count
  would be the ERM defect in kinetic dress: a number that looks like a measurement and is actually a
  latch. It is never a ``BondCount`` here.

So the count is ``basis="per_filament"``, ``value = 2H``, resolved against the MEASURED minifilament
population from ``build/nmii.py`` — never a typed number. **PI decision 6 (2026-08-21) ratified
``n_heads_per_side = 30`` and ``F_stall_head = 2.0 pN``.**

⚠⚠ **THE BAND IS NOT RATIFIED — 30 IS A TEST POINT AND NOTHING IN THIS MODULE MAY SAY OTHERWISE.**
The PI's condition is recorded verbatim: *"band length must be rethought"*. The proposed [10, 30] was
NOT ratified, because its low end comes from a draft row whose own citations field declares them
absent from the corpus, and its high end from a composition review read out of scope. This module
therefore carries an axis with a TEST POINT and no band; :func:`crossbridge_count` refuses to be given
a band and :data:`SPEC` says so in its source string.

WHAT IS REFUSED — every magnitude, and one structure
-----------------------------------------------------
* ``k_xb`` — the crossbridge stiffness, ``params_i0b3.yaml``'s own "MASTER force knob", listed under
  ``still_gap_under_nm2b``. Required argument, no default.
* ``r0_xb`` — the crossbridge REST LENGTH, required for a reason that is not just sourcing. Two
  readings are defensible and they are not close: (a) the as-built head-to-site separation, making the
  family force-free at t0 — the ERM convention; (b) ``~0``, because "a bound head sits on the actin
  site" (``hand.NMIIParams``). Under (b) the built separation becomes strain, and that separation is
  set by the cortex DISCRETISATION — which is a mesh number acquiring a force, the same class as the
  ERM count. Both readings are named in the refusal and neither is chosen here.
* the four Pereverzev constants, ``k_on`` and the capture radius — all I0-B3 PI-GAPs, refused by
  :class:`~aleph.laws.crossbridge_kmc.CrossbridgeKinetics`, which has no defaults either.
* **per-strand cortex polarity** — not a missing number. PI decision 4 ratified it and
  ``StrandPopulation`` carries ``polarity_per_strand`` / :meth:`antiparallel` /
  ``polarity_is_mixed``, and ``build_strand_population`` gained ``polarity_mix`` / ``polarity_basis``
  at ``832e4313``. So this module CALLS that API and lets it fail, rather than routing around it — a
  cortex built WITHOUT ``polarity_mix`` is still uniform, and for that one every station is a spring
  wearing a motor's name, so it raises :class:`~aleph.world.families.ConnectorGapError` instead of
  placing bonds.

THE STRADDLE IS THE PHYSICS, AND BOTH DEGENERATE CASES ARE REJECTED
-------------------------------------------------------------------
A bipolar minifilament pulls two anti-parallel filaments TOWARD each other. Two arrangements look like
a motor station and are not:

* **both sides on the SAME filament** — the minifilament then pulls one filament against itself. That
  is an internal spring wearing a motor's name.
* **two filaments of the SAME polarity** — the two heads walk the same way, so the station TRANSLATES
  the pair instead of contracting it. ``build/nmii.py`` says the same thing from the geometry side:
  "a same-side placement would build a filament that translates rather than contracts".

Both are rejected by :func:`assert_bipolar_station`, and a minifilament with no legal partner within
reach is REPORTED, never dropped — ``STATE.md`` (e) 5 records the precedent: SF geometry refused a
motor station for a MODEL reason and the standing instruction is *do not force it*. Silently dropping
unstationed minifilaments would lower the motor count without lowering the declared density, which is
the ERM defect exactly.

⚠ THE DECLARED CHEMISTRY CARD NAMES TWO LAWS THAT ARE BOTH RETIRED FOR THIS BOND.
``engine/contracts.py:879-902`` gives all four NMII motor edges ``nmii_head_actin_hill_bell``.
**Hill** was superseded by the PI-ratified linear force-velocity on 2026-07-07, and PI decision 7 of
2026-08-21 re-scopes ``KB-DRAFT-7-05`` on the grounds that Hill 1938 is muscle. **Bell** was
superseded for this head by the 2026-07-23 fidelity correction to Kovacs 2007 catch-slip. A card is a
family's IDENTITY under this package's ``DUPLICATE_CRITERION``, so an identity naming two retired laws
is not cosmetic. :data:`CHEMISTRY_CARD` is the corrected one and the disagreement is reported rather
than silently reconciled — ``engine/`` is frozen port source and is not this session's to edit.

DUPLICATES: THIS IS THE CANONICAL SIDE OF FOUR, AND THE ARGUMENT IS EXHAUSTIVE
---------------------------------------------------------------------------
Recorded here because ``families/__init__.py`` requires it: *"The argument for each lives in the
owning module's docstring, which is where it must stay."* Lead merges :data:`DUPLICATE_OF_ENTRY`; the
reasoning is this module's to supply, and ``test_nmii_crossbridge.py`` re-derives every count below
from ``engine/contracts.py`` rather than trusting this prose.

**Counted, not sampled.** ``engine/contracts.py`` declares **five** ``ConnectorFamily.MOTOR`` edges.
**Four** of them carry ``nmii_head_actin_hill_bell`` — the card appears exactly four times in the
whole file and on no other family.

**The four are identical in every declared field except one.** Not "similar": the same
``ConnectorFamily.MOTOR``, the same source population ``"nmii"``, the same ``kinetics=True,
commit_on_accept=True``, the same ``endpoint_role_a="individual NMII head crossbridge"``, the same
``endpoint_role_b="live polar actin material coordinate"``, the same card and the same
``generation_required / remap_on_accept / blocks_sleep_refine``. **The only difference is the target
population** — ``sf_arc`` / ``cortex`` / ``lamellipodium`` / ``filopodium``. And under ``bond.py``
the component pair is not part of a family's identity at all: it is DERIVED from the arena's ID
ranges at query time, which is why a bond stores two global node indices and nothing else.

**THE CONTROL, and it is what makes this a card verdict rather than a name verdict.** The fifth
MOTOR edge is ``mt_cortex_capture``, and it carries a DIFFERENT card —
``microtubule_cortical_dynein``. So folding by chemistry does NOT fold everything in
``ConnectorFamily.MOTOR``; it separates a myosin crossbridge from a dynein capture site, which is
the same shape as D13's control (``filopodium_cortex_root`` sitting in ``TRANSIENT_ACTIN`` under a
different card) and the reason that verdict was trusted.

**The physics, stated so it can be disagreed with.** A myosin head does not know which actin
structure it caught. Binding is head-to-F-actin: the same crossbridge chemistry, the same catch-slip
detachment, the same stall force. What differs across the four is the ARCHITECTURE the actin is
arranged in — and architecture is what decides whether a station is legal (rectified bundle versus
mixed-polarity cortex, PI decision 4), not what the bond is made of. One family, four scopes.

⚠ **AND THE CONTRACT ITSELF ALREADY SAID WHAT ``bond.py`` COULD NOT EXPRESS.** All four declare
``endpoint_role_b="live polar actin material coordinate"`` — *live*, and a *coordinate* rather than a
node. That is the dynamic partner of item 3 above, written into the contract before this family was
built and invisible in ``BondFamily`` until ``f99400a6``. It is independent corroboration that item 3
is a property of the physics rather than of this module's implementation.

⚠ **WHAT A DUPLICATE VERDICT DOES NOT SETTLE.** ``DUPLICATE_CRITERION``'s own warning applies here:
folding answers "is this one family or two" and never "how many". Each of the three siblings still
needs its own straddle over ITS population, and ``sf_arc`` has no builder at all (PI decision 8 — it
is a transverse/dorsal arc, and what PHASE 1 built is a ventral FA-to-FA bundle, a different
structure). Merging them removes three chemistry questions, not three build questions.

WHAT CLOSING THE GAPS WOULD BUY — WRITTEN BEFORE THE VALUES EXIST
-----------------------------------------------------------------
PI queue item 11 will be read with the question *"what does this buy"*, and an answer written after
the constants arrive can be shaped by them. So it is written now, and :data:`MEASUREMENT_TIERS` is
the machine-readable form a queue can enumerate. ``CLAUDE.md``'s rule that a gate is a contract
written BEFORE the run is the same rule; this is that rule applied to a decision rather than a run.

⚠ **CORRECTION, AND IT IS THE FIRST THING WRITING THIS DOWN PRODUCED.** This module's own session
told Lead the runnable minimum was *five* — ``k_xb`` plus the four catch-slip constants. **That is
wrong, and so is six.** Tension needs ``r0_xb`` as well as ``k_xb``, which makes six; and
:class:`~aleph.laws.crossbridge_kmc.CrossbridgeKinetics` refuses without ``k_on`` and
``capture_radius_um``, which makes eight. **There is no partial set that runs.**

**AND THE TEMPTING INTERMEDIATE TIER IS NOT LEGITIMATE, WHICH IS THE REAL FINDING.** A detach-only
experiment — the six tension-and-detachment constants, skipping attachment — looks like a cheap way
to see the catch branch in a populated system. It is not available, and not because of an API: a
detachment measurement needs an INITIAL BOUND SET, :class:`CrossbridgeState` is born all-free on
purpose because a bound ``t0`` is an imposed duty ratio, and the only non-imposed source of a bound
set is attachment. So skipping ``k_on`` does not buy a cheaper experiment; it buys an experiment
whose initial condition is the very thing ``params_i0b3.yaml`` forbids imposing. **The tiers are
zero and eight.**

**TIER 0 — ZERO constants, and it is the measurement with the most at stake.** Nothing in the eight
is needed to ask whether the built cell admits a bipolar station AT ALL:
:func:`plan_crossbridge_stations` needs positions, per-strand polarity and a declared reach. It
answers *how many of the native minifilaments have a legal anti-parallel partner within reach, and
how many have none* — a MODEL question, refutable, and with no magnitude anywhere in it.
⚠ **The likely answer is a refusal**, and a refusal is the result: ``STATE.md`` (e) 5 records SF
geometry returning ``N_stations: 0`` for a model reason, with the standing instruction *do not force
it*. If the cortex's polarity assignment and the minifilament placement do not produce anti-parallel
partners, **no value of any of the eight constants makes this family contractile**, and knowing that
before the PI rules is worth more than the ruling. It needs a driver, which is ``scripts/`` and not
this lane's.

**TIER 1 — ALL EIGHT: the attach/detach cycle, and exactly three things it would buy.**

1. *The catch branch in a populated system, not in a closed form.* ``laws/crossbridge_kmc``'s sign
   claim is currently checked on an analytic curve. Under load from a real cortex the question is
   whether the load a bound head actually sees lands near ``F*`` at all, or two decades away — in
   which case the catch-slip law is doing nothing the retired Bell law would not have done, and the
   2026-07-23 correction, while right, would be unobservable at this operating point.
2. *``F*`` against the ratified ``F_stall_head`` = 2.0 pN.* The oracle is
   :attr:`CrossbridgeKinetics.peak_force_pn`, written before the constants exist, which is the only
   order in which it is a contract rather than a fit.
3. *The emergent bound fraction against the NM2B duty band 0.2-0.3 (Nagy 2013).* A **CROSS-CHECK,
   never a target** — ``params_i0b3.yaml`` is explicit that duty is "the unloaded consistency check
   phi_b(f=0) ~ duty, never a runtime input". A run that misses the band is evidence about the
   constants; a run tuned to hit it is nothing at all.

**WHAT NEITHER TIER BUYS — the half that decides what may be quoted.**

* **No cortical tension γ, and not by an oversight.** This family binds and releases; it does not
  STEP. No powerstroke, no ``active_step``, no abscissa advance — ``laws/crossbridge_kmc`` writes an
  off-rate and nothing that moves a head along actin. A bound crossbridge here is a spring with a
  lifetime. γ from NMII needs the stepping law, which is a separate PI-GAP set (``v0``,
  ``kappa_hill``) and a separate module.
* **No magnitude, at any tier, until acceptance is defined.** PI decision 18 leaves the predicate
  undefined and ``world/step.py`` returns ``ACCEPTANCE_UNDEFINED``, so kinetics that never ran inside
  an accepted step produce a sequence of states, not a trajectory. Tier 1 with all eight constants
  and no acceptance predicate buys the three items above and **still no physics magnitude**.
* ⚠ **A residual-based predicate would be LOOSENED by the events it must judge.** PHASE 4 measured
  the balance gate's tolerance growing 86x when 8.2M force terms were added while the residual did
  not move. Every attach adds a force term, so in this family the term count is not constant between
  steps: more binding raises the tolerance, and the binding is the thing under test. A tolerance for
  a kinetic family has to be set against fixed CAPACITY, never against in-step term count — occupancy
  is emergent and must not become a gate parameter.
* **Nothing transfers to the three sibling populations.** The duplicate verdict above folds four
  chemistry questions into one; it folds zero build questions. ``sf_arc`` has no builder at all
  (PI decision 8).
* **Every number would stand on a DECLARED axis, not a sourced one.** ``hand.NMIICatchSlipParams``
  says of these values *"NEVER chosen to lift a force/gamma band"*. A PI ruling makes them
  declarable, not measured, so any Tier-1 result carries the same expiry ERM's density axis carries
  — the PI's own word for that class is TEMPORARY.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — the count is ``per_filament`` [1] resolved against a minifilament COUNT [filaments];
    ``reach_um``/``r0_xb_um`` [µm]; ``k_xb`` [pN/µm]. ``BondCount`` rejects a support carrying any
    other unit by name, so an areal density passed here fails rather than resolving against a count.
  * boundary — a count that does not resolve to exactly one site per head raises with BOTH numbers
    (heads with no site and sites with no head are both errors, not roundings); a minifilament with no
    anti-parallel partner inside ``reach_um`` raises with the tally; a head whose assigned filament
    has no node inside ``reach_um`` raises rather than being dropped.
  * conservation/invariant — every endpoint is a GLOBAL arena node id inside the live NODE prefix
    (asserted by ``build_bond_family``); each head appears exactly once, so no head carries two
    molecular states; the two sides of one minifilament land on DISJOINT filaments, asserted per
    station. Bonds are never welds: the head and the actin node stay in their own populations.
  * CFL/precision — no integration here; float64 throughout. ``k_xb`` enters the CFL bound through the
    row sum at its endpoints, which is why ``bond.py`` stores stiffness per bond.
  * sign sense — the station contracts rather than translates, which is exactly the anti-parallel
    requirement above and is the one thing this module asserts about force. **No magnitude is claimed
    anywhere**: existence, sign, and that a transition commits only on an accepted step.
  * measurement protocol — host-side construction from downloaded positions; no device is touched and
    nothing is uploaded. Positions are passed in rather than read from the arena, so this module is
    CPU-importable and testable without a card. The KMC kernel it hands to Lead is codegen-checked,
    never executed here.

engine units: length µm, stiffness pN/µm, force pN.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.laws.crossbridge_kmc import NMII_CROSSBRIDGE_GAPS, CrossbridgeKinetics
from aleph.world.bond import BondCount, BondFamily, SourceClass, build_bond_family
from aleph.world.build.nmii import _anchor_beads
from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = [
    "CONNECTOR", "CHEMISTRY_CARD", "DECLARED_CHEMISTRY_CARD", "DUPLICATE_OF_ENTRY",
    "F_STALL_HEAD_PN", "N_HEADS_PER_SIDE_TEST_POINT", "BLOCKED_BY", "SPEC",
    "MEASUREMENT_TIERS", "NOT_BOUGHT_BY_ANY_TIER",
    "CrossbridgeState", "assert_bipolar_station", "crossbridge_count",
    "station_census", "plan_crossbridge_stations", "build_nmii_cortex_crossbridge",
]

#: The declared census name (``engine/contracts.py:883``). The file is named for the physics.
CONNECTOR = "nmii_cortex_motor"

#: The corrected card. See the module docstring: the declared one names Hill (retired 2026-07-07,
#: re-scoped by PI decision 7 as muscle-only) and Bell (retired for this head 2026-07-23).
CHEMISTRY_CARD = "nmii_head_actin_catchslip_kovacs2007"

#: What ``engine/contracts.py`` still declares for all four NMII motor edges. Reported, not edited —
#: ``engine/`` is frozen PORT SOURCE (``STATE.md`` (a)) and a card is an identity, so changing it is a
#: contract change and belongs to the PI.
DECLARED_CHEMISTRY_CARD = "nmii_head_actin_hill_bell"

#: The three sibling edges that are the SAME CARD and therefore the same family. A myosin head does
#: not know which actin structure it caught, and the component pair is derived from ID ranges. Lead
#: merges this into ``families.DUPLICATE_OF``; this module does not write that dict.
DUPLICATE_OF_ENTRY: dict[str, str] = {
    "nmii_sf_motor": CONNECTOR,
    "nmii_lamellipodium_motor": CONNECTOR,
    "nmii_filopodium_motor": CONNECTOR,
}

#: Per-head isometric stall force [pN]. PI decision 6, 2026-08-21; ``laws/myosin_linear.F_HEAD_PN``
#: already carries 2.0 (KB-3.18). Recorded so the ratified scale travels with the family — this module
#: computes NO force and claims no magnitude; it is here to be compared against
#: :attr:`CrossbridgeKinetics.peak_force_pn` once the PI closes the catch-slip constants.
F_STALL_HEAD_PN = 2.0

#: Heads per anti-parallel side. **A TEST POINT, NOT A BAND** — PI decision 6, 2026-08-21, with the
#: PI's condition recorded: *"band length must be rethought"*. The proposed [10, 30] was NOT ratified.
N_HEADS_PER_SIDE_TEST_POINT = 30

#: The magnitudes this family cannot source. Every one has no default anywhere in the path.
BLOCKED_BY: tuple[str, ...] = (
    "k_xb_pn_per_um",
    "r0_xb_um",
    *NMII_CROSSBRIDGE_GAPS,
)

#: What closing the gaps would buy, enumerable by Lead's queue. See the module docstring for the
#: argument; this is the same content in the form a decision queue can read.
#:
#: ⚠ There is no tier between these two, and the tempting one — the six tension-and-detachment
#: constants without ``k_on`` — is not merely unimplemented. A detachment measurement needs an
#: initial bound set; ``CrossbridgeState`` is all-free by construction because a bound t0 is an
#: imposed duty ratio; and attachment is the only non-imposed source of one.
MEASUREMENT_TIERS: dict[str, dict[str, object]] = {
    "tier_0_station_census": {
        "requires": (),
        "buys": "whether the built cell admits a legal bipolar station at all, and how many "
                "minifilaments have no anti-parallel partner within reach. A MODEL question, "
                "refutable, no magnitude. If the answer is zero, no value of any of the eight "
                "constants makes this family contractile.",
        "blocked_on": "a driver, which is scripts/ and not this lane's",
    },
    "tier_1_attach_detach_cycle": {
        "requires": BLOCKED_BY,
        "buys": ("the catch branch in a populated system rather than in a closed form; F* against "
                 "the ratified F_stall_head = 2.0 pN, using an oracle written before the constants "
                 "exist; the emergent bound fraction as a CROSS-CHECK against the NM2B duty band "
                 "0.2-0.3 (Nagy 2013), never as a target."),
        "blocked_on": "PI queue item 11 — all eight constants",
    },
}

#: What NO tier buys, so a reader cannot infer it from what they do. See the module docstring.
NOT_BOUGHT_BY_ANY_TIER: tuple[str, ...] = (
    "cortical tension gamma — this family binds and releases but does not STEP; there is no "
    "powerstroke and no abscissa advance here, so a bound crossbridge is a spring with a lifetime",
    "any physics magnitude — PI decision 18 leaves acceptance undefined, so kinetics that never ran "
    "inside an accepted step produce a sequence of states, not a trajectory",
    "a residual-based acceptance predicate — PHASE 4 measured the tolerance growing 86x with 8.2M "
    "added force terms, and every attach adds a term, so the gate would be loosened by the events "
    "under test; a kinetic tolerance must be set against fixed CAPACITY, never in-step term count",
    "anything about the three sibling populations — the duplicate verdict folds four CHEMISTRY "
    "questions and zero BUILD questions; sf_arc has no builder (PI decision 8)",
    "a SOURCED number — a PI ruling makes these declarable, not measured, so any result carries the "
    "same expiry the ERM density axis carries; the PI's own word for that class is TEMPORARY",
)

_GAP_NOTE_KINETICS = (
    "The head-actin bond is CATCH-SLIP (Kovacs 2007; components/motor/hand.py fidelity correction "
    "2026-07-23), NOT Bell slip — hand.step_detach_kernel stays pure Bell for its ERM/alpha-actinin "
    "consumers and is not this bond's law, and hand_kmc's own NMIIA_MYOSIN preset is catch_slip=False "
    "i.e. the retired law. The FORM and the SIGN are sourced and are written in "
    "laws/crossbridge_kmc.py; every MAGNITUDE is an I0-B3 PI-GAP and hand.NMIICatchSlipParams says of "
    "them 'NEVER chosen to lift a force/gamma band'. k_xb is params_i0b3's own MASTER force knob. "
    "r0_xb is refused for a second reason: as-built (force-free, the ERM convention) and ~0 ('a bound "
    "head sits on the actin site') differ by the whole built separation, and that separation is set "
    "by the cortex discretisation — a mesh number acquiring a force."
)

_GAP_NOTE_POLARITY = (
    "PI decision 4 (2026-08-21) ratified PER-STRAND polarity and StrandPopulation already exposes "
    "polarity_per_strand / antiparallel() / polarity_is_mixed, and build_strand_population gained "
    "polarity_mix / polarity_basis at 832e4313. A cortex built WITHOUT polarity_mix is still uniform "
    "and answers with one scalar for the whole population; in such a population an anti-parallel pair "
    "is not merely unverifiable, it is UNREPRESENTABLE, so every station this module could place "
    "would be a spring wearing a motor's name. This is a BUILD CONFIGURATION refusal, not an open "
    "question, which is why it is a runtime refusal and not a SPEC blocker."
)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("cortex", "nmii"),
    basis="per_filament",
    blocked_by=BLOCKED_BY,
    source=(
        "COUNT answered, EVERY MAGNITUDE absent. Capacity = one crossbridge site per head = 2H per "
        f"bipolar minifilament, resolved against the MEASURED minifilament population from "
        f"build/nmii.py. H = {N_HEADS_PER_SIDE_TEST_POINT} and F_stall_head = {F_STALL_HEAD_PN} pN are "
        "PI decision 6, 2026-08-21 — ratified as a declared axis with a TEST POINT and ⚠ NO BAND: the "
        "proposed [10, 30] was explicitly NOT ratified ('band length must be rethought'), because its "
        "low end is a draft row whose own citations field declares them absent from the corpus and its "
        "high end is a composition review read out of scope. OCCUPANCY is never a count here — the "
        "engaged fraction is EMERGENT from k_on/k_off (params_i0b3.yaml) and declaring it would be a "
        "latch wearing a measurement's label."
    ),
)


# ── the count ───────────────────────────────────────────────────────────────────────────────────
def crossbridge_count(
    n_heads_per_side: int,
    *,
    scope: str,
    source_class: SourceClass = SourceClass.PI_GAP,
    provenance: str | None = None,
) -> BondCount:
    """The CAPACITY count: one crossbridge site per head, ``2H`` per bipolar minifilament.

    Args:
        n_heads_per_side: ``H``. A single value — passing a range is refused, because a band for this
            axis was explicitly NOT ratified (PI decision 6) and a builder that accepts one would be
            asserting admissibility the PI declined to grant.
        scope: the cell line, state, assay and temperature. Required by :class:`BondCount`.
        source_class: defaults to ``PI_GAP``; pass ``SOURCED`` only against a real areal datum.
        provenance: overrides the default citation. The default cites PI decision 6 AND its condition.

    Returns:
        A ``per_filament`` :class:`BondCount` to resolve against a minifilament COUNT.

    Raises:
        TypeError: if given a sequence — see ``n_heads_per_side``.
        ValueError: if ``H < 1``.
    """
    if isinstance(n_heads_per_side, (list, tuple, np.ndarray)):
        raise TypeError(
            "n_heads_per_side takes one value, not a band. PI decision 6 (2026-08-21) ratified the "
            "AXIS and the TEST POINT 30; the band [10, 30] was NOT ratified — the PI's condition is "
            "'band length must be rethought'. A band is a claim about which values are admissible and "
            "neither endpoint supports one."
        )
    h = int(n_heads_per_side)
    if h < 1:
        raise ValueError(f"n_heads_per_side must be >= 1; a side with no heads has no crossbridge. Got {h}")
    return BondCount(
        basis="per_filament", value=float(2 * h), scope=scope, source_class=source_class,
        provenance=provenance or (
            f"one crossbridge SITE per head; 2H = {2 * h} sites per bipolar minifilament. H = {h} is "
            "PI decision 6, 2026-08-21 — a declared axis with a TEST POINT of 30 and ⚠ NO BAND "
            "('band length must be rethought'). This is CAPACITY; occupancy is emergent from "
            "k_on/k_off and is never declared as a count."
        ),
    )


# ── the straddle ────────────────────────────────────────────────────────────────────────────────
def assert_bipolar_station(cortex, strand_plus: int, strand_minus: int) -> None:
    """Refuse the two arrangements that look like a motor station and are not.

    Args:
        cortex: the :class:`~aleph.world.build.cortex.StrandPopulation` the actin comes from.
        strand_plus: filament index the ``+`` head row would bind.
        strand_minus: filament index the ``-`` head row would bind.

    Raises:
        ValueError: if the two are the same filament, or run the same way.
        ConnectorGapError: if the population cannot answer, i.e. it carries no per-strand polarity.
    """
    if int(strand_plus) == int(strand_minus):
        raise ValueError(
            f"{CONNECTOR}: both head rows land on filament {strand_plus}. A bipolar minifilament "
            "straddling ONE filament pulls it against itself — an internal spring wearing a motor's "
            "name, not a contractile station."
        )
    try:
        opposed = cortex.antiparallel(int(strand_plus), int(strand_minus))
    except ValueError as exc:  # the population carries one scalar polarity for all of it
        raise ConnectorGapError(
            CONNECTOR, ["cortex.polarity_per_strand"], f"{exc} {_GAP_NOTE_POLARITY}") from exc
    if not opposed:
        raise ValueError(
            f"{CONNECTOR}: filaments {strand_plus} and {strand_minus} run the SAME way. The two head "
            "rows would walk the same direction, so the station TRANSLATES the pair instead of "
            "contracting it — build/nmii.py names the same failure from the geometry side. A bipolar "
            "minifilament straddles two ANTI-PARALLEL filaments or it is not a motor."
        )


def _strand_of(node_ids: npt.NDArray[np.int64], cortex) -> npt.NDArray[np.int64]:
    """Filament index of each GLOBAL cortex node id — arithmetic over one contiguous range."""
    return (node_ids - cortex.nodes.lo) // cortex.nodes_per_strand


def _select_stations(*, n_mf, n_per, mf, cor, tree, strand_of_all, polarity, k_plus, k_minus,
                     cortex, reach_um) -> dict:
    """Choose one anti-parallel filament pair per minifilament. The shared core of both entry points.

    Split out of :func:`plan_crossbridge_stations` on 2026-08-21 because the two callers want
    OPPOSITE things from the same computation: a builder must REFUSE if any minifilament has no legal
    station, and a census must REPORT how many do. Sharing the selection rather than the exception is
    what keeps those two answers about the same arrangement.
    """
    station = np.full((n_mf, 2), -1, np.int64)
    unstationed: list[int] = []
    n_no_cortex = 0
    n_no_partner = 0
    n_first_choice = 0

    for m in range(n_mf):
        base = m * n_per
        c_plus = mf[base + k_plus].mean(axis=0)
        c_minus = mf[base + k_minus].mean(axis=0)

        near_plus = np.asarray(tree.query_ball_point(c_plus, reach_um), np.int64)
        near_minus = np.asarray(tree.query_ball_point(c_minus, reach_um), np.int64)
        if near_plus.size == 0 or near_minus.size == 0:
            unstationed.append(m)
            n_no_cortex += 1
            continue

        # The + row takes its nearest filament; the - row must then take an ANTI-PARALLEL one. The
        # search is over the - row's own reach, so a station is never manufactured by reaching further
        # than the geometry allows.
        s_plus = int(strand_of_all[near_plus[np.argmin(np.linalg.norm(cor[near_plus] - c_plus, axis=1))]])
        s_minus_all = strand_of_all[near_minus]
        legal = (s_minus_all != s_plus) & (polarity[s_minus_all] != polarity[s_plus])
        if not legal.any():
            unstationed.append(m)
            n_no_partner += 1
            continue
        cand = near_minus[legal]
        best = cand[np.argmin(np.linalg.norm(cor[cand] - c_minus, axis=1))]
        s_minus = int(strand_of_all[best])
        n_first_choice += int(s_minus == int(strand_of_all[near_minus[
            np.argmin(np.linalg.norm(cor[near_minus] - c_minus, axis=1))]]))

        assert_bipolar_station(cortex, s_plus, s_minus)
        station[m] = (s_plus, s_minus)

    return {
        "station": station, "unstationed": unstationed, "_unstationed_set": set(unstationed),
        "n_no_cortex_in_reach": n_no_cortex, "n_no_antiparallel_partner": n_no_partner,
        "n_minus_row_first_choice": n_first_choice,
    }


def station_census(
    *,
    nmii_inventory: dict,
    nmii_pos_um: npt.ArrayLike,
    cortex,
    cortex_pos_um: npt.ArrayLike,
    reach_um: float,
) -> dict[str, object]:
    """TIER 0 — does the built cell admit a legal bipolar station at all, and how many?

    **Needs NONE of the eight PI-GAPs**: positions, per-strand polarity and a declared reach. No
    stiffness, no rate, no force, so nothing it returns is a magnitude and nothing in it can be
    tuned. See :data:`MEASUREMENT_TIERS`; a test asserts this signature stays disjoint from
    :data:`BLOCKED_BY` so the tier cannot quietly acquire a constant.

    ⚠ **THE ANSWER MAY BE ZERO, AND ZERO IS THE RESULT.** ``STATE.md`` (e) 5 records the precedent —
    SF geometry returned ``N_stations: 0`` for a MODEL reason, with the standing instruction *do not
    force it*. If the built cell admits no legal station, **no value of any of the eight constants
    makes this family contractile**, and that is worth more than the ruling it precedes.

    ⚠ **REACH IS A DECLARED AXIS AND A SINGLE POINT IS NOT A RESULT.** The count depends on
    ``reach_um`` monotonically, so one value produces a number whose axis is invisible — and picking
    the value that makes the count acceptable is exactly what ``CLAUDE.md`` forbids ("nothing is
    chosen to make a gate pass"). **SWEEP it and report the curve.** The physically meaningful scale
    to bracket is the built head offset (``nmii_inventory["head_offset_um"]``) and the cortical shell
    thickness; the KINETIC capture radius is a different quantity and is a PI-GAP.

    ⚠ **THE TWO FAILURE REASONS ARE DIFFERENT PROBLEMS AND ARE COUNTED SEPARATELY.**
    ``n_no_cortex_in_reach`` says the minifilament sits too far from any actin — a PLACEMENT
    problem, owned by ``build/nmii.py``. ``n_no_antiparallel_partner`` says it found actin and none
    of it opposed — a POLARITY problem, owned by ``build/cortex.py``'s ``polarity_mix`` (PI decision
    4). Collapsing them into one "unstationed" count would send the fix to the wrong builder.

    Args:
        nmii_inventory / nmii_pos_um / cortex / cortex_pos_um / reach_um: see
            :func:`plan_crossbridge_stations`. ``cortex`` MUST carry per-strand polarity — build it
            with ``polarity_mix=`` and ``polarity_basis=``, or this refuses.

    Returns:
        The census. ``station`` is ``(n_mf, 2)`` with ``-1`` in both columns for an unstationed
        minifilament, so a reader cannot mistake a missing station for filament 0.

    Raises:
        ConnectorGapError: if the cortex carries no per-strand polarity, or is uniform.
        ValueError: on a malformed array or a non-positive reach.

    Example — the whole TIER 0 driver body, for a caller in ``scripts/``::

        from aleph.world.families.nmii_cortex_crossbridge import station_census

        rows = []
        for reach_um in REACH_SWEEP_UM:          # declared before the run, never narrowed after
            rows.append(station_census(
                nmii_inventory=nmii,             # the dict build_nmii() returned
                nmii_pos_um=pos[nmii_lo:nmii_hi],     # host copy, BETWEEN steps
                cortex=cortex_pop,               # built with polarity_mix=/polarity_basis=
                cortex_pos_um=pos[cortex_lo:cortex_hi],
                reach_um=reach_um))
    """
    from scipy.spatial import cKDTree

    if not np.isfinite(reach_um) or reach_um <= 0.0:
        raise ValueError(f"reach_um must be finite and positive; got {reach_um!r}")
    if getattr(cortex, "polarity_per_strand", None) is None:
        raise ConnectorGapError(CONNECTOR, ["cortex.polarity_per_strand"], _GAP_NOTE_POLARITY)
    if not cortex.polarity_is_mixed:
        raise ConnectorGapError(
            CONNECTOR, ["cortex.polarity_per_strand"],
            "the population carries a per-strand array but every filament runs the SAME way, so it "
            f"holds no anti-parallel pair at all. {_GAP_NOTE_POLARITY}")

    mf = np.asarray(nmii_pos_um, np.float64)
    cor = np.asarray(cortex_pos_um, np.float64)
    if mf.ndim != 2 or mf.shape[1] != 3 or cor.ndim != 2 or cor.shape[1] != 3:
        raise ValueError("nmii_pos_um and cortex_pos_um must both be (N, 3)")

    n_mf = int(nmii_inventory["n_minifilaments"])
    n_bb = int(nmii_inventory["n_bb"])
    h = int(nmii_inventory["n_heads_per_side"])
    n_per = int(nmii_inventory["nodes_per_minifilament"])
    if mf.shape[0] != n_mf * n_per:
        raise ValueError(
            f"{CONNECTOR}: nmii_pos_um has {mf.shape[0]} rows for {n_mf} x {n_per} = {n_mf * n_per} "
            "NMII nodes. The positions must be the whole block in arena order.")

    polarity = np.asarray(cortex.polarity_per_strand, np.int8)
    tree = cKDTree(cor)

    # ⚠ HOW FAR IS THE ACTIN, ALWAYS — not only when the answer is "further than this reach".
    #
    # Added 2026-08-21 AFTER the first sweep, and the reason is recorded rather than smoothed over:
    # that sweep reported `n_no_cortex_in_reach` without ever reporting HOW far, so the distance had
    # to be recovered by sweeping nine reaches and reading where the count moved. One query answers
    # it directly. This changes no recorded number and re-scopes no question — it makes the next run
    # self-diagnosing instead of requiring a ladder to locate a distance the instrument already had
    # in hand.
    _c = np.stack([mf[np.arange(n_mf)[:, None] * n_per + (n_bb + np.arange(h))[None, :]].mean(axis=1),
                   mf[np.arange(n_mf)[:, None] * n_per
                      + (n_bb + h + np.arange(h))[None, :]].mean(axis=1)], axis=1)
    _d = tree.query(_c.reshape(-1, 3), k=1)[0].reshape(n_mf, 2)
    nearest_actin = {
        "note": "distance from each head-row centroid to the NEAREST cortex NODE [um], over ALL "
                "minifilaments and independent of reach_um. A floor far above zero means the two "
                "populations do not interpenetrate, which no reach can fix and no rate constant can "
                "reach across.",
        "min": float(_d.min()), "median": float(np.median(_d)), "max": float(_d.max()),
        "p05": float(np.percentile(_d, 5)), "p95": float(np.percentile(_d, 95)),
    }

    sel = _select_stations(
        n_mf=n_mf, n_per=n_per, mf=mf, cor=cor, tree=tree,
        strand_of_all=_strand_of(cortex.nodes.lo + np.arange(cor.shape[0], dtype=np.int64), cortex),
        polarity=polarity, k_plus=n_bb + np.arange(h), k_minus=n_bb + h + np.arange(h),
        cortex=cortex, reach_um=float(reach_um))

    n_unstationed = len(sel["unstationed"])
    return {
        "kind": "nmii_crossbridge_station_census",
        "connector": CONNECTOR,
        "reach_um": float(reach_um),
        "n_minifilaments": n_mf,
        "n_stationed": n_mf - n_unstationed,
        "n_unstationed": n_unstationed,
        # ⚠ Two different problems, two different owners. Never sum these into one number.
        "n_no_cortex_in_reach": sel["n_no_cortex_in_reach"],
        "n_no_antiparallel_partner": sel["n_no_antiparallel_partner"],
        # How often the anti-parallel requirement did NOT cost anything — the - row's nearest
        # filament was already legal. A low value says the constraint is doing real work; a value of
        # n_stationed says polarity is effectively unconstrained at this reach.
        "n_minus_row_first_choice": sel["n_minus_row_first_choice"],
        "station": sel["station"],
        "unstationed": list(sel["unstationed"]),
        "n_heads_per_side": h,
        "crossbridge_capacity_if_all_stationed": n_mf * 2 * h,
        "cortex_polarity": {
            "n_plus": int((polarity > 0).sum()), "n_minus": int((polarity < 0).sum()),
            "basis": getattr(cortex, "polarity_basis", None), "n_strands": int(cortex.n_strands),
        },
        "nearest_cortex_node_um": nearest_actin,
        "no_pi_gap_constant_used": True,
        "magnitude_claimed": None,
    }


def plan_crossbridge_stations(
    *,
    nmii_inventory: dict,
    nmii_pos_um: npt.ArrayLike,
    cortex,
    cortex_pos_um: npt.ArrayLike,
    reach_um: float,
) -> dict[str, object]:
    """Assign each minifilament a straddle over two ANTI-PARALLEL filaments, and each head a site.

    Existence and sign only — no force is computed and no magnitude is produced.

    Args:
        nmii_inventory: the census from :func:`aleph.world.build.nmii.build_nmii` (or ``plan_nmii``
            plus a ``claims`` entry). Read, never recomputed: the head index arithmetic here mirrors
            the kernel's through ``build.nmii._anchor_beads`` rather than restating it.
        nmii_pos_um: ``(n_nodes, 3)`` host positions of the whole NMII node block, in arena order.
        cortex: the :class:`~aleph.world.build.cortex.StrandPopulation` supplying the actin.
        cortex_pos_um: ``(N, 3)`` host cortex node positions [µm], in arena order.
        reach_um: the GEOMETRIC pairing reach [µm] — what could ever be paired at build time.
            ⚠ **Not** the kinetic capture radius ``epsilon``, which decides what binds this tick and is
            a PI-GAP on :class:`~aleph.laws.crossbridge_kmc.CrossbridgeKinetics`. Conflating them
            would let a build detail set a rate.

    Returns:
        ``pairs`` ``(n_heads, 2)`` global ids ordered ``(head, actin)``, ``separation_um``, the
        per-minifilament ``station`` assignment, and the straddle census.

    Raises:
        ConnectorGapError: if the cortex carries no per-strand polarity, or is uniform.
        ValueError: on a malformed array, a non-positive reach, a minifilament with no anti-parallel
            partner in reach, or a head with no site on its assigned filament in reach.
    """
    from scipy.spatial import cKDTree

    if not np.isfinite(reach_um) or reach_um <= 0.0:
        raise ValueError(f"reach_um must be finite and positive; got {reach_um!r}")
    if getattr(cortex, "polarity_per_strand", None) is None:
        raise ConnectorGapError(CONNECTOR, ["cortex.polarity_per_strand"], _GAP_NOTE_POLARITY)
    if not cortex.polarity_is_mixed:
        raise ConnectorGapError(
            CONNECTOR, ["cortex.polarity_per_strand"],
            "the population carries a per-strand array but every filament runs the SAME way, so it "
            f"holds no anti-parallel pair at all. {_GAP_NOTE_POLARITY}")

    mf = np.asarray(nmii_pos_um, np.float64)
    cor = np.asarray(cortex_pos_um, np.float64)
    if mf.ndim != 2 or mf.shape[1] != 3 or cor.ndim != 2 or cor.shape[1] != 3:
        raise ValueError("nmii_pos_um and cortex_pos_um must both be (N, 3)")

    n_mf = int(nmii_inventory["n_minifilaments"])
    n_bb = int(nmii_inventory["n_bb"])
    h = int(nmii_inventory["n_heads_per_side"])
    n_per = int(nmii_inventory["nodes_per_minifilament"])
    nmii_lo = int(nmii_inventory["claims"]["node"][0])
    if mf.shape[0] != n_mf * n_per:
        raise ValueError(
            f"{CONNECTOR}: nmii_pos_um has {mf.shape[0]} rows for {n_mf} x {n_per} = {n_mf * n_per} "
            "NMII nodes. The positions must be the whole block in arena order.")

    # The head rows, from the builder's own arithmetic. `_anchor_beads` is the host mirror of the
    # kernel's anchor split and is used rather than restated, so the two cannot drift apart.
    anchors_plus, anchors_minus = _anchor_beads(n_bb, h)
    if max(anchors_plus) >= min(anchors_minus):
        raise ValueError(
            f"{CONNECTOR}: the builder's anchor split does not separate the two halves of the "
            f"backbone (+ up to bead {max(anchors_plus)}, - from bead {min(anchors_minus)}). Without "
            "an axial split the two head rows are not on opposite ends and the filament is not "
            "bipolar.")
    k_plus = n_bb + np.arange(h)
    k_minus = n_bb + h + np.arange(h)

    tree = cKDTree(cor)
    strand_of_all = _strand_of(cortex.nodes.lo + np.arange(cor.shape[0], dtype=np.int64), cortex)
    polarity = np.asarray(cortex.polarity_per_strand, np.int8)

    census = _select_stations(
        n_mf=n_mf, n_per=n_per, mf=mf, cor=cor, tree=tree, strand_of_all=strand_of_all,
        polarity=polarity, k_plus=k_plus, k_minus=k_minus, cortex=cortex, reach_um=float(reach_um))
    station = census["station"]
    unstationed = census["unstationed"]

    pairs = np.empty((n_mf * 2 * h, 2), np.int64)
    sep = np.empty(n_mf * 2 * h, np.float64)

    for m in range(n_mf):
        if m in census["_unstationed_set"]:
            continue
        base = m * n_per
        s_plus, s_minus = int(station[m, 0]), int(station[m, 1])
        for row, (ks, s) in enumerate(((k_plus, s_plus), (k_minus, s_minus))):
            sites = cortex.nodes.lo + s * cortex.nodes_per_strand + np.arange(
                cortex.nodes_per_strand, dtype=np.int64)
            d = np.linalg.norm(mf[base + ks][:, None, :] - cor[sites - cortex.nodes.lo][None, :, :],
                               axis=2)
            pick = np.argmin(d, axis=1)
            dist = d[np.arange(h), pick]
            if np.any(dist > reach_um):
                raise ValueError(
                    f"{CONNECTOR}: minifilament {m}, {'+-'[row]} row — {int((dist > reach_um).sum())} "
                    f"of {h} heads have no site on filament {s} within the {reach_um:g} um reach (max "
                    f"{dist.max():.4g} um). Refine the cortex resolution or re-place the minifilament "
                    "rather than dropping heads: a dropped head lowers the motor count while the "
                    "declared density stays put, which is the ERM defect exactly.")
            slot = (m * 2 + row) * h + np.arange(h)
            pairs[slot, 0] = nmii_lo + base + ks
            pairs[slot, 1] = sites[pick]
            sep[slot] = dist

    if unstationed:
        raise ValueError(
            f"{CONNECTOR}: {len(unstationed)} of {n_mf} minifilaments have NO legal bipolar station "
            f"within the {reach_um:g} um reach (first: {unstationed[:5]}; "
            f"{census['n_no_cortex_in_reach']} found no cortex node at all, "
            f"{census['n_no_antiparallel_partner']} found cortex but no ANTI-PARALLEL partner). "
            "Reported rather than dropped — STATE.md (e) 5 records the precedent, an SF geometry "
            "that refused a motor station for a MODEL reason with the standing instruction 'do not "
            "force it'. Either the polarity assignment or the placement is what must change, and "
            "both are builder decisions, not this connector's. Call station_census() for the counts "
            "without the exception."
        )

    return {
        "pairs": pairs, "separation_um": sep, "station": station,
        "n_stations": n_mf, "n_heads": n_mf * 2 * h, "n_heads_per_side": h,
        "reach_um": float(reach_um),
        "n_minus_row_first_choice": census["n_minus_row_first_choice"],
        "anchor_beads": {"plus": anchors_plus, "minus": anchors_minus},
    }


# ── the kinetic state ───────────────────────────────────────────────────────────────────────────
@dataclass(slots=True)
class CrossbridgeState:
    """The bound/free state of one crossbridge family, and the only path by which it changes.

    ⚠ **There is no free list and there is no snapshot twin, and both absences are measured rather
    than deferred.** A head owns exactly one crossbridge slot for the whole run, so detaching frees no
    storage; and :func:`~aleph.laws.crossbridge_kmc.propose_crossbridge_transitions_kernel` writes a
    SEPARATE proposal array, so a rejected step is handled by not copying and there is nothing to roll
    back. See the module docstring — these are this family's answers to what ``bond.py`` deferred.

    Attributes:
        bound: ``(B,)`` int8, 1 if the crossbridge is attached. Built ALL-FREE: a bound state at t0
            would be an imposed duty ratio, and the engaged fraction is emergent.
        n_commits: how many accepted steps this state has taken. Recorded so an artifact can show
            whether the kinetics ever ran.
    """

    bound: npt.NDArray[np.int8] = field(repr=False)
    n_commits: int = 0

    @classmethod
    def all_free(cls, n_bonds: int) -> "CrossbridgeState":
        """An all-free population. The only legal initial state — see :attr:`bound`."""
        return cls(bound=np.zeros(int(n_bonds), np.int8))

    @property
    def bound_fraction(self) -> float:
        """Fraction currently attached. MEASURED from the state, never imposed as a duty ratio."""
        return float(self.bound.mean()) if self.bound.size else 0.0

    def commit(self, proposal: npt.ArrayLike, *, accepted: bool) -> int:
        """Apply a KMC proposal, and ONLY on an accepted physical step.

        ``CLAUDE.md`` §Architectural principle: *"a kinetic connector commits ONLY on an accepted
        physical step under one device-resident transaction."*

        Args:
            proposal: ``(B,)`` proposed next bound state, from the KMC kernel.
            accepted: whether the physical step that produced it was ACCEPTED. There is no default:
                the caller must state it, because the whole failure mode is a commit that happened
                because nobody said not to.

        Returns:
            The number of bonds whose state changed.

        Raises:
            ValueError: on a length mismatch or a proposal that is not 0/1.
            RuntimeError: if ``accepted`` is False — the proposal is DISCARDED and the caller is told
                so by name rather than by a silent no-op.
        """
        p = np.asarray(proposal, np.int8).reshape(-1)
        if p.size != self.bound.size:
            raise ValueError(f"proposal has {p.size} entries for {self.bound.size} crossbridges")
        if p.size and not np.isin(p, (0, 1)).all():
            raise ValueError("a crossbridge is bound or free; a proposal must be 0 or 1")
        if not accepted:
            raise RuntimeError(
                f"{CONNECTOR}: refusing to commit a kinetic transition on a step that was not "
                "accepted. The proposal is discarded. ⚠ PI decision 18 (2026-08-21) leaves the "
                "acceptance predicate UNDEFINED — balance_ok is a Higham rounding bound and cannot "
                "fail for physics, descent_ratio reads a post-rollback value and cannot succeed — so "
                "the caller must pass what its own step machinery decided and this family will not "
                "guess it. ⚠ CONSEQUENCE, STATED RATHER THAN WORKED AROUND: world/step.py's default "
                "predicate is UndefinedAcceptance, whose verdict is ACCEPTANCE_UNDEFINED and whose "
                "is_accepted is False on purpose — so on today's step path this family commits "
                "NOTHING, and its kinetics have never run inside an accepted step."
            )
        changed = int(np.count_nonzero(p != self.bound))
        self.bound[:] = p
        self.n_commits += 1
        return changed


# ── the family ──────────────────────────────────────────────────────────────────────────────────
def build_nmii_cortex_crossbridge(
    arena,
    *,
    nmii_inventory: dict,
    nmii_pos_um: npt.ArrayLike,
    cortex,
    cortex_pos_um: npt.ArrayLike,
    reach_um: float,
    count: BondCount | None = None,
    k_xb_pn_per_um: float | None = None,
    r0_xb_um: float | None = None,
    kinetics: CrossbridgeKinetics | None = None,
) -> tuple[BondFamily, CrossbridgeState, dict[str, object]]:
    """Build the crossbridge family, its kinetic state and its census — or refuse by name.

    Args:
        arena: the :class:`~aleph.world.arena.WorldArena` to claim the BOND range from.
        nmii_inventory / nmii_pos_um / cortex / cortex_pos_um / reach_um: see
            :func:`plan_crossbridge_stations`.
        count: the ``per_filament`` :class:`BondCount`, from :func:`crossbridge_count`. ``None``
            raises rather than falling back — there is no fallback in this module to fall back to.
        k_xb_pn_per_um: crossbridge stiffness [pN/µm]. **REQUIRED**, PI-GAP.
        r0_xb_um: crossbridge rest length [µm]. **REQUIRED**, PI-GAP for two reasons — see the module
            docstring.
        kinetics: the catch-slip law. **REQUIRED**: a bond family with no rate law is a static spring
            wearing a motor's name, and this connector's whole point is that it is the FIRST kinetic
            one.

    Returns:
        ``(family, state, census)``. The census carries ``partner_is_static: True`` — the interface
        finding this module exists to produce.

    Raises:
        ConnectorGapError: while any magnitude is undeclared, or the cortex cannot say which
            filaments oppose.
        ValueError: if ``count`` is not ``per_filament``, if it does not resolve to exactly one site
            per head, or on any geometric build error.
    """
    missing = [name for name, value in (
        ("count", count), ("k_xb_pn_per_um", k_xb_pn_per_um),
        ("r0_xb_um", r0_xb_um), ("kinetics", kinetics)) if value is None]
    if missing:
        raise ConnectorGapError(CONNECTOR, missing, _GAP_NOTE_KINETICS)
    if count.basis != "per_filament":
        raise ValueError(
            f"{CONNECTOR}: a crossbridge population is counted PER MINIFILAMENT — one site per head, "
            f"2H per bipolar minifilament — not per unit area or volume. Got basis {count.basis!r}.")
    for name, v in (("k_xb_pn_per_um", k_xb_pn_per_um), ("r0_xb_um", r0_xb_um)):
        if not np.isfinite(float(v)) or float(v) < 0.0:
            raise ValueError(f"{name} must be finite and nonnegative; got {v!r}")

    plan = plan_crossbridge_stations(
        nmii_inventory=nmii_inventory, nmii_pos_um=nmii_pos_um, cortex=cortex,
        cortex_pos_um=cortex_pos_um, reach_um=reach_um)

    n_mf = int(nmii_inventory["n_minifilaments"])
    resolved = count.resolve(float(n_mf))
    if resolved != plan["n_heads"]:
        raise ValueError(
            f"{CONNECTOR}: the count resolves to {resolved} crossbridge sites over {n_mf} "
            f"minifilaments, but the built population carries {plan['n_heads']} heads. One site per "
            "head, exactly: a shortfall leaves heads that can never bind and an excess leaves sites "
            "with no head. Resolve the disagreement in the head count (build/nmii.py "
            "n_heads_per_side), never by rounding it away here.")

    # ⚠ A CROSSBRIDGE THAT CAN NEVER ATTACH IS NOT A CROSSBRIDGE.
    #
    # This tests the SAME condition the kernel's attach branch tests — `L <= capture_radius` in
    # laws/crossbridge_kmc.propose_crossbridge_transitions_kernel — against the BUILT configuration,
    # and it uses the caller's own declared capture radius rather than a number invented here.
    #
    # It is decisive at t0 precisely because CrossbridgeState is born all-free: nothing can move
    # until something binds, and nothing can bind beyond its capture radius. So if every bond is
    # out of range in the built configuration, the population is inert for any value of the eight
    # constants, and building it would produce 26,520 rows that read as a motor.
    #
    # Written after the 2026-08-21 TIER 0 result found exactly this: NMII on shell [6.850, 7.050]
    # and cortex on [7.300, 7.500] are DISJOINT by 0.250 um, against a capture-radius proxy of
    # 0.210 um. The build succeeded and reported n_outside: 0. This is that defect made
    # unrepresentable rather than left to be rediscovered.
    sep_built = np.asarray(plan["separation_um"], np.float64)
    n_beyond = int(np.count_nonzero(sep_built > float(kinetics.capture_radius_um)))
    if n_beyond == sep_built.size and sep_built.size:
        raise ValueError(
            f"{CONNECTOR}: all {sep_built.size} crossbridges are built beyond the declared capture "
            f"radius {kinetics.capture_radius_um:g} um (as-built separation min "
            f"{sep_built.min():.4g}, median {np.median(sep_built):.4g} um). The kernel's attach "
            "branch tests L <= capture_radius, and the state is born all-free, so nothing binds, "
            "nothing moves, and nothing ever binds — the population is inert for EVERY value of the "
            "eight PI-GAPs. This is a PLACEMENT result, not a kinetics one: check whether the NMII "
            "and cortex shells interpenetrate before asking what a rate constant should be. "
            "station_census().nearest_cortex_node_um reports the separation without building "
            "anything."
        )

    family = build_bond_family(
        arena, CONNECTOR, chemistry_card=CHEMISTRY_CARD, count=count, support=float(n_mf),
        pairs=plan["pairs"], rest_um=float(r0_xb_um),
        stiffness_pn_per_um=float(k_xb_pn_per_um),
        # ⚠ node_j is the AS-BUILT partner only. Declared through the builder since `f99400a6`
        # forwarded the pair; the owner is named because BondFamily refuses one without the other,
        # so this family cannot show node_j in an artifact as the partner in force at that step.
        # `__name__` rather than a literal path, so the citation follows the file if it moves.
        partner_is_dynamic=True,
        live_partner_owner=f"{__name__}.CrossbridgeState.bound")
    state = CrossbridgeState.all_free(family.n_bonds)

    sep = np.asarray(plan["separation_um"], np.float64)
    census = {
        "connector": CONNECTOR, "chemistry_card": CHEMISTRY_CARD,
        "declared_chemistry_card": DECLARED_CHEMISTRY_CARD,
        "kinetic": True, "first_kinetic_family": True,
        "n_bonds": family.n_bonds, "n_stations": plan["n_stations"],
        "n_heads_per_side": plan["n_heads_per_side"],
        "n_heads_per_side_is_a_test_point_not_a_band": True,
        "f_stall_head_pn_reference": F_STALL_HEAD_PN,
        "reach_um": plan["reach_um"],
        # The as-built separation is REPORTED beside the declared rest length, so whatever strain the
        # caller's r0_xb implies is visible in the artifact rather than inferable from it.
        "as_built_separation_um": {"min": float(sep.min()), "max": float(sep.max()),
                                   "mean": float(sep.mean())},
        "r0_xb_um": float(r0_xb_um),
        # How many bonds start outside their own capture radius. Zero is the only value that means
        # the whole population can participate; anything else is a placement report, not a rate one.
        "n_beyond_capture_at_t0": n_beyond,
        "initial_bound_fraction": state.bound_fraction,
        # ⚠ Stated, not worked around. world/step.py's default predicate is UndefinedAcceptance:
        # verdict ACCEPTANCE_UNDEFINED, is_accepted False by design (PI decision 18). So the kinetics
        # below have never run inside an accepted step, and this row says so rather than leaving a
        # reader to infer it from a bound fraction that never moved.
        "commit_reachable_today": False,
        "commit_reachable_today_note": (
            "world/step.py defaults to UndefinedAcceptance -> ACCEPTANCE_UNDEFINED -> is_accepted "
            "False (PI decision 18, 2026-08-21: the acceptance predicate is undecidable until "
            "cortical tension is emitted on the resting path). CrossbridgeState.commit refuses an "
            "unaccepted step, so nothing this family proposes can be applied on the current path."),
        # ⛔ THE INTERFACE FINDING — now DECLARED, not described. Raised by this family 2026-08-21
        # and answered in bond.py at `f99400a6`, so the census reads the interface's own fields
        # rather than spelling the same fact a second way. A second spelling is how a typo passes
        # silently, which is the argument SourceClass already makes about PROVENANCE.
        "partner_is_dynamic": family.partner_is_dynamic,
        "live_partner_owner": family.live_partner_owner,
        "n_bonds_is": family.provenance_row()["n_bonds_is"],
        "free_list": "unnecessary — the bond population IS the head population; a slot is a head's "
                     "identity, so detaching frees nothing",
        "snapshot_twin": "unnecessary — the KMC kernel writes a separate proposal array and cannot "
                         "touch committed state, so a rejected step is a non-copy, not a rollback",
        "provenance": family.provenance_row(),
        "kinetics": kinetics.record(),
    }
    return family, state, census


def _demo() -> None:
    """Self-check: the refusals are the product, and the straddle machinery behind them is correct."""
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.build.cortex import StrandPopulation

    assert not SPEC.buildable and SPEC.blocked_by == BLOCKED_BY
    assert SPEC.populations == ("cortex", "nmii") and SPEC.basis == "per_filament"
    assert CHEMISTRY_CARD != DECLARED_CHEMISTRY_CARD, "the declared card names two retired laws"
    assert CONNECTOR not in DUPLICATE_OF_ENTRY, "a connector cannot be a duplicate of itself"
    assert set(DUPLICATE_OF_ENTRY.values()) == {CONNECTOR} and len(DUPLICATE_OF_ENTRY) == 3

    # A BAND is refused. The axis has a test point and PI decision 6 ratified no band.
    for band in ([10, 30], (10, 30), np.array([10, 30])):
        try:
            crossbridge_count(band, scope="x")
        except TypeError as exc:
            assert "was NOT ratified" in str(exc) and "rethought" in str(exc)
        else:  # pragma: no cover
            raise AssertionError("a band must be refused: PI decision 6 ratified a test point only")

    scope = "MCF7 interphase adherent 37C — cortical NMII"
    count = crossbridge_count(N_HEADS_PER_SIDE_TEST_POINT, scope=scope)
    assert count.basis == "per_filament" and count.value == 60.0
    assert count.support_unit == "filaments" and "NO BAND" in count.provenance
    assert count.resolve(442.0) == 26_520, "442 minifilaments x 2 x 30 heads, one site each"

    # A small fixture: 6 cortex filaments of 5 nodes on two parallel lines, alternating polarity, and
    # 2 minifilaments straddling them. Everything is host geometry; no device is touched.
    n_str, n_per_str, h, n_bb = 6, 5, 3, 4
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.SEGMENT: 400, Kind.ANGLE3: 400, Kind.BOND: 400})
    cnodes = arena.claim("cortex", Kind.NODE, n_str * n_per_str)
    cor_pos = np.zeros((n_str * n_per_str, 3))
    for s in range(n_str):
        cor_pos[s * n_per_str:(s + 1) * n_per_str, 0] = np.arange(n_per_str) * 0.1
        cor_pos[s * n_per_str:(s + 1) * n_per_str, 1] = 0.0 if s % 2 == 0 else 0.4
        cor_pos[s * n_per_str:(s + 1) * n_per_str, 2] = (s // 2) * 1.0

    def _cortex(polarity_per_strand):
        return StrandPopulation(
            population="cortex", nodes=cnodes, segments=cnodes, angles=cnodes,
            n_strands=n_str, nodes_per_strand=n_per_str,
            seg_node=None, seg_rest_um=None, seg_arc_um=None, angle_idx=None,
            radius_um=7.4, thickness_um=0.2, centre_um=(0.0, 0.0, 0.0),
            areal_density_um2=100.0, density_provenance="self-check fixture",
            contour_um_requested=0.4, contour_um_realised=0.4, seg_um_requested=0.1,
            seg_um_realised=0.1, seed=0, polarity=+1, topology_bytes=0,
            polarity_per_strand=polarity_per_strand)

    n_per = n_bb + 2 * h
    n_mf = 2
    nmii_nodes = arena.claim("nmii", Kind.NODE, n_mf * n_per)
    mf_pos = np.zeros((n_mf * n_per, 3))
    for m in range(n_mf):
        b = m * n_per
        mf_pos[b:b + n_bb, 0] = 0.2
        mf_pos[b:b + n_bb, 1] = 0.2
        mf_pos[b:b + n_bb, 2] = m * 1.0
        mf_pos[b + n_bb:b + n_bb + h, 1] = 0.05          # + row, near the y=0.0 filament
        mf_pos[b + n_bb + h:b + n_per, 1] = 0.35         # - row, near the y=0.4 filament
        mf_pos[b + n_bb:b + n_per, 0] = 0.2
        mf_pos[b + n_bb:b + n_per, 2] = m * 1.0
    inventory = {"n_minifilaments": n_mf, "n_bb": n_bb, "n_heads_per_side": h,
                 "nodes_per_minifilament": n_per,
                 "claims": {"node": (nmii_nodes.lo, nmii_nodes.hi)}}

    mixed = np.array([+1, -1] * (n_str // 2), np.int8)
    uniform = np.ones(n_str, np.int8)

    # THE FIRST PRODUCT: a cortex that cannot say which filaments oppose refuses by TYPE, and the
    # refusal names PI decision 4 rather than reading as a defect in this file.
    for pol, why in ((None, "carries a single scalar polarity"), (uniform, "every filament runs")):
        try:
            plan_crossbridge_stations(nmii_inventory=inventory, nmii_pos_um=mf_pos,
                                      cortex=_cortex(pol), cortex_pos_um=cor_pos, reach_um=0.3)
        except ConnectorGapError as exc:
            assert exc.missing == ["cortex.polarity_per_strand"], exc.missing
            assert "UNREPRESENTABLE" in str(exc) and "PI decision 4" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a cortex that {why} cannot place a bipolar station")

    cortex = _cortex(mixed)
    # THE SECOND PRODUCT: no magnitude, no family — and the arena is untouched by the refusal.
    try:
        build_nmii_cortex_crossbridge(arena, nmii_inventory=inventory, nmii_pos_um=mf_pos,
                                      cortex=cortex, cortex_pos_um=cor_pos, reach_um=0.3, count=count)
    except ConnectorGapError as exc:
        assert exc.connector == CONNECTOR
        assert exc.missing == ["k_xb_pn_per_um", "r0_xb_um", "kinetics"], exc.missing
        assert "MASTER force knob" in str(exc) and "CATCH-SLIP" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("the crossbridge must refuse: every magnitude on this bond is a PI-GAP")
    assert arena.n_live(Kind.BOND) == 0, "a refused family must claim nothing"

    # BOTH degenerate stations are rejected, and for different reasons.
    try:
        assert_bipolar_station(cortex, 0, 0)
    except ValueError as exc:
        assert "pulls it against itself" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a minifilament straddling ONE filament is not a station")
    try:
        assert_bipolar_station(cortex, 0, 2)      # both +1
    except ValueError as exc:
        assert "TRANSLATES the pair" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("two same-polarity filaments are not a station")
    assert_bipolar_station(cortex, 0, 1)          # +1 against -1: the only legal arrangement

    # Behind the refusals the mechanism is exercised with EXPLICITLY LABELLED PI-GAP values, which is
    # the only source class this connector can carry today.
    provisional = CrossbridgeKinetics(
        k_catch0=0.35, x_catch_um=1.0e-3, k_slip0=0.35, x_slip_um=0.6e-3, k_on=50.0,
        capture_radius_um=0.210, scope="PROVISIONAL — self-check fixture, NOT a run configuration",
        provenance="ac_gate_b_cortex_motor_native.CATCH_SLIP_PROVISIONAL; every value an I0-B3 PI-GAP")
    fixture_count = crossbridge_count(
        h, scope="self-check fixture", provenance="fixture head count, not a physiological axis")
    family, state, census = build_nmii_cortex_crossbridge(
        arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex, cortex_pos_um=cor_pos,
        reach_um=0.3, count=fixture_count, k_xb_pn_per_um=1.0, r0_xb_um=0.0, kinetics=provisional)

    assert family.n_bonds == n_mf * 2 * h == family.bonds.count
    assert family.component_pairs(arena) == {("cortex", "nmii"): family.n_bonds}
    # Every head appears EXACTLY once: a head carries one molecular state, never two.
    assert np.unique(family.node_i).size == family.n_bonds
    # The two rows of each minifilament land on DISJOINT, anti-parallel filaments.
    stations = plan_crossbridge_stations(
        nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex, cortex_pos_um=cor_pos,
        reach_um=0.3)["station"]
    for a, b in stations:
        assert a != b and cortex.antiparallel(int(a), int(b)), (a, b)

    # A count that is not one site per head is refused with BOTH numbers.
    try:
        build_nmii_cortex_crossbridge(
            arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex,
            cortex_pos_um=cor_pos, reach_um=0.3, count=crossbridge_count(h + 1, scope="wrong"),
            k_xb_pn_per_um=1.0, r0_xb_um=0.0, kinetics=provisional)
    except ValueError as exc:
        assert "One site per head, exactly" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a capacity that is not the head count must refuse")

    # An areal density is a category error for a per-minifilament population.
    try:
        build_nmii_cortex_crossbridge(
            arena, nmii_inventory=inventory, nmii_pos_um=mf_pos, cortex=cortex,
            cortex_pos_um=cor_pos, reach_um=0.3, k_xb_pn_per_um=1.0, r0_xb_um=0.0,
            kinetics=provisional,
            count=BondCount(basis="areal", value=1.0, scope="x", source_class=SourceClass.PI_GAP,
                            provenance="wrong basis on purpose"))
    except ValueError as exc:
        assert "counted PER MINIFILAMENT" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an areal crossbridge density must be refused by name")

    # THE THIRD PRODUCT: the state is born ALL FREE — a bound t0 would be an imposed duty ratio — and
    # a transition commits ONLY on an accepted step.
    assert state.bound_fraction == 0.0 and state.n_commits == 0
    proposal = np.ones(family.n_bonds, np.int8)
    try:
        state.commit(proposal, accepted=False)
    except RuntimeError as exc:
        assert "not accepted" in str(exc) and "PI decision 18" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a kinetic transition must not commit on an unaccepted step")
    assert state.bound_fraction == 0.0 and state.n_commits == 0, "the refusal must change nothing"
    assert state.commit(proposal, accepted=True) == family.n_bonds
    assert state.bound_fraction == 1.0 and state.n_commits == 1
    for bad, expect in ((np.ones(3, np.int8), "entries for"), (np.full(family.n_bonds, 2, np.int8),
                                                               "bound or free")):
        try:
            state.commit(bad, accepted=True)
        except ValueError as exc:
            assert expect in str(exc)
        else:  # pragma: no cover
            raise AssertionError("a malformed proposal must refuse")

    # The census carries the interface finding and the provenance, so neither travels separately.
    assert census["partner_is_dynamic"] is True, "the interface field, not a second spelling"
    assert census["live_partner_owner"].endswith("CrossbridgeState.bound")
    assert census["n_bonds_is"].startswith("CAPACITY")
    assert family.partner_is_dynamic and family.provenance_row()["node_j_is"].startswith("the AS-BUILT")
    assert census["n_heads_per_side_is_a_test_point_not_a_band"] is True
    assert census["provenance"]["source_class"] == "PI_GAP"
    assert census["provenance"]["support_unit"] == "filaments"
    assert "PROVISIONAL" in census["kinetics"]["scope"]
    assert census["as_built_separation_um"]["max"] <= 0.3

    print(f"{CONNECTOR} self-check OK — BLOCKED on {len(BLOCKED_BY)}: {', '.join(BLOCKED_BY)}; "
          f"straddle + kinetic state verified on a {family.n_bonds}-crossbridge PI-GAP fixture; "
          f"bond.py interface finding ANSWERED at f99400a6: partner_is_dynamic declared")


if __name__ == "__main__":
    _demo()
