r"""``mt_nucleus_linc`` — the microtubule end against the nuclear envelope. IT DOES NOT BUILD.

This module's product is a refusal. Nothing here resolves a count, nothing here chooses a stiffness,
and the builder raises :class:`~aleph.world.families.ConnectorGapError` naming six unanswered questions.
That is the whole deliverable; the reasoning below is what makes it checkable rather than a shrug.

WHAT WAS PORTED FROM, AND WHAT WAS THERE.  ``aleph/engine/linc_connector.py`` declares the edge
(``MT_NUCLEUS_LINC``), ``microtubule_rig.py:434`` requires it ``ConnectorFamily.LINC``, bidirectional,
kinetic and accepted-step-committed, and the committed census
(``outputs/ac/connector_devicerun/native_record.json``, entry 8) records the result: verdict
``NOT_BOUND``, ``n_calls: 0``, ``"no runtime object in this composed world"``. There is no
implementation to port. There is a declaration, a force law shared with two siblings, and a constant
that is a PI-GAP.

──────────────────────────────────────────────────────────────────────────────────────────────────
IS IT A DUPLICATE?  NO — AND THE INCUMBENT'S CODE IS NOT EVIDENCE THAT IT IS.
──────────────────────────────────────────────────────────────────────────────────────────────────

The three LINC edges (``actin_cap_linc``, ``mt_nucleus_linc``, ``if_nucleus_linc``) run through ONE
adapter class, ONE kernel (``_linc_pair_force_kernel``), ONE law
(``f = k_linc·Δ·(1 + stiffening·Δ²)``, tension-only) and ONE constant
(``components/incumbent/compartments.py:113``, ``k_linc = 1.0e2``). Read as an implementation that is
three names over one population, and the NMII precedent in this package's ``__init__`` would fold them.

**It is the opposite.** The NMII case folds because one motor lands on two populations and the
population is derivable; here the CYTOSKELETAL SIDE IS A DIFFERENT MOLECULE IN EACH CASE:

* ``actin_cap_linc`` — nesprin-1/2 giant, a calponin-homology domain gripping F-actin directly. A
  load-bearing tether.
* ``if_nucleus_linc`` — nesprin-3 → plectin → vimentin/keratin. A crosslinker in series with the
  tether, so the compliance is not the same object's.
* ``mt_nucleus_linc`` — nesprin-4 / KASH5, which recruits **kinesin-1 and dynein**. The MT is not
  gripped; a motor anchored in the outer nuclear membrane WALKS on its lattice. Its own chemistry card
  in the incumbent says so: ``microtubule_motor_linc``.

The incumbent declares three distinct ``chemistry_card`` values and then gives all three the same
number. ``bond.py`` makes the card the family identity precisely so that cannot happen silently — so
what the shared constant records is not that the three are one physics, but that **the difference was
declared in names and erased in values.** That is the ERM defect in its second form: ERM hid a
population behind one declaration; this hides three parameterisations behind one literal. Folding the
three families here would ratify the erasure and make it unrecoverable, because a folded family has one
:class:`~aleph.world.bond.BondCount` and the three counts are not the same quantity (below).

⚠ **The fold is nevertheless a legitimate PI option and is filed as one, not performed here.** If the
PI rules that all three are modelled as one passive nesprin–SUN tether with a single ``k_linc`` and a
single perinuclear areal density — which is what the incumbent silently does — then they ARE one family
and ``DUPLICATE_OF`` gains two entries. That ruling is a modelling decision with a physiological
consequence (it deletes the motor), so it belongs to the PI. :data:`DUPLICATE_OF` is empty until then.

──────────────────────────────────────────────────────────────────────────────────────────────────
THE FOUR QUESTIONS, AND WHY NONE OF THEM HAS AN ANSWER
──────────────────────────────────────────────────────────────────────────────────────────────────

**How many.** LINC areal density is ``NOT FOUND``. ``PI_GAP_EVIDENCE_CARDS_2026-07-25.md`` card C2:
*"Zero LINC/nesprin claim rows in KB. PI source/choose."* ``PI_GAP_LITERATURE_SOURCING_2026-07-23.md``
§7 reaches the same wall independently.

**Against what.** The basis is not merely unmeasured, it is UNDECIDED, and the two candidates are not
convertible without the answer to the first question:

* ``areal`` over the nuclear envelope — resolved against ``Surface.area0_total_um2``, the MEASURED
  triangle area, never the analytic sphere (the property that keeps a density invariant when the
  continuum resolution changes);
* ``per_filament`` over MT plus-ends — which is what a KASH5/dynein site actually is, one per captured
  end, and is bounded above by ``n_mt`` rather than by envelope area.

``areal`` is recorded in :data:`SPEC` because the census needs a slot, **not** because it was chosen.

**Which cell.** The target is MCF7. The only sourced LINC number anywhere in the corpus is a resting
**tension** ≈ 8 pN (Déjardin 2020, ``10.1083/jcb.201908036``) measured in another system. A sourced
value out of scope is a proxy — ``SourceClass.UNRATIFIED_PROXY``, the class the PI's 2026-08-15 ruling
exists for — so even a borrowed density would not be ``SOURCED`` here.

**Whose authority.** Nobody's. There is no Contract-Graph row to cite, and ``ValidationGate`` /
``ModelContract`` rows are PI-authored by charter.

⚠ **AND THE STIFFNESS.** ``k_linc = 1.0e2 pN/µm`` carries its own refutation in the comment beside it:
*"8 pN is a TENSION not a stiffness"*. Card C1: no literature stiffness exists; spectrin-repeat
unfolding at 25–35 pN (Rief 1999) is a rupture scale, not a resting slope. A ``BondFamily`` requires
``stiffness_pn_per_um``, so this blocks the build independently of the count.

⚠ **THE TRAP IS ALREADY ON THE RECORD, PRE-REGISTERED.**
``COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md`` gate N9 defines the weld parameter
``W = k_linc·N_linc/k_shell`` and states the negative control: *"leave LINC density unset → the reported
density changes with subdivision and the manifest gate must FAIL"*. It also names the far side: ``W ≫ 1``
is auto-labelled RIGID-WELD REGIME and cannot support a transmission claim. So the exact failure this
module refuses to commit — a mesh number wearing a physiological label — was predicted for LINC by name,
a month before ERM demonstrated it. And ``laws/intermediate_filaments.py:133`` shows the sibling edge
already taking the fallback: one LINC bond per IF spoke, ``inner[s] → nearest nucleus bead``, count set
by ``n_fil``, itself unsourced.

──────────────────────────────────────────────────────────────────────────────────────────────────

WHAT IS DELIBERATELY ABSENT.  No kinetic law — ``bond.py`` defers attach/detach to the first family
that has kinetics, and if the motor reading survives PI review this edge IS that family, at which point
it is a dynein/kinesin force-velocity law and belongs in ``laws/``, not here. No fallback count, no
"provisional" stiffness, no default anywhere. A default is how the ERM count became a mesh artefact.

POPULATIONS.  The arena names are ``microtubule`` and ``nuclear_envelope``. The incumbent's endpoint is
``"nucleus"``, one component covering envelope + lamina + nucleoplasm; the bond lands on the ENVELOPE —
``NuclearSurfaceSocket`` requires ``ElementKind.TRIANGLE`` on a nuclear surface face. Recorded for the
census only: the runtime pair is DERIVED from the arena's ID ranges at query time and is never stored.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — nothing is computed. The gate that WOULD apply: a density is [1/µm²] against
    [µm²] or dimensionless against a filament count; ``rest_um`` [µm]; ``k_linc`` [pN/µm]. The basis
    ambiguity above is exactly a unit ambiguity, which is why it is a blocker and not a preference.
  * boundary — the refusal is total: there is no input that makes this module build, so there is no
    boundary case to bracket. ``SPEC.buildable`` is ``False`` unconditionally.
  * conservation/invariant — no bond is created, so no bond can address an unclaimed node. The
    invariant preserved is negative and it is the one that matters: **the arena's claim ledger and
    ``STATE.md`` gain nothing from this module**, so no count can later be quoted back from it.
  * CFL/precision — no integration, no array, no float. Host only.
  * sign sense — n/a. Noting for the unblocked case: the incumbent's law is TENSION-ONLY
    (``Δ ≥ 0``, ``linc_analytic.py:53`` — "a tether does not push"), and if the motor reading wins,
    sign sense stops being a property of the bond at all.
  * measurement protocol — no device is touched, nothing is uploaded, nothing is read back. The
    self-check runs on the host and reads one committed JSON to prove the connector name is the
    census's name and not a fresh typo.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable, and it
imports nothing from ``aleph.engine`` or ``aleph.components`` — those were READ, never bound.
"""

from __future__ import annotations

from typing import NoReturn

from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = ["CONNECTOR", "SPEC", "DUPLICATE_OF", "build_mt_nucleus_linc"]

CONNECTOR = "mt_nucleus_linc"

#: The unanswered questions, in ``BondCount``'s order plus the two this edge adds. Named, not valued:
#: a name Lead can queue is the product; a value would be an invention.
BLOCKED_BY: tuple[str, ...] = (
    "how-many: LINC/nesprin count per nucleus — NOT FOUND (PI_GAP_EVIDENCE_CARDS C2; zero LINC "
    "claim rows in the Contract-Graph)",
    "against-what: basis UNDECIDED — areal over measured envelope area, or per_filament over MT "
    "plus-ends (KASH5/dynein is one site per captured end); the two are not inter-convertible "
    "without the count",
    "which-cell: no LINC datum is MCF7 — Dejardin 2020 (10.1083/jcb.201908036) is a ~8 pN tension in "
    "another system, so any borrowed density is UNRATIFIED_PROXY, not SOURCED",
    "whose-authority: PI-authored value required; no ValidationGate/ModelContract row exists and one "
    "may not be auto-created",
    "stiffness: k_linc = 1.0e2 pN/um is a PI-GAP TEST value whose own comment refutes it — "
    "'8 pN is a TENSION not a stiffness' (card C1: no literature stiffness exists)",
    "is-it-a-bond-at-all: nesprin-4/KASH5 recruits kinesin-1 and dynein, so this coupling may be a "
    "MOTOR rather than a Hookean tether. ConnectorFamily.LINC vs MOTOR is a PI modelling decision, "
    "and if MOTOR wins this is the first family with kinetics — which bond.py defers by design",
)

SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("microtubule", "nuclear_envelope"),
    basis="areal",
    blocked_by=BLOCKED_BY,
    source="none — see BLOCKED_BY. The incumbent supplies a shared k_linc literal and no count at all.",
)

#: Empty on purpose, and the docstring argues why: the three LINC edges share one kernel and one
#: constant but NOT one molecule, and folding them would ratify an erasure rather than remove a
#: redundancy. The fold remains a PI option; performing it here would make it unrecoverable.
DUPLICATE_OF: dict[str, str] = {}


def build_mt_nucleus_linc(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse to build, naming what is missing.

    Takes any signature and honours none of it: the signature this family will have is not decidable
    until the basis is, because an ``areal`` count resolves against the envelope's measured area and a
    ``per_filament`` count against the MT census. Committing to arguments now would smuggle in the
    choice this module exists to escalate.

    Raises:
        ConnectorGapError: always. Carries :data:`BLOCKED_BY` verbatim.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(BLOCKED_BY),
        note=(
            "NOT a duplicate of actin_cap_linc / if_nucleus_linc — different molecule on the "
            "cytoskeletal side; see this module's docstring. Closing this needs a PI decision, not a "
            "search: card C1 (k_linc stiffness) + card C2 (LINC density) are both NOT FOUND, and gate "
            "N9 already pre-registered the failure that inventing them would cause."
        ),
    )


def _demo() -> None:
    """Self-check: the refusal is total, complete, and about the connector the census names."""
    import json
    import pathlib

    assert not SPEC.buildable, "a spec with unanswered questions must never report buildable"
    assert SPEC.connector == CONNECTOR
    assert len(SPEC.blocked_by) == 6, SPEC.blocked_by
    # Every question is NAMED, so Lead can queue it by its name rather than by its prose.
    names = [q.split(":", 1)[0] for q in SPEC.blocked_by]
    assert len(set(names)) == len(names), f"duplicate question names: {names}"
    assert all(name and " " not in name for name in names), names

    # The builder refuses whatever it is handed, and loses nothing on the way out.
    for args, kwargs in (((), {}), ((object(),), {"count": 235.0, "support": 314.0})):
        try:
            build_mt_nucleus_linc(*args, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert exc.missing == list(BLOCKED_BY)
            assert "NOT a duplicate" in exc.note
        else:  # pragma: no cover - the builder has no success path
            raise AssertionError("build_mt_nucleus_linc must refuse")

    # Not a duplicate, and it does not quietly defer to a sibling either.
    assert DUPLICATE_OF == {}, "the fold is the PI's to make; this module records it, never performs it"

    # The name is the CENSUS's name. A typo here would invent a connector that nothing declared, and
    # the census is the only committed record of which thirty exist.
    census = pathlib.Path(__file__).resolve().parents[2] / "outputs/ac/connector_devicerun/native_record.json"
    if census.is_file():
        record = json.loads(census.read_text())
        entry = next(c for c in record["connectors"] if c["connector"] == CONNECTOR)
        assert entry["verdict"] == "NOT_BOUND" and entry["n_calls"] == 0, entry
        # The census does not fold the three either — it records no sharing.
        assert entry["shared_with"] == [], entry["shared_with"]
        assert CONNECTOR in record["summary"]["unbound_connectors"]

    print(f"{CONNECTOR}: BLOCKED on {len(SPEC.blocked_by)} questions, builds nothing —", names)


if __name__ == "__main__":
    _demo()
