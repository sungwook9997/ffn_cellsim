"""What the production cell does NOT contain, pinned so an absence cannot stay invisible.

⚠ **WHY THIS EXISTS.** The PHASE 4 cell the 2026-08-21 tau runs were measured on stands eleven
populations and `build/nmii.py` is not among the builders that made them: **there is no motor in that
cell.** `build/cytosol.py` never ran either, which is why `census.live.grid_cell` is 0. Neither fact
appeared in any run record, and `aleph/docs/PLAN_PROTEINS_2026-08-22.md` was written -- an ordering
for adding myosin and integrin -- without either being stated as already missing from production.

The finding was already written down, in `test_observe_gamma.py`'s docstring for a test defending
something else entirely: `gamma_source` must be `None` and not `0.0`, because *"the motors contributed
nothing"* and *"there are no motors"* are different findings. **The distinction stood; the thing it
was distinguishing never reached a record or a test.**

⚠ **This is a RATCHET, not a target.** The pinned set is what is missing TODAY. Standing NMII or the
cytosol makes this test FAIL, and that is the intent -- the set may only shrink, and shrinking it is a
deliberate edit by whoever changes the cell, not a silent drift.
"""

from __future__ import annotations

import json
from pathlib import Path

from aleph.world.build import builders_not_standing

#: Builders that contribute NOTHING to the production cell. May only SHRINK, and shrinking it is a
#: deliberate edit. ⚠ `nmii` left this set 2026-08-24 by PI ruling: the minifilaments now stand, and
#: their crossbridges deliberately do NOT bind — a crossbridge that can never release is a permanent
#: weld, which the charter forbids, and with `UndefinedAcceptance` shipping nothing may commit an
#: attach. So the cell gains mass, drag and sterics and no contractility.
MISSING_TODAY = ("cytosol",)

#: What was missing from the CELL THE TAU SERIES RAN ON. ⚠ A separate constant from `MISSING_TODAY`
#: on purpose: one is a historical fact about a committed file and can never change, the other tracks
#: the current builders and shrinks as they are wired. Sharing one constant made wiring NMII fail the
#: record test, which was the first sign that the two were being conflated.
MISSING_FROM_THE_TAU_CELL = ("cytosol", "nmii")

#: The record the tau series -- and every stationarity verdict taken from it -- was measured on.
TAU_RECORD = Path(__file__).resolve().parents[3] / "aleph/outputs/ac/world_phase4/tau3_seed1.json"


def test_the_helper_reports_an_absent_builder_rather_than_an_empty_population() -> None:
    """A population the cell never stood must be REPORTED, not inferred from a list's length."""
    whole = ("chromatin", "cortex", "cytosol", "filopodium", "intermediate_filament",
             "lamellipodium", "lamina", "membrane", "microtubule", "nmii", "nuclear_envelope",
             "sf_arc", "stress_fiber")
    assert builders_not_standing(whole) == ()
    assert builders_not_standing([p for p in whole if p != "nmii"]) == ("nmii",)
    # ⚠ The one alias, and the reason it is here: `build/envelope.py` stands `nuclear_envelope`. Left
    # unmapped, EVERY cell would report `envelope` missing and bury the real absences under a naming
    # argument that is an OPEN PI disagreement (see test_families_census.py, xfail strict).
    assert "envelope" not in builders_not_standing(whole)


def test_the_tau_cell_has_no_motor_in_it() -> None:
    """The cell three tau seeds and every gamma verdict were measured on, named rather than implied.

    ⚠ This asserts a property of a COMMITTED RECORD, so it cannot be made true by changing code. If
    the record is ever re-generated from a cell that stands NMII, this fails and the tau verdicts
    taken on the old one must be re-examined rather than carried across.
    """
    record = json.loads(TAU_RECORD.read_text())
    assert "nmii" not in record["populations"], record["populations"]
    assert builders_not_standing(record["populations"]) == MISSING_FROM_THE_TAU_CELL
    assert record["census"]["live"]["grid_cell"] == 0, "cytosol absent but grid_cell claimed"
    assert record["census"]["live"]["bond"] == 0, "no motor, so nothing may hold a crossbridge"


def _populations_the_builders_stand() -> set[str]:
    """What the CURRENT builders stand, read from source without a device.

    ⚠ Source inspection, not a build: `build_all` and `build_remaining_populations` need CUDA, and a
    completeness ratchet that only runs on the GPU host is a ratchet that runs almost never.
    """
    import inspect
    import re

    import aleph.world.build as B
    import aleph.world.populations as P

    stood = set(re.findall(r'out\["([a-z_]+)"\]',
                           inspect.getsource(P.build_remaining_populations)))
    # `build_all` stands these three; it assigns them to a dataclass rather than a dict.
    stood |= {"cortex", "membrane", "nuclear_envelope"}
    assert "cortex" in stood and len(stood) > 5, stood
    del B
    return stood


def test_the_missing_set_may_only_shrink() -> None:
    """The ratchet on the CURRENT BUILDERS. ⚠ It read the frozen record instead, and so never fired.

    The first version of this test asserted `builders_not_standing(record["populations"])` — a
    property of `tau3_seed1.json`, a file. Wiring NMII into `populations.py` on 2026-08-24 left it
    GREEN, because the record still said eleven. A ratchet that reads a frozen artifact ratchets the
    artifact, not the code, and the two had already diverged by the time anyone looked.

    Derived from source on both sides, so a builder added tomorrow is covered with no edit here.
    """
    now = builders_not_standing(_populations_the_builders_stand())
    grew = set(now) - set(MISSING_TODAY)
    assert not grew, (
        f"builders that stand in NO production cell grew by {sorted(grew)}. A new builder the "
        f"production path never calls is the defect this file exists for -- wire it into "
        f"world/build/build_all or world/populations.py, or state here why it may not be.")
    shrank = set(MISSING_TODAY) - set(now)
    assert not shrank, (
        f"{sorted(shrank)} now stand(s) in production. That is the intended direction -- update "
        f"MISSING_TODAY to {tuple(sorted(now))} deliberately, so the change is a decision in a diff "
        f"rather than a threshold that moved on its own.")
