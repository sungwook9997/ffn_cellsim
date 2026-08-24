"""Where each PHASE 1 population is allowed to be, DERIVED from the cell that was built.

⚠ **Moved out of `aleph/scripts/` on 2026-08-21, and the move is a defect report.** It was written
there because both callers are drivers, and the Lead then told a family session to import it — which
`world/` may not do. `tests/architecture/test_layer_directions.py` caught it on the first commit that
tried: *"`world/` reached upward. Aleph's runtime/vertical cycle reached 53 edges one good reason at a
time; this is what stops the first one."* A placement rule the arena's own builders have to obey belongs
in `world/`, not beside the scripts that happened to need it first.

**Why this file exists.** The PHASE 1 driver and the renderer each carried their own copy of the
placement arguments — ``origin=(0,0,-7.0)``, ``length_um=10.0``, ``pitch_um=0.5``, ``width_um=8.0``.
Those numbers were typed into a call site and carried no provenance at all. Rendered on 2026-08-20 and
then measured, they put **95.5% of the stress-fibre population and 92.7% of the lamellipodium outside
the cell**: 50 fibres at a 0.5 µm pitch need a 24.5 µm lattice, and the basal contact disc of this cell
is 5.39 µm across — 4.5x too small.

⚠ **This is the ERM defect class again, in its worst costume so far.** The COUNTS went through
``BondCount`` and came out labelled ``PI_GAP``, which at least announces itself as a blank. The
GEOMETRY went through nothing. It was not a wrong physiological value, it was not even a placeholder —
it was an unlabelled number in a call site, and the population it placed did not fit inside the cell
while ``assert_partitioned`` reported success. Partition is an ID-range property. **It says nothing
about where a node is.**

**The rule, PI 2026-08-20: every node of every population lies inside the built membrane.**
:func:`assert_inside_membrane` enforces it after the build rather than trusting the arithmetic here,
because an invariant that is only argued is not an invariant.

⚠ **What this rule costs, stated rather than hidden.** A filopodium is a membrane PROTRUSION — the
membrane wraps around it and follows it out. This cell's membrane is a rigid icosphere that cannot
deform, so "inside the membrane" and "protruding" cannot both hold. The filopodia here are therefore
rooted so their tips stop a clearance SHORT of the membrane: geometrically they satisfy the rule,
biologically they are now cortical bundles rather than protrusions. **A protruding filopodium needs a deformable
membrane, which PHASE 1 does not have.** That is a PHASE 2 property, and it is recorded here so the
picture is not read as a cell whose filopodia have been quietly retracted.

engine units: length µm.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["CellFootprint", "assert_inside_membrane", "footprint"]


class CellFootprint:
    """The placement envelope every population is derived from.

    Attributes:
        r_cell_um: the built membrane radius [µm]. Read from the cell, never assumed.
        z_basal_um: the substrate plane [µm].
        r_footprint_um: radius of the basal contact disc, ``sqrt(R^2 - z^2)`` [µm]. This is the number
            that was missing: it is what a ventral stress fibre and a lamellipodium have to fit inside,
            and it is much smaller than a reader expects from ``R_cell`` alone.
        sf_pitch_um: lateral spacing between ventral stress fibres [µm]. Kept at the physiological
            value and NOT squeezed to fit — the count yields instead, which is the direction that keeps
            a sourced quantity sourced.
        n_stress_fibres: fibres the INSCRIBED footprint admits at that pitch. **DERIVED, so this is no
            longer a PI-GAP.** ``floor(2 half_side / pitch) + 1``, the same shape as the cortex's
            density x area.
        sf_length_um: fibre length [µm] — a side of the inscribed square, so both FA ends land on the
            contact disc and so do the corners.
        lamellipodium_width_um / lamellipodium_depth_um: the leading-edge sheet, sized to the same disc.
        filopodium_contour_um / filopodium_root_r_um: rooted so the TIP stops inside the membrane. See
            the module warning about what that costs.
    """

    __slots__ = ("r_cell_um", "z_basal_um", "r_footprint_um", "sf_pitch_um", "n_stress_fibres",
                 "sf_length_um", "lamellipodium_width_um", "lamellipodium_depth_um",
                 "filopodium_contour_um", "filopodium_root_r_um", "clearance", "tip_clearance_um",
                 "sf_pitch_source_class", "n_stress_fibres_source_class",
                 "sf_origin_um", "lamellipodium_origin_um", "lamellipodium_contour_um")

    def __init__(self, r_cell_um: float, z_basal_um: float, sf_pitch_um: float,
                 filopodium_contour_um: float, sf_pitch_source_class: str = "PI_GAP") -> None:
        if not (0.0 < abs(z_basal_um) < r_cell_um):
            raise ValueError(
                f"the basal plane z={z_basal_um} does not cut the cell of radius {r_cell_um}: there is "
                "no contact disc, so there is nowhere for a ventral fibre or a lamellipodium to be."
            )
        self.r_cell_um = float(r_cell_um)
        self.z_basal_um = float(z_basal_um)
        self.r_footprint_um = math.sqrt(r_cell_um * r_cell_um - z_basal_um * z_basal_um)
        # ⚠ THE CONTACT REGION IS A DISC AND EVERY BUILDER LAYS DOWN A RECTANGLE. Sizing that
        # rectangle to the disc DIAMETER is what the second attempt got wrong: a corner at
        # (r_f, r_f, z) sits at sqrt(2 r_f^2 + z^2), which is OUTSIDE the sphere even though every
        # edge midpoint is inside. The rectangle has to be INSCRIBED — half-side r_f/sqrt(2) — and
        # then its corner lands at exactly sqrt(r_f^2 + z^2) = R, i.e. ON the membrane.
        #
        # `clearance` pulls it strictly inside and also absorbs what the builders add on top of the
        # footprint they are given: a stress fibre is a BUNDLE (bundle_count x spacing_um) and so is a
        # filopodium, so the outermost filament of the outermost bundle sits further out than the
        # centreline this envelope places. The margin is not tuned to make a check pass — it is the
        # bundle half-width plus a round number, and the check is run afterwards either way.
        self.clearance = 0.85
        half_side = self.clearance * self.r_footprint_um / math.sqrt(2.0)
        self.sf_pitch_um = float(sf_pitch_um)
        self.n_stress_fibres = int(2.0 * half_side / sf_pitch_um) + 1
        # ⚠ A DERIVED count is only as sourced as the thing it was derived FROM, and this class first
        # labelled n_stress_fibres `DERIVED` outright. Session D's sf_arc builder got it right and the
        # pattern is the general one: **a derivation must not launder its input.** The count inherits
        # the pitch's source class, so the day a sourced lamellar spacing is registered the same call
        # starts reporting DERIVED without anyone editing this line — and until then it reports what
        # the pitch actually is.
        self.sf_pitch_source_class = str(sf_pitch_source_class)
        self.n_stress_fibres_source_class = self.sf_pitch_source_class
        self.sf_length_um = 2.0 * half_side
        self.lamellipodium_width_um = 2.0 * half_side
        self.lamellipodium_depth_um = half_side
        # ⚠ `origin` is where the rectangle STARTS, not where it is centred — `stress_fiber.py`'s own
        # docstring says "the basal-plane point the centre fibre starts at, its first FA end". Passing
        # (0, 0, z) therefore pushed the whole sheet into +x and put its far corners outside, which is
        # what the third attempt measured. These origins put the SPAN across the disc centre.
        self.sf_origin_um = (-self.sf_length_um / 2.0, 0.0, self.z_basal_um)
        # The lamellipodium's filaments radiate `contour_um` beyond the sheet edge at `mode_deg`, so the
        # sheet itself has to stop that much short of the envelope. Subtracted, not absorbed by margin.
        self.lamellipodium_contour_um = 1.0
        self.lamellipodium_depth_um = max(half_side - self.lamellipodium_contour_um, 0.2)
        self.lamellipodium_width_um = max(2.0 * (half_side - self.lamellipodium_contour_um), 0.4)
        self.lamellipodium_origin_um = (-self.lamellipodium_depth_um / 2.0, 0.0, self.z_basal_um)
        # A filopodium's TIP must land inside too, so it is rooted a clearance short of the membrane
        # rather than exactly on it — a tip at exactly R is on the boundary, and in float32 half of
        # those nodes read as outside. That is not a tolerance question to be widened away: a
        # structure placed ON a boundary has no margin for the bundle spread around it.
        self.tip_clearance_um = 0.25
        self.filopodium_contour_um = float(filopodium_contour_um)
        self.filopodium_root_r_um = self.r_cell_um - self.tip_clearance_um - self.filopodium_contour_um
        # ⚠ THE NMII SHELL IS NOT DERIVED HERE ANY MORE, AND THAT IS THE POINT.
        #
        # This class used to carry `nmii_radius_um = r_cell - tip_clearance - 0.30` and
        # `nmii_thickness_um = 0.20` — a SECOND formula for the cortical shell, which
        # `build/__init__.py` already derives from the native counts. The two drifted by 0.45 µm, and
        # the consequence was not cosmetic: the NMII population stood clear of the cortex it exists to
        # pull, and 0 of 442 minifilaments could station (`STATE.md` (c) 20). PI queue 14 named the
        # cause exactly — *nothing owned the relationship between the two shells.*
        #
        # ⚠ The 0.30 also charged a TANGENTIAL head arm as a radial one: `build/nmii.py` builds it as
        # `ey = cross(c, ex)` and computes its true radial cost as offset²/2r ≈ 2.7 nm, about seventy
        # times less than reserved. But repairing the constant would not have been enough — with
        # `tip_clearance` alone the formula reaches 7.25 against a shell at 7.40, so no subtraction
        # gets there. The SHAPE was wrong: the NMII shell is not the membrane pulled in by a
        # clearance, it is the cortical shell.
        #
        # Callers take it from `aleph.world.build.cortex_shell()`, which is the one derivation.
        if self.filopodium_root_r_um <= 0.0:
            raise ValueError(
                f"a filopodium of contour {filopodium_contour_um} µm cannot be rooted inside a cell of "
                f"radius {r_cell_um} µm and still end on the membrane."
            )

    def as_record(self) -> dict[str, object]:
        """The row that travels with the run, so the picture states what sized it."""
        return {
            "r_cell_um": self.r_cell_um,
            "z_basal_um": self.z_basal_um,
            "r_footprint_um": round(self.r_footprint_um, 4),
            "footprint_basis": "sqrt(R_cell^2 - z_basal^2) — the basal contact disc, DERIVED from the "
                               "built membrane radius, not a call-site constant",
            "sf_pitch_um": self.sf_pitch_um,
            "n_stress_fibres": self.n_stress_fibres,
            "sf_pitch_source_class": self.sf_pitch_source_class,
            "n_stress_fibres_source_class": self.n_stress_fibres_source_class,
            "n_stress_fibres_inheritance": "the count INHERITS the pitch's source class — a derivation "
                                           "must not launder its input. Register a sourced pitch and "
                                           "this reports DERIVED with no edit here.",
            "n_stress_fibres_basis": "floor(2 half_side / pitch) + 1. The pitch is held and "
                                     "the count yields, so the sourced quantity stays sourced.",
            "sf_length_um": round(self.sf_length_um, 4),
            "lamellipodium_width_um": round(self.lamellipodium_width_um, 4),
            "lamellipodium_depth_um": round(self.lamellipodium_depth_um, 4),
            "filopodium_root_r_um": round(self.filopodium_root_r_um, 4),
            "sf_origin_um": [round(x, 4) for x in self.sf_origin_um],
            "lamellipodium_origin_um": [round(x, 4) for x in self.lamellipodium_origin_um],
            "origin_note": "`origin` is where a builder STARTS its rectangle, not where it centres it; "
                           "passing the disc centre pushed both sheets into +x and their far corners "
                           "outside the membrane",
            "inscribed_square_clearance": self.clearance,
            "inscribed_square_note": "the contact region is a DISC and every builder lays down a "
                                     "RECTANGLE; sizing it to the disc diameter puts the corners "
                                     "outside the sphere while every edge midpoint is inside",
            "filopodium_note": "rooted so the tip REACHES the membrane and stops. A real filopodium "
                               "protrudes and the membrane follows it; this membrane is a rigid "
                               "icosphere, so protrusion and 'inside the membrane' cannot both hold. "
                               "PHASE 2 property, recorded rather than silently traded away.",
        }


def footprint(r_cell_um: float, z_basal_um: float = -7.0, *, sf_pitch_um: float = 0.5,
              filopodium_contour_um: float = 3.0,
              sf_pitch_source_class: str = "PI_GAP") -> CellFootprint:
    """Build the placement envelope for a cell of radius ``r_cell_um``.

    Args:
        sf_pitch_source_class: what the PITCH is, in ``bond.SourceClass`` terms. The derived fibre
            count inherits it — see :class:`CellFootprint`. Defaults to ``PI_GAP`` because no sourced
            ventral-fibre lamellar spacing is registered in the Contract-Graph today.
    """
    return CellFootprint(r_cell_um, z_basal_um, sf_pitch_um, filopodium_contour_um,
                         sf_pitch_source_class)


def assert_inside_membrane(pos: np.ndarray, pops: dict, r_cell_um: float, *,
                           tol_um: float = 1e-6) -> dict[str, dict[str, float]]:
    """Check the PI rule against the BUILT positions, and report per population.

    Args:
        pos: the arena's host node array.
        pops: population name -> object exposing either ``claims["node"]`` or ``.nodes``.
        r_cell_um: the built membrane radius.
        tol_um: slack for the membrane's own vertices, which sit AT the radius and land either side of
            it in float32. Nothing else needs it.

    Returns:
        Per population: node count, max radius, and how many nodes broke the rule.

    Raises:
        AssertionError: naming every population that has a node outside, with its worst radius. The
            message is the finding — a driver that placed a population outside the cell should say so
            in the run that did it, not in a later audit.
    """
    # ⚠ An empty population dict passes vacuously, and that is the same defect class this function
    # exists to catch: `assert_partitioned` reported success while 95.5% of a population sat outside
    # the cell, because it was checking a property that did not cover the question. A caller that
    # built nothing and got back a clean "all inside" would be in exactly that position. Found on a
    # fresh read 2026-08-21, alongside the identical bug in world/step.py's as_record.
    if not pops:
        raise ValueError(
            "refused: nothing to check. An empty population set is not 'every node is inside' — it "
            "is a caller that built no cell, and returning a clean report for one is how a check "
            "certifies a question it was never asked."
        )

    report: dict[str, dict[str, float]] = {}
    broke: list[str] = []
    for name, obj in pops.items():
        lo, hi = (obj["claims"]["node"] if isinstance(obj, dict) else (obj.nodes.lo, obj.nodes.hi))
        if hi <= lo:
            # Said, not skipped and not crashed. numpy raises "zero-size array to reduction" here,
            # which tells a reader nothing about which population is empty or why that matters.
            raise ValueError(
                f"population {name!r} claims no nodes ([{lo}, {hi})). A population that built nothing "
                "cannot be checked for being inside anything, and passing it silently would put an "
                "empty row in a report that reads as a pass."
            )
        r = np.linalg.norm(pos[lo:hi], axis=1)
        n_out = int((r > r_cell_um + tol_um).sum())
        report[name] = {"n_nodes": int(hi - lo), "r_max_um": float(r.max()), "n_outside": n_out}
        if n_out:
            broke.append(f"{name}: {n_out:,}/{hi - lo:,} nodes out, worst r = {r.max():.3f} µm")
    assert not broke, (
        f"PI rule 2026-08-20 — every node lies inside the built membrane (R = {r_cell_um} µm) — is "
        f"violated by: {'; '.join(broke)}. `assert_partitioned` passing does NOT cover this: partition "
        "is a property of ID ranges and says nothing about where a node is."
    )
    return report


def _demo() -> None:
    """Self-check: the derivation, the rule, and the case that motivated both."""
    fp = footprint(7.5, -7.0)

    # 1. The number that was missing, and how far off the old call site was.
    assert abs(fp.r_footprint_um - 2.6926) < 1e-3, fp.r_footprint_um
    expected = int(2.0 * fp.clearance * fp.r_footprint_um / math.sqrt(2.0) / 0.5) + 1
    assert fp.n_stress_fibres == expected, (fp.n_stress_fibres, expected)   # not the 50 typed in
    # the count is exactly as sourced as the pitch, and says so
    assert fp.n_stress_fibres_source_class == "PI_GAP"
    assert footprint(7.5, sf_pitch_source_class="SOURCED").n_stress_fibres_source_class == "SOURCED"
    old_lattice_um = (50 - 1) * 0.5
    assert old_lattice_um / (2 * fp.r_footprint_um) > 4.0, "the regression this file exists for"

    # 2. Everything the envelope sizes must fit the disc it was derived from.
    # the CORNER, not the edge — this is the check the second attempt did not have
    for o, along, across, pad in (
            (fp.sf_origin_um, fp.sf_length_um, (fp.n_stress_fibres - 1) * fp.sf_pitch_um, 0.0),
            (fp.lamellipodium_origin_um, fp.lamellipodium_depth_um, fp.lamellipodium_width_um,
             fp.lamellipodium_contour_um)):
        # the two ends are o[0] and o[0] + along; the far one is whichever has the larger |x|. Writing
        # this as |o[0]| + along double-counts the offset — caught by this very check on 2026-08-20.
        far_x = max(abs(o[0]), abs(o[0] + along)) + pad
        corner = math.sqrt(far_x ** 2 + (across / 2 + pad) ** 2 + fp.z_basal_um ** 2)
        assert corner < fp.r_cell_um, (o, along, across, corner)
    assert abs(fp.filopodium_root_r_um + fp.filopodium_contour_um + fp.tip_clearance_um
               - fp.r_cell_um) < 1e-9
    # ⚠ The NMII-shell assertion that stood here is GONE WITH THE FIELD IT CHECKED. It read
    #     `nmii_radius + nmii_thickness/2 + 0.20 < r_cell`
    # and it PASSED throughout the period when the motors could not reach the actin, because it
    # checked the shell against the MEMBRANE and nothing checked it against the CORTEX. An assertion
    # that holds while the thing it is about is broken is worse than none. The relationship now has
    # one owner (`build.cortex_shell`), and what enforces it is `build/nmii.py`'s inset plus
    # `assert_inside_membrane` — a fit check and a containment check, each in its own scope.

    # 3. A basal plane that misses the cell is refused rather than yielding an imaginary disc.
    for z in (-7.5, -9.0, 0.0):
        try:
            footprint(7.5, z)
        except ValueError:
            pass
        else:
            raise AssertionError(f"z={z} must be refused")

    # 4. The checker fires, and its message names the population. A guard that cannot fail certifies
    #    nothing — this engine has been bitten by exactly that (`balance_ok`, `descent_ratio`).
    class _P:
        def __init__(self, lo, hi):
            self.nodes = type("R", (), {"lo": lo, "hi": hi})()

    pos = np.zeros((4, 3), dtype=np.float64)
    pos[2] = (0.0, 0.0, 9.0)                                     # one node 1.5 µm outside
    try:
        assert_inside_membrane(pos, {"good": _P(0, 2), "bad": _P(2, 4)}, 7.5)
    except AssertionError as exc:
        assert "bad" in str(exc) and "good" not in str(exc), exc
        assert "9.0" in str(exc) or "9.000" in str(exc), exc
    else:
        raise AssertionError("the rule checker did not fire on a node outside the membrane")

    rep = assert_inside_membrane(pos[:2], {"good": _P(0, 2)}, 7.5)
    assert rep["good"]["n_outside"] == 0

    # 5. The two vacuous cases, both found on a fresh read rather than by a failure. An empty set is
    #    not "everything is inside", and an empty population cannot be checked for being inside.
    for bad, needle in (({}, "built no cell"), ({"empty": _P(1, 1)}, "claims no nodes")):
        try:
            assert_inside_membrane(pos, bad, 7.5)
        except ValueError as exc:
            assert needle in str(exc), (bad, exc)
        else:
            raise AssertionError(f"{bad} must be refused, not reported as a pass")
    print("world_cell_geometry self-check OK — "
          f"r_footprint {fp.r_footprint_um:.3f} µm, {fp.n_stress_fibres} stress fibres derived")


if __name__ == "__main__":
    _demo()
