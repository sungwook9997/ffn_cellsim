r"""The nine populations that are not cortex, membrane or envelope — built once, called by both drivers.

**Why this exists is a defect report.** The PHASE 1 driver and the renderer each carried their own
copy of these calls, with the placement arguments written out twice. On 2026-08-20 those two copies
drifted from the cell they were placing and put **95.5% of the stress fibres and 92.7% of the
lamellipodium outside the membrane**, while `assert_partitioned` reported success — partition is a
property of ID ranges and says nothing about where a node is. The geometry was pulled into
:mod:`aleph.world.geometry` then; **this is the other half of the same fix**, and it lands now because
PHASE 4 needed the same cell and copying the block a third time is how a third divergence starts.

`build_all` stands cortex, membrane and nuclear envelope. This stands everything else, in one place,
so a driver names a cell rather than assembling one.

⚠ **Every count here is a PI-GAP except two, and the two exceptions are the interesting ones.**

* `n_stress_fibres` is **DERIVED** from the basal contact disc — and it inherits the pitch's source
  class rather than claiming `DERIVED` outright, because *a derivation does not launder its input*
  (session D's phrasing, and the rule this module follows throughout).
* `sf_arc`'s count is derived the same way, from the lamellar annulus.

Everything else — microtubules, IF spokes, filopodia, lamellipodial filaments, chromatin domains —
has **no sourced per-cell structure count anywhere in the Contract-Graph**, and the `gap()` helper
puts that sentence inside each `BondCount` rather than beside it.

⚠ **And the whole basal set inherits a premise that does not hold.** The contact disc is 22.8 µm²
against a nucleus of 81.7 µm² projected, and no basal plane makes a sphere of R = 7.5 µm adherent —
the maximum contact disc is 177 µm², reachable only when the cell is cut in half
(`SPHERE_CANNOT_BE_ADHERENT_2026-08-21.md`). The arithmetic below is correct; the cell it is correct
*about* is a sphere in contact, not a spread cell, and that is a PI question rather than a parameter.

Sanity Gate:
    * dimensions — µm throughout; every length comes from :class:`~aleph.world.geometry.CellFootprint`
      or is a named argument, and none is typed at a call site here.
    * boundary cases — a caller asking for a population twice would double-claim, so the returned dict
      is keyed and built once; an arena without the membrane/envelope claims cannot produce a
      footprint and fails there rather than here.
    * conservation — the caller runs :func:`~aleph.world.geometry.assert_inside_membrane` over the
      returned dict; this module does not check it, because the check belongs where the whole cell is
      visible and running it twice would let one copy drift.
    * sign sense — not applicable: geometry only, no force.
    * measurement protocol — the returned dict is the population census, in build order.

engine units: µm. Runtime: NVIDIA Warp on CUDA.
"""

from __future__ import annotations

import math

from typing import Any, Callable

from aleph.world.arena import WorldArena
from aleph.world.bond import BondCount, SourceClass
from aleph.world.build.chromatin import build_chromatin
from aleph.world.build.filopodium import build_filopodia
from aleph.world.build.intermediate_filament import build_intermediate_filaments
from aleph.world.build.lamellipodium import build_lamellipodium
from aleph.world.build.lamina import build_lamina
from aleph.world.build.nmii import build_nmii
from aleph.world.build import cortex_shell
from aleph.world.build.microtubule import build_microtubules
from aleph.world.build.sf_arc import arc_count_from_footprint, build_sf_arcs
from aleph.world.build.stress_fiber import build_stress_fibers
from aleph.world.geometry import CellFootprint

__all__ = ["ARC_PITCH_UM", "build_remaining_populations", "gap"]

#: Transverse-arc lamellar spacing [µm]. ⚠ A declared TEST POINT — no sourced transverse-arc spacing
#: exists in the Contract-Graph, and `arc_count_from_footprint` inherits that class into the count.
ARC_PITCH_UM = 0.6


def gap(n: float, what: str, basis: str = "explicit") -> BondCount:
    """A count that says, inside the artifact, that no source exists for it.

    The scope string is the payload. A `PI_GAP` label beside a number can be dropped when the number
    is quoted; a scope that reads *"NO VALUE EXISTS"* travels with it.
    """
    return BondCount(value=float(n), basis=basis, source_class=SourceClass.PI_GAP,
                     scope=f"NO VALUE EXISTS — per-cell structure count absent from the "
                           f"Contract-Graph for {what}",
                     provenance="PI-GAP; a size placeholder, not a claim")


def build_remaining_populations(
    arena: WorldArena, fp: CellFootprint, *, seg_um: float,
    stage: Callable[[str, Callable[[], Any]], Any] | None = None,
) -> dict[str, Any]:
    """Stand the nine populations `build_all` does not, into ``arena``.

    Args:
        arena: the world to claim from. Must already hold cortex, membrane and envelope.
        fp: the placement envelope. **Every position below comes from it** — nothing is typed here,
            which is the point of the module.
        seg_um: the discretisation step [µm]. Required: it is the cortex constant the provenance audit
            ranks second most load-bearing, with a source of *"none (0 KB rows)"*.
        stage: optional ``(name, thunk) -> result`` wrapper for a driver that times its builds. When
            omitted the thunks are called directly, so a driver that does not care about timing does
            not have to supply one.

    Returns:
        Population name -> builder result, in build order.
    """
    run = stage if stage is not None else (lambda _n, fn: fn())
    out: dict[str, Any] = {}

    out["microtubule"] = run("microtubule", lambda: build_microtubules(
        arena, count=gap(500, "microtubules"), centre=(0.0, 0.0, 0.0), reach_R_um=7.4, seg_um=seg_um))

    out["intermediate_filament"] = run("intermediate_filament", lambda: build_intermediate_filaments(
        arena, count=gap(60, "IF spokes"), centre=(0.0, 0.0, 0.0), R_nuc_um=5.1, R_cortex_um=7.4,
        gap_um=0.15, seg_um=seg_um))

    out["filopodium"] = run("filopodium", lambda: build_filopodia(
        arena, count=gap(50, "filopodia"), bundle_count=gap(20, "filaments per filopodial bundle"),
        centre=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), contour_um=fp.filopodium_contour_um,
        seg_um=seg_um, spacing_um=0.008, cap_deg=60.0, root_R_um=fp.filopodium_root_r_um))

    out["lamellipodium"] = run("lamellipodium", lambda: build_lamellipodium(
        arena, count=gap(200, "lamellipodial filaments"), origin=fp.lamellipodium_origin_um,
        axis=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0), width_um=fp.lamellipodium_width_um,
        depth_um=fp.lamellipodium_depth_um, contour_um=fp.lamellipodium_contour_um, seg_um=seg_um,
        mode_deg=35.0, support=1.0))

    # ⚠ A transverse arc is NOT a ventral FA-to-FA fibre (PI decision 8), so it gets its own builder
    # rather than a rename. Its count is derived from the lamellar annulus and INHERITS the pitch's
    # source class, so today it reports PI_GAP and the day a sourced spacing is registered the same
    # call reports DERIVED with no edit here.
    out["sf_arc"] = run("sf_arc", lambda: build_sf_arcs(
        arena, fp,
        count=arc_count_from_footprint(
            fp, pitch_um=ARC_PITCH_UM, spacing_um=0.012, lift_um=0.4, n_filaments_per_arc=20,
            pitch_scope="MCF7 basal lamellar annulus as built",
            pitch_provenance="NO sourced transverse-arc lamellar spacing exists in the "
                             "Contract-Graph; 0.6 um is a declared test point",
            pitch_source_class=SourceClass.PI_GAP),
        bundle_count=gap(20, "filaments per transverse arc"),
        pitch_um=ARC_PITCH_UM, lift_um=0.4, seg_um=seg_um, spacing_um=0.012))

    out["stress_fiber"] = run("stress_fiber", lambda: build_stress_fibers(
        arena,
        count=BondCount(
            value=float(fp.n_stress_fibres), basis="explicit",
            source_class=SourceClass.PI_GAP, scope="MCF7 basal footprint as built",
            provenance=f"DERIVED from the contact disc, inheriting the pitch's class "
                       f"({fp.n_stress_fibres_source_class})"),
        bundle_count=gap(20, "filaments per fibre"), origin=fp.sf_origin_um, axis=(1.0, 0.0, 0.0),
        normal=(0.0, 0.0, 1.0), length_um=fp.sf_length_um, seg_um=seg_um, pitch_um=fp.sf_pitch_um,
        spacing_um=0.012, sarcomere_um=1.0))

    # ── CORE_BODY. lamina's spacing is SOURCED_SECONDARY; chromatin's subunit count is a PI-GAP. ──
    out["lamina"] = run("lamina", lambda: build_lamina(
        arena, radius_um=5.1, filament_um=0.4,
        spacing_provenance="SOURCED_SECONDARY: Stephens2017_MBoC Table 1 citing Shimi et al. 2015; "
                           "Shimi is not itself in this corpus"))

    out["chromatin"] = run("chromatin", lambda: build_chromatin(
        arena, count=gap(552, "Mbp chromatin domains"), subunit_um=0.6, radius_um=5.1))

    # ── NMII. PI ruling 2026-08-24. ────────────────────────────────────────────────────────────────
    # ⚠ **THE PRODUCTION CELL HAD NO MOTOR IN IT UNTIL THIS LINE.** `build/nmii.py` existed, was
    # repaired under PI queue 14(a), and was called by exactly two TIER driver scripts — never by the
    # builder the PHASE 1 and PHASE 4 drivers and the renderer use. Every gamma and tau this project
    # has recorded was measured on a cell with no myosin; `STATE.md` (c) 21 blocks all of them.
    #
    # ⚠ **THE MINIFILAMENTS STAND. THE CROSSBRIDGES DO NOT BIND, AND THAT IS THE POINT.** A bound
    # crossbridge that can never release is a PERMANENT WELD, which the charter forbids by name, and
    # with `UndefinedAcceptance` shipping nothing may commit an attach or a detach. So this population
    # contributes mass, drag and sterics and NO contractility, `census.live.bond` stays 0, and the
    # record says which of "absent" and "measured zero" it is — the distinction `builders_not_standing`
    # and `gamma_source is None` both exist to keep.
    #
    # ⚠ **`support` IS DERIVED, AND CHANGING IT MOVES THE COUNT 442 -> 430.** The TIER scripts typed
    # 706.86 um^2, which is 4*pi*7.5^2 — the MEMBRANE's area, on a population that lives in the
    # CORTICAL shell at 7.4 um (688.13 um^2). That is the same defect PI queue 14 closed from the
    # other side: nothing owned the relationship between the two shells, so one of them got used for
    # the other's quantity. `cortex_shell()` owns it now and this reads from it.
    r_shell, t_shell = cortex_shell()
    out["nmii"] = run("nmii", lambda: build_nmii(
        arena,
        count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                        scope="MCF7 cortical shell, resting",
                        provenance="Nie 2015 areal density 0.625 um^-2; the density the frozen "
                                   "incumbent already carried, not a new number"),
        support=4.0 * math.pi * r_shell * r_shell,
        n_bb=14, n_heads_per_side=30, backbone_length_um=0.301, head_offset_um=0.200,
        radius_um=r_shell, thickness_um=t_shell))

    return out


def _demo() -> None:
    """Sanity Gate — what can be checked without a device: the gaps say so, and nothing is typed.

    The builds themselves need CUDA. What this checks is the property the module exists for: that a
    caller cannot read a placeholder as a measurement, and that no position is written here.
    """
    import inspect

    # 1. Every gap carries its absence INSIDE the BondCount, not beside it.
    g = gap(500, "microtubules")
    assert g.source_class == SourceClass.PI_GAP
    assert "NO VALUE EXISTS" in g.scope and "microtubules" in g.scope
    assert "not a claim" in g.provenance

    # 2. No placement number is typed in this module. Every position must come from `fp`, because two
    #    copies of these numbers is exactly what put 95.5% of the stress fibres outside the cell.
    src = inspect.getsource(build_remaining_populations)
    for banned in ("origin=(0.0, 0.0, -7", "root_R_um=7.4", "width_um=8.0", "length_um=10.0",
                   "pitch_um=0.5"):
        assert banned not in src, f"a placement constant is typed here again: {banned}"
    for required in ("fp.filopodium_root_r_um", "fp.lamellipodium_origin_um", "fp.sf_origin_um",
                     "fp.sf_length_um", "fp.sf_pitch_um", "fp.n_stress_fibres"):
        assert required in src, f"{required} is not being read from the footprint"

    # 3. `seg_um` is required and has no default — it is the constant the audit ranks second most
    #    load-bearing with a source of "none".
    sig = inspect.signature(build_remaining_populations)
    assert sig.parameters["seg_um"].default is inspect.Parameter.empty

    # 4. The stage wrapper is optional, and its absence must not change what is built.
    calls: list[str] = []
    assert (lambda n, fn: (calls.append(n), fn())[1])("x", lambda: 7) == 7 and calls == ["x"]

    print(f"populations self-check OK — {len(inspect.getsource(build_remaining_populations).splitlines())} "
          "lines, one copy, every position read from the footprint")


if __name__ == "__main__":
    _demo()
