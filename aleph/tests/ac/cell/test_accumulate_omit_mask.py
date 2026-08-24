"""``_accumulate_all(omit=)`` — the parameter four validation tracks were blocked behind.

WHY IT EXISTS.  Until 2026-07-29 there was no way to hand ONE component to the engine while the rest of
the cell kept the incumbent force-assembly path.  ``ac/engine/cortex_state.py`` gives the cortex its own
``position_d``/``force_d`` and was verified at full native (node-by-node parity 2.80e-16), but binding it
into a production run would have added the cortex force TWICE — once on the component's arrays, once by
``_accumulate_all`` over the global ones.  The only pre-existing "omit" was build-time
(``--no-myosin``/``--no-nucleus``/``--no-membrane``), which REMOVES the physics rather than relocating
it.  PI decision **`COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md` §6 D2** (approved 2026-07-29, option (a):
additive, bit-identical default, per component).  Name the document: `AC_EXECUTION_PLAN_2026-07-25.md` §9 has a
DIFFERENT D2 (membrane area-stiffness, already landed), and a bare "D2" cannot be resolved to either.

WHAT THESE TESTS PIN, and why each is the failure that would otherwise be invisible:

* **The default changes nothing.** A mask parameter whose default perturbed the assembly would
  invalidate every result measured before it existed, and the perturbation would be a reordering nobody
  would think to look for.
* **A typo raises.** ``omit={"steric"}`` instead of ``"steric_wca"`` must not be read as "omit nothing":
  that reads as success at the call site and doubles the force the caller believes it moved.
* **Omitting is per-channel and additive** — omitting one leaves every other launch exactly as it was.
* **The vocabulary is the manifest's.** Two spellings of a channel name would make the typo guard
  unenforceable, so the omittable set is checked against
  ``ac/engine/forces_manifest`` rather than maintained beside it.

No device is involved at all: `wp.launch` is replaced by a recorder and the cell's ``device`` attribute
is a sentinel that is never resolved, because the CONTRACT under test is WHICH CHANNELS FIRE — not what
they compute, which is the native gate's job.
"""

from __future__ import annotations

import pytest

from aleph.components.incumbent.driver import OMITTABLE_CHANNELS


class _Recorder:
    """Stands in for the composed sub-objects, recording that its ``accumulate`` was reached."""

    def __init__(self, log: list[str], name: str) -> None:
        self._log = log
        self._name = name

    def accumulate(self, *args: object) -> None:
        self._log.append(self._name)


class _FakeCell:
    """The minimum ``AssembledCell`` surface ``_accumulate_all`` touches, with every channel present."""

    #: NOT a Warp device. `_accumulate_all` forwards this to `wp.launch(device=...)`, which these tests
    #: replace with a recorder, so nothing ever resolves it.  The I0-A static contract forbids a
    #: CPU-device literal anywhere under `tests/ac/` regardless of intent — correctly, because neither a
    #: reader nor a grep can tell a passed-through sentinel from a real CPU launch.  It caught this file
    #: twice: once for the attribute, and once for a comment that quoted the forbidden literal while
    #: explaining why it is forbidden.
    _NOT_A_DEVICE = "<recording-launcher, never resolved>"

    def __init__(self, log: list[str]) -> None:
        self.device = self._NOT_A_DEVICE
        self.n_total = 8
        self.n_tri = 3
        self.n_xl = 4
        self.n_branch = 2
        self.tri_d = self.alpha_d = self.xl_d = self.kxl_d = self.r0xl_d = object()
        self.branch_triples_d = self.branch_active_d = object()
        self.state = object()
        self.myosin = _Recorder(log, "nmii_minifilament_force")
        self.nucleus = _Recorder(log, "nucleus_compartment_force")
        self.membrane = _Recorder(log, "membrane_compartment_force")
        self.membrane_pressure = _Recorder(log, "membrane_pressure_traction")
        self.steric = _Recorder(log, "steric_wca")
        self.pressure = _Recorder(log, "biot_pressure_coupling")


def _launched(monkeypatch: pytest.MonkeyPatch, **kwargs: object) -> list[str]:
    """Run ``_accumulate_all`` against recording doubles and return the channels that fired, in order."""
    from aleph.components.incumbent import driver as drv

    log: list[str] = []
    kernel_names = {
        id(drv.cytosim_bending_kernel): "actin_bending_cytosim",
        id(drv.link_spring_kernel): "actin_crosslink_link_spring",
        id(drv.branch_angle_kernel): "arp23_branch_angle",
    }

    def fake_launch(kernel: object, **_: object) -> None:
        name = kernel_names.get(id(kernel))
        if name is not None:
            log.append(name)

    monkeypatch.setattr(drv.wp, "launch", fake_launch)
    drv._accumulate_all(_FakeCell(log), object(), object(), **kwargs)  # type: ignore[arg-type]
    return log


ALL_CHANNELS = [
    "actin_bending_cytosim",
    "actin_crosslink_link_spring",
    "arp23_branch_angle",
    "nmii_minifilament_force",
    "nucleus_compartment_force",
    "membrane_compartment_force",
    "membrane_pressure_traction",
    "steric_wca",
    "biot_pressure_coupling",
]


def test_the_default_launches_every_channel_in_the_incumbent_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bit-identical default: omitting nothing must reach every launch, in the order it always did.

    Order is asserted, not just membership — the sum is order-independent but float64 atomic ARRIVAL is
    not, so a reordering would perturb the last bits of every historical comparison.
    """
    assert _launched(monkeypatch) == ALL_CHANNELS


def test_an_empty_mask_is_the_same_as_no_mask(monkeypatch: pytest.MonkeyPatch) -> None:
    """Passing the parameter explicitly must not differ from not passing it."""
    assert _launched(monkeypatch, omit=frozenset()) == _launched(monkeypatch)


@pytest.mark.parametrize("channel", ALL_CHANNELS)
def test_omitting_one_channel_removes_exactly_that_channel(
    monkeypatch: pytest.MonkeyPatch, channel: str
) -> None:
    """Additive and per-channel: every other launch survives untouched."""
    fired = _launched(monkeypatch, omit={channel})
    assert channel not in fired
    assert fired == [c for c in ALL_CHANNELS if c != channel]


def test_omitting_the_cortex_pair_leaves_the_rest_of_the_cell(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The T2 case: hand the cortex's own channels to the engine, keep every other compartment here."""
    fired = _launched(monkeypatch, omit={"actin_bending_cytosim", "actin_crosslink_link_spring"})
    assert "actin_bending_cytosim" not in fired and "actin_crosslink_link_spring" not in fired
    assert "nmii_minifilament_force" in fired and "membrane_compartment_force" in fired


def test_an_unknown_channel_raises_rather_than_omitting_nothing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A typo must not read as success — that is the exact path to a silently doubled force.

    ``"steric"`` for ``"steric_wca"`` is the realistic slip: plausible, close, and under a permissive
    implementation it would omit nothing while the caller's component adds the same force again.
    """
    with pytest.raises(ValueError, match="unknown force channel"):
        _launched(monkeypatch, omit={"steric"})
    with pytest.raises(ValueError, match="omit nothing"):
        _launched(monkeypatch, omit={"cortex"})


def test_the_omittable_set_is_the_manifest_vocabulary() -> None:
    """One spelling of each channel, or the typo guard above guards nothing.

    The manifest is the source: if a channel is added there and not here it cannot be handed to a
    component, and if it is spelled differently here the guard passes a name the manifest never knew.
    """
    from aleph.engine.forces_manifest import _FORCE_TERMS

    manifest_names = {term[0] for term in _FORCE_TERMS}
    assert OMITTABLE_CHANNELS <= manifest_names, (
        f"omittable names absent from the manifest: {sorted(OMITTABLE_CHANNELS - manifest_names)}"
    )
    assert set(ALL_CHANNELS) == OMITTABLE_CHANNELS
