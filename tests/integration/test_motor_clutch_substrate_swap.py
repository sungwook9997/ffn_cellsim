"""Cross-unit integration: motor-clutch ↔ substrate adapter swap pattern.

Phase 1 ships a single substrate model
(``acs_kb.bridge.substrate_stub.LinearElasticSubstrate``). Unit 2.2+
will introduce ``acs_kb.bridge.ecm_adapter`` that projects an ECM
fibre network into a local effective stiffness. The motor-clutch
module reads only the substrate's ``stiffness`` and
``compute_displacement``; this test pins that contract so the future
swap is mechanical, not philosophical.

Phase 1 self-consistency: two ``LinearElasticSubstrate`` instances
with the same KU-1.21 parameters must produce identical
``stiffness`` and ``compute_displacement`` outputs (±10 %
tolerance, which is the long-term swap target for the ECM adapter).

KU tags: KU-1.21 (linear-elastic substrate), KU-2.4 (motor-clutch),
         KU-2.12 (traction).
"""

from __future__ import annotations

import importlib
import inspect
from typing import Protocol, runtime_checkable

import numpy as np
import pytest


def _try_import(modpath: str):
    try:
        return importlib.import_module(modpath)
    except Exception as e:  # noqa: BLE001
        return e


substrate_stub = _try_import("acs_kb.bridge.substrate_stub")
motor_clutch = _try_import("acs_kb.bridge.motor_clutch")
bridge_types = _try_import("acs_kb.bridge.types")

skip_if_bridge_missing = pytest.mark.skipif(
    any(isinstance(m, Exception) for m in (substrate_stub, motor_clutch, bridge_types)),
    reason=(
        f"Worker B bridge stack not importable: "
        f"substrate_stub={substrate_stub!r}, motor_clutch={motor_clutch!r}, "
        f"types={bridge_types!r}"
    ),
)


@runtime_checkable
class _SubstrateAdapterProtocol(Protocol):
    """The week-5 frozen substrate contract.

    Any future ECM adapter must implement these two members. The test
    suite checks the LinearElasticSubstrate against this Protocol so
    that a silent rename / signature drift fails loudly.
    """

    @property
    def stiffness(self) -> float: ...
    def compute_displacement(self, position: np.ndarray, force):  # noqa: ANN001
        ...


@skip_if_bridge_missing
def test_substrate_stub_conforms_to_swap_protocol():
    """KU-1.21: LinearElasticSubstrate satisfies the swap contract."""
    sub = substrate_stub.LinearElasticSubstrate()
    assert isinstance(sub, _SubstrateAdapterProtocol), (
        "LinearElasticSubstrate fails the swap protocol; the future "
        "ECM adapter must conform to the same two-member API."
    )


@skip_if_bridge_missing
def test_two_stub_instances_agree_within_tolerance():
    """KU-1.21 self-consistency: identical params ⇒ identical outputs."""
    a = substrate_stub.LinearElasticSubstrate(
        young_modulus=5.0e3, poisson_ratio=0.45, contact_radius=1.0e-6,
    )
    b = substrate_stub.LinearElasticSubstrate(
        young_modulus=5.0e3, poisson_ratio=0.45, contact_radius=1.0e-6,
    )
    # stiffness identity
    assert a.stiffness == pytest.approx(b.stiffness, rel=1e-12)
    # compute_displacement identity over a force sweep
    pos = np.array([0.0, 0.0])
    for F in (0.0, 1e-12, 1e-10, 5e-9, 1e-8):
        u_a = float(a.compute_displacement(pos, F))
        u_b = float(b.compute_displacement(pos, F))
        # 10 % swap tolerance — the long-term target for stub ↔ ECM
        if u_a == 0.0:
            assert u_b == 0.0
        else:
            assert abs(u_a - u_b) / abs(u_a) <= 0.10


@skip_if_bridge_missing
def test_motor_clutch_only_uses_swap_contract():
    """KU-2.4: MotorClutchFA must not access private substrate fields.

    Inspect MotorClutchFA's source for attribute access into ``substrate``;
    only ``.stiffness`` and ``.compute_displacement`` are allowed (plus
    Python dunders). A failure here means a future ECM adapter would
    have to expose more than the swap contract promises.
    """
    src = inspect.getsource(motor_clutch.MotorClutchFA)
    # Tokenise the source looking for `self.substrate.<attr>` references.
    import re
    attrs_used = set(re.findall(r"\.substrate\.([A-Za-z_]\w*)", src))
    # Allow the two contracted names plus any dunder.
    contract = {"stiffness", "compute_displacement"}
    forbidden = {a for a in attrs_used if not a.startswith("__") and a not in contract}
    assert not forbidden, (
        f"MotorClutchFA reaches into substrate attributes outside the "
        f"swap contract: {sorted(forbidden)}. Extend the swap Protocol "
        f"before merging."
    )


@skip_if_bridge_missing
def test_motor_clutch_step_with_two_identical_substrates_agrees():
    """KU-2.4 + KU-1.21: identical substrate stubs ⇒ identical traction
    trajectory, within 10 % over a 200-step run."""
    rng_a = np.random.default_rng(7)
    rng_b = np.random.default_rng(7)
    fa_a = bridge_types.make_focal_adhesion(position=np.array([0.0, 0.0]))
    fa_b = bridge_types.make_focal_adhesion(position=np.array([0.0, 0.0]))
    sub_a = substrate_stub.LinearElasticSubstrate()
    sub_b = substrate_stub.LinearElasticSubstrate()
    mc_a = motor_clutch.MotorClutchFA(fa_a, sub_a)
    mc_b = motor_clutch.MotorClutchFA(fa_b, sub_b)
    dt = 1.0e-3
    n = 200
    Fs_a, Fs_b = [], []
    for _ in range(n):
        d_a = mc_a.step(dt, rng_a)
        d_b = mc_b.step(dt, rng_b)
        Fs_a.append(float(d_a["force_total"]))
        Fs_b.append(float(d_b["force_total"]))
    # Identical RNG seeds + identical stubs ⇒ trajectories must be
    # bit-equal. The 10 % tolerance is the long-term swap target;
    # for self-consistency we assert exact equality.
    np.testing.assert_array_equal(Fs_a, Fs_b)
