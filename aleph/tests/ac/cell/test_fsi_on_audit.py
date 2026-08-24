"""Oracle A — FSI-on audit: two-way poroelastic coupling is WIRED and at the physiological setpoint.

Rationale (BLEBBING_VALIDATION_PLAN_2026-07-20.md §3). The physiological-baseline HARD rule warns that an
`additive/default-off` module silently becomes de-facto production while running passive physics. The Biot
FSI back-coupling was historically stubbed (`div_vs` left ZEROS -> a deforming skeleton did NOT drive the
pore fluid). This audit is the standing regression that the two-way coupling is actually invoked by the
production driver AND that its constants sit at the sourced poroelastic operating point — NOT a physics run
(the div_vs sign, conservation, and Green's-function response are covered by the I1a fluid oracles).

Everything here is pure-Python / AST + constant reads, so it runs in the default (CPU) pytest on the dev Mac;
no Warp launch, no CUDA. What it certifies:

  * A1  solid -> fluid is wired: the driver constructs ``SolidDilatationCoupling`` and fills ``div_vs`` from
        the live solid displacement every accepted step (not the ZEROS stub).
  * A2  fluid -> solid is wired: the driver adds the ``-alpha grad p`` body force via ``PressureCoupling``.
  * A3  the Biot source consumes ``div_vs``: the scheduler feeds it into the p/mass update.
  * A4  physiological setpoint: M, c_v, alpha, Pi_0 equal the c_v-anchored / PI-ratified values.
  * A5  no yaml drift: the runtime constants equal the (unit-converted) ``params_i0b1.yaml`` ledger for the
        reconciled parameters, and the KNOWN open L_p discrepancy (PI action item 1) is tracked so it cannot
        silently change on either side without tripping this test.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from aleph.components.incumbent import assemble

ROOT = Path(__file__).resolve().parents[4]
CELL = ROOT / "aleph" / "components" / "incumbent"
FLUID = ROOT / "aleph" / "components" / "fluid"
PARAMS_YAML = FLUID / "params_i0b1.yaml"

# 1 Pa == 1 pN/um^2 exactly (engine convention), and 1 m^2 == 1e12 um^2, 1 m == 1e6 um.
_M2_TO_UM2 = 1.0e12          # mobility k/mu:  m^2/(Pa*s) -> um^4/(pN*s)  (Pa->pN/um^2 cancels, m^2->um^2)
_M_TO_UM = 1.0e6             # L_p:            m/(s*Pa)   -> um/(s*Pa)    (Pa->pN/um^2 cancels)
# storage S = 1/M and pressure are numerically invariant under Pa <-> pN/um^2 (factor 1).


def _params():
    yaml = pytest.importorskip("yaml")
    doc = yaml.safe_load(PARAMS_YAML.read_text(encoding="utf-8"))
    return doc["parameters"]


# ── A1/A2/A3 — the two-way coupling is invoked by the production runtime (not stubbed) ────────────────

def test_solid_to_fluid_dilatation_is_wired_into_the_driver() -> None:
    """A1: the driver builds SolidDilatationCoupling and fills div_vs from live displacement each step."""
    driver = CELL / "driver.py"
    text = driver.read_text(encoding="utf-8")
    assert "from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling" in text
    assert "SolidDilatationCoupling(cell.grid)" in text, "solid->fluid coupling must be constructed"
    # the dilatation source must actually be refreshed from the moving skeleton, not left at ZEROS
    assert "solid_coupling.update_from_displacement(" in text, "div_vs must be filled (anti-stub)"


def test_fluid_to_solid_pressure_body_force_is_wired_into_the_driver() -> None:
    """A2: the driver adds the -alpha grad p fluid->solid body force via PressureCoupling."""
    text = (CELL / "driver.py").read_text(encoding="utf-8")
    assert "cell.pressure.accumulate(" in text, "fluid->solid -alpha p I coupling must be accumulated"


def test_biot_update_consumes_div_vs() -> None:
    """A3: the scheduler feeds the solid dilatation into the p/mass update (the source is live)."""
    text = (FLUID / "scheduler.py").read_text(encoding="utf-8")
    # div_vs is passed into the p/mass substep launch alongside the Biot-Willis alpha: the -alpha*div(v_s)
    # source term. Both present => the dilatation source is genuinely consumed, not ignored.
    assert "g.div_vs" in text, "the Biot p/mass update must receive the solid dilatation source"
    assert "sub.alpha" in text, "the -alpha*div(v_s) coupling term must enter the p/mass update"


# ── A4 — physiological operating point (the c_v-anchored / PI-ratified closure) ───────────────────────

def test_biot_constants_are_at_the_physiological_setpoint() -> None:
    """A4: M, c_v, alpha, Pi_0 are the sourced poroelastic operating point, not water/zero defaults."""
    # storage S = 1/M -> M in Pa (== pN/um^2). c_v = mobility / S.
    M = 1.0 / assemble.BIOT_STORAGE_S
    c_v = assemble.BIOT_MOBILITY / assemble.BIOT_STORAGE_S
    assert M == pytest.approx(1.0e4, rel=1e-9), "M must be the c_v-anchored ~10 kPa, not the stale 1 kPa"
    assert c_v == pytest.approx(50.0, rel=1e-9), "mobility/S must equal the Moeendarbary c_v anchor (40-60)"
    assert 40.0 <= c_v <= 60.0
    assert assemble.BIOT_ALPHA == 1.0, "Biot-Willis alpha is PI-ratified 1.0 (incompressible-constituent)"
    assert assemble.PI_0_PA == pytest.approx(40.0), "resting turgor must be the 40 Pa baseline, never zero"


# ── A5 — no drift between the runtime constants and the params ledger; L_p gap tracked ────────────────

def test_reconciled_constants_match_the_params_ledger() -> None:
    """A5: mobility/S/alpha/dP_hyd equal the unit-converted params_i0b1.yaml values (drift guard).

    The resting-turgor ledger key was renamed ``Pi_0`` -> ``delta_p_hydrostatic_rest`` (symbol dP_hyd) in the
    2026-07-22 provenance correction (it is the net hydrostatic ΔP the cortex bears, NOT an osmotic pressure);
    the 40 Pa value and the ``assemble.PI_0_PA`` constant are unchanged. Renaming ``PI_0_PA`` -> ``DP_HYD_PA``
    in code + splitting the driver ``osmotic_difference`` remains an open TODO (see params_i0b1.yaml).
    """
    p = _params()
    mobility_engine = float(p["mobility"]["value"]) * _M2_TO_UM2      # 5e-15 m^2/(Pa*s) -> 5e-3
    storage_engine = float(p["M_biot"]["storativity_S"])             # /Pa == um^2/pN (factor 1)
    assert assemble.BIOT_MOBILITY == pytest.approx(mobility_engine, rel=1e-6)
    assert assemble.BIOT_STORAGE_S == pytest.approx(storage_engine, rel=1e-6)
    assert assemble.BIOT_ALPHA == pytest.approx(float(p["alpha_biot"]["value"]))
    assert assemble.PI_0_PA == pytest.approx(float(p["delta_p_hydrostatic_rest"]["value"]))


def test_Lp_open_discrepancy_is_tracked_not_silent() -> None:
    """A5(L_p): the code L_p and the Jung-2011 ledger value differ by ~62x (PI action item 1).

    This is a KNOWN-OPEN discrepancy, not a passing reconciliation: the code value is anchored to the
    gamma_floor, the ledger to the Jung-2011 MCF7/AQP5 draft. Encoding the ratio as an invariant means the
    test fails loudly the moment either side changes, forcing a PI-ledger update rather than a silent fix.
    """
    p = _params()
    lp_ledger_engine = float(p["L_p"]["value"]) * _M_TO_UM             # 1e-12 m/(s*Pa) -> 1e-6 um/(s*Pa)
    ratio = lp_ledger_engine / assemble.L_P
    assert ratio == pytest.approx(62.5, rel=0.05), (
        "L_p code-vs-ledger ratio changed; resolve PI action item 1 and update the ledger + this invariant"
    )
