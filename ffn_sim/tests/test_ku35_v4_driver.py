"""KU-3.5 v4 FA-anchored driver — light CI tests.

Covers (CI-fast):
  * import-clean (the driver module + its v3 measurement-protocol re-imports
    resolve without error);
  * argparse: the CLI parses, including the ``--smoke`` preset and the
    ``--fa-capture-radius`` override;
  * a tiny in-process build of the v4 cell (cortex + myosin + xlinks + FA +
    enclosed-volume + turnover) yields ``n_fa_integrins > 0`` — i.e. the FA
    layer is actually wired in by the driver's resolver stack.

The full simulation (warm-up + equilibrate + production sampling) is NOT run
here — it is the driver's SMOKE entry-point, exercised opt-in via
``RUN_KU35_V4_SIM=1``. That keeps this module CI-fast.
"""
from __future__ import annotations

import os
import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
import hoomd

PKG = Path(__file__).resolve().parents[1]
H3_CFG = PKG / "configs" / "phase1_h3.yaml"
H4_CFG = PKG / "configs" / "phase1_h4.yaml"


def test_driver_imports_clean():
    """The v4 driver module imports without error (incl. v3 re-imports)."""
    import ffn_sim.scripts.h3_ku35_v4_fa as drv

    assert hasattr(drv, "run")
    assert hasattr(drv, "main")
    assert hasattr(drv, "_run_fingerprint")
    # v3 measurement protocol must have been re-imported (not re-implemented).
    assert callable(drv._tension_method_of_planes)
    assert callable(drv._tension_method_of_planes_rigid)


def test_cli_parses_smoke_and_override():
    """argparse accepts the smoke preset + fa-capture-radius override.

    Drive ``main``'s parser by monkeypatching argv, but stop before the
    (expensive) ``run`` call by patching it to a recorder.
    """
    import ffn_sim.scripts.h3_ku35_v4_fa as drv

    captured = {}

    def _fake_run(n_fil, seed, **kw):
        captured["n_fil"] = n_fil
        captured["seed"] = seed
        captured.update(kw)
        return {}

    orig_run = drv.run
    orig_argv = sys.argv
    try:
        drv.run = _fake_run  # type: ignore[assignment]
        sys.argv = [
            "h3_ku35_v4_fa",
            "--smoke",
            "--fa-capture-radius", "3e-6",
            "--out", "/tmp/_ku35_v4_test.json",
        ]
        drv.main()
    finally:
        drv.run = orig_run  # type: ignore[assignment]
        sys.argv = orig_argv

    # Smoke preset must have shrunk the defaults.
    assert captured["n_fil"] == 60
    assert captured["n_sample"] == 3
    assert captured["interval"] == 500
    assert captured["equilibrate_steps"] == 500
    assert captured["n_warmup"] == 2_000
    assert captured["fa_capture_radius"] == pytest.approx(3e-6)


def test_fingerprint_changes_with_capture_radius():
    """The resume fingerprint must distinguish a capture-radius override."""
    import ffn_sim.scripts.h3_ku35_v4_fa as drv

    base = dict(
        n_fil=60, seed=1, dt_factor=0.001, equilibrate_steps=500,
        n_warmup=2000, n_sample=3, interval=500,
    )
    fp_phys = drv._run_fingerprint(fa_capture_radius=None, **base)
    fp_over = drv._run_fingerprint(fa_capture_radius=3e-6, **base)
    assert fp_phys != fp_over


def test_v4_cell_builds_with_fa_integrins():
    """A tiny in-process v4 build wires the FA layer (n_fa_integrins > 0).

    Mirrors the driver's resolver stack but builds the cheaper unconstrained
    sim (no equilibration, no production) so it stays CI-fast: the assertion
    is purely that the FA integrins are present in the constructed state.
    """
    from ffn_sim.cortex.cortex import resolve_h3_derived
    from ffn_sim.cortex.myosin import resolve_cortex_myosin
    from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
    from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume
    from ffn_sim.cortex.turnover import resolve_turnover
    from ffn_sim.bridge.fa import resolve_h4
    from ffn_sim.cell.cell import build_cortex_full_simulation

    cfg = yaml.safe_load(open(H3_CFG))
    cfg["cortex"]["n_filaments"] = 60
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    nca = p.n_filaments * p.beads_per_filament
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend

    p_myo = resolve_cortex_myosin(cfg, dt=dtc)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    p_ev = resolve_enclosed_volume(cfg, R_cell=p.R_cell)
    p_to = resolve_turnover(cfg, dt=dtc, rest_length=p.rest_length)
    p_fa = resolve_h4(yaml.safe_load(open(H4_CFG)))

    h = build_cortex_full_simulation(
        p,
        p_xlinks=p_xl,
        p_myosin=p_myo,
        p_fa=p_fa,
        p_enclosed_volume=p_ev,
        p_turnover=p_to,
        device=hoomd.device.CPU(notice_level=0),
        constrained=False,
        rng=np.random.default_rng(0),
    )
    assert int(h["n_fa_integrins"]) > 0
    assert int(h["n_fa_clutch_bonds"]) >= 0  # may be 0 pre-equilibration


@pytest.mark.skipif(
    os.environ.get("RUN_KU35_V4_SIM") != "1",
    reason="opt-in full smoke sim (set RUN_KU35_V4_SIM=1)",
)
def test_v4_smoke_runs():
    """Opt-in: the actual smoke sim builds + runs without NaN/crash."""
    import ffn_sim.scripts.h3_ku35_v4_fa as drv

    res = drv.run(
        60, 1, dt_factor=0.001, n_warmup=500, n_sample=2, interval=300,
        equilibrate_steps=300, fa_capture_radius=3e-6, device="cpu",
        out=None, skip_integrity=True, resume=False,
    )
    assert res["n_fa_integrins"] > 0
    assert np.isfinite(res["tension_plateau_mN_per_m"])
