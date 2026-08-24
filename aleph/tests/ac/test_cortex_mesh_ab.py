"""CPU smoke for the cortex fine-vs-coarse mesh A/B wiring (no CUDA — arg-parse + CellConfig threading only).

Guards the additive change readying the mesh-physics verification:
  * ``CellConfig`` accepts and stores ``cortex_seg_um`` / ``cortex_density_per_fil`` / ``cortex_length_um`` so
    ``--cortex-seg-um 0.075 --cortex-density 40`` builds the fine cortex, and the defaults reproduce today's
    committed baseline bit-identically;
  * the ``ac_cortex_mesh_ab`` orchestrator threads those flags into the two existing native drivers' commands;
  * both native drivers expose the new flags via ``--help`` (which returns BEFORE any CUDA init).

No ``build_cell`` / CUDA here (I0-A — the native runs are gbook-only); this only exercises the host wiring.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from aleph.components.incumbent.assemble import CellConfig
from aleph.scripts import ac_cortex_mesh_ab as ab

REPO_ROOT = Path(__file__).resolve().parents[3]
FFN_SIM = REPO_ROOT / "aleph"


def test_cellconfig_defaults_are_todays_baseline() -> None:
    """Bit-identical guard: the default cortex config is the committed 0.5 / 20 / 3.0 baseline."""
    c = CellConfig(n_filaments=70686)
    assert c.cortex_seg_um == 0.5
    assert c.cortex_density_per_fil == 20.0
    assert c.cortex_length_um == 3.0


def test_cellconfig_threads_fine_mesh() -> None:
    """The fine derived config (75 nm mesh / density 40) is accepted and stored on CellConfig."""
    c = CellConfig(n_filaments=70686, cortex_seg_um=0.075, cortex_density_per_fil=40, cortex_length_um=3.0)
    assert c.cortex_seg_um == 0.075
    assert c.cortex_density_per_fil == 40
    assert c.cortex_length_um == 3.0


def test_node_count_helper_matches_geometry() -> None:
    """nodes/filament = round(L/ℓ₀)+1 → coarse 494,802, fine 2,898,126 (geometry-verified counts)."""
    assert ab._nodes(0.5, 3.0, 70686) == 494_802
    assert ab._nodes(0.075, 3.0, 70686) == 2_898_126


def _args(**over) -> argparse.Namespace:
    base = dict(native=True, outer=40, seed=0, myosin_fraction=None, membrane_subdiv=None,
                steps=40, dt=0.01, catch_slip=False, filaments=70686, python=sys.executable)
    base.update(over)
    return argparse.Namespace(**base)


def test_gate_a_command_threads_cortex_flags() -> None:
    coarse = ab._gate_a_cmd(sys.executable, ab.COARSE, _args())
    fine = ab._gate_a_cmd(sys.executable, ab.FINE, _args())
    assert "--cortex-seg-um" in coarse and "0.5" in coarse and "20.0" in coarse
    assert "--cortex-seg-um" in fine and "0.075" in fine and "40.0" in fine
    assert "--native" in fine


def test_gate_b_command_threads_cortex_flags() -> None:
    coarse = ab._gate_b_cmd(sys.executable, ab.COARSE, _args())
    fine = ab._gate_b_cmd(sys.executable, ab.FINE, _args())
    assert "--cortex-seg-um" in coarse and "0.5" in coarse
    assert "0.075" in fine and "40.0" in fine
    assert "--steps" in fine and "--catch-slip" not in fine


def _help(script: str) -> str:
    import os

    env = dict(os.environ)
    env["PYTHONPATH"] = f"{REPO_ROOT}:{FFN_SIM}"
    proc = subprocess.run(
        [sys.executable, str(FFN_SIM / "scripts" / script), "--help"],
        capture_output=True, text=True, env=env,
    )
    return proc.stdout + proc.stderr


def test_gate_a_help_exposes_flags() -> None:
    out = _help("ac_gate_a_fq_coarse_test.py")
    assert "--cortex-seg-um" in out and "--cortex-density" in out and "--cortex-length-um" in out


def test_gate_b_help_exposes_flags() -> None:
    out = _help("ac_gate_b_cortex_motor_native.py")
    assert "--cortex-seg-um" in out and "--cortex-density" in out and "--cortex-length-um" in out
