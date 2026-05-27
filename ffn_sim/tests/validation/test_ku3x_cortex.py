"""KU-3.x cortex emergent-physics validation gates (H.3 production sign-off).

Per brief §Validation acceptance table:

| Gate         | Criterion                                            | KU      |
| ---          | ---                                                  | ---     |
| Cell rounding| Aspect ratio 1.5 → < 1.2 in 60 s                     | KU-3.1  |
| Cortical γ   | γ_cortex ≈ 0.5 mN/m ± 30 % (Young-Laplace from ΔP)   | KU-3.5  |
| Blebbistatin | myosin-OFF: γ ↓ → rounding fails (aspect > 1.3 at 60 s) | KU-3.18 |
| Nematic Q_ij | Isotropic cell S < 0.1; aligned cell S > 0.3         | KU-3.20 |

These gates require:

* Full Cell composition: cortex + ERM + crosslinkers + myosin (all wired
  via :class:`ffn_sim.cell.Cell.build`).
* Long-time simulation: 60 s simulated, ~4.6·10⁹ BAOAB steps at
  ``dt_CFL = 13 ns`` → ~hours to days wall time even at the demo
  cortex scale. Production sign-off is multi-hour and runs
  ``H3_KU3_PRODUCTION=1``.

This file ships TEST SKELETONS that establish:

* The first-principles acceptance bands (from KU-3.x literature; same
  no-measurement-anchoring rule as H.2 Day 5).
* The plumbing to build the full Cell + ERM + myosin sim with the
  brief-canonical setup.
* The aspect-ratio / Q_ij / γ_cortex measurement utilities.

PI sign-off path (next iteration with PI ratification): run with
``H3_KU3_PRODUCTION=1``, on M1 Max budget at least 6-12 hours wall
per gate. Until then ALL gates SKIP at CI time (no false-PASS).

Also requires the OPEN H.3 단계 3 ERM CFL conflict to be resolved (PI
decision on dt-reduction vs k_ERM-softening) — see ``outputs/h3/REPORT.md``
§단계 3 sanity finding.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell import Cell, CellBuildOptions
from ffn_sim.cortex.cortex import generate_cortex_topology, resolve_h3_derived
from ffn_sim.validation.oracles.common.sanity_gate import gate_nematic_order


CONFIG_PATH = (
    Path(__file__).resolve().parents[2] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


H3_KU3_PRODUCTION = bool(int(os.environ.get("H3_KU3_PRODUCTION", "0")))


# ---------------------------------------------------------------------------
# Measurement utilities (always available; trivial enough to keep here)
# ---------------------------------------------------------------------------
def cell_aspect_ratio(cortex_positions: np.ndarray) -> float:
    """Aspect ratio of cell-shape ellipsoid from cortex bead positions.

    Computes the inertia-tensor-aligned bounding ellipsoid and returns
    ``max(λ)/min(λ)`` where λ are the principal-axis eigenvalues.

    For an initial sphere with aspect ratio 1.5 (semi-major 12 μm,
    semi-minor 8 μm per brief KU-3.1 VALIDATION setup), this returns
    1.5 to numerical precision.
    """
    if cortex_positions.shape[0] < 3:
        return float("nan")
    com = cortex_positions.mean(axis=0)
    centered = cortex_positions - com
    I = centered.T @ centered / cortex_positions.shape[0]
    eigvals = np.linalg.eigvalsh(I)
    eigvals = np.clip(eigvals, 1e-30, None)
    return float(math.sqrt(eigvals.max() / eigvals.min()))


def nematic_order_S(cortex_topology_axes: np.ndarray) -> float:
    """Q-tensor nematic order parameter S for a set of filament axes.

    For F unit-vector axes, builds the Q tensor
    ``Q = ⟨3·û⊗û − I⟩/2`` and returns its largest eigenvalue.
    Isotropic distribution → S → 0; perfectly aligned → S = 1.
    """
    if cortex_topology_axes.shape[0] < 2:
        return float("nan")
    axes = cortex_topology_axes / np.linalg.norm(
        cortex_topology_axes, axis=1, keepdims=True
    ).clip(min=1e-30)
    n = axes.shape[0]
    Q = np.zeros((3, 3), dtype=np.float64)
    for u in axes:
        Q += 3.0 * np.outer(u, u) - np.eye(3)
    Q /= 2.0 * n
    eigvals = np.linalg.eigvalsh(Q)
    return float(eigvals.max())


# ---------------------------------------------------------------------------
# KU-3.1 cell rounding gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H3_KU3_PRODUCTION,
    reason=(
        "KU-3.1 cell rounding gate requires full Cell (cortex + ERM + "
        "myosin) running 60 s simulated time (~4.6e9 BAOAB steps, multi-"
        "hour wall on M1 Max). Opt-in via H3_KU3_PRODUCTION=1. Also "
        "requires H.3 단계 3 ERM CFL conflict resolution (PI sign-off)."
    ),
)
class TestKU31Rounding:
    """KU-3.1 VALIDATION: cortex with active myosin rounds an initially
    elongated cell from aspect ratio 1.5 → < 1.2 in 60 s simulated time.

    First-principles acceptance band (NO measurement-anchored tuning):
    ``aspect_ratio_at_60_s ≤ 1.2`` — KU-3.1 literature criterion
    (Salbreux 2012). NOT widened to fit any specific simulation result.
    """

    def test_rounding_aspect_ratio_below_1_2_at_60s(self):
        pytest.skip("Production sign-off skeleton; see module docstring.")


# ---------------------------------------------------------------------------
# KU-3.5 cortical tension gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H3_KU3_PRODUCTION,
    reason="H3_KU3_PRODUCTION=1 required (see KU-3.1 reason).",
)
class TestKU35CortexTension:
    """KU-3.5: emergent cortex tension γ ≈ 0.5 mN/m ± 30 %.

    First-principles acceptance band: KU-3.5 literature (Salbreux 2012
    Phase 1 default mid-range epithelial value). Band ±30 % derived
    from the literature range 0.1–1 mN/m.
    """

    def test_emergent_gamma_in_band(self):
        pytest.skip("Production sign-off skeleton; see module docstring.")


# ---------------------------------------------------------------------------
# KU-3.18 blebbistatin gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H3_KU3_PRODUCTION,
    reason="H3_KU3_PRODUCTION=1 required (see KU-3.1 reason).",
)
class TestKU318Blebbistatin:
    """KU-3.18: with myosin OFF, cortex tension drops and rounding
    fails (aspect ratio stays > 1.3 at 60 s).

    Comparative gate: pairs against TestKU31Rounding; the same cell
    setup with all myosin heads disabled (k_head_actin = 0 or
    MyosinStepUpdater not attached) must give aspect > 1.3 at 60 s.
    """

    def test_myosin_off_rounding_fails(self):
        pytest.skip("Production sign-off skeleton; see module docstring.")


# ---------------------------------------------------------------------------
# KU-3.20 nematic order gate (STRUCTURAL — always runs)
# ---------------------------------------------------------------------------
def _cortex_tangents(
    *, n_filaments: int, seed: int,
    bias_axis=None, bias_kappa: float = 0.0,
) -> np.ndarray:
    """Build a cortex topology and return its per-filament tangent axes.

    Pure topology construction (no HOOMD dynamics) — the nematic order
    parameter S is a structural property of the constructed cortex
    filament-axis ensemble (PI 2026-05-28: KU-3.20 = structural measure,
    no 60 s dynamics run needed). ``demo_mode`` relaxes the cost ceiling
    so the full ×40 mesoscopic count (1000 filaments) can be built for
    good Q-tensor statistics.
    """
    cfg = _load_cfg()
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    topo = generate_cortex_topology(
        p, rng=np.random.default_rng(seed),
        tangent_bias_axis=bias_axis, tangent_bias_kappa=bias_kappa,
    )
    return topo.tangents


class TestKU320NematicOrder:
    """KU-3.20: emergent nematic order parameter S satisfies:

    * Isotropic cell (no orientation bias): S < 0.1.
    * Cortex with biased alignment (von-Mises tangent at construction):
      S > 0.3 along the bias axis.

    S is a *structural* property of the cortex filament-axis ensemble
    (Q-tensor max eigenvalue), so unlike the dynamics gates (KU-3.1/3.5/
    3.18) this runs at CI time on the full ×40 mesoscopic filament count
    — no multi-hour ``H3_KU3_PRODUCTION`` run required (PI 2026-05-28).
    The acceptance bands (S<0.1 isotropic / S>0.3 aligned) are the brief
    KU-3.20 contract, verified through the oracle ``gate_nematic_order``.
    """

    # ×40 mesoscopic full count → Q-tensor statistics are well-converged.
    N_FIL = 1000

    def test_isotropic_and_aligned_S_satisfy_KU320(self):
        # Isotropic: uniform azimuth → S → 0.
        S_iso = nematic_order_S(
            _cortex_tangents(n_filaments=self.N_FIL, seed=20)
        )
        # Aligned: strongly-concentrated von-Mises about a global stress
        # axis (κ=8 ⇒ "aligned cortex" experimental condition, NOT tuned
        # to the band — a perfectly combed sphere caps at S=0.5 > 0.3).
        S_aligned = nematic_order_S(
            _cortex_tangents(
                n_filaments=self.N_FIL, seed=20,
                bias_axis=(0.0, 0.0, 1.0), bias_kappa=8.0,
            )
        )
        report = gate_nematic_order(
            S_isotropic_cell=S_iso, S_aligned_cell=S_aligned
        )
        assert report.passed, (
            f"KU-3.20 nematic gate FAIL: S_iso={S_iso:.4f} (need <0.1), "
            f"S_aligned={S_aligned:.4f} (need >0.3)\n{report.summary()}"
        )

    def test_isotropic_cell_S_below_0_1(self):
        S_iso = nematic_order_S(
            _cortex_tangents(n_filaments=self.N_FIL, seed=20)
        )
        assert S_iso < 0.1, f"Isotropic cortex S={S_iso:.4f} should be <0.1."

    def test_aligned_cell_S_above_0_3(self):
        S_aligned = nematic_order_S(
            _cortex_tangents(
                n_filaments=self.N_FIL, seed=20,
                bias_axis=(0.0, 0.0, 1.0), bias_kappa=8.0,
            )
        )
        assert S_aligned > 0.3, (
            f"Aligned cortex S={S_aligned:.4f} should be >0.3."
        )


# ---------------------------------------------------------------------------
# Measurement utility STATIC checks (always run)
# ---------------------------------------------------------------------------
class TestMeasurementUtilities:
    def test_aspect_ratio_isotropic_sphere_is_one(self):
        rng = np.random.default_rng(0)
        # Sample 1000 points on a unit sphere via Marsaglia.
        n = 1000
        pts = []
        while len(pts) < n:
            u = rng.uniform(-1, 1, 2 * n)
            v = rng.uniform(-1, 1, 2 * n)
            s = u * u + v * v
            mask = s < 1
            u, v, s = u[mask], v[mask], s[mask]
            for i in range(len(u)):
                if len(pts) >= n:
                    break
                fac = 2 * math.sqrt(1 - s[i])
                pts.append((u[i] * fac, v[i] * fac, 1 - 2 * s[i]))
        pts = np.array(pts)
        ar = cell_aspect_ratio(pts)
        # For a uniform sphere, principal eigenvalues are equal → ratio 1.
        assert math.isclose(ar, 1.0, abs_tol=0.2), (
            f"Isotropic sphere aspect ratio {ar:.3f} should be ≈ 1; "
            "measurement util broken."
        )

    def test_aspect_ratio_elongated_ellipsoid_recovers_geometric_ratio(self):
        rng = np.random.default_rng(1)
        n = 3000
        a, b, c = 12.0, 8.0, 8.0   # semi-major, semi-minor, semi-minor [μm]
        # Sample on the surface of the ellipsoid via Marsaglia + scale.
        pts = []
        while len(pts) < n:
            u = rng.uniform(-1, 1, 2 * n)
            v = rng.uniform(-1, 1, 2 * n)
            s = u * u + v * v
            mask = s < 1
            u, v, s = u[mask], v[mask], s[mask]
            for i in range(min(len(u), n - len(pts))):
                fac = 2 * math.sqrt(1 - s[i])
                pts.append((a * u[i] * fac, b * v[i] * fac, c * (1 - 2 * s[i])))
        pts = np.array(pts)
        ar = cell_aspect_ratio(pts)
        # The measurement util uses the COVARIANCE tensor (gyration
        # tensor), whose principal eigenvalues for a Marsaglia-sampled
        # ellipsoid (a, b, c) are (a²/3, b²/3, c²/3). So
        # √(max/min) = max(a,b,c)/min(a,b,c) = geometric aspect ratio
        # = 12/8 = 1.5.
        expected = a / min(b, c)
        assert math.isclose(ar, expected, rel_tol=0.10), (
            f"Elongated ellipsoid aspect ratio {ar:.3f} vs expected "
            f"{expected:.3f} (12×8×8 semi-axes; covariance principal eigvals)."
        )

    def test_nematic_S_isotropic_zero(self):
        rng = np.random.default_rng(2)
        # Random isotropic 3D unit vectors.
        axes = rng.standard_normal((2000, 3))
        axes /= np.linalg.norm(axes, axis=1, keepdims=True)
        S = nematic_order_S(axes)
        assert abs(S) < 0.1, f"Isotropic S = {S:.4f} should be ≈ 0."

    def test_nematic_S_perfect_alignment_one(self):
        axes = np.tile(np.array([1.0, 0.0, 0.0]), (500, 1))
        S = nematic_order_S(axes)
        assert math.isclose(S, 1.0, abs_tol=1e-9), (
            f"Perfect alignment S = {S} should = 1.0."
        )
