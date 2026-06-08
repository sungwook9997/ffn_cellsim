"""Sanity gates for the H.7 Gate-B relaxed-constraint (unilateral M-SHAKE) mode.

The rigid bond constraint ``|s_a| = ℓ₀`` becomes the one-sided ``|s_a| ≤ ℓ₀`` when
``compression_release=True``: a bond is projected back to ℓ₀ only while in TENSION
(stretched); a COMPRESSED bond is released (λ_a = 0, free to shorten), permitting
buckling/condensation (contract §3 — the transmission lever under test).

Sanity Gate Protocol (CLAUDE.md) for ``shake_project_chains(..., compression_release)``:

1. **Boundary / superset parity (control §5.1)** — with the flag OFF the projection
   is the exact rigid M-SHAKE; with the flag ON but every bond already in tension the
   result is identical (relaxed is a clean superset of rigid).
2. **Sign-sense (control §5.2/§5.3)** — a stretched bond is pulled in (both modes); a
   compressed bond is pulled OUT in rigid mode but LEFT compressed in relaxed mode.
3. **Tension-side inextensibility preserved** — active (tension) bonds converge to ℓ₀
   to tolerance even when a neighbour is released.
4. **Chain decomposition correctness** — a released mid-chain bond must decouple its
   two active neighbours so each active run solves independently; the relaxed solve of
   a 3-bond chain with a compressed middle must equal the SHAKE of the two outer bonds
   solved as independent dimers (this is the non-trivial tridiagonal-masking claim).
5. **Conservation** — equal-mobility SHAKE corrections move the per-bond COM rigidly;
   released bonds apply no impulse, so the mobility-weighted COM of each independent
   active run is preserved.
6. **Condensation engages** — a straight chain under end-to-end compression shortens
   (end-to-end < contour) in relaxed mode but is held rigid (end-to-end = contour) in
   rigid mode (the §5.3 observable).

Both the vectorised uniform fast path (stacked ``(F, m+1)`` chains) and the ragged
per-chain fallback are exercised.
"""
from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.integrator.constrained_baoab import shake_project_chains

L0 = 0.5e-6  # 500 nm segment (physiological ℓ₀); units are immaterial to the projection
BIG_BOX = np.array([1.0e3 * L0, 1.0e3 * L0, 1.0e3 * L0])  # no min-image wrap
TOL = 1.0e-12


def _bond_lengths(pos, chain):
    p = np.asarray(chain)
    d = pos[p[:-1]] - pos[p[1:]]
    return np.linalg.norm(d, axis=1)


def _collinear_chain(offsets):
    """N beads on the x-axis at the given cumulative x-offsets (×L0)."""
    pos = np.zeros((len(offsets), 3), dtype=np.float64)
    pos[:, 0] = np.array(offsets) * L0
    return pos


# ---------------------------------------------------------------------------
# Gate 1 — superset parity
# ---------------------------------------------------------------------------
def test_flag_off_is_rigid_all_bonds_to_l0():
    """compression_release=False: every bond → ℓ₀ (exact rigid M-SHAKE)."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])  # stretch / compress / stretch
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    out = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                               tol=TOL, max_iter=200, compression_release=False)
    lens = _bond_lengths(out, chains[0])
    assert np.allclose(lens, L0, rtol=1e-9), lens / L0


def test_flag_on_all_tension_equals_rigid():
    """Relaxed is a clean superset: when no bond is compressed, identical to rigid."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.3, 2.6, 3.9])  # all bonds stretched
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    rigid = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                 tol=TOL, max_iter=200, compression_release=False)
    relaxed = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=200, compression_release=True)
    assert np.allclose(rigid, relaxed, atol=1e-14), np.abs(rigid - relaxed).max()


# ---------------------------------------------------------------------------
# Gate 2 — sign-sense: tension pulled in, compression released
# ---------------------------------------------------------------------------
def test_single_compressed_bond_released():
    """A lone compressed bond (|s|<ℓ₀): rigid pushes it to ℓ₀, relaxed leaves it short."""
    ref = _collinear_chain([0, 1])
    pred = _collinear_chain([0, 0.4])  # compressed to 0.4 ℓ₀
    chains = [np.array([0, 1])]  # 1-bond chain → ragged path (stacked needs ≥3 beads)
    inv_mass = np.ones(2)
    rigid = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                 tol=TOL, max_iter=200, compression_release=False)
    relaxed = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=200, compression_release=True)
    assert np.isclose(_bond_lengths(rigid, [0, 1])[0], L0, rtol=1e-9)      # rigid restores
    assert np.allclose(relaxed, pred, atol=1e-18)                         # relaxed untouched
    assert _bond_lengths(relaxed, [0, 1])[0] < L0


def test_single_tension_bond_enforced_both_modes():
    """A stretched bond is pulled to ℓ₀ in BOTH modes (tension side inextensible)."""
    ref = _collinear_chain([0, 1])
    pred = _collinear_chain([0, 1.5])
    chains = [np.array([0, 1])]  # 1-bond chain → ragged path
    inv_mass = np.ones(2)
    for cr in (False, True):
        out = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=200, compression_release=cr)
        assert np.isclose(_bond_lengths(out, [0, 1])[0], L0, rtol=1e-9), cr


# ---------------------------------------------------------------------------
# Gate 3 + 4 — tension-side preserved next to a release; chain decomposition
# ---------------------------------------------------------------------------
def test_released_middle_decouples_outer_bonds():
    """3-bond chain, compressed middle: outer (tension) bonds must each reach ℓ₀, and
    the result must equal the two outer bonds solved as INDEPENDENT dimers — the
    tridiagonal-masking decomposition claim (§4)."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])  # b0 stretch, b1 compress, b2 stretch
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    relaxed = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=300, compression_release=True)
    lens = _bond_lengths(relaxed, [0, 1, 2, 3])
    assert np.isclose(lens[0], L0, rtol=1e-8), lens / L0   # outer bond 0 active → ℓ₀
    assert np.isclose(lens[2], L0, rtol=1e-8), lens / L0   # outer bond 2 active → ℓ₀
    assert lens[1] < L0                                    # middle stays released

    # Ground truth: solve the two outer bonds as independent dimers (ragged path).
    d0 = shake_project_chains(pred, ref, [np.array([0, 1])], L0, inv_mass, BIG_BOX,
                              tol=TOL, max_iter=200, compression_release=False)
    d2 = shake_project_chains(pred, ref, [np.array([2, 3])], L0, inv_mass, BIG_BOX,
                              tol=TOL, max_iter=200, compression_release=False)
    assert np.allclose(relaxed[[0, 1]], d0[[0, 1]], atol=1e-12), "bond0 ≠ independent dimer"
    assert np.allclose(relaxed[[2, 3]], d2[[2, 3]], atol=1e-12), "bond2 ≠ independent dimer"


def test_released_run_preserves_mobility_com():
    """Each independent active run conserves its equal-mobility COM (internal forces)."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    relaxed = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=300, compression_release=True)
    # Active run {0,1} and {2,3}; equal mobility → mean position preserved per run.
    assert np.allclose(relaxed[[0, 1]].mean(0), pred[[0, 1]].mean(0), atol=1e-12)
    assert np.allclose(relaxed[[2, 3]].mean(0), pred[[2, 3]].mean(0), atol=1e-12)


# ---------------------------------------------------------------------------
# Gate 5 — condensation engages (the §5.3 observable)
# ---------------------------------------------------------------------------
def test_chain_condenses_under_compression_relaxed_only():
    """A straight chain whose ends are pushed inward shortens (end-to-end < contour)
    in relaxed mode, but is held at full contour by rigid M-SHAKE."""
    n = 7
    ref = _collinear_chain(list(range(n)))                  # contour = (n-1) ℓ₀
    # Compress every bond uniformly to 0.7 ℓ₀ (end-to-end pushed in).
    pred = _collinear_chain([0.7 * k for k in range(n)])
    chains = np.array([list(range(n))])
    inv_mass = np.ones(n)
    contour = (n - 1) * L0

    rigid = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                 tol=TOL, max_iter=300, compression_release=False)
    relaxed = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                   tol=TOL, max_iter=300, compression_release=True)
    e2e_rigid = np.linalg.norm(rigid[n - 1] - rigid[0])
    e2e_relaxed = np.linalg.norm(relaxed[n - 1] - relaxed[0])

    assert np.isclose(e2e_rigid, contour, rtol=1e-7)        # rigid: held at contour
    assert e2e_relaxed < 0.9 * contour                      # relaxed: condensed
    # every bond left compressed (≤ ℓ₀) in relaxed mode
    assert np.all(_bond_lengths(relaxed, chains[0]) <= L0 * (1 + 1e-7))


# ---------------------------------------------------------------------------
# Ragged per-chain fallback parity (mixed-length chains)
# ---------------------------------------------------------------------------
def test_ragged_path_matches_stacked_for_uniform_chains():
    """The ragged (list-of-arrays) fallback and the stacked fast path agree in
    relaxed mode for uniform chains."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])
    inv_mass = np.ones(4)
    stacked = shake_project_chains(pred, ref, np.array([[0, 1, 2, 3]]), L0, inv_mass,
                                   BIG_BOX, tol=TOL, max_iter=300, compression_release=True)
    ragged = shake_project_chains(pred, ref, [np.array([0, 1, 2, 3])], L0, inv_mass,
                                  BIG_BOX, tol=TOL, max_iter=300, compression_release=True)
    assert np.allclose(stacked, ragged, atol=1e-12)


def test_ragged_mixed_length_relaxed():
    """Two chains of different lengths (forces the ragged path): each condenses under
    compression with its outer/active bonds respected."""
    # chain A: 2 beads compressed; chain B: 4 beads with compressed middle.
    pos_ref = np.zeros((6, 3)); pos_ref[:, 0] = np.array([0, 1, 0, 1, 2, 3]) * L0
    pos_pred = np.zeros((6, 3))
    pos_pred[:, 0] = np.array([0, 0.4, 0, 1.2, 1.5, 2.7]) * L0
    chains = [np.array([0, 1]), np.array([2, 3, 4, 5])]
    inv_mass = np.ones(6)
    out = shake_project_chains(pos_pred, pos_ref, chains, L0, inv_mass, BIG_BOX,
                               tol=TOL, max_iter=300, compression_release=True)
    assert _bond_lengths(out, [0, 1])[0] < L0                       # A released
    lensB = _bond_lengths(out, [2, 3, 4, 5])
    assert np.isclose(lensB[0], L0, rtol=1e-8)                      # B outer active
    assert np.isclose(lensB[2], L0, rtol=1e-8)
    assert lensB[1] < L0                                            # B middle released


# ---------------------------------------------------------------------------
# Euler-threshold (load-based) release — contract §8: buckle only above F_crit
# ---------------------------------------------------------------------------
def test_euler_threshold_gates_release():
    """A compressed bond releases only if its compressive constraint force exceeds
    F_crit (release_load_crit). Below threshold it stays rigid (the rod holds);
    above threshold it buckles (releases)."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])  # only the MIDDLE bond is compressed
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    dt = 1.0e-7

    # rigid λ → the compressive force the middle bond would carry if held rigid
    _, lam_rigid = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                        tol=TOL, max_iter=200, return_lambdas=True,
                                        compression_release=False)
    lam_mid = float(np.asarray(lam_rigid)[0, 1])
    assert lam_mid < 0.0, "middle bond should be compressive (λ<0)"
    T_mid = abs(lam_mid) * L0 / dt  # compressive constraint force [arb force units]

    # threshold ABOVE the load → middle not buckle-eligible → stays rigid (= ℓ₀)
    hi = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                              tol=TOL, max_iter=300, compression_release=True,
                              release_load_crit=2.0 * T_mid, dt=dt)
    assert np.allclose(_bond_lengths(hi, [0, 1, 2, 3]), L0, rtol=1e-7), \
        "below-threshold compression must stay rigid"

    # threshold BELOW the load → middle eligible → releases (stays compressed)
    lo = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                              tol=TOL, max_iter=300, compression_release=True,
                              release_load_crit=0.5 * T_mid, dt=dt)
    lens = _bond_lengths(lo, [0, 1, 2, 3])
    assert lens[1] < L0, "above-threshold compression must release"
    assert np.isclose(lens[0], L0, rtol=1e-7) and np.isclose(lens[2], L0, rtol=1e-7), \
        "tension bonds stay inextensible"


def test_euler_release_requires_dt():
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])
    chains = np.array([[0, 1, 2, 3]])
    with pytest.raises(ValueError, match="dt"):
        shake_project_chains(pred, ref, chains, L0, np.ones(4), BIG_BOX,
                             compression_release=True, release_load_crit=1.0, dt=None)


# ---------------------------------------------------------------------------
# return_lambdas: released bonds carry zero multiplier
# ---------------------------------------------------------------------------
def test_released_bond_lambda_is_zero():
    """The accumulated Lagrange multiplier on a released bond is exactly 0 (it applies
    no constraint impulse → contributes nothing to the rigid-MOP tension channel)."""
    ref = _collinear_chain([0, 1, 2, 3])
    pred = _collinear_chain([0, 1.2, 1.5, 2.7])
    chains = np.array([[0, 1, 2, 3]])
    inv_mass = np.ones(4)
    out, lam = shake_project_chains(pred, ref, chains, L0, inv_mass, BIG_BOX,
                                    tol=TOL, max_iter=300, return_lambdas=True,
                                    compression_release=True)
    lam = np.asarray(lam)
    assert lam.shape == (1, 3)
    assert lam[0, 1] == 0.0                 # middle (released) bond → λ ≡ 0
    assert lam[0, 0] != 0.0 and lam[0, 2] != 0.0  # active outer bonds carry impulse
