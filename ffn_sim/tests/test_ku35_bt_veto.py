"""TRACK-3 — tests for the ``_bipolar_accepts`` veto-anatomy assay.

Validates (HOOMD-free) that the instrumented replica in
``scripts/h3_ku35_bt_veto`` (a) reproduces the core accept/reject decision,
(b) attributes each decision to the correct clause (polarity vs
different-filament), and (c) encodes the key physics fact that two beads of ONE
actin filament can never be antiparallel — so the different-filament clause never
masks a would-be same-filament antiparallel pair.

These are the contracts that make the veto-breakdown number trustworthy: if the
replica diverged from the core, the per-clause attribution (and hence the TRACK-3
verdict) would be meaningless.
"""
from __future__ import annotations

import types
from collections import Counter

import numpy as np

from ffn_sim.scripts.h3_ku35_bt_veto import (
    _FakeAction,
    _make_instrumented_accepts,
    self_test,
)


def test_self_test_passes():
    """The script's HOOMD-free self-test (clause attribution) must pass."""
    assert self_test(verbose=False) is True


def _two_filament_fixture():
    """Two 5-bead filaments of OPPOSITE polarity in û = +x.

    Filament 0 beads ascend in x → minus-ward tangent points −x → m̂·û < 0 (a
    '− side' filament). Filament 1 beads descend in x → m̂·û > 0 (a '+ side'
    filament). Returns (pos, nb, axes).
    """
    H, nb = 1, 5
    axes = [[1.0, 0.0, 0.0]]
    ell = 100e-9
    n_act = 2 * nb
    pos = np.zeros((n_act + 2, 3))
    for j in range(nb):
        pos[j] = [j * ell, ell, 0.0]
    for j in range(nb):
        pos[nb + j] = [(nb - 1 - j) * ell, -ell, 0.0]
    return pos, nb, axes, H


def test_polarity_clause_rejects_wrong_sign():
    pos, nb, axes, H = _two_filament_fixture()
    act = _FakeAction(H=H, nb=nb, axes=axes, head_bound_filament=[-1, -1])
    cnt = Counter()
    f = types.MethodType(_make_instrumented_accepts(cnt), act)
    # + head (head_local 0) onto filament 0 (m̂·û < 0) → polarity reject.
    assert f(0, 3, pos) is False
    assert cnt["polarity_plus"] == 1
    assert cnt["different_filament"] == 0


def test_polarity_clause_accepts_right_sign():
    pos, nb, axes, H = _two_filament_fixture()
    act = _FakeAction(H=H, nb=nb, axes=axes, head_bound_filament=[-1, -1])
    cnt = Counter()
    f = types.MethodType(_make_instrumented_accepts(cnt), act)
    # + head onto filament 1 (m̂·û > 0), other side free → accept.
    assert f(0, nb + 3, pos) is True
    assert cnt["accept"] == 1


def test_different_filament_clause_only_with_matching_polarity():
    """Different-filament fires ONLY when polarity already passed and the other
    side grips this filament."""
    pos, nb, axes, H = _two_filament_fixture()
    # − head (idx 1) already grips filament 1 (the m̂·û>0 one).
    act = _FakeAction(H=H, nb=nb, axes=axes, head_bound_filament=[-1, 1])
    cnt = Counter()
    f = types.MethodType(_make_instrumented_accepts(cnt), act)
    # + head onto another bead of filament 1: polarity OK (m̂·û>0) → diff-fil reject.
    assert f(0, nb + 2, pos) is False
    assert cnt["different_filament"] == 1
    assert cnt["polarity_plus"] == 0


def test_relaxation_flips_diff_fil_to_accept():
    pos, nb, axes, H = _two_filament_fixture()
    act = _FakeAction(H=H, nb=nb, axes=axes, head_bound_filament=[-1, 1])
    cnt = Counter()
    f = types.MethodType(
        _make_instrumented_accepts(cnt, relax_diff_fil=True), act)
    assert f(0, nb + 2, pos) is True
    assert cnt["relaxed_diff_fil_accept"] == 1


def test_single_filament_never_antiparallel():
    """KEY PHYSICS: every bead of one filament shares one polarity, so a + head
    is polarity-rejected at EVERY bead of a − side filament — there is no
    'antiparallel region farther along' the different-filament clause could mask.
    """
    pos, nb, axes, H = _two_filament_fixture()
    act = _FakeAction(H=H, nb=nb, axes=axes, head_bound_filament=[-1, -1])
    cnt = Counter()
    f = types.MethodType(_make_instrumented_accepts(cnt), act)
    # filament 0 is a − side filament (m̂·û < 0): a + head must be rejected at
    # EVERY one of its beads, all via polarity, never via different-filament.
    decisions = [f(0, j, pos) for j in range(nb)]
    assert all(d is False for d in decisions)
    assert cnt["polarity_plus"] == nb
    assert cnt["different_filament"] == 0
