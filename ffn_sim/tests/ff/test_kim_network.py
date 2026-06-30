"""FF actin network vs Taeyoon Kim's BD model (ff/kim_network) — Stage 6i peer-model validation.

The single-cell FF network's correct oracle is a peer fine-grained model (Kim 2007 MIT thesis), NOT
the multi-cell spheroid law. Anchors: the semiflexible mesh ξ(151µM)≈64nm (Kim/actin lit ~50-100nm);
the FF box network reproduces the ξ ∝ C_A^(−1/2) scaling (Schmidt 1989 / MacKintosh / Kim); and the
network connectivity scales as z = 2R with the crosslinker ratio.
"""

import numpy as np
import pytest

from ffn_sim.ff.kim_network import (
    build_box_network,
    connectivity_z,
    length_density_per_um2,
    measure_mesh_size_um,
    mesh_size_theory_um,
)


def test_mesh_theory_matches_actin_literature():
    xi = mesh_size_theory_um(151.0) * 1e3      # nm
    assert 50.0 < xi < 100.0                    # Kim cortex mesh ~50-100 nm; analytic ~64 nm
    # ρ_L scales linearly in C_A → ξ = ρ_L^(-1/2) ∝ C_A^(-1/2)
    assert length_density_per_um2(302.0) == pytest.approx(2 * length_density_per_um2(151.0))


def test_ff_box_reproduces_mesh_scaling():
    """The FF 3D box actin network reproduces ξ ∝ C_A^(−1/2) (Kim / Schmidt-MacKintosh)."""
    Cs = [75.0, 150.0, 300.0, 600.0, 1200.0]
    xi = []
    for C in Cs:
        net, xl, meta = build_box_network(C, box_um=1.0, rng=np.random.default_rng(0))
        xi.append(measure_mesh_size_um(net, box_um=1.0))
    slope = np.polyfit(np.log(Cs), np.log(xi), 1)[0]
    assert slope == pytest.approx(-0.5, abs=0.06)   # measured exponent ≈ −0.5


def test_connectivity_scales_with_crosslinker_ratio():
    """Mean crosslinks/filament z = 2R (connectivity linear in the ACP ratio R)."""
    for R in (0.1, 0.3, 0.5, 1.0):
        net, xl, meta = build_box_network(151.0, box_um=1.0, R_acp=R, rng=np.random.default_rng(0))
        assert connectivity_z(net, xl) == pytest.approx(2 * R, abs=0.05)


def test_shear_modulus_rises_with_connectivity():
    """The FF cross-linked network stiffens monotonically with connectivity (floppy→rigid trend,
    Head/Levine/MacKintosh 2003 / Kim 2007). NOTE: with FF-native INEXTENSIBLE actin (reshape) +
    bending, there is no sharp G≈0 floppy phase — bending gives a finite small modulus below the
    central-force threshold (the bending-dominated regime), so we assert the monotone RISE, not G≈0."""
    from ffn_sim.ff.kim_network import shear_modulus
    G_lo, z_lo, _ = shear_modulus(C_A_uM=300.0, R_acp=0.2, k_xl=10.0, n_steps=4000)
    G_mid, z_mid, _ = shear_modulus(C_A_uM=300.0, R_acp=1.0, k_xl=10.0, n_steps=4000)
    G_hi, z_hi, _ = shear_modulus(C_A_uM=300.0, R_acp=2.0, k_xl=10.0, n_steps=4000)
    assert z_lo < z_mid < z_hi
    assert G_lo < G_mid < G_hi                      # stiffens with crosslinking
    assert G_hi > 3.0 * G_lo                        # substantial rise across the transition window


def test_shear_modulus_crosslink_limited_linear():
    """In the crosslink-limited regime (k_xl ≪ EA/L_seg, so actin is the rigid reshape backbone) the
    shear modulus is LINEAR in the crosslink junction stiffness k_xl, and matches the analytic affine
    form G ≈ (non-affine factor)·k_xl·ρ_L·ℓc with a STABLE non-affine factor < 1. This is the robust,
    no-tuning Kim/analytic comparison (assert the linearity + the analytic FORM, not a calibrated
    absolute — the absolute at sourced α-actinin stiffness, where k_xl ~ EA/L_seg, needs finite-EA)."""
    import numpy as np

    from ffn_sim.ff.kim_network import shear_modulus
    res = [shear_modulus(C_A_uM=300.0, R_acp=1.5, k_xl=k, n_steps=6000, seed=0) for k in (3.0, 30.0)]
    (G1, _, m1), (G2, _, m2) = res
    assert 8.0 < G2 / G1 < 12.0                                  # ×10 in k_xl → ~×10 in G (linear)
    # analytic affine form is an UPPER bound; the non-affine factor is stable across k_xl
    f1, f2 = G1 / m1["G_analytic_crosslink"], G2 / m2["G_analytic_crosslink"]
    assert 0.0 < f1 < 1.0 and 0.0 < f2 < 1.0                     # G < affine bound (non-affine softening)
    assert abs(f1 - f2) / f1 < 0.1                               # factor stable → the FORM G∝k_xl·ρ_L·ℓc holds
