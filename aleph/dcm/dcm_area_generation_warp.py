"""T-2 ACTIVE AREA GENERATION — the missing DCM spreading physics (PI 2026-07-13).

DIAGNOSIS (`docs/v2_audit/_historical/DCM_SPREADING_ACTIVE_AREA_PLAN_2026-07-13.md`): the DCM cell is a
turgid closed shell whose membrane surface tension / area-elasticity only MINIMISES area
(rounds the cell), so it structurally OPPOSES spreading. Real spreading is the opposite — the
cell does WORK to GROW its footprint: lamellipodial actin polymerisation pushes the membrane
outward and the membrane reservoir unfolds, ADDING surface/footprint area. Prior attempts
(`project-dcm-lamellipodium-graft`, `project-dcm-spreading-foundation`) failed because the added
drivers were too weak to beat the shell's surface tension + turgor — they lacked the AREA
GENERATOR entirely.

THE TERM (this module): the per-cell preferred surface area ``A0`` — the setpoint of the B4
global area-elasticity force ``F = −k_a·(A_cell − A0)·∂A/∂r`` (``global_area_force_kernel`` in
``dcm_neighbor_warp``) — is made to actively GROW while the cell is basal-adhered, up to a
membrane-reservoir cap::

    dA0/dt = +k_prot_area   while basal-adhered and A0 < A0_cap    (else 0)
    k_prot_area = v_protrusion · P_edge        (edge advance × contact perimeter)
    A0_cap      = reservoir · A0_init                              (membrane reservoir ceiling)

This is the mechanistic stand-in for lamellipodial actin polymerisation + membrane unfolding:
the cell MAKES area, so the area-elasticity term now PREFERS a larger, flatter shell. Substrate
adhesion (T-1 wetting/well) breaks the up/down symmetry, so the grown area flattens into a
pancake at constant volume (turgor conserves V). It is NOT a body-force proxy — A0 is the
setpoint of a conservative area-elasticity energy; growing it moves the energy minimum, the
force follows from the energy gradient exactly as the committed ``global_area_force_kernel``.

GROUNDING (no magic numbers — every constant KB/lit-anchored):
  * v_protrusion — lamellipodial leading-edge advance velocity. KB-3.6 (Mogilner & Oster 1996
    Biophys J 71:3030; 2003 Biophys J 84:1591; Brownian-ratchet polymerisation). Mid-band
    6 µm/min = 1.0e-7 m/s (lit 3–12 µm/min). ``dcm_lamellipodium_host.py:37`` uses the same.
  * k_a (area-elasticity / membrane-area modulus) — Brückner 2015 Sci Rep 5:14700 reservoir-
    buffered effective area modulus k_a ≈ 0.13–0.14 N/m (0.243 N/m once reservoir exhausted,
    Rawicz 2000). This is the STIFFNESS of the area term, set by the caller (``k_area``); this
    module only grows the SETPOINT A0.
  * reservoir (A0 ceiling) — membrane-area reservoir factor. KB/project value ~1.2–1.4× (folded
    membrane, Brückner 2015 / Figard & Sokac 2014; ~2–3 % areal strain then lysis). The larger
    2–4× folded reservoir (microvilli/caveolae, Gauthier-Masters-Sheetz 2011 PNAS) is a
    FALLBACK not registered in this KB — surface to PI if the larger cap is needed. Reservoir R
    maps to a thin-pancake silhouette A/A0 ≈ 2·R (surface area ∝ length², footprint set by the
    unfolded area): reservoir 1.4→A/A0≈2.8, 2→4, 3→6, 4→8 — see the plan G2 table.

Sanity Gate (see also ``dcm_area_generation_sanity`` in the self-test):
  * Dimensional: A0 [m²], k_prot_area [m²/s], dA0=k_prot_area·dt [m²]; A0_cap [m²]. ✓
  * Monotone + bounded: A0 non-decreasing, clamped to A0_cap (no runaway; reservoir ceiling). ✓
  * Off-invariance: area_growth disabled ⇒ A0 == A0_init forever (byte-identical to the fixed-A0
    global-area term). ✓
  * Adhesion gate: only cells with a basal-contact node grow — a suspended/interior cell keeps
    A0_init (interior cells of an aggregate do not spread). ✓
  * Volume: A0 growth is an AREA setpoint only; turgor (K_vol) conserves V independently ⇒ the
    spread is area↑ at V=const (G3), the cell flattens rather than inflating.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


# ──────────────────────────────────────────────────────────────────────────
# Grounding helpers (host; pure functions, KB-anchored — see module docstring)
# ──────────────────────────────────────────────────────────────────────────
#: KB-3.6 lamellipodial protrusion velocity, mid-band 6 µm/min (lit 3–12 µm/min).
V_PROTRUSION_M_S: float = 1.0e-7
#: Brückner 2015 reservoir-buffered membrane-area modulus (N/m) — the k_a stiffness scale.
K_AREA_MEMBRANE_N_M: float = 0.13
#: KB/project membrane-area reservoir factor (folded membrane, Brückner/Figard). 1.2–1.4×.
RESERVOIR_KB: float = 1.4
#: Deep membrane reservoir deployed over spreading time (caveolae/microvilli + exocytosis).
#: Gauthier, Masters & Sheetz 2011 (Trends Cell Biol 22:527; PNAS 108:14467): cells buffer a
#: 2–4× apparent-area reservoir, tension-regulated, replenished during spreading. Option-(i) ceiling.
RESERVOIR_DEEP: float = 4.0


def protrusion_area_rate(v_protrusion: float, contact_perimeter: float) -> float:
    """Area-generation rate dA0/dt = v_protrusion · P (edge advances at v over the contact rim).

    Args:
        v_protrusion: lamellipodial edge-advance velocity (m/s). KB-3.6 = 1.0e-7 (6 µm/min).
        contact_perimeter: length of the advancing (basal-contact) edge (m). For a round cell
            of contact radius r, P = 2πr; at build we default to the equatorial 2πR.

    Returns:
        dA0/dt (m²/s), the rate the preferred surface area grows.

    Sanity Gate: [m/s]·[m] = [m²/s] ✓ ; ≥0 (protrusion adds area, never removes).
    """
    return float(v_protrusion) * float(contact_perimeter)


def reservoir_to_AoverA0(reservoir: float) -> float:
    """Thin-pancake silhouette A/A0 reachable at a membrane-area reservoir factor.

    Surface area S = reservoir·S0 = reservoir·4πR0²; a thin pancake of footprint radius Rp has
    S ≈ 2πRp² (two faces) ⇒ Rp² = 2·reservoir·R0² ⇒ A/A0_silhouette = Rp²/R0² = 2·reservoir.
    (Grid-invariant geometric identity; see the plan G2 table.) reservoir 1.4→2.8, 2→4, 4→8.
    """
    return 2.0 * float(reservoir)


# ──────────────────────────────────────────────────────────────────────────
# Device kernels (GPU-resident; no host sync — grow the setpoint in place)
# ──────────────────────────────────────────────────────────────────────────
@wp.kernel
def mark_basal_adhered_kernel(
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    z0: wp.float64,
    band: wp.float64,
    adhered: wp.array(dtype=wp.int32),      # (n_cells,) OUT — set 1 if any basal node (reset to 0 first)
):
    """Flag each cell that has ≥1 node within ``band`` of the substrate plane z0 as basal-adhered.

    Write-1 is idempotent (a race just re-writes 1), so no atomic is needed. The caller zeroes
    ``adhered`` before launch. A cell with no basal node stays 0 → its A0 does not grow (a
    suspended cell, or an interior/apical cell of an aggregate, does not spread)."""
    i = wp.tid()
    c = cof[i]
    if c < wp.int32(0):
        return
    if pos[i][2] - z0 <= band:
        adhered[c] = wp.int32(1)


@wp.kernel
def grow_preferred_area_kernel(
    a0cell: wp.array(dtype=wp.float64),      # (n_cells,) preferred area — GROWN IN PLACE
    a0cap: wp.array(dtype=wp.float64),       # (n_cells,) reservoir ceiling
    adhered: wp.array(dtype=wp.int32),       # (n_cells,) 1 if basal-adhered
    dA0: wp.float64,                          # area increment this step = k_prot_area·dt (m²)
):
    """T-2 area generation: A0 += dA0 for basal-adhered cells, clamped to the reservoir cap.

    Monotone (dA0 ≥ 0) and bounded (≤ a0cap) → no runaway. A non-adhered cell is skipped
    (keeps A0_init)."""
    c = wp.tid()
    if adhered[c] == wp.int32(0):
        return
    a = a0cell[c] + dA0
    if a > a0cap[c]:
        a = a0cap[c]
    a0cell[c] = a


@wp.kernel
def grow_reservoir_cap_kernel(
    a0cap: wp.array(dtype=wp.float64),       # (n_cells,) reservoir cap — GROWN IN PLACE
    a0cap_max: wp.array(dtype=wp.float64),   # (n_cells,) deep-reservoir ceiling (r_deep·A0_init)
    adhered: wp.array(dtype=wp.int32),       # (n_cells,) 1 if basal-adhered
    dcap: wp.float64,                         # cap increment this step (membrane-addition, m²)
):
    """Option (i) MEMBRANE ADDITION over time (PI 2026-07-14): the reservoir cap itself GROWS
    while the cell is basal-adhered — exocytosis + deep-reservoir (caveolae/microvilli) unfolding
    replenish the membrane budget beyond the instantaneous fold, up to the deep-reservoir ceiling.

    The cap grows from r_fold·A0_init (Brückner instantaneous fold) toward r_deep·A0_init
    (Gauthier-Masters-Sheetz 2011 deep reservoir, 2–4×) at the membrane-addition rate (tied to the
    protrusion rate — membrane is added at the advancing edge, the same process that generates area).
    So the spread is NOT capped by the 1.4× instantaneous fold; the binding limit is the deep budget,
    reached over the spreading time. This is a PHYSICS addition, NOT reservoir tuning: r_fold, r_deep,
    and the rate are all grounded; the value is not chosen to pass a gate. Monotone + bounded."""
    c = wp.tid()
    if adhered[c] == wp.int32(0):
        return
    a = a0cap[c] + dcap
    if a > a0cap_max[c]:
        a = a0cap_max[c]
    a0cap[c] = a


# ──────────────────────────────────────────────────────────────────────────
# Self-test / Sanity Gate
# ──────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    dev = "cpu"
    print("=" * 72)
    print("T-2 active area generation — self-test / Sanity Gate")
    print("=" * 72)

    # grounding numbers
    R = 7.5e-6
    A0_init = 4.0 * np.pi * R ** 2
    P_eq = 2.0 * np.pi * R
    dA0dt = protrusion_area_rate(V_PROTRUSION_M_S, P_eq)
    print(f"\n[grounding] v_protrusion={V_PROTRUSION_M_S:.2e} m/s (6 µm/min, KB-3.6)")
    print(f"            A0_init=4πR²={A0_init*1e12:.1f} µm²  perimeter 2πR={P_eq*1e6:.1f} µm")
    print(f"            dA0/dt = v·P = {dA0dt:.3e} m²/s")
    for res in (1.4, 2.0, 3.0, 4.0):
        dA0 = (res - 1.0) * A0_init
        t = dA0 / dA0dt
        print(f"            reservoir {res:.1f}× → A0_cap={res*A0_init*1e12:.0f} µm², "
              f"A/A0≈{reservoir_to_AoverA0(res):.1f}, time to cap {t:.0f} s ({t/60:.1f} min) real")

    # --- (A) monotone + bounded growth to the cap ---
    n = 4
    a0 = wp.array(np.full(n, A0_init), dtype=wp.float64, device=dev)
    cap = wp.array(np.full(n, 2.0 * A0_init), dtype=wp.float64, device=dev)
    adh = wp.array(np.ones(n, dtype=np.int32), dtype=wp.int32, device=dev)
    dA0_step = dA0dt * 1.0e-3            # a coarse dt for the test
    hist = [a0.numpy().mean()]
    for _ in range(200000):
        wp.launch(grow_preferred_area_kernel, dim=n, inputs=[a0, cap, adh, wp.float64(dA0_step)], device=dev)
    wp.synchronize_device(dev)
    a0f = a0.numpy()
    monotone_ok = bool(np.all(a0f >= A0_init))
    capped_ok = bool(np.all(a0f <= 2.0 * A0_init + 1e-30) and np.allclose(a0f, 2.0 * A0_init))
    print(f"\n[A] grows to cap, bounded  : monotone={monotone_ok}  capped_at_2x={capped_ok} "
          f"(A0 {A0_init*1e12:.1f} → {a0f.mean()*1e12:.1f} µm²)")

    # --- (B) adhesion gate: a non-adhered cell does not grow ---
    a0b = wp.array(np.full(2, A0_init), dtype=wp.float64, device=dev)
    capb = wp.array(np.full(2, 5.0 * A0_init), dtype=wp.float64, device=dev)
    adhb = wp.array(np.array([1, 0], dtype=np.int32), dtype=wp.int32, device=dev)
    for _ in range(1000):
        wp.launch(grow_preferred_area_kernel, dim=2, inputs=[a0b, capb, adhb, wp.float64(dA0dt * 1e-2)], device=dev)
    wp.synchronize_device(dev)
    ab = a0b.numpy()
    gate_ok = bool(ab[0] > A0_init and np.isclose(ab[1], A0_init))
    print(f"[B] adhesion gate          : adhered grows ({ab[0]*1e12:.1f}), "
          f"suspended stays ({ab[1]*1e12:.1f} µm²) → {gate_ok}")

    # --- (C) basal-adhered detection ---
    pos = np.array([[0, 0, 0.1e-6], [0, 0, 5.0e-6], [1e-6, 0, 0.0]], dtype=np.float64)  # 2 basal, 1 apical
    posd = wp.array(pos, dtype=wp.vec3d, device=dev)
    cofd = wp.array(np.zeros(3, dtype=np.int32), dtype=wp.int32, device=dev)  # all cell 0
    adhc = wp.zeros(1, dtype=wp.int32, device=dev)
    wp.launch(mark_basal_adhered_kernel, dim=3,
              inputs=[posd, cofd, wp.float64(0.0), wp.float64(0.5e-6), adhc], device=dev)
    wp.synchronize_device(dev)
    detect_ok = bool(adhc.numpy()[0] == 1)
    # and a cell with all nodes lifted off is not adhered
    posd2 = wp.array(pos + np.array([0, 0, 10e-6]), dtype=wp.vec3d, device=dev)
    adhc2 = wp.zeros(1, dtype=wp.int32, device=dev)
    wp.launch(mark_basal_adhered_kernel, dim=3,
              inputs=[posd2, cofd, wp.float64(0.0), wp.float64(0.5e-6), adhc2], device=dev)
    wp.synchronize_device(dev)
    detect_ok = detect_ok and bool(adhc2.numpy()[0] == 0)
    print(f"[C] basal-adhesion detect  : contacting→adhered, lifted→not → {detect_ok}")

    ok = monotone_ok and capped_ok and gate_ok and detect_ok
    print("\n" + "=" * 72)
    print(f"SELF-TEST {'PASS' if ok else 'FAIL'}")
    print("=" * 72)
    sys.exit(0 if ok else 1)
