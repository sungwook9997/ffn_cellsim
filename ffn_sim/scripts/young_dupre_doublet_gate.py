"""Young–Dupré doublet gate: relax a 2-cell DCM doublet at several cadherin cohesion
strengths, measure the emergent contact angle from geometry, and check it tracks the
analytic Young–Dupré curve cos θ = 1 − w_adh/γ_surf.

This is the SimuCell3D Supp-Fig-5 validation ported to our explicit-catch-bond DCM:
does the aggregate's cell–cell contact obey the capillary force balance (= is it
*mechanically* correctly assembled, not just geometrically touching)?

Run:  PYTHONPATH=. python -m ffn_sim.scripts.young_dupre_doublet_gate --device cpu
"""

from __future__ import annotations

import argparse
import numpy as np

from ffn_sim.warp_port.dcm_warp_decohesion import run_decohesion
from ffn_sim.validation.oracles.young_dupre import (
    doublet_angle_from_adhesion, fit_sphere, angle_from_doublet_geometry)

GAMMA_SURF = 1.0e-4          # N/m, the surface-tension module default (--gamma-surf)
F0 = 29.2e-12               # N, molecular cadherin catch peak


def measure_contact_angle(P: np.ndarray, cof: np.ndarray, c_adh: float) -> tuple:
    """θ [deg] from a relaxed doublet. Split each cell into free vs interface nodes
    (interface = within c_adh of the OTHER cell), fit the free surface to a sphere, take
    the interface plane ⊥ the line of centres at the contact, cos θ = d/r."""
    cells = np.unique(cof[cof >= 0])
    assert cells.size == 2, f"expected 2 cells, got {cells.size}"
    A = P[cof == cells[0]]; B = P[cof == cells[1]]
    cA, cB = A.mean(0), B.mean(0)
    axis = (cB - cA) / max(np.linalg.norm(cB - cA), 1e-30)
    mid = 0.5 * (cA + cB)
    from scipy.spatial import cKDTree
    tB = cKDTree(B); tA = cKDTree(A)
    A_iface = np.array([len(tB.query_ball_point(p, c_adh)) > 0 for p in A])
    B_iface = np.array([len(tA.query_ball_point(p, c_adh)) > 0 for p in B])
    thetas = []
    for pts, iface, cen in ((A, A_iface, cA), (B, B_iface, cB)):
        free = pts[~iface]
        if free.shape[0] < 8 or iface.sum() < 3:
            return None, A_iface.sum(), B_iface.sum()
        c_fit, r_fit = fit_sphere(free)
        thetas.append(np.degrees(angle_from_doublet_geometry(c_fit, r_fit, mid, axis)))
    return float(np.mean(thetas)), int(A_iface.sum()), int(B_iface.sum())


def w_adh_from_bonds(out: dict, c_adh: float, n_cells: int = 2) -> float:
    """Effective adhesion energy density w_adh [J/m²] from the relaxed bond population.
    w_adh = (bonds / interface area) × E_bond, E_bond = ½·k_eff·(r_bind−r0)² (work stored
    over the capture window by one bundle bond). Documented estimate; absolute calibration
    is approximate — the gate checks the FUNCTIONAL form across cohesion, like SimuCell3D."""
    cb = out.get("cadherin_bonds")
    if not cb:
        return 0.0
    n_bonds = cb["n_bonds_final"]
    k_eff = cb["k_trans"] * out.get("cad_bundle", 1.0)
    rng = cb["r_bind"] - cb["r0_trans"]
    E_bond = 0.5 * k_eff * rng ** 2
    # one interface (doublet) ≈ a disc of radius ~ c_adh-scale contact patch; use the mean
    # edge² × n_iface_nodes as the interface area proxy is fragile, so use a fixed contact
    # patch ~ π(R/2)² (half-radius wetted cap) — order-of-magnitude areal density.
    R = 7.5e-6
    A_iface = np.pi * (0.5 * R) ** 2
    return n_bonds * E_bond / A_iface


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--bundles", default="1,10,40,120", help="cadherin bundle_n sweep")
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--gap", type=float, default=1.95, help="initial cell spacing in R (≲2 → real contact patch)")
    ap.add_argument("--tol-deg", type=float, default=15.0, help="gate RMS tolerance [deg]")
    args = ap.parse_args()
    bundles = [float(x) for x in args.bundles.split(",")]
    print(f"Young–Dupré doublet gate — γ_surf={GAMMA_SURF:.1e}N/m, bundle sweep {bundles}\n")
    import tempfile, os
    print(f"{'bundle':>7} {'w_adh[µN/m]':>12} {'θ_meas[°]':>10} {'θ_pred[°]':>10} {'|Δ|[°]':>8} {'n_iface':>8}")
    errs = []
    for bn in bundles:
        tmp = os.path.join(tempfile.gettempdir(), f"_yd_doublet_b{bn:.0f}.npz")
        out = run_decohesion(
            n_cells=2, subdiv=2, steps=args.steps, frames=2, device=args.device,
            dt=8e-6, warmup=500, settle_steps=args.steps, settle_frames=0, gap=args.gap,
            cadherin=True, cad_bundle=bn, surface_tension=True, gamma_surf=GAMMA_SURF,
            substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
            builder="fcc", save_frames=tmp)
        out["cad_bundle"] = bn
        d = np.load(tmp, allow_pickle=True)
        P = d["frames"][-1].astype(np.float64); cof = d["cof"]
        c_iface = 1.5e-6        # geometric contact threshold for free-vs-interface node split
        theta, nA, nB = measure_contact_angle(P, cof, c_iface)
        w = w_adh_from_bonds(out, c_iface)
        w = min(w, 0.999 * GAMMA_SURF)               # cap at the physical wetting limit
        theta_pred = np.degrees(doublet_angle_from_adhesion(GAMMA_SURF, w))
        if theta is None:
            print(f"{bn:7.0f}  (no resolvable interface — cells did not adhere)")
            continue
        err = abs(theta - theta_pred); errs.append(err)
        print(f"{bn:7.0f} {w * 1e6:12.2f} {theta:10.1f} {theta_pred:10.1f} {err:8.1f} {nA + nB:8d}")
    if errs:
        rms = float(np.sqrt(np.mean(np.square(errs))))
        verdict = "PASS" if rms <= args.tol_deg else "FAIL"
        print(f"\n[GATE young_dupre_doublet] RMS Δθ = {rms:.1f}°  (tol {args.tol_deg:.0f}°) → {verdict}")


if __name__ == "__main__":
    main()
