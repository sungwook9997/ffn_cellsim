"""Emit the committed Warp-vs-HOOMD parity-result JSON the results-gate checks.

The pytest gate (``tests/warp_port/test_*.py``) is the real verification — it
asserts each Warp piece matches the committed HOOMD reference fixture within the
plan's thresholds. THIS produces the lightweight, committed JSON summary that
``outputs/tag_kb/verify_runs.py`` (the results-integrity gate) reads on disk, so
a verified Warp piece is enforced by the same gate that guards every other
headline result (the C12 lesson: claim a committed lightweight artifact, never
an uncommitted bulk file).

Reproducible on CPU (Warp CPU backend): re-run to regenerate the JSON.

    python ffn_sim/warp_port/parity_report.py
"""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
OUT = os.path.join(FIX, "warp_parity_results.json")


def _baoab_parity() -> dict:
    """Run the B1 BAOAB Warp kernel against both committed fixtures."""
    from ffn_sim.warp_port.baoab_warp import run_baoab_warp

    results = {}
    for label in ("kt0", "ktpos"):
        path = os.path.join(FIX, f"baoab_ref_{label}.npz")
        fx = dict(np.load(path))
        got = run_baoab_warp(
            pos0=fx["pos0"], image0=fx["image0"], force=fx["force"],
            gamma=fx["gamma"], kT=float(fx["kT"]), dt=float(fx["dt"]),
            noise=fx["noise"], box=tuple(float(x) for x in fx["box"]),
            device="cpu",
        )
        results[label] = {
            "kT": float(fx["kT"]),
            "N": int(fx["pos0"].shape[0]),
            "K": int(fx["K"]),
            "max_abs_pos_diff": float(np.abs(got["pos"] - fx["ref_pos"]).max()),
            "max_abs_image_diff": int(np.abs(got["image"] - fx["ref_image"]).max()),
            "max_abs_prv_diff": float(np.abs(got["prv"] - fx["ref_prv"]).max()),
            "max_abs_image_flag": int(np.abs(fx["ref_image"]).max()),  # wrap exercised
        }
    return results


def _radial_shell_parity() -> dict:
    """Run the B2 radial-shell Warp kernel against the 3 committed law fixtures.

    Reports both the host-reduced force-LAW parity (gated < 1e-12) and the full
    Warp atomic-reduction parity (diagnostic, reduction-order < 1e-8)."""
    from ffn_sim.warp_port.radial_shell_warp import run_radial_shell_warp

    out = {}
    for name in ("nucleus", "membrane", "turgor"):
        fx = dict(np.load(os.path.join(FIX, f"radial_ref_{name}.npz")))
        kw = dict(
            pos=fx["pos"], tag=fx["tag"],
            tag_range=(int(fx["t0"]), int(fx["t1"])),
            law=int(fx["law"]), R0=float(fx["R0"]), pa=float(fx["pa"]),
            pb=float(fx["pb"]), pc=float(fx["pc"]), pd=float(fx["pd"]),
            device="cpu",
        )
        rec = {"law": int(fx["law"]), "N": int(fx["N"])}
        for mode in ("host", "warp"):
            got = run_radial_shell_warp(reduce=mode, **kw)
            dF = float(np.abs(got["force"] - fx["ref_force"]).max())
            sF = float(np.abs(fx["ref_force"]).max()) + 1e-300
            dU = float(np.abs(got["energy"] - fx["ref_energy"]).max())
            sU = float(np.abs(fx["ref_energy"]).max()) + 1e-300
            rec[f"{mode}_force_rel"] = dF / sF
            rec[f"{mode}_energy_rel"] = dU / sU
        out[name] = rec
    return out


def _shake_parity() -> dict:
    """M-SHAKE: Warp chain-constraint projection vs the committed Python reference."""
    from ffn_sim.warp_port.shake_warp import run_shake_warp

    fx = dict(np.load(os.path.join(FIX, "shake_ref.npz")))
    got = run_shake_warp(
        pred_pos=fx["pred"], ref_pos=fx["ref"], chains=fx["chains"],
        inv_mass=fx["inv_mass"], rest_length=float(fx["r0"]),
        box_L=(float(fx["L"]),) * 3, tol=float(fx["tol"]),
        max_iter=int(fx["max_iter"]), device="cpu",
    )
    pscale = float(np.abs(fx["ref_proj"]).max())
    lscale = float(np.abs(fx["ref_lam"]).max()) + 1e-30
    s = got["pos"][fx["chains"][:, :-1]] - got["pos"][fx["chains"][:, 1:]]
    bond = np.sqrt((s * s).sum(axis=2))
    return {
        "F": int(fx["F"]), "m": int(fx["m"]),
        "nonconverged": got["nonconverged"],
        "pos_rel": float(np.abs(got["pos"] - fx["ref_proj"]).max()) / pscale,
        "lam_rel": float(np.abs(got["lam"] - fx["ref_lam"]).max()) / lscale,
        "bond_rel_err": float(np.abs(bond - float(fx["r0"])).max() / float(fx["r0"])),
    }


def _bond_parity() -> dict:
    """Harmonic bond (cortex axial spring) vs HOOMD md.bond.Harmonic."""
    from ffn_sim.warp_port.network_warp import run_harmonic_bond_warp

    fx = dict(np.load(os.path.join(FIX, "network_bond_ref.npz")))
    got = run_harmonic_bond_warp(
        pos=fx["pos"], bonds=fx["bonds"], k=float(fx["k"]), r0=float(fx["r0"]),
        box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = float(np.abs(fx["ref_energy"]).max()) + 1e-30
    return {
        "N": int(fx["N"]), "B": int(fx["B"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
        "energy_rel": float(np.abs(got["energy"] - fx["ref_energy"]).max()) / Uscale,
    }


def _angle_parity() -> dict:
    """Harmonic angle (cortex bending) vs HOOMD md.angle.Harmonic."""
    from ffn_sim.warp_port.network_warp import run_harmonic_angle_warp

    fx = dict(np.load(os.path.join(FIX, "network_angle_ref.npz")))
    got = run_harmonic_angle_warp(
        pos=fx["pos"], angles=fx["angles"], k=float(fx["k"]), t0=float(fx["t0"]),
        box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = float(np.abs(fx["ref_energy"]).max()) + 1e-30
    return {
        "N": int(fx["N"]), "A": int(fx["A"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
        "energy_rel": float(np.abs(got["energy"] - fx["ref_energy"]).max()) / Uscale,
    }


def _lj_parity() -> dict:
    """LJ excluded volume (WCA) vs HOOMD md.pair.LJ (mode=shift, r_cut=2^(1/6)σ)."""
    from ffn_sim.warp_port.network_warp import run_wca_pair_warp

    fx = dict(np.load(os.path.join(FIX, "network_lj_ref.npz")))
    got = run_wca_pair_warp(
        pos=fx["pos"], epsilon=float(fx["epsilon"]), sigma=float(fx["sigma"]),
        r_cut=float(fx["r_cut"]), box_L=(float(fx["L"]),) * 3, shift=True, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = float(np.abs(fx["ref_energy"]).max()) + 1e-30
    return {
        "N": int(fx["N"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
        "energy_rel": float(np.abs(got["energy"] - fx["ref_energy"]).max()) / Uscale,
    }


def _fixman_parity() -> dict:
    """Fixman: Warp metric pseudo-force vs the committed Python LAPACK reference."""
    from ffn_sim.warp_port.fixman_warp import run_fixman_warp

    fx = dict(np.load(os.path.join(FIX, "fixman_ref.npz")))
    got = run_fixman_warp(
        pos=fx["pos"], chains=fx["chains"], inv_gamma=fx["inv_gamma"],
        kT=float(fx["kT"]), box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = abs(float(fx["ref_U"])) + 1e-30
    return {
        "F": int(fx["F"]), "m": int(fx["m"]),
        "bad_sign": got["bad_sign"],
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
        "U_rel": abs(got["U_F"] - float(fx["ref_U"])) / Uscale,
    }


def _dcm_contact_parity() -> dict:
    """DCM node-face contact: Warp vs the numpy brute-force ground truth."""
    from ffn_sim.warp_port.dcm_contact_warp import run_node_face_contact_warp

    fx = dict(np.load(os.path.join(FIX, "dcm_contact_ref.npz")))
    got = run_node_face_contact_warp(
        pos=fx["pos"], cell_of_node=fx["cell_of_node"], faces=fx["faces"],
        face_cell=fx["face_cell"], rep_strength=float(fx["rep_strength"]),
        adh_strength=float(fx["adh_strength"]), c_rep=float(fx["c_rep"]),
        c_adh=float(fx["c_adh"]), cad_mult=fx["cad_mult"], device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    return {
        "N": int(fx["N"]), "M": int(fx["M"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
    }


def _dcm_turgor_parity() -> dict:
    """DCM exact-volume turgor force: Warp vs the committed HOOMD DcmTurgorForce."""
    from ffn_sim.warp_port.dcm_turgor_warp import run_dcm_turgor_warp

    fx = dict(np.load(os.path.join(FIX, "dcm_turgor_ref.npz")))
    kw = dict(
        pos=fx["pos"], faces=fx["faces"], face_cell=fx["face_cell"],
        n_cells=int(fx["n_cells"]), V0=float(fx["V0"]),
        turgor_dP0=float(fx["turgor_dP0"]), K_vol=float(fx["K_vol"]), device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    rec = {"N": int(fx["N"]), "n_cells": int(fx["n_cells"])}
    for mode in ("host", "warp"):
        got = run_dcm_turgor_warp(reduce=mode, **kw)
        rec[f"{mode}_force_rel"] = float(
            np.abs(got["force"] - fx["ref_force"]).max()) / Fscale
    return rec


def _dcm_cohesion_parity() -> dict:
    """DCM node-node cohesion: Warp vs the committed HOOMD DcmTentContact."""
    from ffn_sim.warp_port.dcm_cohesion_warp import run_dcm_cohesion_warp

    fx = dict(np.load(os.path.join(FIX, "dcm_cohesion_ref.npz")))
    got = run_dcm_cohesion_warp(
        pos=fx["pos"], cell_of_node=fx["cell_of_node"],
        r_contact=float(fx["r_contact"]), c_adh=float(fx["c_adh"]),
        rep_strength=float(fx["rep"]), adh_strength=float(fx["adh"]),
        patch_area=float(fx["patch_area"]), force_cap=float(fx["force_cap"]),
        cad_mult=fx["cad_mult"], device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    return {
        "N": int(fx["N"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
    }


def _lamellipodium_parity() -> dict:
    """Lamellipodial traction-tether: Warp vs committed HOOMD LamellipodialTractionTether."""
    from ffn_sim.warp_port.lamellipodium_warp import run_lamellipodium_tether_warp

    fx = dict(np.load(os.path.join(FIX, "lamellipodium_ref.npz")))
    got = run_lamellipodium_tether_warp(
        pos=fx["pos"], typeid=fx["typeid"], cell_of_node=fx["cell_of_node"],
        rim_cells=fx["rim_cells"], actin_typeid=int(fx["actin_typeid"]),
        mem_typeid=int(fx["mem_typeid"]), z_basal=float(fx["z_basal"]),
        basal_band=float(fx["basal_band"]), k_tether=float(fx["k_tether"]),
        force_cap=float(fx["force_cap"]), tether_radius=float(fx["tether_radius"]),
        lead_frac=float(fx["lead_frac"]), device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    return {
        "N": int(fx["N"]), "n_tethered": int(got["n_tethered"]),
        "force_rel": float(np.abs(got["force"] - fx["ref_force"]).max()) / Fscale,
    }


def _b4_differentiability() -> dict:
    """B4: Warp reverse-mode autodiff through the membrane force law w.r.t. γ_mem,
    vs closed-form analytic + finite-difference."""
    from ffn_sim.warp_port.differentiability_b4 import run_check

    return run_check(device="cpu")


def main() -> None:
    report = {
        "_about": (
            "Warp-CPU vs HOOMD-numpy bit-parity per ported piece. Reference = "
            "committed HOOMD fixtures (frozen Action / production forces). Gate "
            "thresholds: B1 kT=0 < 1e-9 / kT>0 < 1e-7; B2 host force-law < 1e-12, "
            "full warp-reduce < 1e-8 (reduction-order)."
        ),
        "B1_baoab": _baoab_parity(),
        "B2_radial_shell": _radial_shell_parity(),
        "MSHAKE_chain_constraint": _shake_parity(),
        "Fixman_metric_force": _fixman_parity(),
        "compartment_harmonic_bond": _bond_parity(),
        "compartment_harmonic_angle": _angle_parity(),
        "compartment_lj_wca": _lj_parity(),
        "DCM_node_face_contact": _dcm_contact_parity(),
        "DCM_turgor_force": _dcm_turgor_parity(),
        "DCM_cohesion_force": _dcm_cohesion_parity(),
        "lamellipodium_tether": _lamellipodium_parity(),
        "B4_differentiability": _b4_differentiability(),
    }
    with open(OUT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"wrote {OUT}")
    print(json.dumps({k: report[k] for k in
                      ("B1_baoab", "B2_radial_shell", "MSHAKE_chain_constraint",
                       "Fixman_metric_force", "compartment_harmonic_bond",
                       "compartment_harmonic_angle", "compartment_lj_wca",
                       "DCM_node_face_contact", "DCM_turgor_force",
                       "DCM_cohesion_force", "lamellipodium_tether",
                       "B4_differentiability")}, indent=2))


if __name__ == "__main__":
    main()
