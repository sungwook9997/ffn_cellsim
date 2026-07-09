"""Figure (S2): cadherin junction MATURATION now ENGAGES — the fix for the DCM_CADHERIN_MATURATION
'never engages at lit tau' finding. A sustained apposed contact, driven through CadherinBondHost with
the load-sharing cluster + per-node contact-age maturation, matures over tau_mature: its capacity
N_b(contact_age) climbs nascent->mature and its turnover collapses (LOCKS). The single-molecule model
(pre-cluster) can never reach this — it turns over at ~1/k_off forever.

Left:  N_b(t) and contact_age(t) for a sustained contact (nascent -> matured/locked).
Right: cumulative junction breaks (turnover) vs time — maturing cluster plateaus (locks); single
       molecule keeps breaking linearly. tau_mature scaled to 2 s for a legible plot (shape-identical).

Integrity: real units on both axes, tau_mature marked, both mechanisms overlaid.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ffn_sim.dcm.dcm_cadherin_host import CadherinBondHost, CadherinParams
from ffn_sim.validation.cadherin_sliding_rebinding import effective_k_off

K_TRANS = 5.84e-5
R0 = 0.5e-6
K_ON = 27.96
OUT = "ffn_sim/outputs/h_dcm_two_stage/figs"
os.makedirs(OUT, exist_ok=True)


def trace(mature, tau, steps, dt, n_nascent=4, n_mature=25, single_molecule=False, seed=1):
    P = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R0]], dtype=np.float64)
    cof = np.array([0, 1], dtype=np.int64)
    nb = 1 if single_molecule else n_mature
    p = CadherinParams(cluster=not single_molecule, mature=mature, bundle_n=float(nb),
                       n_nascent=n_nascent, tau_mature=tau, seed=seed)
    h = CadherinBondHost(cof=cof, n_cells=2, dt=dt, params=p)
    h.bonds = np.array([[0, 1]], dtype=np.int64)
    if not single_molecule:
        h.m = np.array([n_nascent if mature else n_mature], dtype=np.int64)
    t, nbt, caget, broke = [], [], [], []
    for s in range(steps):
        h._tick(P, dt)
        t.append(s * dt)
        caget.append(float(h.contact_age[0]) if mature else 0.0)
        nbt.append(int(h._nb_of_age(np.array([h.contact_age[0]]))[0]) if (mature and not single_molecule)
                   else (1 if single_molecule else n_mature))
        broke.append(h.n_broken)
    return np.array(t), np.array(nbt), np.array(caget), np.array(broke)


def main():
    tau = 2.0
    dt = 0.05 / K_ON
    steps = int(5.0 * tau / dt)
    t, nb, cage, broke_mat = trace(True, tau, steps, dt)
    _, _, _, broke_single = trace(True, tau, steps, dt, single_molecule=True)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.0))

    ax2 = axL.twinx()
    axL.plot(t, nb, color="tab:blue", lw=2.2, label="N_b(contact_age)  [engaged capacity]")
    ax2.plot(t, cage, color="tab:orange", lw=1.6, ls="--", label="contact_age [s]")
    axL.axvline(tau, color="gray", ls=":", lw=1.2)
    axL.text(tau * 1.03, 6, "τ_mature", color="gray", fontsize=9)
    axL.axhline(25, color="tab:blue", alpha=0.3, ls=":")
    axL.set_xlabel("time [s]")
    axL.set_ylabel("cluster capacity  N_b", color="tab:blue")
    ax2.set_ylabel("contact_age [s]", color="tab:orange")
    axL.set_title("A sustained contact MATURES: nascent (4) → locked (25)\nover τ_mature")
    axL.set_ylim(0, 27)
    axL.legend(loc="center right", fontsize=8)

    axR.plot(t, broke_mat, color="tab:blue", lw=2.2, label="load-sharing cluster + maturation")
    axR.plot(t, broke_single, color="crimson", lw=2.0, label="single-molecule (pre-cluster model)")
    axR.axvline(tau, color="gray", ls=":", lw=1.2)
    axR.text(tau * 1.03, broke_single[-1] * 0.5, "τ_mature", color="gray", fontsize=9)
    axR.set_xlabel("time [s]")
    axR.set_ylabel("cumulative junction breaks (turnover)")
    axR.set_title("Turnover COLLAPSES as it matures (locks);\nthe single molecule breaks forever")
    axR.legend(loc="upper left", fontsize=9)

    fig.tight_layout()
    path = f"{OUT}/dcm_cadherin_maturation_engages.png"
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")
    print(f"  maturing: N_b {nb[0]}→{nb[-1]}, contact_age→{cage[-1]:.1f}s, total breaks={broke_mat[-1]}")
    print(f"  single-molecule total breaks over {t[-1]:.0f}s = {broke_single[-1]} (never locks)")


if __name__ == "__main__":
    main()
