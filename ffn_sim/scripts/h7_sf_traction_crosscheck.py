"""loop24b/24d seed-instability cross-check — independent estimator vs coherent FA-traction.

QUESTION (boot autonomous step 2a). The per-SF coherent FA-traction differential
(``h7_ventral_sf_traction.py:anchor_traction``, ON−OFF) is seed-unstable: the
single-seed "+131 pN" did not reproduce (loop24d: paired Δ +62±156 pN, NS; per-sample
sign-flips within one run). Is that instability **physics** (the realised bond-tension
STATE genuinely fluctuates — the myosin signal is below the per-realisation noise floor)
or an **estimator artefact** (the specific left/right anchor-bond split fabricates it)?

METHOD. Re-run the SAME same-seed paired ON/OFF design (the SF states were never
checkpointed — only scalar observables are in mband_ab/*.json), and at every matched
sample compute, on the SAME states, THREE readouts:

  1. ``coherent``  — the existing observable (anchor bonds, left/right inward split).
  2. ``midplane``  — INDEPENDENT: net axial tension Σ T·|û·axis| carried across the
                     bundle mid cross-section, over ALL backbone bonds (no left/right
                     split, different bonds, different spatial location). By force
                     balance this is the SAME physical tension the anchors feel — so if
                     it tracks ``coherent`` sample-by-sample, the instability lives in
                     the STATE, not the coherent aggregation.
  3. ``sigma(r)``  — the new Slater radial stress estimator
                     (``ffn_sim.cortex.radial_stress``) applied to the SF bundle (≈1D
                     along the axis → spherical shells sample the axial tension profile).

VERDICT LOGIC.
  - corr(coherent_diff, midplane_diff) ≈ 1 AND both ON−OFF differentials are
    noise-dominated (|mean| < 2·sem) → instability is PHYSICAL (state divergence under a
    weak contractile perturbation); a better estimator cannot fix it → the boot HALT
    "protocol fix / GPU seed-ensemble" is the right call. loop24b "generation/
    engagement-limited" stands, now corroborated by an independent estimator.
  - midplane stable while coherent flips → the coherent split is the artefact.

This is a Mac-dev cross-check; scale is reduced vs the loop24d production (flagged in the
REPORT). It tests the ESTIMATOR question, not the magnitude band.

Run (smoke):
    PYTHONPATH=. python ffn_sim/scripts/h7_sf_traction_crosscheck.py \
        --seeds 1 --equilibrate 1000 --contract 2000 --n-samples 4
Run (reduced replication):
    PYTHONPATH=. python ffn_sim/scripts/h7_sf_traction_crosscheck.py \
        --seeds 1 2 3 --equilibrate 20000 --contract 40000 --n-samples 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ffn_sim.cortex.radial_stress import fit_power_law, radial_stress_profile
from ffn_sim.scripts.h7_ventral_sf_traction import anchor_traction, build_sf_sim

_PN = 1e12          # N → pN
_UM = 1e6

OUT_DIR = Path("ffn_sim/outputs/h7/production/sf_traction_crosscheck")


def _bond_tension_state(handles):
    """(rA, rB, u, T, axis, sc) for ALL backbone bonds, from live geometry.

    Tension uses the SAME source as ``anchor_traction`` (T = k·(r − r0) with k =
    ``p.bond_k``, r0 = ``ell0``) — no force-object dependency. Positions are RAW
    (un-centred); the radial estimator centres them itself, the mid-plane uses the
    axial centre ``sc``.
    """
    snap = handles["sim"].state.get_snapshot()
    if snap.communicator.rank != 0:
        return None
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    axis = np.asarray(handles["layout"].axis, dtype=np.float64)
    k = float(handles["p"].bond_k)
    r0 = float(handles["ell0"])
    bonds = np.asarray(handles["layout"].backbone_bonds, dtype=np.int64)
    a, b = bonds[:, 0], bonds[:, 1]
    rA, rB = pos[a], pos[b]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    Lsafe = np.where(L > 0.0, L, 1.0)
    u = d / Lsafe[:, None]
    T = k * (L - r0)
    s = pos @ axis
    sc = float(np.mean(s[np.unique(bonds)]))   # bundle centre along the fiber axis
    return rA, rB, u, T, axis, sc


def midplane_axial_tension(handles) -> float:
    """INDEPENDENT estimator: net axial tension Σ T·|û·axis| across the mid cross-section.

    All backbone bonds crossing the plane {s = s_centre} (s = position·axis); the
    ``|û·axis|`` (T keeps its sign) is the same labeling-invariant convention as the
    radial estimator. Force balance ⇒ this equals the anchor inward pull, but measured
    from the fiber MIDDLE over ALL bonds — independent of the coherent left/right split.
    """
    st = _bond_tension_state(handles)
    if st is None:
        return 0.0
    rA, rB, u, T, axis, sc = st
    sA, sB = rA @ axis, rB @ axis
    crossing = (sA - sc) * (sB - sc) < 0.0
    if not crossing.any():
        return 0.0
    return float(np.sum(T[crossing] * np.abs(u[crossing] @ axis)))


def sigma_r(handles, n_radii: int = 24):
    """Slater σ(r) profile on the SF bundle (the new tool, live state)."""
    st = _bond_tension_state(handles)
    if st is None:
        return None
    rA, rB, u, T, _, _ = st
    center = 0.5 * (rA + rB).mean(axis=0)
    rAc, rBc = rA - center, rB - center
    endpoint_r = np.concatenate([np.linalg.norm(rAc, axis=1), np.linalg.norm(rBc, axis=1)])
    rmax = float(endpoint_r.max())
    radii = np.logspace(np.log10(0.05 * rmax), np.log10(rmax), n_radii)
    stress = radial_stress_profile(rAc, rBc, u, T, radii, geometry="sphere")
    return radii, stress


def run_seed(seed: int, args) -> dict:
    import hoomd

    def _dev():
        return (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
                else hoomd.device.CPU(notice_level=0))

    common = dict(
        n_fil=args.n_filaments, fiber_length=args.fiber_length_um * 1e-6,
        bundle_radius=args.bundle_radius_nm * 1e-9, seed=seed,
        with_myosin=True, n_motors=args.n_motors,
        mband=False, sarcomeric=True, n_sarcomeres=args.n_sarcomeres,
        overlap_frac=0.3,
    )
    h_off = build_sf_sim(device=_dev(), myosin_force_scale=0.0, **common)
    h_on = build_sf_sim(device=_dev(), myosin_force_scale=1.0, **common)
    h_off["sim"].run(args.equilibrate)
    h_on["sim"].run(args.equilibrate)
    chunk = max(1, args.contract // args.n_samples)

    coh_diff, mid_diff, coh_on, coh_off, mid_on, mid_off = [], [], [], [], [], []
    for _ in range(args.n_samples):
        h_off["sim"].run(chunk)
        h_on["sim"].run(chunk)
        _, _, c_off = anchor_traction(h_off)
        _, _, c_on = anchor_traction(h_on)
        m_off = midplane_axial_tension(h_off)
        m_on = midplane_axial_tension(h_on)
        coh_on.append(c_on * _PN); coh_off.append(c_off * _PN)
        mid_on.append(m_on * _PN); mid_off.append(m_off * _PN)
        coh_diff.append((c_on - c_off) * _PN)
        mid_diff.append((m_on - m_off) * _PN)

    coh_diff = np.array(coh_diff); mid_diff = np.array(mid_diff)
    # σ(r) profile (ON & OFF) from the final state — the new-tool panel.
    sig_on = sigma_r(h_on); sig_off = sigma_r(h_off)
    return {
        "seed": seed,
        "coherent_diff_pN": coh_diff.tolist(),
        "midplane_diff_pN": mid_diff.tolist(),
        "coherent_on_pN": coh_on, "coherent_off_pN": coh_off,
        "midplane_on_pN": mid_on, "midplane_off_pN": mid_off,
        "coherent_diff_mean_pN": float(coh_diff.mean()),
        "coherent_diff_sem_pN": float(coh_diff.std() / max(1, len(coh_diff) ** 0.5)),
        "midplane_diff_mean_pN": float(mid_diff.mean()),
        "midplane_diff_sem_pN": float(mid_diff.std() / max(1, len(mid_diff) ** 0.5)),
        "engaged_mean": float(np.mean([h_on["myosin_action"].n_engaged])),
        "_sigma_on": sig_on, "_sigma_off": sig_off,
    }


def _verdict(results: list[dict]) -> dict:
    coh = np.concatenate([np.asarray(r["coherent_diff_pN"]) for r in results])
    mid = np.concatenate([np.asarray(r["midplane_diff_pN"]) for r in results])
    corr = float(np.corrcoef(coh, mid)[0, 1]) if coh.size > 1 else float("nan")
    coh_seed = np.array([r["coherent_diff_mean_pN"] for r in results])
    mid_seed = np.array([r["midplane_diff_mean_pN"] for r in results])
    coh_ens_mean = float(coh_seed.mean())
    coh_ens_sem = float(coh_seed.std() / max(1, len(coh_seed) ** 0.5))
    mid_ens_mean = float(mid_seed.mean())
    mid_ens_sem = float(mid_seed.std() / max(1, len(mid_seed) ** 0.5))
    coh_noise = abs(coh_ens_mean) < 2 * coh_ens_sem
    mid_noise = abs(mid_ens_mean) < 2 * mid_ens_sem
    # PRIMARY criterion = is the ON−OFF differential noise-dominated in the INDEPENDENT
    # estimator too? If both are, no estimator yields a robust per-SF traction → the boot
    # HALT (protocol fix / seed-ensemble) stands, independently corroborated. Correlation
    # is a SECONDARY diagnostic of WHY (shared state vs both-near-zero-noise).
    cs = "high" if corr > 0.6 else "weak" if corr > 0.2 else "absent"
    if coh_noise and mid_noise:
        verdict = ("BOOT HALT CORROBORATED — the ON−OFF differential is noise-dominated "
                   "(|mean|<2·sem across seeds) for BOTH the coherent FA-traction AND the "
                   "independent mid-plane estimator. No robust per-SF traction exists by "
                   "either readout → it is NOT fixable by a better estimator; the boot HALT "
                   "(protocol fix / GPU seed-ensemble) is the right call, and loop24b "
                   "'generation/engagement-limited' stands on an independent estimator. "
                   "Per-sample corr={:.2f} ({} — {}).".format(
                       corr, cs,
                       "same fluctuating state" if corr > 0.6
                       else "both ~0 ± noise, weakly coupled"))
    elif (not mid_noise) and coh_noise:
        verdict = ("DIVERGENT — the independent mid-plane differential is robustly non-zero "
                   "while the coherent one is noise-dominated (corr={:.2f}); the coherent "
                   "left/right anchor split may UNDER-report or destabilise the signal. "
                   "Worth a protocol look before quoting either.".format(corr))
    elif (not coh_noise) and (not mid_noise):
        verdict = ("BOTH ROBUST — both estimators give a non-noise differential (corr={:.2f}); "
                   "a per-SF traction number may be recoverable. Re-examine vs loop24d "
                   "(this is reduced scale).".format(corr))
    else:
        verdict = ("MIXED — coherent noise-dominated={}, midplane noise-dominated={}, "
                   "corr={:.2f}. See per-seed table.".format(coh_noise, mid_noise, corr))
    return {
        "corr_persample_coherent_vs_midplane": corr,
        "coherent_ensemble_mean_pN": coh_ens_mean, "coherent_ensemble_sem_pN": coh_ens_sem,
        "midplane_ensemble_mean_pN": mid_ens_mean, "midplane_ensemble_sem_pN": mid_ens_sem,
        "coherent_noise_dominated": bool(coh_noise),
        "midplane_noise_dominated": bool(mid_noise),
        "verdict": verdict,
    }


def _figure(results: list[dict], verdict: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(results)))

    # Panel 1: per-sample coherent vs midplane differential (should lie on a line).
    allc, allm = [], []
    for r, c in zip(results, colors):
        cd = np.asarray(r["coherent_diff_pN"]); md = np.asarray(r["midplane_diff_pN"])
        allc.append(cd); allm.append(md)
        ax[0].scatter(cd, md, s=22, color=c, alpha=0.8, label=f"seed {r['seed']}")
    allc = np.concatenate(allc); allm = np.concatenate(allm)
    lim = max(np.abs(np.r_[allc, allm]).max(), 1.0) * 1.1
    ax[0].plot([-lim, lim], [-lim, lim], "k--", lw=1, alpha=0.6, label="identity")
    ax[0].axhline(0, color="grey", lw=0.6); ax[0].axvline(0, color="grey", lw=0.6)
    ax[0].set_xlim(-lim, lim); ax[0].set_ylim(-lim, lim)
    ax[0].set_xlabel("coherent FA-traction differential ON−OFF [pN]")
    ax[0].set_ylabel("INDEPENDENT mid-plane axial tension ON−OFF [pN]")
    ax[0].set_title(f"independent vs coherent (per sample)\nPearson r = "
                    f"{verdict['corr_persample_coherent_vs_midplane']:.3f}")
    ax[0].legend(fontsize=7); ax[0].grid(True, alpha=0.25)

    # Panel 2: per-seed differential — both estimators, with SEM.
    seeds = [r["seed"] for r in results]
    x = np.arange(len(seeds)); w = 0.36
    ax[1].bar(x - w / 2, [r["coherent_diff_mean_pN"] for r in results], w,
              yerr=[r["coherent_diff_sem_pN"] for r in results], capsize=3,
              label="coherent", color="#1f77b4")
    ax[1].bar(x + w / 2, [r["midplane_diff_mean_pN"] for r in results], w,
              yerr=[r["midplane_diff_sem_pN"] for r in results], capsize=3,
              label="midplane (indep.)", color="#d62728")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xticks(x); ax[1].set_xticklabels([f"s{s}" for s in seeds])
    ax[1].set_xlabel("seed"); ax[1].set_ylabel("ON−OFF differential [pN]")
    ax[1].set_title("per-seed ON−OFF differential (coherent vs independent)\n"
                    f"coherent {verdict['coherent_ensemble_mean_pN']:+.0f}±"
                    f"{verdict['coherent_ensemble_sem_pN']:.0f}, "
                    f"midplane {verdict['midplane_ensemble_mean_pN']:+.0f}±"
                    f"{verdict['midplane_ensemble_sem_pN']:.0f} pN")
    ax[1].legend(fontsize=7); ax[1].grid(True, axis="y", alpha=0.25)

    # Panel 3: the new σ(r) tool on the SF bundle (seed 1, ON vs OFF).
    r0 = results[0]
    if r0.get("_sigma_on") is not None:
        rr, son = r0["_sigma_on"]; _, soff = r0["_sigma_off"]
        ax[2].loglog(rr * _UM, np.abs(son), "o-", ms=3, color="#1f77b4", label="|σ(r)| ON")
        ax[2].loglog(rr * _UM, np.abs(soff), "s-", ms=3, color="#999999", label="|σ(r)| OFF")
        ax[2].set_xlabel("shell radius r [µm]"); ax[2].set_ylabel("|σ(r)| [Pa]")
        ax[2].set_title(f"Slater σ(r) on SF bundle (seed {r0['seed']})\n"
                        "new radial_stress tool, live state")
        ax[2].legend(fontsize=7); ax[2].grid(True, which="both", alpha=0.25)
    fig.suptitle("H.7 seed-instability cross-check — coherent FA-traction vs independent "
                 "estimators", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    # n_filaments = filaments PER CROSS-SECTION; total = n_sarcomeres·2·n_fil. loop24d
    # (mband_ab/run_ab.sh) used 6 → 36 total filaments / 12 FA anchors. Matching it is
    # REQUIRED for a valid replication (36 here would give a 6× bundle, different physics).
    ap.add_argument("--n-filaments", type=int, default=6)
    ap.add_argument("--n-motors", type=int, default=16)
    ap.add_argument("--n-sarcomeres", type=int, default=3)
    ap.add_argument("--fiber-length-um", type=float, default=6.0)
    ap.add_argument("--bundle-radius-nm", type=float, default=400.0)
    ap.add_argument("--equilibrate", type=int, default=20000)
    ap.add_argument("--contract", type=int, default=40000)
    ap.add_argument("--n-samples", type=int, default=20)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--tag", type=str, default="reduced")
    args = ap.parse_args()

    print(f"H.7 SF traction cross-check — seeds {args.seeds}, sarcomeric n_fil={args.n_filaments} "
          f"n_motors={args.n_motors}, equil={args.equilibrate} contract={args.contract} "
          f"nsamp={args.n_samples} [{args.device}]", flush=True)
    results = []
    for sd in args.seeds:
        r = run_seed(sd, args)
        results.append(r)
        print(f"  seed {sd}: coherent ΔON−OFF = {r['coherent_diff_mean_pN']:+.1f}±"
              f"{r['coherent_diff_sem_pN']:.1f} pN | midplane = "
              f"{r['midplane_diff_mean_pN']:+.1f}±{r['midplane_diff_sem_pN']:.1f} pN", flush=True)

    verdict = _verdict(results)
    print(f"\n  corr(coherent, midplane) per-sample = "
          f"{verdict['corr_persample_coherent_vs_midplane']:.3f}", flush=True)
    print(f"  VERDICT: {verdict['verdict']}", flush=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = OUT_DIR / f"sf_traction_crosscheck_{args.tag}.png"
    _figure(results, verdict, fig_path)
    # serialise σ(r) profiles too (figure fully regenerable from JSON without re-running)
    serial = []
    for r in results:
        rr = {k: v for k, v in r.items() if not k.startswith("_sigma")}
        son, soff = r.get("_sigma_on"), r.get("_sigma_off")
        if son is not None:
            rr["sigma_radii_um"] = (son[0] * _UM).tolist()
            rr["sigma_on_Pa"] = son[1].tolist()
            rr["sigma_off_Pa"] = soff[1].tolist()
        serial.append(rr)
    out = {
        "config": {
            "seeds": args.seeds, "n_filaments": args.n_filaments, "n_motors": args.n_motors,
            "n_sarcomeres": args.n_sarcomeres, "equilibrate": args.equilibrate,
            "contract": args.contract, "n_samples": args.n_samples, "device": args.device,
            "scale_note": "REDUCED vs loop24d production (equil/contract 120k, nsamp 20); "
                          "Mac-dev estimator cross-check, not a magnitude-band run.",
        },
        "per_seed": serial,
        "verdict": verdict,
    }
    (OUT_DIR / f"crosscheck_{args.tag}.json").write_text(json.dumps(out, indent=2))
    print(f"  json → {OUT_DIR / f'crosscheck_{args.tag}.json'}\n  fig  → {fig_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
