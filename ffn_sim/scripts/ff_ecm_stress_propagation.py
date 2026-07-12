"""FF ECM — long-range stress propagation from a contractile inclusion (ROADMAP MID, KB-1.10).

A contractile cell embedded in a fibrous matrix transmits stress far beyond its own radius — the basis of
long-range mechanical cell–cell communication (KB-1.10). We embed a spherical CONTRACTILE INCLUSION (a sphere
of collagen nodes displaced radially inward by a prescribed contraction strain ε and frozen — an Eshelby-type
inclusion, NO force tuning: ε is the physical input) at the centre of a collagen-I network, pin the far outer
boundary, relax the free interior on its own bending + segment + crosslink elasticity, and measure the radial
stress profile σ_rr(r) via the shell virial (`ecm_mechanics.stress_field_radial`). Fitting |σ|~r^−n gives the
decay exponent:
  * a linear-elastic continuum point source gives n≈2–3 (short range);
  * a sub-isostatic FIBROUS network propagates LONGER-range (smaller n, toward ~1) because tension localizes
    into fiber chains (Notbohm 2015, Rosakis, Han 2018, Wang 2014) — the emergent novelty.
We also contrast ISOTROPIC vs ALIGNED collagen: aligned fibers should channel stress even further along the
director (an anisotropic, longer-range profile). Everything unfitted; ε and the microstructure are the inputs,
the exponent EMERGES. Compares to the elastic references as acceptance oracles (per PI-2026-06-30, not KB-registered).

Run (COARSE local smoke — NON-AUTHORITATIVE):
    python -m ffn_sim.scripts.ff_ecm_stress_propagation --box 40 --steps 400 --device cpu --tag sp_smoke
Run (NATIVE gbook A5000):
    ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.ff_ecm_stress_propagation --box 80 --steps 3000 \
        --device cuda:0 --tag sp_native
Out: ffn_sim/outputs/ff/ecm_lib/{stress_propagation_<tag>.json, figs/stress_propagation_<tag>.png}
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import warp as wp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.ff import ecm_library as L
from ffn_sim.ff.ecm_mechanics import stress_field_radial, stress_field_directional
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.network_warp import _zero, link_spring_kernel
from ffn_sim.ff.motility_warp import axpy_physical_kernel

OUT = "ffn_sim/outputs/ff/ecm_lib"
CFL_SAFETY = 0.05


def _relax_inclusion(ecm, center, R_incl_um, eps, *, steps, device):
    """Contract a central spherical inclusion (radial strain ε, frozen) + pin the outer BC; relax the interior.

    Returns the relaxed node positions. The inclusion nodes are displaced inward by ε toward ``center`` and
    frozen (huge drag); the outermost shell is pinned (far-field zero displacement); every other node relaxes on
    the network's own bending + segment + crosslink springs under overdamped explicit dynamics."""
    d = device
    pos0 = np.ascontiguousarray(ecm.net.pos, np.float64).copy()
    c = np.asarray(center, float)
    rad = np.linalg.norm(pos0 - c, axis=1)
    box = np.asarray(ecm.box_hi, float) - np.asarray(ecm.box_lo, float)
    R_outer = 0.46 * float(box.min())                             # SPHERICAL far-field BC (radial symmetry —
    incl = rad < R_incl_um                                        # avoids the cube-corner truncation artifact)
    outer = rad > R_outer
    frozen = incl | outer
    pos = pos0.copy()
    pos[incl] = c + (1.0 - eps) * (pos0[incl] - c)                 # prescribed inward contraction
    N = pos.shape[0]
    gam = np.where(frozen, 1.0e18, 1.0).astype(np.float64)
    tri = np.ascontiguousarray(ecm.net.bend_triples, np.int32); nT = tri.shape[0]
    al = np.ascontiguousarray(_per_triple_alpha(ecm.net), np.float64) if nT else np.zeros(0)
    Eb = np.ascontiguousarray(ecm.links, np.int32); nB = Eb.shape[0]
    Ebk = np.ascontiguousarray(ecm.link_k, np.float64)
    Ebr = np.ascontiguousarray(ecm.link_rest, np.float64)
    seg_e = float(ecm.net.seg_rest.mean()) if ecm.net.seg_rest.size else 1.0
    bend_term = (float(ecm.net.kappa.max()) / seg_e ** 3) if nT else 1.0
    kmax = max(bend_term, float(Ebk.max()) if Ebk.size else 1.0)
    dt = CFL_SAFETY / kmax
    pos_d = wp.array(pos, dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d) if nT else None
    al_d = wp.array(al, dtype=wp.float64, device=d) if nT else None
    Eb_d = wp.array(Eb, dtype=wp.int32, ndim=2, device=d) if nB else None
    Ebk_d = wp.array(Ebk, dtype=wp.float64, device=d) if nB else None
    Ebr_d = wp.array(Ebr, dtype=wp.float64, device=d) if nB else None
    gam_d = wp.array(gam, dtype=wp.float64, device=d)
    for _ in range(steps):
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, al_d, f_d], device=d)
        if nB:
            wp.launch(link_spring_kernel, dim=nB, inputs=[pos_d, Eb_d, Ebk_d, Ebr_d, f_d], device=d)
        wp.launch(axpy_physical_kernel, dim=N, inputs=[pos_d, wp.float64(dt), gam_d, f_d], device=d)
    wp.synchronize_device(d)
    return pos_d.numpy(), dict(n_incl=int(incl.sum()), n_outer=int(outer.sum()), dt=float(dt),
                               R_incl_um=float(R_incl_um), eps=float(eps))


def run_case(S, *, box, conc, R_incl_um, eps, steps, n_shells, device, seed=7, material="collagen_I", dim=3, composite_with=None):
    """Build a fibrillar matrix at alignment S, contract the inclusion, relax, and measure σ_rr(r)."""
    spec = L.get_spec(material)
    lo, hi = [0.0, 0.0, 0.0], [box, box, box]
    rng = np.random.default_rng(seed)
    if composite_with:                                            # interpenetrating composite: aligned fibrillar + a 2nd material
        ecm = L.build_composite([dict(material=material, concentration=conc, alignment_S=float(S),
                                      director=(1.0, 0.0, 0.0), target_z=3.2), dict(material=composite_with)],
                                lo, hi, dim=dim, interlink_um=0.75, interlink_k=100.0, pin_faces=(), rng=rng)
        ecm.director = np.array([1.0, 0.0, 0.0])
    else:
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=conc, dim=dim, alignment_S=float(S),
                                    director=(1.0, 0.0, 0.0), pin_faces=(), target_z=3.2, rng=rng)
    center = 0.5 * (np.asarray(lo) + np.asarray(hi))
    t0 = time.time()
    pos, meta = _relax_inclusion(ecm, center, R_incl_um, eps, steps=steps, device=device)
    prof = stress_field_radial(ecm, pos, center, n_shells=n_shells,
                               r_min=R_incl_um * 1.3, r_max=0.32 * box)   # buffer from both BCs → clean power-law
    # DIRECTIONAL: does stress channel FARTHER along the fiber director? (n_∥ < n_⊥ for aligned)
    dvec = np.asarray(getattr(ecm, "director", (1.0, 0.0, 0.0)), float)
    dprof = stress_field_directional(ecm, pos, center, dvec, r_min=R_incl_um * 1.3, r_max=0.32 * box)
    return dict(S=float(S), S_measured=float(getattr(ecm, "S_measured", S)), n_exp=float(prof["n_exp"]),
                n_par=float(dprof["n_par"]), n_perp=float(dprof["n_perp"]), n_aniso=float(dprof["anisotropy"]),
                r=[float(x) for x in prof["r"]], sigma_rr=[float(x) for x in prof["sigma_rr"]],
                sigma_tt=[float(x) for x in prof["sigma_tt"]], pressure=[float(x) for x in prof["pressure"]],
                r_fit=[prof["r_fit_lo"], prof["r_fit_hi"]], n_nodes=int(ecm.meta["n_nodes"]),
                wall_s=time.time() - t0, **meta)


def _by_S(cases):
    """Group cases by S_target → {S: [cases]} preserving order."""
    out = {}
    for c in cases:
        out.setdefault(round(float(c["S"]), 3), []).append(c)
    return out


def plot(cases, meta, path):
    """|σ_rr|(r) log-log: per-S ensemble-mean curve + fitted exponent (mean±sd over seeds) + elastic refs → PNG."""
    coarse = meta["box"] < 60
    note = f"COARSE box={meta['box']}µm — DEV SMOKE (non-authoritative)" if coarse else f"NATIVE box={meta['box']}µm"
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    cols = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    grp = _by_S(cases)
    first_curve = None
    for i, (S, cs) in enumerate(grp.items()):
        r = np.array(cs[0]["r"])
        s_all = np.array([np.abs(c["sigma_rr"]) for c in cs])       # (seeds, shells)
        s = s_all.mean(0)
        ns = np.array([c["n_exp"] for c in cs], float); ns = ns[np.isfinite(ns)]
        Sm = float(np.mean([c["S_measured"] for c in cs]))
        ok = s > 0
        lbl = (f"S={S:g} (meas {Sm:.2f}) — n={ns.mean():.2f}" + (f"±{ns.std():.2f}" if ns.size > 1 else ""))
        ax.plot(r[ok], s[ok], "-o", color=cols[i % len(cols)], lw=2, ms=6, label=lbl)
        if s_all.shape[0] > 1:                                       # per-seed thin lines (stochastic spread)
            for sc in s_all:
                ok2 = sc > 0
                ax.plot(r[ok2], sc[ok2], "-", color=cols[i % len(cols)], lw=0.6, alpha=0.35)
        if first_curve is None and ok.any():
            first_curve = (r[ok][0], s[ok][0], r[ok][-1])
    if first_curve:
        r0, s0, r1 = first_curve
        rr = np.linspace(r0, r1, 40)
        for n_ref, ls, lbl in [(1.0, ":", "elastic ~r⁻¹ (KB-1.10 long-range)"),
                               (2.0, "--", "linear-elastic ~r⁻² (short range)")]:
            ax.plot(rr, s0 * (rr / r0) ** (-n_ref), ls, color="0.5", lw=1.4, label=lbl)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("distance from inclusion centre  r  [µm]")
    ax.set_ylabel("radial stress  |σ_rr(r)|  [Pa]")
    ax.set_title(f"FF ECM — stress propagation from a contractile inclusion  ·  {note}\n"
                 f"{meta.get('material', 'collagen_I')} ε={meta['eps']:.0%} contraction  ·  |σ_rr|~r⁻ⁿ decay "
                 f"exponent EMERGES vs the elastic references (aligned vs isotropic)  ·  tag={meta['tag']}", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140); plt.close(fig)


def plot_directional(cases, meta, path):
    """n∥(S) vs n⊥(S): the DIRECT channeling signature — along-fiber decay stays shallow (n∥ ~const) while
    across-fiber decay steepens with alignment (n⊥↑), so the stress becomes a fiber-guided waveguide."""
    coarse = meta["box"] < 60
    note = f"COARSE box={meta['box']}µm — DEV SMOKE (non-authoritative)" if coarse else f"NATIVE box={meta['box']}µm"
    grp = _by_S(cases)
    S = sorted(grp)
    npar = [np.nanmean([c["n_par"] for c in grp[s]]) for s in S]
    npar_sd = [np.nanstd([c["n_par"] for c in grp[s]]) for s in S]
    nperp = [np.nanmean([c["n_perp"] for c in grp[s]]) for s in S]
    nperp_sd = [np.nanstd([c["n_perp"] for c in grp[s]]) for s in S]
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.errorbar(S, npar, yerr=npar_sd, fmt="-o", color="#2ca02c", lw=2.2, ms=8, capsize=4,
                label="n∥  (decay ALONG the fiber director)")
    ax.errorbar(S, nperp, yerr=nperp_sd, fmt="-s", color="#d62728", lw=2.2, ms=8, capsize=4,
                label="n⊥  (decay ACROSS the fibers)")
    ax.axhline(3.0, color="0.5", ls="--", lw=1.2, label="linear-elastic n=3 (KB-1.10)")
    ax.set_xlabel("nematic alignment order  S"); ax.set_ylabel("stress-decay exponent n  (|σ_rr|~r⁻ⁿ)")
    mat = meta.get("material", "collagen_I")
    cw = meta.get("composite_with")
    gap_hi = (nperp[-1] - npar[-1]) if (S and np.isfinite(nperp[-1]) and np.isfinite(npar[-1])) else float("nan")
    if cw:                                                        # composite: the continuum short-circuits the waveguide
        sub = (f"{mat}+{cw} COMPOSITE: n∥≈n⊥≈{np.nanmean(npar+nperp):.1f} FLAT (alignment-independent) — the {cw} "
               f"continuum SHORT-CIRCUITS the fiber waveguide → isotropic (~continuum n=3)")
    elif np.isfinite(gap_hi) and gap_hi > 3.0:                    # pure fibrillar with a real waveguide
        sub = f"{mat}: stress propagates FAR along fibers (n∥ low) but dies ACROSS them (n⊥↑ with S) — a fiber WAVEGUIDE"
    else:
        sub = f"{mat}: n∥ vs n⊥ (directional stress-decay exponents)"
    ax.set_title(f"FF ECM — DIRECTIONAL stress channeling  ·  {note}\n{sub}", fontsize=9)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--box", type=float, default=40.0, help="box side [µm] (native ≥60)")
    ap.add_argument("--material", default="collagen_I", help="fibrillar ecm_library key (collagen_I / fibrin)")
    ap.add_argument("--dim", type=int, default=3, choices=(2, 3), help="2D planar sheet vs 3D bulk")
    ap.add_argument("--composite-with", default=None, help="interpenetrate the aligned fibrillar matrix with a 2nd material (e.g. matrigel)")
    ap.add_argument("--conc", type=float, default=1.5, help="fibrillar concentration [mg/mL]")
    ap.add_argument("--R-incl", type=float, default=None, help="inclusion radius [µm] (default box/8)")
    ap.add_argument("--eps", type=float, default=0.2, help="inclusion contraction strain (physical input)")
    ap.add_argument("--S", default="0.0,0.83", help="comma-separated alignment orders (isotropic vs aligned)")
    ap.add_argument("--seeds", type=int, default=1, help="matrix realizations per S (ensemble; ≥3 for a clean exponent)")
    ap.add_argument("--steps", type=int, default=400, help="relaxation steps")
    ap.add_argument("--n-shells", type=int, default=16)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--tag", default="sp")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    R_incl = a.R_incl if a.R_incl is not None else a.box / 8.0
    S_list = [float(x) for x in a.S.split(",") if x.strip()]
    wp.init(); t0 = time.time()
    coarse = a.box < 60
    if coarse:
        print(f"[!] COARSE box={a.box}µm — DEV SMOKE, NON-AUTHORITATIVE; reconfirm at NATIVE box≥80 --device cuda:0.")
    cases = []
    for S in S_list:
        for sd in range(a.seeds):
            c = run_case(S, box=a.box, conc=a.conc, R_incl_um=R_incl, eps=a.eps, steps=a.steps,
                         n_shells=a.n_shells, device=a.device, seed=7 + 13 * sd, material=a.material, dim=a.dim, composite_with=a.composite_with)
            c["seed"] = int(sd)
            print(f"[S={S:>4.2f} seed={sd}] n_iso={c['n_exp']:.2f}  n∥={c['n_par']:.2f}  n⊥={c['n_perp']:.2f}  "
                  f"(Δ={c['n_aniso']:+.2f} → {'channels ∥' if c['n_aniso'] > 0 else 'no ∥ channeling'})  "
                  f"incl={c['n_incl']} nodes={c['n_nodes']}  ({c['wall_s']:.0f}s)", flush=True)
            cases.append(c)
    n_by_S = {}
    for S in S_list:
        ns = np.array([c["n_exp"] for c in cases if round(c["S"], 3) == round(S, 3) and np.isfinite(c["n_exp"])])
        n_by_S[f"S={S:g}"] = dict(mean=float(ns.mean()) if ns.size else float("nan"),
                                  sd=float(ns.std()) if ns.size > 1 else 0.0, n_seeds=int(ns.size))
    print(f"[exponents] " + "  ".join(f"S={S:g}:n={n_by_S[f'S={S:g}']['mean']:.2f}±{n_by_S[f'S={S:g}']['sd']:.2f}"
                                      for S in S_list))
    dir_by_S = {}
    for S in S_list:
        npar = np.array([c["n_par"] for c in cases if round(c["S"], 3) == round(S, 3) and np.isfinite(c["n_par"])])
        nperp = np.array([c["n_perp"] for c in cases if round(c["S"], 3) == round(S, 3) and np.isfinite(c["n_perp"])])
        dir_by_S[f"S={S:g}"] = dict(n_par=float(npar.mean()) if npar.size else float("nan"),
                                    n_perp=float(nperp.mean()) if nperp.size else float("nan"))
    print(f"[directional] " + "  ".join(f"S={S:g}:n∥={dir_by_S[f'S={S:g}']['n_par']:.1f}/n⊥="
                                        f"{dir_by_S[f'S={S:g}']['n_perp']:.1f}" for S in S_list)
          + "  (n∥<n⊥ ⇒ stress channels along the fibers)")
    meta = dict(tag=a.tag, material=a.material, composite_with=a.composite_with, dim=a.dim, box=a.box, conc=a.conc, R_incl_um=R_incl, eps=a.eps,
                S_list=S_list, seeds=a.seeds, steps=a.steps, device=a.device,
                coarse_nonauthoritative=bool(coarse), n_exp=n_by_S, directional=dir_by_S,
                reference="linear-elastic point source n≈2–3; fibrous network longer-range n→~1 (KB-1.10; "
                          "Notbohm2015/Han2018 oracle-overlay)")
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/stress_propagation_{a.tag}.json", "w") as f:
        json.dump(dict(meta=meta, cases=cases), f, indent=2, default=float)
    plot(cases, meta, f"{a.out}/figs/stress_propagation_{a.tag}.png")
    plot_directional(cases, meta, f"{a.out}/figs/stress_propagation_directional_{a.tag}.png")
    print(f"[done] wrote {a.out}/stress_propagation_{a.tag}.json + figs/stress_propagation_{a.tag}.png "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
