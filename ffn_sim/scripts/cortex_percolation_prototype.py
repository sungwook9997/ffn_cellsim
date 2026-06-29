"""Cortex percolation PROTOTYPE — pure geometry + graph, no HOOMD.

Find the construction regime (bimodal length distribution, filament count,
crosslinker reach, n_xl/n_fil) that turns the fragmented mesh (z=1.3, giant 7%,
56% same-filament staples — outputs/h3/production/cortex_network.json) into a
CONNECTED PERCOLATED mesh: target z=3.0-3.5, giant>=0.9, L/lc>=5.9.

Mirrors the real construction geometry (Marsaglia sphere shell + tangent-plane
isotropic filaments + bimodal-exponential length: Arp2/3-short ~90% count +
formin-long ~10% = the connecting backbone, Fritzsche 2016/2017) and the real
crosslink rule (BRIDGE DIFFERENT FILAMENTS — forbid same-filament staples).
The bridges are SEEDED at construction (physiological/adhered baseline: the
cortex STARTS connected, not emergently), then the giant component is measured.

This is a tuning harness only — the production code lands in cortex.py +
crosslinkers.py once a regime is confirmed here.

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.cortex_percolation_prototype
"""
from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# Reuse the EXACT placement primitives the real cortex uses.
from ffn_sim.archive.hoomd_legacy.cortex.cortex import _sample_sphere_surface, _tangent_plane_basis


def place_bimodal_filaments(
    *, n_fil, R_cell, ell0, formin_fraction, L_short, L_long,
    L_min_beads=2, rng,
):
    """Place n_fil disordered-isotropic filaments on the R_cell shell with a
    bimodal-EXPONENTIAL length distribution.

    - Each filament is formin (long) with prob ``formin_fraction``, else
      Arp2/3 (short).  Length within each subpopulation ~ Exponential(mean),
      mean = L_long (formin) or L_short (Arp2/3).  Quantized to ell0 multiples
      -> n_beads = round(L/ell0)+1, clamped >= L_min_beads.
    - Center uniform on sphere (Marsaglia); tangent in local plane, uniform
      azimuth (isotropic, S->0).  Beads laid at ell0 spacing about the CoM.

    Returns flat bead positions (n_beads_total,3), per-bead filament index,
    per-filament bead count, per-filament length, is_formin mask.
    """
    is_formin = rng.random(n_fil) < formin_fraction
    means = np.where(is_formin, L_long, L_short)
    L = rng.exponential(means)
    n_beads = np.clip(np.round(L / ell0).astype(int) + 1, L_min_beads, None)
    L_real = (n_beads - 1) * ell0

    centers = _sample_sphere_surface(rng, n_fil, R_cell)
    normals = centers / R_cell
    e1, e2 = _tangent_plane_basis(normals)
    phi = rng.uniform(0.0, 2.0 * math.pi, n_fil)
    tang = np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
    tang /= np.linalg.norm(tang, axis=1, keepdims=True).clip(1e-30)

    starts = np.zeros(n_fil, dtype=np.int64)
    starts[1:] = np.cumsum(n_beads[:-1])
    n_total = int(n_beads.sum())
    pos = np.empty((n_total, 3), dtype=np.float64)
    fil_idx = np.empty(n_total, dtype=np.int64)
    for f in range(n_fil):
        N = int(n_beads[f]); s = int(starts[f])
        off = (np.arange(N) - 0.5 * (N - 1)) * ell0
        pos[s:s + N] = centers[f] + off[:, None] * tang[f]
        fil_idx[s:s + N] = f
    # Project every bead onto the 200 nm cortex shell band so even the long
    # formin backbone curves ALONG the membrane (radial drift stays < 200 nm
    # for any length; KU-3.17).  Radius drawn uniformly in [R-200nm, R].
    cortex_thickness = 200e-9
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    r_target = R_cell - rng.uniform(0.0, cortex_thickness, (n_total, 1))
    pos = pos / r.clip(1e-30) * r_target
    return pos, fil_idx, n_beads, L_real, is_formin


def seed_bridges_random_anchor(pos, fil_idx, *, n_xl, reach, rng):
    """[v1] Random-anchor seeding (mirrors the prior fragmented build)."""
    tree = cKDTree(pos)
    n_beads = pos.shape[0]
    bridges = []
    n_homeless = 0
    for _ in range(n_xl):
        a = int(rng.integers(0, n_beads))
        fa = int(fil_idx[a])
        nbr = tree.query_ball_point(pos[a], r=reach)
        cand = [j for j in nbr if fil_idx[j] != fa]
        if not cand:
            n_homeless += 1
            continue
        cand = np.asarray(cand)
        d = np.linalg.norm(pos[cand] - pos[a], axis=1)
        b = int(cand[np.argmin(d)])
        bridges.append((fa, int(fil_idx[b])))
    return bridges, n_homeless


def seed_connected_mesh(pos, fil_idx, n_fil, *, z_struct, bundle_mult, reach, rng):
    """Build the connected mesh DIRECTLY at the literature connectivity set-point.

    Two-stage, both BRIDGE-DIFFERENT-FILAMENT:
      Stage 1 (structural percolation): give every filament a distinct-neighbour
      degree ~ z_struct (target 3.0-3.5, the rigidity-percolation coordination,
      Kim 2007 / Kadzik-Munro 2026) by linking to its z_struct nearest distinct
      different-filament partners (capped so the mean DISTINCT degree ~ z_struct,
      not 2*z_struct).
      Stage 2 (bundling): add (bundle_mult-1) extra parallel crosslinks per
      connected pair so crosslinks-per-filament L/lc = z_struct*bundle_mult >= 5.9
      (Flormann 2024 stiffness-from-bundling; Head 2003 L/lc).

    Returns (bridges_all, degree_struct array, homeless count).  bridges_all
    counts every crosslink (incl. bundling); the structural graph dedups pairs.
    """
    tree = cKDTree(pos)
    beads_of = [np.flatnonzero(fil_idx == f) for f in range(n_fil)]
    deg = np.zeros(n_fil, dtype=np.int64)
    # Per-filament coordination cap: fractional z_struct realised by randomising
    # each cap between floor and ceil so the ENSEMBLE mean coordination = z_struct
    # (the Kim 2007 / Kadzik-Munro physiological z), not pinned to an integer.
    lo = int(math.floor(z_struct))
    cap = lo + (rng.random(n_fil) < (z_struct - lo)).astype(np.int64)
    edges = set()
    bridges = []
    order = rng.permutation(n_fil)
    for f in order:
        if deg[f] >= cap[f]:
            continue
        my_beads = beads_of[f]
        nbr_lists = tree.query_ball_point(pos[my_beads], r=reach)
        best = {}
        for k, nbrs in enumerate(nbr_lists):
            for j in nbrs:
                fj = int(fil_idx[j])
                if fj == f:
                    continue
                d = float(np.linalg.norm(pos[j] - pos[my_beads[k]]))
                if fj not in best or d < best[fj]:
                    best[fj] = d
        for fp in sorted(best, key=best.get):
            if deg[f] >= cap[f]:
                break
            if deg[fp] >= cap[fp]:          # HARD cap both endpoints -> mean degree ~ z_struct
                continue
            key = (min(f, fp), max(f, fp))
            if key in edges:
                continue
            edges.add(key)
            bridges.append((f, fp))
            deg[f] += 1
            deg[fp] += 1
    n_homeless = int((deg == 0).sum())
    # Stage 2: bundling — duplicate each structural edge (bundle_mult-1) times.
    if bundle_mult > 1:
        struct_edges = list(edges)
        for _ in range(int(bundle_mult) - 1):
            for (a, b) in struct_edges:
                bridges.append((a, b))
    return bridges, deg, n_homeless


def seed_bridges_per_filament(pos, fil_idx, n_fil, *, z_target, reach, rng):
    """Seed bridges PER FILAMENT so the connected mesh is built directly
    (physiological adhered baseline: cortex STARTS percolated).

    For every filament f: gather different-filament beads within ``reach`` of
    any of f's beads, rank candidate PARTNER FILAMENTS by nearest bead-distance,
    and add bridges to the ``z_target`` nearest DISTINCT partners that f is not
    already bridged to.  BRIDGE-DIFFERENT-FILAMENT is enforced (partners != f);
    a filament with no partner in reach stays homeless (singleton).

    This guarantees degree ~ z_target for every non-isolated filament, so the
    giant component spans (giant -> 1) at exactly the literature connectivity
    set-point (Kim 2007; Kadzik-Munro 2026) instead of leaving the random tail
    isolated.

    Returns list of (f_a,f_b) bridges + homeless (isolated-filament) count.
    """
    tree = cKDTree(pos)
    bead_starts = {}
    for f in range(n_fil):
        bead_starts.setdefault(f, [])
    beads_of = [np.flatnonzero(fil_idx == f) for f in range(n_fil)]
    edges = set()
    bridges = []
    n_homeless = 0
    for f in range(n_fil):
        my_beads = beads_of[f]
        nbr_lists = tree.query_ball_point(pos[my_beads], r=reach)
        # nearest bead-distance to each candidate partner filament
        best = {}
        for k, nbrs in enumerate(nbr_lists):
            for j in nbrs:
                fj = int(fil_idx[j])
                if fj == f:
                    continue
                d = float(np.linalg.norm(pos[j] - pos[my_beads[k]]))
                if fj not in best or d < best[fj]:
                    best[fj] = d
        if not best:
            n_homeless += 1
            continue
        partners = sorted(best, key=best.get)
        added = 0
        for fp in partners:
            if added >= z_target:
                break
            key = (min(f, fp), max(f, fp))
            if key in edges:
                continue
            edges.add(key)
            bridges.append((f, fp))
            added += 1
    return bridges, n_homeless


def measure(n_fil, n_beads, bridges):
    """z (mean graph degree = distinct neighbors), giant fraction, L/lc."""
    if not bridges:
        return dict(z=0.0, z_xl=0.0, giant=0.0, n_bridges=0, n_components=n_fil,
                    L_over_lc=0.0)
    rows = np.array([b[0] for b in bridges] + [b[1] for b in bridges])
    cols = np.array([b[1] for b in bridges] + [b[0] for b in bridges])
    # distinct-neighbor degree: dedup edges for graph degree
    g = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n_fil, n_fil))
    g.data[:] = 1.0
    g.sum_duplicates()
    n_comp, labels = connected_components(g, directed=False)
    sizes = np.bincount(labels, minlength=n_fil)
    giant = int(sizes.max())
    # degree counting DISTINCT neighbors
    deg_distinct = np.asarray((g > 0).sum(axis=1)).ravel()
    z = float(deg_distinct.mean())
    # crosslinks-per-filament (counts every bridge end, incl. bundling)
    ends = np.zeros(n_fil)
    for a, b in bridges:
        ends[a] += 1; ends[b] += 1
    z_xl = float(ends.mean())
    # L/lc: crosslinks-per-filament (lc = mean spacing -> L/lc = #crosslinks)
    # averaged over filaments that carry >=1 crosslink (in-mesh filaments).
    in_mesh = ends > 0
    L_over_lc = float(ends[in_mesh].mean()) if in_mesh.any() else 0.0
    return dict(z=z, z_xl=z_xl, giant=giant / n_fil, n_bridges=len(bridges),
                n_components=int((sizes > 0).sum()), L_over_lc=L_over_lc,
                homeless_frac=None)


def run_one(*, n_fil, R_cell, ell0, formin_fraction, L_short, L_long,
            z_target, reach, seed):
    rng = np.random.default_rng(seed)
    pos, fil_idx, n_beads, L_real, is_formin = place_bimodal_filaments(
        n_fil=n_fil, R_cell=R_cell, ell0=ell0, formin_fraction=formin_fraction,
        L_short=L_short, L_long=L_long, rng=rng)
    bridges, n_homeless = seed_bridges_per_filament(
        pos, fil_idx, n_fil, z_target=z_target, reach=reach, rng=rng)
    res = measure(n_fil, n_beads, bridges)
    # implied mesoscale mesh size: sqrt(area / total contour) ~ pore length
    total_contour = float(L_real.sum())
    area = 4 * math.pi * R_cell ** 2
    res.update(
        n_fil=n_fil, n_xl=len(bridges), n_homeless=n_homeless,
        homeless_frac=n_homeless / max(1, n_fil),
        mean_beads=float(n_beads.mean()), max_beads=int(n_beads.max()),
        mean_L_um=float(L_real.mean() * 1e6),
        formin_frac_real=float(is_formin.mean()),
        areal_contour_per_nm=total_contour / area * 1e-9,
        mesh_xi_nm=math.sqrt(area / total_contour) * 1e9,
    )
    return res


def run_lock(*, n_fil, R_cell, ell0, formin_fraction, L_short, L_long,
             z_struct, bundle_mult, reach, seed):
    rng = np.random.default_rng(seed)
    pos, fil_idx, n_beads, L_real, is_formin = place_bimodal_filaments(
        n_fil=n_fil, R_cell=R_cell, ell0=ell0, formin_fraction=formin_fraction,
        L_short=L_short, L_long=L_long, rng=rng)
    bridges, deg, n_homeless = seed_connected_mesh(
        pos, fil_idx, n_fil, z_struct=z_struct, bundle_mult=bundle_mult,
        reach=reach, rng=rng)
    res = measure(n_fil, n_beads, bridges)  # z=distinct degree, L/lc=crosslinks/fil
    res.update(n_homeless=n_homeless, homeless_frac=n_homeless / n_fil,
               cortex_beads=int(n_beads.sum()),
               max_drift_nm=float((R_cell - np.linalg.norm(pos, axis=1)).max() * 1e9))
    return res


def make_viz(cfg, seed, out_png):
    """Render the LOCKED-recipe connected mesh vs the prior fragmented baseline."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    rng = np.random.default_rng(seed)
    pos, fil_idx, n_beads, L_real, is_formin = place_bimodal_filaments(
        n_fil=cfg["n_fil"], R_cell=cfg["R_cell"], ell0=cfg["ell0"],
        formin_fraction=cfg["formin_fraction"], L_short=cfg["L_short"],
        L_long=cfg["L_long"], rng=rng)
    bridges, deg, n_homeless = seed_connected_mesh(
        pos, fil_idx, cfg["n_fil"], z_struct=cfg["z_struct"],
        bundle_mult=cfg["bundle_mult"], reach=cfg["reach"], rng=rng)
    res = measure(cfg["n_fil"], n_beads, bridges)
    n_fil = cfg["n_fil"]

    # filament graph -> connected components -> giant label
    rows = np.array([b[0] for b in bridges] + [b[1] for b in bridges])
    cols = np.array([b[1] for b in bridges] + [b[0] for b in bridges])
    g = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n_fil, n_fil))
    g.data[:] = 1.0; g.sum_duplicates()
    _, labels = connected_components(g, directed=False)
    sizes = np.bincount(labels, minlength=n_fil)
    giant_label = int(np.argmax(sizes))

    bp = pos * 1e6
    in_giant = labels[fil_idx] == giant_label
    has_bridge = sizes[labels[fil_idx]] > 1
    is_formin_bead = is_formin[fil_idx]

    fig = plt.figure(figsize=(20, 6.8))
    # Panel A: connectivity (giant vs clusters vs isolated)
    axA = fig.add_subplot(1, 3, 1, projection="3d")
    axA.scatter(bp[~has_bridge, 0], bp[~has_bridge, 1], bp[~has_bridge, 2],
                s=2, c="lightgray", label="isolated filaments")
    axA.scatter(bp[has_bridge & ~in_giant, 0], bp[has_bridge & ~in_giant, 1],
                bp[has_bridge & ~in_giant, 2], s=4, c="tab:orange",
                label="small clusters")
    axA.scatter(bp[in_giant, 0], bp[in_giant, 1], bp[in_giant, 2], s=3,
                c="tab:blue",
                label=f"giant component ({res['giant']*100:.0f}%)")
    axA.set_title(f"A  filaments by connectivity\n"
                  f"z(distinct)={res['z']:.2f}  giant={res['giant']*100:.0f}%  "
                  f"L/lc={res['L_over_lc']:.1f}")
    axA.legend(fontsize=7, loc="upper left")
    for ax in (axA,):
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)"); ax.set_zlabel("z (µm)")

    # Panel B: bimodal architecture (Arp2/3 short vs formin long backbone)
    axB = fig.add_subplot(1, 3, 2, projection="3d")
    axB.scatter(bp[~is_formin_bead, 0], bp[~is_formin_bead, 1],
                bp[~is_formin_bead, 2], s=2, c="silver", alpha=0.5,
                label=f"Arp2/3 short ({(1-is_formin.mean())*100:.0f}%)")
    axB.scatter(bp[is_formin_bead, 0], bp[is_formin_bead, 1],
                bp[is_formin_bead, 2], s=6, c="crimson",
                label=f"formin backbone ({is_formin.mean()*100:.0f}%)")
    axB.set_title(f"B  bimodal length architecture\n"
                  f"mean L={L_real.mean()*1e6:.2f}µm  max L={L_real.max()*1e6:.1f}µm  "
                  f"(formin = spanning backbone)")
    axB.legend(fontsize=7, loc="upper left")
    axB.set_xlabel("x (µm)"); axB.set_ylabel("y (µm)"); axB.set_zlabel("z (µm)")

    # Panel C: the bridges (mesh edges), drawn at filament centroids
    axC = fig.add_subplot(1, 3, 3, projection="3d")
    axC.scatter(bp[:, 0], bp[:, 1], bp[:, 2], s=1, c="lightgray", alpha=0.25)
    cent = np.array([bp[fil_idx == f].mean(axis=0) for f in range(n_fil)])
    # draw a representative subsample of structural edges to keep it readable
    struct = list({(min(a, b), max(a, b)) for a, b in bridges})
    rng2 = np.random.default_rng(0)
    sub = [struct[i] for i in rng2.choice(len(struct),
                                          size=min(2500, len(struct)), replace=False)]
    segs = [[cent[a].tolist(), cent[b].tolist()] for a, b in sub]
    axC.add_collection3d(Line3DCollection(segs, colors="seagreen",
                                          linewidths=0.4, alpha=0.5))
    axC.set_title(f"C  crosslink bridges (mesh edges)\n"
                  f"{len(struct)} distinct filament-pairs bridged  "
                  f"(+bundling -> {res['n_bridges']} crosslinks)")
    axC.set_xlabel("x (µm)"); axC.set_ylabel("y (µm)"); axC.set_zlabel("z (µm)")

    verdict = "CONNECTED SPANNING MESH" if res["giant"] >= 0.9 else "FRAGMENTED"
    fig.suptitle(
        f"Cortex construction REBUILD — {verdict}    "
        f"[BEFORE: z=1.3, giant 7%, 56% same-filament staples  →  "
        f"AFTER: z={res['z']:.2f}, giant {res['giant']*100:.0f}%, L/lc={res['L_over_lc']:.1f}]\n"
        f"bimodal length + bridge-different-filament + bundling + adhered-baseline "
        f"seeding   (n_fil={n_fil}, prototype geometry)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return res, verdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R-cell", type=float, default=7.5e-6)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--lock", action="store_true",
                    help="confirm the LOCKED recipe (distinct z_struct + bundling)")
    ap.add_argument("--viz", action="store_true",
                    help="render the locked-recipe connected mesh to a PNG")
    args = ap.parse_args()

    if args.viz:
        from pathlib import Path
        cfg = dict(n_fil=2000, R_cell=args.R_cell, ell0=0.5e-6,
                   formin_fraction=0.12, L_short=0.5e-6, L_long=5.0e-6,
                   z_struct=3.7, bundle_mult=2, reach=594e-9)
        out = (Path(__file__).resolve().parents[1] / "outputs" / "h3" / "figs"
               / "cortex_connected_mesh_prototype.png")
        out.parent.mkdir(parents=True, exist_ok=True)
        res, verdict = make_viz(cfg, seed=1, out_png=out)
        print(f"[viz] z(distinct)={res['z']:.2f} giant={res['giant']*100:.0f}% "
              f"L/lc={res['L_over_lc']:.1f} => {verdict}")
        print(f"[viz] {out}")
        return 0

    if args.lock:
        print("LOCKED recipe: bimodal(formin 0.12, L_short 0.5um, L_long 5um), "
              "z_struct=3.7 frac-cap + bundle x2, reach=sqrt(A/n_fil), "
              "shell-projected.\nTargets: z(distinct)=3.0-3.5, giant>=0.9, "
              "L/lc>=5.9.   budget = cortex_beads + 2*n_xl <= 15000\n")
        area = 4 * math.pi * args.R_cell ** 2
        for n_fil in (1000, 1200, 1500, 2000):
            reach = math.sqrt(area / n_fil)
            cfg = dict(n_fil=n_fil, R_cell=args.R_cell, ell0=0.5e-6,
                       formin_fraction=0.12, L_short=0.5e-6, L_long=5.0e-6,
                       z_struct=3.7, bundle_mult=2, reach=reach)
            zs, gis, lls, nxl, cb = [], [], [], [], []
            for s in range(max(3, args.seeds)):
                r = run_lock(seed=1 + s, **cfg)
                zs.append(r["z"]); gis.append(r["giant"]); lls.append(r["L_over_lc"])
                nxl.append(r["n_bridges"]); cb.append(r.get("cortex_beads", 0))
            budget = int(np.mean(cb)) + 2 * int(np.mean(nxl))
            ok = "FIT" if budget <= 15000 else "OVER"
            print(f"  n_fil={n_fil:>4} reach={reach*1e9:>3.0f}nm | z={np.mean(zs):.2f} "
                  f"giant={np.mean(gis)*100:>4.1f}% L/lc={np.mean(lls):.1f} | "
                  f"n_xl={int(np.mean(nxl))} cortex_beads={int(np.mean(cb))} "
                  f"=> particles={budget} [{ok}]")
        return 0

    print("Flormann physiological areal-contour density: 0.06-0.10 nm^-1.  "
          "Target: z=3.0-3.5, giant>=0.9, L/lc>=5.9.")
    print("Per-filament target-degree seeding (build the connected state "
          "directly = adhered baseline).  Sweep reach to find where every "
          "filament finds partners.\n")

    grids = []
    for n_fil in (1500, 2000, 3000):
        for L_long in (3.0e-6, 5.0e-6):
            for z_target in (3, 4):
                for reach in (380e-9, 600e-9, 900e-9, 1200e-9):
                    grids.append(dict(
                        n_fil=n_fil, formin_fraction=0.12,
                        L_short=0.5e-6, L_long=L_long,
                        z_target=z_target, reach=reach))
    hdr = (f"{'n_fil':>5} {'L_lng':>6} {'zT':>2} {'reach':>5} | "
           f"{'z':>4} {'giant':>5} {'L/lc':>4} {'home%':>5} {'xi_nm':>6} "
           f"{'n_xl/f':>6} {'mbeads':>6}")
    print(hdr); print("-" * len(hdr))
    for g in grids:
        agg = {k: [] for k in ("z", "giant", "L_over_lc", "homeless_frac",
                               "mesh_xi_nm", "n_xl", "mean_beads")}
        for s in range(args.seeds):
            r = run_one(R_cell=args.R_cell, ell0=0.5e-6, seed=1 + s, **g)
            for k in agg:
                agg[k].append(r[k])
        z = np.mean(agg["z"]); gi = np.mean(agg["giant"])
        ll = np.mean(agg["L_over_lc"])
        nxlf = np.mean(agg["n_xl"]) / g["n_fil"]
        passed = (3.0 <= z <= 3.6 and gi >= 0.9 and ll >= 5.9)
        flag = "  <== PASS" if passed else ("  ~giant" if gi >= 0.9 else "")
        print(f"{g['n_fil']:>5} {g['L_long']*1e6:>6.1f} {g['z_target']:>2} "
              f"{g['reach']*1e9:>5.0f} | {z:>4.2f} {gi:>5.2f} {ll:>4.1f} "
              f"{np.mean(agg['homeless_frac'])*100:>5.1f} "
              f"{np.mean(agg['mesh_xi_nm']):>6.0f} {nxlf:>6.2f} "
              f"{np.mean(agg['mean_beads']):>6.1f}{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
