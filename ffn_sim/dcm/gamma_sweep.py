"""Cortical-tension (γ) controlled-variable sweep for the Warp DCM aggregate.

PI 2026-06-29 (DCM autonomous loop): cortical surface tension γ is treated as a
**controlled independent variable** — we vary γ across a ladder and REPORT the
faceting / morphology that emerges, with the SimuCell3D faceting band overlaid.
We do NOT tune γ to hit any target morphology (that would be the magic-number sin
the PI forbade); γ is the knob, morphology is the read-out.

Physics basis (SimuCell3D, Runser, Vetter & Iber 2024): a deformable-cell
aggregate facets when the dimensionless surface tension

    γ̃ = γ / (K · ℓ),     ℓ = V₀^(1/3)  (cell length scale), K = bulk modulus

enters the band γ̃ ∈ [0.02, 0.10]. Below it cells stay rounded marbles (turgor
dominates, no flat junction faces); inside it cohesive cells flatten into
foam-like polygonal contact faces (tissue faceting); far above it the cortex
crushes the aggregate. This sweep walks γ through and across that band so the
faceting transition is *observed*, not assumed.

Mechanics: each γ point is one run of the parity-gated engine
(``dcm_warp_decohesion``) as an isolated **free aggregate** (substrate well +
wetting OFF, no division, no spreading drivers) so γ is the only thing that
changes. Runs launch as concurrent subprocesses (``--jobs``) — on the gbook
A5000 they share the GPU (slower but parallel, PI-sanctioned); on CPU they share
cores. Each worker saves its per-frame mesh (``--save-frames``); this harness
reads the final frame, computes morphology, and assembles a table + a
faceting=f(γ) figure.

Run (CPU smoke — FORM only):
    PYTHONPATH=. python -m ffn_sim.dcm.gamma_sweep --device cpu --n-cells 12 \
        --subdiv 1 --steps 4000 --gammas 0,3e-4,1e-3,2.5e-3 --jobs 2

Run (gbook production):
    PYTHONPATH=. python -m ffn_sim.dcm.gamma_sweep --device cuda:0 --n-cells 200 \
        --subdiv 2 --steps 40000 --gammas 0,3e-4,5e-4,1e-3,1.5e-3,2.5e-3,5e-3 --jobs 4
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

#: SimuCell3D faceting band in the dimensionless surface tension γ̃ = γ/(K·ℓ).
FACETING_BAND = (0.02, 0.10)
#: bulk modulus the Warp driver uses (``run_decohesion`` k_vol default) — γ̃ denominator.
K_VOL_DEFAULT = 7.73e5
UM2 = 1.0e12


# ---------------------------------------------------------------------------
# morphology read-out (computed from a saved final mesh — never tuned)
# ---------------------------------------------------------------------------
def _cell_volume(verts: np.ndarray, tris: np.ndarray) -> float:
    """Signed-tetrahedron enclosed volume of one closed triangle mesh [m³]."""
    t = verts[tris]
    return abs(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()) / 6.0


def _per_cell_asphericity(pos: np.ndarray, cof: np.ndarray) -> float:
    """Mean per-cell shape asphericity from the node-cloud gyration tensor.

    0 = perfect sphere (rounded marble); grows as a cell flattens into faceted /
    elongated polygonal contact. Averaged over live cells (cof ≥ 0).
    """
    vals = []
    for c in np.unique(cof[cof >= 0]):
        p = pos[cof == c]
        if p.shape[0] < 4:
            continue
        d = p - p.mean(0)
        ev = np.sort(np.linalg.eigvalsh(d.T @ d))[::-1]
        s = ev.sum()
        if s <= 0:
            continue
        vals.append((ev[0] - 0.5 * (ev[1] + ev[2])) / s)
    return float(np.mean(vals)) if vals else 0.0


def _contact_area_fraction(pos: np.ndarray, faces: np.ndarray, cof: np.ndarray,
                           d_contact: float) -> float:
    """Fraction of total surface area on **cell–cell contact faces** — the direct faceting signal.

    A free cell under turgor+cortical tension is a sphere (no contact faces). As γ̃ enters the
    SimuCell3D faceting band, cohesive cells flatten into foam-like polygonal junctions, so the
    contact-apposed area fraction rises. A face is "in contact" if its centroid lies within
    ``d_contact`` of a node belonging to a *different* cell. Complements per-cell asphericity
    (which measures whole-cell roundness, not junction flattening). All variables but γ are held
    fixed in the sweep, so the *change* with γ isolates the faceting response.
    """
    live = cof >= 0
    fcell = cof[faces[:, 0]]                       # owner cell of each face (node-0)
    surf = faces[fcell >= 0]
    if surf.shape[0] == 0:
        return 0.0
    cent = pos[surf].mean(axis=1)                  # face centroids
    # face areas
    e1 = pos[surf[:, 1]] - pos[surf[:, 0]]
    e2 = pos[surf[:, 2]] - pos[surf[:, 0]]
    area = 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1)
    own = cof[surf[:, 0]]
    contact = np.zeros(surf.shape[0], dtype=bool)
    pts = pos[live]
    ptc = cof[live]
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(pts)
        for fi in range(surf.shape[0]):
            idx = tree.query_ball_point(cent[fi], d_contact)
            if any(ptc[j] != own[fi] for j in idx):
                contact[fi] = True
    except Exception:  # noqa: BLE001 — brute fallback
        for fi in range(surf.shape[0]):
            dd = np.linalg.norm(pts - cent[fi], axis=1)
            if np.any((dd < d_contact) & (ptc != own[fi])):
                contact[fi] = True
    tot = float(area.sum())
    return float(area[contact].sum() / tot) if tot > 0 else 0.0


def _mean_edge(pos: np.ndarray, faces: np.ndarray) -> float:
    """Mean triangle-edge length [m] (contact-band scale)."""
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
    return float(np.linalg.norm(pos[e[:, 0]] - pos[e[:, 1]], axis=1).mean())


def _cluster_asphericity(centroids: np.ndarray) -> float:
    """Asphericity of the cell-centroid cloud (whole-aggregate roundness)."""
    if centroids.shape[0] < 4:
        return 0.0
    d = centroids - centroids.mean(0)
    ev = np.sort(np.linalg.eigvalsh(d.T @ d))[::-1]
    s = ev.sum()
    return float((ev[0] - 0.5 * (ev[1] + ev[2])) / s) if s > 0 else 0.0


@dataclass
class GammaPoint:
    gamma: float           # cortical surface tension [N/m]
    gamma_tilde: float     # γ̃ = γ/(K·ℓ) dimensionless
    in_band: bool          # γ̃ inside the SimuCell3D faceting band?
    n_cells: int
    vv0_mean: float        # turgor: mean V/V0 (1 = held; <1 = compressed)
    cell_asphericity: float  # mean per-cell shape (0 sphere → faceting)
    cluster_asphericity: float
    contact_area_frac: float  # fraction of surface on cell-cell contact faces (direct faceting)
    pen_frac: float        # interpenetration depth / mean_edge (contact fidelity)
    topdown_um2: float
    wall_s: float
    ok: bool


def _measure(npz_path: Path, out_json: Path, k_vol: float) -> GammaPoint | None:
    """Load a worker's final mesh + out-json; compute the γ-point morphology."""
    if not npz_path.exists() or not out_json.exists():
        return None
    out = json.loads(out_json.read_text())
    d = np.load(npz_path, allow_pickle=True)
    pos = d["frames"][-1].astype(np.float64)            # final-frame node positions [m]
    # topology: under no-remesh the static faces/cof apply to every frame.
    faces = d["faces"].astype(np.int64)
    cof = d["cof"].astype(np.int64)
    live = cof >= 0
    cells = np.unique(cof[live])
    # per-cell volume vs rest volume V0 (rest = first-frame volume of the same cell)
    pos0 = d["frames"][0].astype(np.float64)
    vrs = []
    cents = []
    for c in cells:
        fmask = (cof[faces[:, 0]] == c)
        ftri = faces[fmask]
        if ftri.shape[0] < 4:
            continue
        v_now = _cell_volume(pos, ftri)
        v_rest = _cell_volume(pos0, ftri)
        if v_rest > 0:
            vrs.append(v_now / v_rest)
        cents.append(pos[cof == c].mean(0))
    gamma = float(out.get("gamma_surf") or 0.0)
    k_run = float(out.get("k_vol") or k_vol)            # K the run actually used (faceting denom)
    ell = (float(np.mean([_cell_volume(pos0, faces[(cof[faces[:, 0]] == c)])
                          for c in cells]))) ** (1.0 / 3.0)
    gt = gamma / (k_run * ell) if (k_run * ell) > 0 else 0.0
    xy = pos[live][:, :2]
    try:
        from scipy.spatial import ConvexHull
        topd = float(ConvexHull(xy).volume) * UM2 if xy.shape[0] >= 3 else 0.0
    except Exception:  # noqa: BLE001
        topd = 0.0
    return GammaPoint(
        gamma=gamma, gamma_tilde=gt,
        in_band=(FACETING_BAND[0] <= gt <= FACETING_BAND[1]),
        n_cells=int(cells.size),
        vv0_mean=float(np.mean(vrs)) if vrs else 0.0,
        cell_asphericity=_per_cell_asphericity(pos, cof),
        cluster_asphericity=_cluster_asphericity(np.array(cents)) if cents else 0.0,
        contact_area_frac=_contact_area_fraction(
            pos, faces, cof, d_contact=0.6 * _mean_edge(pos0, faces)),
        pen_frac=float(out.get("pen_frac_final", 0.0)),
        topdown_um2=topd,
        wall_s=float(out.get("steps", 0)) / max(float(out.get("steps_per_s", 1.0)), 1e-9),
        ok=bool(out.get("truncated_at") in (None, 0)),
    )


# ---------------------------------------------------------------------------
# worker launch
# ---------------------------------------------------------------------------
def _worker_cmd(gamma: float, *, device: str, n_cells: int, subdiv: int, steps: int,
                gap: float, frames: int, k_vol: float, npz_path: Path) -> list[str]:
    """One free-aggregate morphology run at a fixed γ (substrate + division OFF)."""
    cmd = [sys.executable, "-m", "ffn_sim.dcm.dcm_warp_decohesion",
           "--device", device, "--n-cells", str(n_cells), "--subdiv", str(subdiv),
           "--steps", str(steps), "--frames", str(frames), "--gap", str(gap),
           "--k-vol", repr(k_vol), "--builder", "sphere", "--no-well", "--no-wetting",
           "--save-frames", str(npz_path)]
    if gamma > 0.0:
        cmd += ["--surface-tension", "--gamma-surf", repr(gamma)]
    return cmd


def _run_one(gamma: float, outdir: Path, args) -> tuple[float, int, Path, Path, Path]:
    tag = f"g{gamma:.2e}".replace("-", "m").replace("+", "").replace(".", "p")
    npz = outdir / f"agg_{tag}.npz"
    log = outdir / f"agg_{tag}.log"
    jout = outdir / f"agg_{tag}.out.json"
    cmd = _worker_cmd(gamma, device=args.device, n_cells=args.n_cells,
                      subdiv=args.subdiv, steps=args.steps, gap=args.gap,
                      frames=args.frames, k_vol=args.k_vol, npz_path=npz)
    t0 = time.perf_counter()
    with log.open("w") as f:
        # the driver prints its out-dict json to stdout tail; capture it to jout via a wrapper
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    # the driver's stdout has the json out-dict at the end — extract the last {...} block
    _extract_out_json(log, jout)
    dt = time.perf_counter() - t0
    print(f"  [γ={gamma:.2e}] done rc={proc.returncode} in {dt:.0f}s -> {npz.name}", flush=True)
    return gamma, proc.returncode, npz, log, jout


def _extract_out_json(log: Path, jout: Path) -> None:
    """Pull the trailing JSON out-dict the driver prints to stdout into ``jout``."""
    txt = log.read_text()
    start = txt.rfind("\n{")
    if start < 0:
        return
    blob = txt[start:].strip()
    # find the matching closing brace by bracket counting
    depth = 0
    end = -1
    for i, ch in enumerate(blob):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end > 0:
        try:
            json.loads(blob[:end])
            jout.write_text(blob[:end])
        except json.JSONDecodeError:
            pass


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n-cells", type=int, default=12)
    ap.add_argument("--subdiv", type=int, default=1)
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--frames", type=int, default=8)
    ap.add_argument("--gap", type=float, default=2.05)
    ap.add_argument("--gammas", default="0,3e-4,5e-4,1e-3,1.5e-3,2.5e-3",
                    help="comma-separated γ ladder [N/m]; 0 = cortical tension OFF (control)")
    ap.add_argument("--k-vol", type=float, default=K_VOL_DEFAULT,
                    help="bulk modulus for the γ̃=γ/(K·ℓ) denominator (driver k_vol)")
    ap.add_argument("--jobs", type=int, default=2, help="concurrent worker subprocesses")
    ap.add_argument("--out", default="ffn_sim/outputs/h_dcm_two_stage/gamma_sweep")
    args = ap.parse_args()

    gammas = [float(g) for g in args.gammas.split(",")]
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"[γ-sweep] {len(gammas)} points {gammas} · device={args.device} · "
          f"N={args.n_cells} subdiv={args.subdiv} steps={args.steps} jobs={args.jobs}", flush=True)

    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        results = list(ex.map(lambda g: _run_one(g, outdir, args), gammas))

    points: list[GammaPoint] = []
    for gamma, rc, npz, log, jout in results:
        gp = _measure(npz, jout, args.k_vol)
        if gp is not None:
            points.append(gp)
    points.sort(key=lambda p: p.gamma)

    summary = {
        "device": args.device, "n_cells": args.n_cells, "subdiv": args.subdiv,
        "steps": args.steps, "gap": args.gap, "k_vol": args.k_vol,
        "faceting_band": FACETING_BAND, "wall_s_total": time.perf_counter() - t0,
        "controlled_variable": "gamma_surf",
        "note": "γ is the controlled independent variable; morphology is the read-out. "
                "NOT tuned to any target (no magic-number).",
        "points": [asdict(p) for p in points],
    }
    sj = outdir / "gamma_sweep_summary.json"
    sj.write_text(json.dumps(summary, indent=2))
    print(f"[γ-sweep] {len(points)}/{len(gammas)} ok -> {sj}", flush=True)
    for p in points:
        print(f"  γ={p.gamma:.2e}  γ̃={p.gamma_tilde:.4f} "
              f"{'IN-BAND' if p.in_band else '       '}  "
              f"V/V0={p.vv0_mean:.3f}  cell_asph={p.cell_asphericity:.3f}  "
              f"contact_area={p.contact_area_frac:.3f}  pen={p.pen_frac:.2f}", flush=True)

    try:
        _plot(points, outdir / "gamma_sweep_faceting.png", args)
    except Exception as e:  # noqa: BLE001
        print(f"  [plot skipped] {e}", flush=True)


def _plot(points: list[GammaPoint], path: Path, args) -> None:
    """Faceting/morphology = f(γ̃), SimuCell3D band overlaid. No axis truncation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not points:
        return
    gt = np.array([p.gamma_tilde for p in points])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    series = [("contact_area_frac", "cell–cell contact-area fraction\n(↑ = foam-like faceting)"),
              ("cell_asphericity", "per-cell asphericity\n(0 = round sphere)"),
              ("vv0_mean", "mean V/V₀\n(turgor held=1)")]
    for ax, (attr, label) in zip(axes, series):
        y = np.array([getattr(p, attr) for p in points])
        ax.axvspan(*FACETING_BAND, color="gold", alpha=0.25,
                   label=f"SimuCell3D faceting band {FACETING_BAND}")
        ax.plot(gt, y, "o-", color="tab:blue")
        for p in points:
            ax.annotate(f"{p.gamma:.0e}", (p.gamma_tilde, getattr(p, attr)),
                        fontsize=7, alpha=0.7, xytext=(0, 5), textcoords="offset points")
        ax.set_xlabel("γ̃ = γ / (K·ℓ)")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="best")
    fig.suptitle(f"DCM cortical-tension sweep · N={args.n_cells} subdiv={args.subdiv} "
                 f"· γ = controlled variable (morphology read-out, NOT tuned)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"  wrote {path}", flush=True)


if __name__ == "__main__":
    main()
