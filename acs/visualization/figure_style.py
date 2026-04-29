"""acs.visualization.figure_style — paper-grade figure helpers.

Per `CODEX_FIGURE_GUIDE.md` (project rule, 2026-04-30):
- Helvetica typography (TeX Gyre Heros / Arial / DejaVu Sans fallback)
- SVG draft + 600 dpi transparent PNG dual save
- `Figures_for_draft/` + `Figures_for_PPT/` parallel directories with same stem
- JSON manifest under `Logs/_captions/<step>/<stem>.json`
- No Markdown sidecar files
- Real Greek symbols in publication labels
"""

from __future__ import annotations

import json
import platform
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

import matplotlib as mpl
import matplotlib.font_manager as fm

# Helvetica setup — prefer Helvetica, then TeX Gyre Heros, then Arial,
# then DejaVu Sans. On Linux/WSL, register TeX Gyre Heros if available.
def install_helvetica_style() -> None:
    """Configure Matplotlib rcParams per CODEX_FIGURE_GUIDE.md."""
    # Optional: register TeX Gyre Heros otf files if present.
    texgyre_dir = Path.home() / ".local/share/fonts/texgyre"
    if texgyre_dir.exists():
        for otf in texgyre_dir.glob("*.otf"):
            try:
                fm.fontManager.addfont(str(otf))
            except Exception:
                continue
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = [
        "Helvetica",
        "TeX Gyre Heros",
        "Arial",
        "DejaVu Sans",
    ]
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["mathtext.fontset"] = "custom"
    mpl.rcParams["mathtext.it"] = "sans:italic"
    mpl.rcParams["mathtext.rm"] = "sans"
    # Recommended sizes from the figure guide.
    mpl.rcParams["axes.labelsize"] = 7.5
    mpl.rcParams["xtick.labelsize"] = 6.0
    mpl.rcParams["ytick.labelsize"] = 6.0
    mpl.rcParams["legend.fontsize"] = 5.5
    mpl.rcParams["axes.titlesize"] = 7.0


def save_figure_dual(
    fig,
    stem: str,
    *,
    run_dir: Path,
    caption: Optional[str] = None,
    description: Optional[str] = None,
    metadata: Optional[dict] = None,
    step: str = "stage1dc",
) -> dict:
    """Save figure to both SVG (draft, white bg) and PNG (PPT, 600 dpi
    transparent bg) under the run's `figures/Figures_for_draft/` and
    `figures/Figures_for_PPT/` directories, plus a JSON manifest under
    `figures/Logs/_captions/<step>/<stem>.json`.

    Returns dict with 'svg', 'png', 'manifest' Path entries.
    """
    fig_root = Path(run_dir) / "figures"
    draft_dir = fig_root / "Figures_for_draft"
    ppt_dir = fig_root / "Figures_for_PPT"
    log_dir = fig_root / "Logs" / "_captions" / step
    draft_dir.mkdir(parents=True, exist_ok=True)
    ppt_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    svg_path = draft_dir / f"{stem}.svg"
    png_path = ppt_dir / f"{stem}.png"
    manifest_path = log_dir / f"{stem}.json"

    # Draft SVG — white background.
    fig.savefig(svg_path, format="svg", bbox_inches="tight",
                facecolor="white", edgecolor="none")
    # PPT PNG — 600 dpi, transparent background.
    fig.savefig(png_path, format="png", dpi=600, bbox_inches="tight",
                facecolor="none", transparent=True)

    # Manifest. Provenance metadata per the figure guide.
    md = {
        "stem": stem,
        "step": step,
        "run_dir": str(Path(run_dir).resolve()),
        "svg_path": str(svg_path.relative_to(fig_root.parent)),
        "png_path": str(png_path.relative_to(fig_root.parent)),
        "caption": caption,
        "description": description,
        "metadata": metadata or {},
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "matplotlib_version": mpl.__version__,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    # Try to capture git commit (best-effort; figure guide allows).
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(Path(run_dir).parent.parent),
            stderr=subprocess.DEVNULL,
            timeout=5,
        ).decode().strip()
        md["git_commit"] = sha
    except Exception:
        md["git_commit"] = None

    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(md, f, indent=2, ensure_ascii=False)

    return {"svg": svg_path, "png": png_path, "manifest": manifest_path}
