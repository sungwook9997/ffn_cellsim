"""One command from the CUDA run artefacts to a figure of the evidence.

`ALEPH-PORT-3612`.

    python -m aleph.viz.evidence_figure

No arguments a person has to invent. It looks for the runner's artefacts locally, falls back to
fetching them off the GPU host over ssh, and writes an SVG. A refusal exits **2** and says what is missing;
exit **1** is an unexpected fault, so a caller can tell a refusal from a crash. That convention comes
from `scripts/gpu_preflight.py` by way of `aleph/viz/figure.py`, and the reason carries over exactly:
a tool that refuses and exits 0 lets ``make-figure && publish`` treat "nothing was drawn" as success.

What the figure shows that reading the JSON does not
----------------------------------------------------
**One row per channel, not per case.** A wrong kernel that is loud in the energy and *bit-identical to
the true kernel* in both force arrays occupies three rows in three different colours, and the pattern
is visible without reading a number. A per-case verdict prints ``False`` once and the finding is gone.

**A log axis, with exact agreement off it.** The measured ULP spans nine decades. On a linear axis
every real result lands on the first pixel. Exact zero is not on the axis at all -- it is a different
kind of statement and it gets its own lane.

**Both devices on the same row.** The CUDA mark is filled, the warp CPU mark is a ring, and a line
joins them when they differ. 95 of 130 channels differ between devices; the blind set does not.

**The gate and the round-off bound, drawn.** A budget sitting less than a binary order under the
round-off bound its own ledger derived is drawn as thin, because it is.

Units: float32 ULP of each field's own reference scale, throughout.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from aleph.viz.evidence import (
    CHANNEL_STYLE,
    DECADE_HIGH,
    DECADE_LOW,
    DEFAULT_HOST,
    EXACT_LANE,
    THIN_HEADROOM,
    ChannelState,
    EvidenceBundle,
    EvidenceError,
    ParityEvidence,
    load_bundle,
    log_position,
)

__all__ = [
    "EXIT_OK",
    "EXIT_FAULT",
    "EXIT_REFUSED",
    "evidence_figure",
    "main",
    "render_evidence_svg",
]

EXIT_OK = 0
EXIT_FAULT = 1
EXIT_REFUSED = 2

_BACKGROUND = "#fbfbfd"
_INK = "#111"
_MUTED = "#555"
_FAINT = "#d8d8e0"
_GATE = "#37474f"
_GATE_THIN = "#ef6c00"
_ROUNDOFF = "#b0bec5"

_ROW_HEIGHT = 14.0
_FAMILY_GAP = 20.0
_LABEL_MAX_CHARS = 38
_GLYPH_X = 326.0
#: The exact lane sits left of the axis with room for its own label, the axis floor label, and a
#: visible break between them. The three collided when the gap was 30px.
_EXACT_X = 348.0
_AXIS_LEFT = 412.0
_AXIS_WIDTH = 690.0
_PANEL_X = 1136.0
_WIDTH = 1680.0

#: Header layout, written as absolute baselines rather than offsets from the total: the summary band,
#: the two prose lines and the axis labels are four things that must not overlap, and expressing them
#: as deltas from one constant is how they came to.
_SUMMARY_Y = (100.0, 113.0, 134.0, 147.0)
_DECLARED_Y = 168.0
_FAULT_Y = 183.0
_AXIS_TITLE_Y = 204.0
_DECADE_Y = 222.0
_HEADER = 236.0


def _esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def _x_for(ulp: float | None) -> tuple[float | None, bool]:
    """Axis x for a ULP figure. ``None`` means it does not belong on the continuum."""
    if ulp is None:
        return None, False
    if ulp == 0.0:
        return _EXACT_X, False
    try:
        fraction, clamped = log_position(ulp)
    except EvidenceError:
        return None, False
    return _AXIS_LEFT + fraction * _AXIS_WIDTH, clamped


def _fmt(value: float | None, *, digits: int = 4) -> str:
    if value is None:
        return "—"
    if value == 0.0:
        return "0 (exact)"
    if abs(value) >= 1.0e5 or abs(value) < 1.0e-2:
        return f"{value:.{digits - 1}e}"
    return f"{value:.{digits}g}"


def _measure(parity: ParityEvidence) -> tuple[float, list[tuple[str, list]]]:
    """Group rows by family in artefact order, and return the canvas height they need."""
    groups: list[tuple[str, list]] = []
    seen: dict[str, list] = {}
    for row in parity.rows:
        if row.family not in seen:
            seen[row.family] = []
            groups.append((row.family, seen[row.family]))
        seen[row.family].append(row)
    height = _HEADER + sum(_FAMILY_GAP + len(rows) * _ROW_HEIGHT for _, rows in groups) + 56.0
    return max(height, 900.0), groups


def render_evidence_svg(bundle: EvidenceBundle, *, title: str | None = None) -> str:
    """Render the whole bundle as one SVG.

    Args:
        bundle: What was read, including what was not there.
        title: Figure title; a default naming the device and run time is used when omitted.

    Returns:
        SVG text, newline-terminated.

    Raises:
        EvidenceError: If the parity evidence has no rows. An empty matrix is a refusal, never an
            empty figure -- a renderer that draws nothing and reports success is this module's most
            likely silent failure.
    """
    parity = bundle.parity
    if not parity.rows:
        raise EvidenceError(
            "the parity evidence carries no field rows, so there is nothing to draw. An empty "
            "figure that exits 0 is the failure this tool exists to prevent."
        )

    height, groups = _measure(parity)
    device_name = parity.device.get("name", "unknown device")
    out: list[str] = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{_WIDTH:.0f}' height='{height:.0f}' "
        f"viewBox='0 0 {_WIDTH:.0f} {height:.0f}' font-family='ui-monospace,monospace' "
        "font-size='11'>",
        f"<rect width='{_WIDTH:.0f}' height='{height:.0f}' fill='{_BACKGROUND}'/>",
    ]

    heading = title or "Aleph — law-parity evidence, per channel"
    out.append(f"<text x='40' y='30' font-size='17' fill='{_INK}'>{_esc(heading)}</text>")
    source = bundle.source if len(bundle.source) <= 62 else "..." + bundle.source[-59:]
    out.append(
        f"<text x='40' y='48' fill='{_MUTED}'>{_esc(device_name)} · "
        f"{_esc(parity.cuda_key)} vs {_esc(parity.cpu_key)} · positions="
        f"{_esc(parity.position_precision or 'unrecorded')} · run finished "
        f"{_esc(parity.finished_utc or 'unrecorded')} on host {_esc(parity.host or '?')}"
        f"<title>read from {_esc(bundle.source)}</title></text>"
    )
    out.append(
        f"<text x='40' y='64' fill='{_MUTED}'>read from {_esc(source)}</text>"
    )
    out.append(
        f"<text x='40' y='80' fill='{_MUTED}'>evidence rung "
        f"<tspan font-weight='bold'>{_esc(parity.evidence_rung or 'UNDECLARED')}</tspan> · the "
        f"artefact's own verdict is <tspan font-weight='bold'>{_esc(parity.verdict or '—')}</tspan>"
        " — that is the artefact's claim, not this figure's</text>"
    )

    out.extend(_summary_band(parity, bundle))
    out.extend(_axis(height))

    y = _HEADER
    for family, rows in groups:
        y += _FAMILY_GAP
        out.append(
            f"<text x='40' y='{y - 6:.1f}' font-size='12' fill='{_INK}' font-weight='bold'>"
            f"{_esc(family)}</text>"
        )
        out.append(
            f"<line x1='40' y1='{y - 2:.1f}' x2='{_AXIS_LEFT + _AXIS_WIDTH:.1f}' "
            f"y2='{y - 2:.1f}' stroke='{_FAINT}' stroke-width='1'/>"
        )
        for row in rows:
            out.extend(_row(row, y))
            y += _ROW_HEIGHT

    out.extend(_side_panel(bundle, height))
    out.append("</svg>")
    return "\n".join(out) + "\n"


def _summary_band(parity: ParityEvidence, bundle: EvidenceBundle) -> list[str]:
    """The headline counts, so the ten-second read does not depend on scanning 130 rows."""
    gate = parity.gate
    faults = parity.faults()
    cells: list[tuple[str, str, str]] = [
        (str(len(parity.true_kernels)), "true kernels", "#2e7d32"),
        (str(len(parity.wrong_kernels)), "wrong on purpose", _INK),
        (
            f"{gate.get('caught_on_cuda', '?')}/{len(parity.wrong_kernels)}",
            f"caught on {parity.cuda_key}",
            "#c62828",
        ),
        (
            f"{gate.get('caught_on_cpu', '?')}/{len(parity.wrong_kernels)}",
            f"caught on {parity.cpu_key}",
            "#c62828",
        ),
        (
            str(len(gate.get("caught_on_cpu_but_not_cuda") or [])
                + len(gate.get("caught_on_cuda_but_not_cpu") or [])),
            "caught on one device only",
            "#2e7d32",
        ),
        (str(len(parity.blind_rows)), "channels blind to the defect", "#6a1b9a"),
        (str(len(parity.thin_budgets)), f"gates with headroom < {THIN_HEADROOM:g}", _GATE_THIN),
        (str(len(faults)), "faults", "#c62828" if faults else "#2e7d32"),
    ]
    out: list[str] = []
    value_row_one, label_row_one, value_row_two, label_row_two = _SUMMARY_Y
    for index, (value, label, colour) in enumerate(cells):
        column, row = index % 4, index // 4
        x = 40.0 + column * 250.0
        value_y = value_row_one if row == 0 else value_row_two
        label_y = label_row_one if row == 0 else label_row_two
        out.append(
            f"<text x='{x:.0f}' y='{value_y:.0f}' font-size='21' fill='{colour}'>"
            f"{_esc(value)}</text>"
        )
        out.append(
            f"<text x='{x:.0f}' y='{label_y:.0f}' font-size='10' fill='{_MUTED}'>"
            f"{_esc(label)}</text>"
        )
    declared = gate.get("blind_on_both_devices") or []
    if declared:
        out.append(
            f"<text x='40' y='{_DECLARED_Y:.0f}' fill='#6a1b9a'>declared blind on both devices "
            f"({len(declared)}): {_esc(', '.join(str(n) for n in declared))} — "
            "ALEPH-PORT-3604 §9a, an IEEE-754 identity, so no gate can see it and none should"
            "</text>"
        )
    else:
        out.append(
            f"<text x='40' y='{_DECLARED_Y:.0f}' fill='{_MUTED}'>no variant is declared blind on "
            "both devices</text>"
        )
    if faults:
        out.append(
            f"<text x='40' y='{_FAULT_Y:.0f}' fill='#c62828' font-weight='bold'>"
            f"{_esc('; '.join(faults)[:150])}</text>"
        )
    else:
        out.append(
            f"<text x='40' y='{_FAULT_Y:.0f}' fill='#2e7d32'>no undeclared blind variant, no "
            "true-kernel regression, no case dropped — the recount reproduces the artefact's own "
            "gate comparison on every key</text>"
        )
    return out


def _axis(height: float) -> list[str]:
    """The shared log axis, its decade gridlines, and the exact lane that is not on it."""
    out: list[str] = []
    top = _HEADER - 8.0
    bottom = height - 44.0
    for decade in range(DECADE_LOW, DECADE_HIGH + 1):
        fraction = (decade - DECADE_LOW) / float(DECADE_HIGH - DECADE_LOW)
        x = _AXIS_LEFT + fraction * _AXIS_WIDTH
        out.append(
            f"<line x1='{x:.1f}' y1='{top:.1f}' x2='{x:.1f}' y2='{bottom:.1f}' "
            f"stroke='{_FAINT}' stroke-width='1'/>"
        )
        out.append(
            f"<text x='{x:.1f}' y='{_DECADE_Y:.1f}' font-size='10' fill='{_MUTED}' "
            f"text-anchor='middle'>1e{decade}</text>"
        )
    out.append(
        f"<text x='{(_AXIS_LEFT + _AXIS_WIDTH / 2):.1f}' y='{_AXIS_TITLE_Y:.1f}' font-size='11' "
        f"fill='{_INK}' text-anchor='middle'>disagreement with the frozen NumPy law — float32 ULP "
        "(log scale)</text>"
    )
    # The exact lane sits left of the axis with a visible break, because 0.0 ULP is exact agreement
    # and not a small disagreement. See `ALEPH-PORT-3612` §4b.
    out.append(
        f"<line x1='{_EXACT_X:.1f}' y1='{top:.1f}' x2='{_EXACT_X:.1f}' y2='{bottom:.1f}' "
        f"stroke='#0277bd' stroke-width='1' stroke-dasharray='2 3'/>"
    )
    out.append(
        f"<text x='{_EXACT_X:.1f}' y='{_DECADE_Y:.1f}' font-size='10' fill='#0277bd' "
        f"text-anchor='middle'>{_esc(EXACT_LANE)}</text>"
    )
    mid = (_EXACT_X + _AXIS_LEFT) / 2.0
    out.append(
        f"<text x='{mid:.1f}' y='{_DECADE_Y:.1f}' font-size='11' fill='{_MUTED}' "
        "text-anchor='middle'>//</text>"
    )
    return out


def _row(row, y: float) -> list[str]:
    """One channel: its label, its gate, its round-off bound, and both devices' marks."""
    out: list[str] = []
    style = CHANNEL_STYLE[row.state]
    baseline = y + _ROW_HEIGHT * 0.72
    label = f"{row.variant} · {row.field_name}"
    if len(label) > _LABEL_MAX_CHARS:
        label = label[: _LABEL_MAX_CHARS - 3] + "..."
    colour = "#c62828" if row.state is ChannelState.OVER_BUDGET else _INK
    if row.state is ChannelState.BLIND_TO_TRUE:
        colour = "#6a1b9a"
    out.append(
        f"<text x='56' y='{baseline:.1f}' font-size='10' fill='{colour}'>{_esc(label)}"
        f"<title>{_esc(_row_title(row))}</title></text>"
    )
    out.append(
        f"<text x='{_GLYPH_X:.0f}' y='{baseline:.1f}' font-size='10' fill='{_MUTED}' "
        f"text-anchor='end'>{_esc(style.glyph)}</text>"
    )

    if row.state is ChannelState.ABSENT:
        out.append(
            f"<text x='{_AXIS_LEFT + 6:.1f}' y='{baseline:.1f}' font-size='10' fill='#8d6e63'>"
            f"absent — {_esc(row.absent_reason[:96])}</text>"
        )
        return out

    top = y + 1.0
    bottom = y + _ROW_HEIGHT - 2.0

    # The round-off bound the owning ledger derived, then the gate that was actually applied.
    bound_x, _ = _x_for(row.roundoff_bound)
    if bound_x is not None:
        out.append(
            f"<line x1='{bound_x:.1f}' y1='{top:.1f}' x2='{bound_x:.1f}' y2='{bottom:.1f}' "
            f"stroke='{_ROUNDOFF}' stroke-width='1' stroke-dasharray='1 2'/>"
        )
    gate_x, _ = _x_for(row.budget)
    if gate_x is not None:
        gate_colour = _GATE_THIN if row.thin_headroom else _GATE
        out.append(
            f"<line x1='{gate_x:.1f}' y1='{top:.1f}' x2='{gate_x:.1f}' y2='{bottom:.1f}' "
            f"stroke='{gate_colour}' stroke-width='{1.8 if row.thin_headroom else 1.1}'/>"
        )
        if row.thin_headroom and bound_x is not None:
            out.append(
                f"<line x1='{gate_x:.1f}' y1='{(top + bottom) / 2:.1f}' x2='{bound_x:.1f}' "
                f"y2='{(top + bottom) / 2:.1f}' stroke='{_GATE_THIN}' stroke-width='0.8' "
                "stroke-dasharray='2 1'/>"
            )

    cuda_x, cuda_clamped = _x_for(row.ulp_cuda)
    cpu_x, cpu_clamped = _x_for(row.ulp_cpu)
    mid_y = (top + bottom) / 2.0

    # Overshoot: from the gate to the mark, so "how far over" is a length and not a ratio to work out.
    if gate_x is not None and cuda_x is not None and row.state is ChannelState.OVER_BUDGET:
        out.append(
            f"<line x1='{gate_x:.1f}' y1='{mid_y:.1f}' x2='{cuda_x:.1f}' y2='{mid_y:.1f}' "
            "stroke='#c62828' stroke-width='1.6' stroke-opacity='0.55'/>"
        )
    if cpu_x is not None and cuda_x is not None and cpu_x != cuda_x:
        out.append(
            f"<line x1='{cpu_x:.1f}' y1='{mid_y:.1f}' x2='{cuda_x:.1f}' y2='{mid_y:.1f}' "
            f"stroke='{_MUTED}' stroke-width='0.8' stroke-opacity='0.6'/>"
        )
    if cpu_x is not None:
        out.append(
            f"<circle cx='{cpu_x:.1f}' cy='{mid_y:.1f}' r='4.2' fill='none' "
            f"stroke='{style.stroke}' stroke-width='1.2' stroke-opacity='0.85'/>"
        )
    if cuda_x is not None:
        edge = "#000" if cuda_clamped else style.stroke
        out.append(
            f"<circle cx='{cuda_x:.1f}' cy='{mid_y:.1f}' r='3.1' fill='{style.fill}' "
            f"stroke='{edge}' stroke-width='{1.4 if cuda_clamped else 0.7}'>"
            f"<title>{_esc(_row_title(row))}</title></circle>"
        )
        if cuda_clamped:
            out.append(
                f"<text x='{cuda_x - 8:.1f}' y='{baseline:.1f}' font-size='9' fill='#000'>&lt;"
                "</text>"
            )
    if cpu_clamped and cpu_x is not None and not cuda_clamped:
        out.append(
            f"<text x='{cpu_x - 10:.1f}' y='{baseline:.1f}' font-size='9' fill='#000'>&lt;</text>"
        )
    return out


def _row_title(row) -> str:
    """The full record for one channel, recoverable on hover since the label is elided."""
    parts = [
        f"{row.case} · {row.field_name}",
        f"state={row.state.value}",
        f"{row.family} shape={list(row.shape) or 'scalar'}",
        f"ULP cuda={_fmt(row.ulp_cuda)}  cpu={_fmt(row.ulp_cpu)}",
        f"budget={_fmt(row.budget)}  round-off bound={_fmt(row.roundoff_bound)}",
        f"pressure={_fmt(row.pressure)}  headroom={_fmt(row.headroom)}",
    ]
    if row.bit_identical_to_true is None:
        parts.append("bit-identical to true kernel: unknown (no true variant in this family)")
    else:
        parts.append(f"bit-identical to true kernel: {row.bit_identical_to_true}")
    if row.scatter_assembled:
        parts.append("assembled by a scatter — the ORDERED/ATOMIC question applies here")
    if row.thin_headroom:
        parts.append(f"THIN: the gate is under {THIN_HEADROOM:g}x from its own round-off bound")
    return "\n".join(parts)


def _side_panel(bundle: EvidenceBundle, height: float) -> list[str]:
    """Legend, the two other jobs, and what the artefact says it does not establish."""
    parity = bundle.parity
    out: list[str] = []
    x = _PANEL_X
    y = _HEADER + 6.0

    def heading(text: str, colour: str = _INK) -> None:
        nonlocal y
        y += 8.0
        out.append(
            f"<text x='{x:.0f}' y='{y:.0f}' font-size='12' fill='{colour}' font-weight='bold'>"
            f"{_esc(text)}</text>"
        )
        y += 15.0

    def line(text: str, colour: str = "#333", indent: float = 0.0) -> None:
        nonlocal y
        out.append(
            f"<text x='{x + indent:.0f}' y='{y:.0f}' font-size='10' fill='{colour}'>"
            f"{_esc(text[:74])}</text>"
        )
        y += 12.5

    heading("channel states")
    for state in ChannelState:
        style = CHANNEL_STYLE[state]
        out.append(
            f"<circle cx='{x + 6:.0f}' cy='{y - 3.5:.0f}' r='3.6' fill='{style.fill}' "
            f"stroke='{style.stroke}' stroke-width='0.9'/>"
        )
        out.append(
            f"<text x='{x + 16:.0f}' y='{y:.0f}' font-size='10' fill='{style.stroke}'>"
            f"{_esc(style.glyph)}</text>"
        )
        line(style.label, colour="#333", indent=26.0)
    line("filled = CUDA · ring = warp CPU · a line joins them when they differ", _MUTED)
    line("solid vertical = the gate · dotted = its round-off bound", _MUTED)
    line(f"orange gate = headroom < {THIN_HEADROOM:g}x, no room to move", _GATE_THIN)
    line("'<' = clamped to the axis floor, not a measured position", _MUTED)

    heading("cross-device gate comparison")
    gate = parity.gate
    line(f"wrong kernels: {gate.get('wrong_kernels_total', '?')}")
    line(f"caught on {parity.cuda_key}: {gate.get('caught_on_cuda', '?')}")
    line(f"caught on {parity.cpu_key}: {gate.get('caught_on_cpu', '?')}")
    only_cuda = gate.get("caught_on_cuda_but_not_cpu") or []
    only_cpu = gate.get("caught_on_cpu_but_not_cuda") or []
    line(
        f"caught on CUDA and missed on CPU: {len(only_cuda)}",
        "#c62828" if only_cuda else "#2e7d32",
    )
    line(
        f"caught on CPU and missed on CUDA: {len(only_cpu)}",
        "#c62828" if only_cpu else "#2e7d32",
    )
    line("(the second is the finding this job existed to look for)", _MUTED)
    differing = sum(1 for row in parity.rows if row.moved_between_devices)
    line(f"channels whose ULP differs between devices: {differing}/{len(parity.rows)}", _MUTED)
    line(f"channels blind to their defect: {len(parity.blind_rows)}", "#6a1b9a")

    heading(f"thinnest gates (headroom < {THIN_HEADROOM:g})", _GATE_THIN)
    thin = sorted(parity.thin_budgets, key=lambda r: r.headroom or 0.0)
    if not thin:
        line("none — every gate sits at least a binary order under its round-off bound", "#2e7d32")
    for row in thin[:6]:
        line(f"{row.headroom:.2f}x  {row.case}", _GATE_THIN)
        line(f"      {row.field_name}: gate {_fmt(row.budget)} vs bound "
             f"{_fmt(row.roundoff_bound)}", _MUTED)
    if len(thin) > 6:
        line(f"... and {len(thin) - 6} more", _MUTED)

    scatter_lines, y = _scatter_panel(bundle, x, y)
    out.extend(scatter_lines)
    descent_lines, y = _descent_panel(bundle, x, y)
    out.extend(descent_lines)

    heading("what this figure does not establish", "#6d4c41")
    for sentence in parity.not_established:
        words = sentence.split()
        buffer = ""
        for word in words:
            if len(buffer) + len(word) + 1 > 72:
                line(buffer, "#6d4c41")
                buffer = word
            else:
                buffer = f"{buffer} {word}".strip()
        if buffer:
            line(buffer, "#6d4c41")
        y += 3.0
    if bundle.missing:
        heading("artefacts looked for and not read", "#8d6e63")
        for name, reason in bundle.missing:
            line(f"{name}: {reason}", "#8d6e63")

    out.append(
        f"<text x='40' y='{height - 18:.0f}' font-size='10' fill='{_MUTED}'>"
        f"ALEPH-PORT-3612 · {len(parity.rows)} channels over {len(parity.cases)} cases · "
        "a ULP figure is STRUCTURAL evidence about an implementation, not evidence about a cell"
        "</text>"
    )
    return out


def _scatter_panel(bundle: EvidenceBundle, x: float, y: float) -> tuple[list[str], float]:
    """ORDERED versus ATOMIC, or a statement that the artefact was not there."""
    out: list[str] = []
    y += 8.0
    out.append(
        f"<text x='{x:.0f}' y='{y:.0f}' font-size='12' fill='{_INK}' font-weight='bold'>"
        "scatter determinism (j2)</text>"
    )
    y += 15.0

    def line(text: str, colour: str = "#333") -> None:
        nonlocal y
        out.append(
            f"<text x='{x:.0f}' y='{y:.0f}' font-size='10' fill='{colour}'>{_esc(text[:74])}</text>"
        )
        y += 12.5

    scatter = bundle.scatter
    if scatter is None:
        line("absent — j2_scatter_determinism.json was not read", "#8d6e63")
        line("this panel is drawn as absent rather than omitted", "#8d6e63")
        return out, y
    for row in scatter.modes:
        identical = bool(row.get("bit_identical"))
        line(
            f"{str(row.get('mode', '?')).upper()}  x{row.get('repeats', '?')} repeats: "
            f"{'bit-identical' if identical else 'NOT bit-identical'}",
            "#2e7d32" if identical else "#c62828",
        )
        line(
            f"      spread {_fmt(row.get('ulp_spread'))} ULP over "
            f"{row.get('contributions', '?')} contributions",
            _MUTED,
        )
    line("a reduction that is not bit-identical is not wrong —", _MUTED)
    line("it is not reproducible, which is a different property.", _MUTED)
    return out, y


def _descent_panel(bundle: EvidenceBundle, x: float, y: float) -> tuple[list[str], float]:
    """The float32 descent window: a number line of probed multiples, and a disagreement."""
    out: list[str] = []
    y += 8.0
    out.append(
        f"<text x='{x:.0f}' y='{y:.0f}' font-size='12' fill='{_INK}' font-weight='bold'>"
        "float32 descent window (g6)</text>"
    )
    y += 15.0

    def line(text: str, colour: str = "#333") -> None:
        nonlocal y
        out.append(
            f"<text x='{x:.0f}' y='{y:.0f}' font-size='10' fill='{colour}'>{_esc(text[:74])}</text>"
        )
        y += 12.5

    descent = bundle.descent
    if descent is None:
        line("absent — g6_descent_window.json was not read", "#8d6e63")
        line("this panel is drawn as absent rather than omitted", "#8d6e63")
        return out, y

    multiples = descent.multiples
    if multiples:
        step = min(46.0, 420.0 / max(len(multiples), 1))
        base_y = y + 6.0
        for index, multiple in enumerate(multiples):
            mark_x = x + 10.0 + index * step
            converging = [
                label for label, values in descent.per_config.items() if multiple in values
            ]
            total = len(descent.per_config) or 1
            if len(converging) == total:
                fill, edge = "#2e7d32", "#005005"
            elif converging:
                fill, edge = "#ffb74d", "#e65100"
            else:
                fill, edge = "#ffffff", "#c62828"
            detail = (
                f"{multiple:g}*u32 — converged on {len(converging)}/{total} configurations: "
                f"{', '.join(converging) if converging else 'none'}"
            )
            out.append(
                f"<circle cx='{mark_x:.1f}' cy='{base_y:.1f}' r='5' fill='{fill}' "
                f"stroke='{edge}' stroke-width='1.3'>"
                f"<title>{_esc(detail)}</title></circle>"
            )
            out.append(
                f"<text x='{mark_x:.1f}' y='{base_y + 16:.1f}' font-size='9' fill='{_MUTED}' "
                f"text-anchor='middle'>{multiple:g}</text>"
            )
        y = base_y + 30.0
        line("slack, in multiples of u32 — filled = every configuration", _MUTED)

    intersection = descent.intersection
    line(
        f"converges on every configuration: {{{', '.join(f'{v:g}' for v in intersection) or '—'}}}",
        "#2e7d32" if intersection else "#c62828",
    )
    line(f"prior CPU-device window: {{{', '.join(f'{v:g}' for v in descent.prior_cpu) or '—'}}}")
    if not descent.has_margin:
        line("ONE multiple wide — no margin on either side.", "#c62828")
        line("the multiple below stalls, the one above wanders.", "#c62828")
    if descent.disagrees_with_prior:
        line("DISAGREES with the prior CPU window — they share", "#c62828")
        line("no multiple. ALEPH-PORT-3609 §14.6 named three", "#c62828")
        line("repairs and took none; this figure takes none either.", "#c62828")
    if descent.prior_source:
        line(f"prior: {descent.prior_source}", _MUTED)
    return out, y


def evidence_figure(
    directory: Path | str | None = None,
    *,
    host: str | None = DEFAULT_HOST,
    refresh: bool = False,
    title: str | None = None,
) -> tuple[str, EvidenceBundle]:
    """Load the artefacts and render them.

    Args:
        directory: An explicit artefact directory, or ``None`` to discover one.
        host: ssh alias for the fallback fetch; ``None`` disables it.
        refresh: Re-fetch even when a local copy exists.
        title: Figure title.

    Returns:
        ``(svg_text, bundle)``, so a caller can assert on the model rather than on pixels.
    """
    bundle = load_bundle(directory, host=host, refresh=refresh)
    return render_evidence_svg(bundle, title=title), bundle


def _summary_json(bundle: EvidenceBundle, destination: Path) -> dict[str, object]:
    parity = bundle.parity
    return {
        "figure": str(destination),
        "source": bundle.source,
        "device": parity.device.get("name", ""),
        "finished_utc": parity.finished_utc,
        "cases": len(parity.cases),
        "channels": len(parity.rows),
        "gate_comparison_declared": parity.gate,
        "gate_comparison_recounted": parity.recounted,
        "blind_channels": [row.label() for row in parity.blind_rows],
        "thin_gates": [
            {"case": row.case, "field": row.field_name, "headroom": row.headroom}
            for row in sorted(parity.thin_budgets, key=lambda r: r.headroom or 0.0)
        ],
        "faults": parity.faults(),
        "missing_artefacts": [
            {"artefact": name, "reason": reason} for name, reason in bundle.missing
        ],
        "not_established": parity.not_established,
        "quantitative_status": "STRUCTURAL",
        "quantitative_status_reason": (
            "a ULP figure is evidence about an implementation; nothing here is evidence about a cell"
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    """The one command. Returns an exit code; never raises for an unreadable artefact."""
    parser = argparse.ArgumentParser(
        prog="python -m aleph.viz.evidence_figure",
        description=(
            "Render the CUDA law-parity evidence to an SVG figure. With no arguments it finds the "
            "artefacts locally or fetches them from the runner's host over ssh."
        ),
    )
    parser.add_argument(
        "--from", dest="source", default=None,
        help="artefact directory (default: ~/.aleph_runner/artefacts, else fetch over ssh)",
    )
    parser.add_argument(
        "-o", "--out", default="docs/results/2026-08-01-cuda-evidence/law_parity_channels.svg",
        help="output .svg path",
    )
    parser.add_argument(
        "--host", default=DEFAULT_HOST, help="ssh alias to fetch from when nothing is local"
    )
    parser.add_argument("--no-fetch", action="store_true", help="never fetch over ssh")
    parser.add_argument("--refresh", action="store_true", help="re-fetch even if a local copy exists")
    parser.add_argument("--title", default=None, help="figure title")
    parser.add_argument("--json", action="store_true", help="print a JSON summary on stdout")
    args = parser.parse_args(argv)

    try:
        svg, bundle = evidence_figure(
            args.source,
            host=None if args.no_fetch else args.host,
            refresh=args.refresh,
            title=args.title,
        )
    except EvidenceError as error:
        print(f"refused: {error}", file=sys.stderr)
        return EXIT_REFUSED

    destination = Path(args.out)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(svg, encoding="utf-8")

    parity = bundle.parity
    summary = _summary_json(bundle, destination)
    if args.json:
        print(json.dumps(summary, indent=2, default=str))
    else:
        print(
            f"wrote {destination} — {len(parity.rows)} channels over {len(parity.cases)} cases, "
            f"read from {bundle.source}"
        )
        print(
            f"  gates: {parity.gate.get('caught_on_cuda', '?')} caught on {parity.cuda_key}, "
            f"{parity.gate.get('caught_on_cpu', '?')} on {parity.cpu_key}, "
            f"{len(parity.wrong_kernels)} wrong kernels"
        )
        print(
            f"  {len(parity.blind_rows)} channels bit-identical to the true kernel — "
            "a per-case verdict would show none of them"
        )
        print(f"  {len(parity.thin_budgets)} gates with headroom < {THIN_HEADROOM:g}x")
        for name, reason in bundle.missing:
            print(f"  absent  {name}: {reason}")
        for fault in parity.faults():
            print(f"  FAULT  {fault}", file=sys.stderr)

    # A fault means an undeclared blind variant, a true-kernel regression, or a case that could not be
    # turned into a row. Any of those makes the figure a picture of something other than the run.
    return EXIT_REFUSED if parity.faults() else EXIT_OK


if __name__ == "__main__":  # pragma: no cover - exercised through `main` in the tests
    raise SystemExit(main())
