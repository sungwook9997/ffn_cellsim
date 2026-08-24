"""A self-contained interactive render of a cortical shell and its non-affinity field.

What this draws, and why it is not another picture of a cell
------------------------------------------------------------
`docs/results/2026-08-06-overlap-free-native-cortex` established that every areal modulus this
project has reported is the **affine bound**, and that the cortex's real modulus is 35× below it.
The difference between those two numbers is a *field*: for each node,

    n_i = |x_relaxed,i - x_affine,i| / (eps * R)

— how far the node went in a direction the imposed strain never asked for, in units of the affine
displacement itself. `n = 0` is a bundle of parallel springs. `n` of order one is a network
rearranging.

**That field is the finding, so it is what the picture is of.** A render of filaments is a render of
the construction; a render of `n` is a render of the mechanics. The affine/relaxed comparison is
what this lane measured, so this is not a prettier version of a picture that already existed.

(The sentence that stood here named the reference project by name, inside `aleph/`.
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` caught it,
which is the guard working: a package that has learned to say the provider's name in its own
docstrings is one where the boundary has already started to blur.)

Why a canvas and not SVG
------------------------
`cell_figure.py` draws whole-cell scenes as SVG with per-element classes, which is right when the
scene is thousands of elements and the point is that each one is inspectable. This scene is 62,330
filament segments at ρ = 8 and 780,000 at the sourced density; SVG stops being a document and
becomes a memory problem. The canvas path keeps the same painter's-algorithm depth ordering and the
same orthographic projection, and it can animate.

Everything is inlined — geometry as base64 `Float32Array`, no external anything — because the
publishing target forbids network requests and a silent fallback in a figure is worse than no
figure.
"""

from __future__ import annotations

import base64
import json
from typing import Any

import numpy as np

__all__ = ["cortex_field_comparison_page", "cortex_field_page", "encode_field"]


def encode_field(
    *,
    nodes: np.ndarray,
    nonaffinity: np.ndarray,
    offsets: np.ndarray,
    crosslinks: np.ndarray,
    nodes_per_filament: int,
    meta: dict[str, Any],
) -> dict[str, Any]:
    """Pack a shell into the payload the page decodes. Float32 throughout — this is a picture.

    The positions are already centred and divided by the shell radius by the exporter, so the page
    never needs a scale and the camera is the same for every density.
    """

    def b64(values, dtype):
        return base64.b64encode(np.ascontiguousarray(values, dtype=dtype).tobytes()).decode("ascii")

    return {
        "meta": dict(meta),
        "nodes_b64": b64(np.asarray(nodes).ravel(), np.float32),
        "nonaffinity_b64": b64(np.asarray(nonaffinity).ravel(), np.float32),
        "offset_b64": b64(np.asarray(offsets).ravel(), np.float32),
        "crosslinks_b64": b64(np.asarray(crosslinks).ravel(), np.int32),
        "nodes_per_filament": int(nodes_per_filament),
    }


def cortex_field_page(payload: dict[str, Any], *, title: str = "Aleph &mdash; the cortex, relaxed") -> str:
    """One HTML page, self-contained, no external requests.

    `payload` is what :func:`encode_field` produced, or what the exporter wrote to JSON.
    """
    return _render([("", payload)], title=title)


def cortex_field_comparison_page(
    payloads: list[tuple[str, dict[str, Any]]],
    *,
    title: str = "Aleph &mdash; the cortex, before and after",
) -> str:
    """The same scene under two constructions, on **one** camera, switchable in place.

    Two figures side by side are read as two pictures; one figure that *changes* is read as a
    difference — which is what a before/after is for. And the camera, the projection and the colour
    ramp are then provably identical between the two, because they are literally the same objects.

    ``payloads`` is ``[(label, payload), ...]``; the first entry is the control and is what the page
    opens on, so a reader who does not touch the switch sees the *old* behaviour rather than the new.

    Every variant's measurement panel is rendered into the page and shown by class, rather than
    built in JavaScript from the payload. The numbers therefore exist in the document whether or not
    the script runs, which is the property `tests/viz/test_cortex_field_page.py` audits.
    """
    if len(payloads) < 2:
        raise ValueError(f"a comparison needs at least two payloads; got {len(payloads)}")
    return _render(list(payloads), title=title)


def _render(variants: list[tuple[str, dict[str, Any]]], *, title: str) -> str:
    bundle = json.dumps({"variants": [{"label": label, "data": data} for label, data in variants]})
    panels = "\n".join(
        f'<dl class="variant-panel" data-variant="{index}"'
        f'{"" if index == 0 else " hidden"}>{_summary_rows(payload["meta"])}</dl>'
        for index, (_, payload) in enumerate(variants)
    )
    switch = ""
    if len(variants) > 1:
        buttons = "\n".join(
            f'<button class="seg var" data-variant="{index}" '
            f'aria-pressed="{"true" if index == 0 else "false"}">{label}</button>'
            for index, (label, _) in enumerate(variants)
        )
        switch = (
            '<div class="control"><h2>Construction</h2>'
            f'<div class="segs" role="group" aria-label="construction">{buttons}</div>'
            '<p class="note">One camera, one ramp, one projection. Only the shell changes.</p>'
            "</div>"
        )
    return (
        _TEMPLATE.replace("__TITLE__", title)
        .replace("__BUNDLE__", bundle)
        .replace("__SWITCH__", switch)
        .replace("__PANELS__", panels)
    )


def _summary_rows(meta: dict[str, Any]) -> str:
    def fmt(value, digits=4):
        if isinstance(value, bool):
            return "yes" if value else "no"
        if isinstance(value, float):
            return f"{value:,.{digits}g}"
        if isinstance(value, int):
            return f"{value:,}"
        return str(value)

    ratio = meta["k_affine_pn_per_um"] / meta["k_relaxed_pn_per_um"]
    rows = [
        ("filament areal density", f"{fmt(meta['density_per_um2'])} &micro;m&#8315;&sup2;"),
        ("filaments", fmt(meta["filaments"])),
        ("nodes", fmt(meta["nodes"])),
        ("crosslinks", fmt(meta["crosslinks"])),
        ("imposed dilation", f"{fmt(meta['strain'])}"),
        ("K<sub>A</sub> affine", f"{fmt(meta['k_affine_pn_per_um'], 6)} pN &micro;m&#8315;&sup1;"),
        ("K<sub>A</sub> relaxed", f"{fmt(meta['k_relaxed_pn_per_um'], 6)} pN &micro;m&#8315;&sup1;"),
        ("the affine bound is", f"{ratio:,.1f}&times; too stiff"),
        ("descent iterations", fmt(meta["iterations"])),
        ("peak residual force", f"{fmt(meta['peak_residual_pn'], 3)} pN"),
        ("converged", fmt(meta["converged"])),
        ("mesh-sized host copies", fmt(meta["bulk_readbacks_during_descent"])),
        ("device", str(meta["device"])),
        ("relaxation wall clock", f"{fmt(meta['relax_seconds'], 4)} s"),
    ]
    return "\n".join(
        f'<div class="row"><dt>{name}</dt><dd>{value}</dd></div>' for name, value in rows
    )


_TEMPLATE = r"""<title>__TITLE__</title>
<style>
  /* ---------------------------------------------------------------------------------------------
     A laboratory instrument, not a poster. The scene is a dark field because that is the vernacular
     of the subject -- a cortex is seen by fluorescence, light on dark, and a filament reads as an
     emitter rather than as ink. The panel is the opposite: paper, so numbers are read and not
     admired. The two grounds are deliberate and the page does not try to reconcile them.
     ------------------------------------------------------------------------------------------ */
  :root {
    --ink: #14161b;            /* near-black with a blue cast, so it sits under the field's indigo */
    --paper: #f4f3ef;
    --rule: #d9d6cc;
    --muted: #6c6a63;
    --field-0: #2b3a67;        /* affine: quiet, deep indigo                                       */
    --field-1: #d98d2b;        /* rearranging: warm amber. One accent, spent in one place.          */
    --panel-bg: var(--paper);
    --panel-fg: var(--ink);
  }
  :root[data-theme="dark"], :root:not([data-theme="light"]) {
    color-scheme: dark light;
  }
  @media (prefers-color-scheme: dark) {
    :root { --panel-bg: #1b1d23; --panel-fg: #e8e6df; --rule: #33363f; --muted: #9a978d; }
  }
  :root[data-theme="dark"] { --panel-bg: #1b1d23; --panel-fg: #e8e6df; --rule: #33363f; --muted: #9a978d; }
  :root[data-theme="light"] { --panel-bg: #f4f3ef; --panel-fg: #14161b; --rule: #d9d6cc; --muted: #6c6a63; }

  html, body { height: 100%; }
  body {
    margin: 0; background: #0b0c10; color: var(--panel-fg);
    font: 400 15px/1.5 ui-sans-serif, -apple-system, "Segoe UI", system-ui, sans-serif;
    overflow: hidden;
  }
  /* Explicit width/height, not `inset: 0` alone: a canvas without a CSS size falls back to its
     intrinsic 300x150, and clientWidth/clientHeight then report that instead of the viewport. */
  #stage { position: fixed; top: 0; left: 0; width: 100vw; height: 100vh;
           display: block; cursor: grab; touch-action: none; }
  #stage.dragging { cursor: grabbing; }

  .panel {
    position: fixed; top: 0; right: 0; bottom: 0; width: min(360px, 92vw);
    background: var(--panel-bg); color: var(--panel-fg);
    border-left: 1px solid var(--rule);
    display: flex; flex-direction: column; gap: 0;
    overflow-y: auto; overscroll-behavior: contain;
  }
  .panel > * { padding: 18px 22px; }
  .panel > * + * { border-top: 1px solid var(--rule); }

  h1 { margin: 0 0 6px; font-size: 19px; font-weight: 600; letter-spacing: -0.01em; text-wrap: balance; }
  .lede { margin: 0; color: var(--muted); font-size: 13.5px; }
  h2 { margin: 0 0 12px; font-size: 11px; font-weight: 600; letter-spacing: 0.09em;
       text-transform: uppercase; color: var(--muted); }

  dl { margin: 0; display: grid; gap: 5px; }
  /* The `hidden` attribute works through the UA stylesheet's `display: none`, and ANY author
     `display` rule beats it -- `dl { display: grid }` above is exactly such a rule. Without this
     the variant panels do not hide, they stack, and the reader sees variant 1's numbers under
     variant 4's button. Caught by screenshotting the page rather than trusting it. */
  [hidden] { display: none !important; }
  .row { display: grid; grid-template-columns: 1fr auto; gap: 12px; align-items: baseline; }
  dt { color: var(--muted); font-size: 12.5px; }
  dd { margin: 0; font-variant-numeric: tabular-nums; font-size: 13px;
       font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }

  .control { display: grid; gap: 9px; }
  label.check { display: flex; align-items: center; gap: 9px; font-size: 13.5px; cursor: pointer; }
  input[type=range] { width: 100%; accent-color: var(--field-1); }
  .slider-head { display: flex; justify-content: space-between; font-size: 12.5px; color: var(--muted); }
  .slider-head b { font-family: ui-monospace, Menlo, monospace; font-weight: 500; color: var(--panel-fg); }
  button.seg {
    font: inherit; font-size: 12.5px; padding: 5px 10px; border: 1px solid var(--rule);
    background: transparent; color: var(--panel-fg); cursor: pointer; border-radius: 2px;
  }
  button.seg[aria-pressed="true"] { background: var(--field-1); border-color: var(--field-1); color: #14161b; }
  .segs { display: flex; gap: 6px; flex-wrap: wrap; }
  button.seg:focus-visible, label.check:focus-within { outline: 2px solid var(--field-1); outline-offset: 2px; }

  .ramp { height: 9px; border-radius: 1px; margin: 8px 0 4px;
          background: linear-gradient(90deg, #687cb0, #a894a8, #f4a83e); }
  .ramp-ends { display: flex; justify-content: space-between; font-size: 11.5px; color: var(--muted);
               font-variant-numeric: tabular-nums; }
  .note { font-size: 12.5px; color: var(--muted); }
  .note b { color: var(--panel-fg); font-weight: 600; }
  @media (max-width: 720px) {
    .panel { top: auto; width: 100%; max-height: 52vh; border-left: none; border-top: 1px solid var(--rule); }
  }
  @media (prefers-reduced-motion: reduce) { * { animation: none !important; transition: none !important; } }
</style>

<canvas id="stage"></canvas>

<aside class="panel">
  <div>
    <h1>The cortex, relaxed</h1>
    <p class="lede">Colour is <b>non-affinity</b>: how far each node moved in a direction the imposed
    dilation never asked for, in units of the affine displacement. Zero is a bundle of parallel
    springs. Order one is a network rearranging.</p>
  </div>

  <div class="control">
    <h2>Non-affinity</h2>
    <div class="ramp"></div>
    <div class="ramp-ends"><span id="ramp-lo">0</span><span id="ramp-hi">&mdash;</span></div>
    <p class="note">Median <b id="na-med">&mdash;</b>, maximum <b id="na-max">&mdash;</b>.</p>
  </div>

  <div class="control">
    <h2>Rearrangement</h2>
    <div class="slider-head"><span>amplify displacement</span><b id="amp-val">0&times;</b></div>
    <input id="amp" type="range" min="0" max="400" value="0" aria-label="displacement amplification">
    <p class="note">The strain is 10&#8315;&sup3; &mdash; about 5&nbsp;nm on a 5&nbsp;&micro;m cell, invisible at any honest
    scale. Drag to amplify the <em>non-affine part only</em>; the affine part is never drawn, so what
    moves is exactly what the affine bound failed to predict.</p>
  </div>

  <div class="control">
    <h2>Layers</h2>
    <label class="check"><input type="checkbox" id="l-fil" checked> filaments</label>
    <label class="check"><input type="checkbox" id="l-link"> crosslinks</label>
    <label class="check"><input type="checkbox" id="l-shell"> shell horizon</label>
    <div class="segs" role="group" aria-label="colour by">
      <button class="seg" id="c-na" aria-pressed="true">non-affinity</button>
      <button class="seg" id="c-depth" aria-pressed="false">depth</button>
    </div>
  </div>

__SWITCH__

  <div>
    <h2>Measurement</h2>
__PANELS__
  </div>

  <div>
    <h2>Read this before believing the bands</h2>
    <p class="note">The high-non-affinity filaments fall in <b>bands</b>, and the bands are a
    property of the <em>construction</em>, not of the material. Measured at &rho;&nbsp;=&nbsp;100
    over 31,165 filaments: correlation with the golden-angle <em>depth</em> sequence
    <b>+0.048</b> and with the <em>orientation</em> sequence <b>+0.015</b> &mdash; neither explains
    them &mdash; but with <b>latitude, &minus;0.479</b>.</p>
    <p class="note">A uniform dilation of a sphere is perfectly isotropic, so <b>the load has no
    polar axis</b>. Only the Fibonacci lattice does. The golden-angle placement is isotropic in
    filament <em>density</em> and is <b>not</b> isotropic in mechanical <em>response</em>, and this
    field is how that becomes visible. It is a limit on the construction, found by looking at the
    picture.</p>
    <p class="note">Drag to rotate, scroll to zoom. Every number here is produced by
    <code>scripts/_ws_export_cortex_field.py</code>; nothing in this page is drawn from a value that
    was not measured.</p>
  </div>
</aside>

<script>
(() => {
  "use strict";
  const BUNDLE = __BUNDLE__;
  let ACTIVE = 0;
  let PAYLOAD = BUNDLE.variants[ACTIVE].data;

  const f32 = (b64) => {
    const bin = atob(b64), buf = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
    return new Float32Array(buf.buffer);
  };
  const i32 = (b64) => {
    const bin = atob(b64), buf = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
    return new Int32Array(buf.buffer);
  };

  // Rebound by selectVariant(). A comparison swaps the shell and NOTHING else -- the camera, the
  // ramp and the projection are the same objects throughout, so a difference on screen is a
  // difference in the data and cannot be a difference in how it was drawn.
  //
  // The ramp normalises on a high percentile rather than the maximum: the distribution is skewed
  // (median 0.28 against a max of 1.18 at rho = 100), so a max-normalised ramp would put every
  // typical filament in the bottom third and show one bright outlier against a flat field. The
  // panel says "98th pct", so the reader is told what the ramp's top means.
  let P, NA, OFF, LINK, PER, N, FIL, naMax, naHi, proj;

  function selectVariant(index) {
    ACTIVE = index;
    PAYLOAD = BUNDLE.variants[ACTIVE].data;
    P = f32(PAYLOAD.nodes_b64);                      // (N,3) centred, radius-normalised
    NA = f32(PAYLOAD.nonaffinity_b64);               // (N,)
    OFF = PAYLOAD.offset_b64 ? f32(PAYLOAD.offset_b64) : new Float32Array(P.length);
    LINK = i32(PAYLOAD.crosslinks_b64);              // (M,2) node indices
    PER = PAYLOAD.nodes_per_filament;
    N = P.length / 3;
    FIL = N / PER;
    proj = new Float32Array(N * 3);
    const sorted = Float32Array.from(NA).sort();
    naHi = sorted[Math.floor(sorted.length * 0.98)] || 1;
    naMax = PAYLOAD.meta.nonaffinity_max || naHi;
    document.getElementById("ramp-hi").textContent = naHi.toFixed(2) + " (98th pct)";
    document.getElementById("na-med").textContent = (PAYLOAD.meta.nonaffinity_median || 0).toFixed(3);
    document.getElementById("na-max").textContent = naMax.toFixed(3);
    document.querySelectorAll(".variant-panel").forEach((el) => {
      el.hidden = Number(el.dataset.variant) !== ACTIVE;
    });
    document.querySelectorAll("button.var").forEach((b) => {
      b.setAttribute("aria-pressed", String(Number(b.dataset.variant) === ACTIVE));
    });
  }

  // -- the ramp. Two ends and one waypoint; a rainbow would imply an ordering the data does not have.
  // Luminance rises along the whole ramp, because the field is drawn as light on dark and a dark
  // low end is not "quiet", it is invisible. The hue still carries the story -- cool where the
  // material followed the strain, warm where it did not -- and the accent is spent only at the top.
  const C0 = [104, 124, 176], C1 = [168, 148, 168], C2 = [244, 168, 62];
  const ramp = (t) => {
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    const [a, b, u] = t < 0.5 ? [C0, C1, t * 2] : [C1, C2, (t - 0.5) * 2];
    return [Math.round(a[0] + (b[0] - a[0]) * u),
            Math.round(a[1] + (b[1] - a[1]) * u),
            Math.round(a[2] + (b[2] - a[2]) * u)];
  };

  const canvas = document.getElementById("stage");
  const ctx = canvas.getContext("2d", { alpha: false });
  let W = 0, H = 0, dpr = 1;
  const view = { yaw: 0.6, pitch: 0.32, zoom: 1.0 };
  const state = { amp: 0, filaments: true, links: false, shell: false, colour: "na" };

  function resize() {
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    W = canvas.clientWidth; H = canvas.clientHeight;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
  }

  // Orthographic, like `cell_figure.py`: a perspective divide makes two filaments at the same depth
  // different sizes, and this picture is read for depth ordering rather than for distance.
  function project() {
    const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw);
    const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
    const a = state.amp;
    for (let i = 0, j = 0; i < N; i++, j += 3) {
      const x = P[j] + a * OFF[j], y = P[j + 1] + a * OFF[j + 1], z = P[j + 2] + a * OFF[j + 2];
      const x1 = cy * x + sy * z, z1 = -sy * x + cy * z;
      const y1 = cp * y - sp * z1, z2 = sp * y + cp * z1;
      proj[j] = x1; proj[j + 1] = y1; proj[j + 2] = z2;
    }
  }

  // Painter's algorithm on the primitive's own mean depth. Filaments are three nodes, so a
  // per-primitive key is not an approximation here -- it is the exact depth of its centre node.
  let order = null;
  function sortOrder(count, depthOf) {
    const idx = new Int32Array(count);
    for (let i = 0; i < count; i++) idx[i] = i;
    const key = new Float32Array(count);
    for (let i = 0; i < count; i++) key[i] = depthOf(i);
    return Array.prototype.sort.call(idx, (p, q) => key[p] - key[q]);
  }

  function draw() {
    project();
    ctx.fillStyle = "#0b0c10";
    ctx.fillRect(0, 0, W, H);
    // Centre in the space the panel leaves, not in the window: a scene centred behind an opaque
    // panel is a scene half of which cannot be seen.
    const panel = document.querySelector(".panel");
    const wide = window.matchMedia("(min-width: 721px)").matches;
    const free = wide ? W - panel.getBoundingClientRect().width : W;
    const freeH = wide ? H : H - panel.getBoundingClientRect().height;
    const s = Math.min(free, freeH) * 0.42 * view.zoom;
    const ox = free * 0.5, oy = freeH * 0.5;
    const X = (i) => ox + proj[i * 3] * s, Y = (i) => oy - proj[i * 3 + 1] * s;

    if (state.shell) {
      ctx.strokeStyle = "rgba(150,160,190,0.28)"; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(ox, oy, s, 0, Math.PI * 2); ctx.stroke();
    }

    if (state.links) {
      const M = LINK.length / 2;
      ctx.lineWidth = 0.5;
      ctx.strokeStyle = "rgba(120,132,168,0.20)";
      ctx.beginPath();
      for (let m = 0; m < M; m++) {
        const a = LINK[m * 2], b = LINK[m * 2 + 1];
        if (proj[a * 3 + 2] < 0 && proj[b * 3 + 2] < 0) continue;   // far side stays quiet
        ctx.moveTo(X(a), Y(a)); ctx.lineTo(X(b), Y(b));
      }
      ctx.stroke();
    }

    if (state.filaments) {
      order = sortOrder(FIL, (f) => proj[(f * PER + ((PER - 1) >> 1)) * 3 + 2]);
      ctx.lineCap = "round";
      for (let k = 0; k < FIL; k++) {
        const f = order[k], base = f * PER;
        const mid = base + ((PER - 1) >> 1);
        const depth = proj[mid * 3 + 2];
        let t;
        if (state.colour === "depth") t = (depth + 1) * 0.5;
        else t = NA[mid] / naHi;
        const [r, g, b] = ramp(t);
        // Depth cue in alpha, not in hue: the far side recedes without changing what colour means.
        const alpha = 0.42 + 0.55 * ((depth + 1) * 0.5);
        ctx.strokeStyle = `rgba(${r},${g},${b},${alpha.toFixed(3)})`;
        ctx.lineWidth = 1.05 + 1.05 * ((depth + 1) * 0.5);
        ctx.beginPath();
        ctx.moveTo(X(base), Y(base));
        for (let n = 1; n < PER; n++) ctx.lineTo(X(base + n), Y(base + n));
        ctx.stroke();
      }
    }

    // scale bar: the shell radius is 1 in these coordinates
    const barUm = 1.0, R = PAYLOAD.meta.radius_um;
    const px = s * (barUm / R);
    ctx.strokeStyle = "rgba(232,230,223,0.75)"; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(28, freeH - 30); ctx.lineTo(28 + px, freeH - 30); ctx.stroke();
    ctx.fillStyle = "rgba(232,230,223,0.75)";
    ctx.font = "12px ui-monospace, Menlo, monospace";
    ctx.fillText("1 \u00b5m", 28, freeH - 38);
  }

  // -- interaction ------------------------------------------------------------------------------
  let dragging = false, lastX = 0, lastY = 0;
  canvas.addEventListener("pointerdown", (e) => {
    dragging = true; lastX = e.clientX; lastY = e.clientY;
    canvas.classList.add("dragging"); canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!dragging) return;
    view.yaw += (e.clientX - lastX) * 0.006;
    view.pitch = Math.max(-1.5, Math.min(1.5, view.pitch + (e.clientY - lastY) * 0.006));
    lastX = e.clientX; lastY = e.clientY; draw();
  });
  const stop = (e) => { dragging = false; canvas.classList.remove("dragging"); };
  canvas.addEventListener("pointerup", stop);
  canvas.addEventListener("pointercancel", stop);
  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    view.zoom = Math.max(0.35, Math.min(6, view.zoom * (e.deltaY < 0 ? 1.08 : 1 / 1.08)));
    draw();
  }, { passive: false });

  const amp = document.getElementById("amp"), ampVal = document.getElementById("amp-val");
  amp.addEventListener("input", () => {
    state.amp = Number(amp.value) * 0.01;
    ampVal.textContent = (Number(amp.value)).toFixed(0) + "&times;";
    draw();
  });
  const bind = (id, key) => document.getElementById(id).addEventListener("change", (e) => {
    state[key] = e.target.checked; draw();
  });
  bind("l-fil", "filaments"); bind("l-link", "links"); bind("l-shell", "shell");
  const cna = document.getElementById("c-na"), cdep = document.getElementById("c-depth");
  const setColour = (mode) => {
    state.colour = mode;
    cna.setAttribute("aria-pressed", String(mode === "na"));
    cdep.setAttribute("aria-pressed", String(mode === "depth"));
    draw();
  };
  cna.addEventListener("click", () => setColour("na"));
  cdep.addEventListener("click", () => setColour("depth"));

  document.querySelectorAll("button.var").forEach((b) => {
    b.addEventListener("click", () => { selectVariant(Number(b.dataset.variant)); resize(); });
  });

  window.addEventListener("resize", resize);
  selectVariant(0);
  resize();
})();
</script>
"""
