#!/usr/bin/env python3
"""Build the Korean paper-style visual report for the Outer Library.

The report is deliberately generated only from committed result artifacts and a
deterministically generated real-image evidence panel.  It does not execute or
modify Aleph physics.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = ROOT / "output" / "pdf"
TMP = ROOT / "tmp" / "pdfs"
REAL_IMAGE_SOURCE = Path(
    "/Users/sw1/.codex/visualizations/2026/08/04/"
    "019fccbd-4129-78f1-a56a-626a3addcbce/outer-library-real-image-evidence.png"
)
REAL_IMAGE_META = REAL_IMAGE_SOURCE.with_suffix(".json")

PDF_PATH = OUT / "aleph_outer_library_visual_paper_2026-08-07.pdf"
IMAGE_COPY = OUT / "outer_library_real_image_evidence.png"

NAVY = "#14213d"
BLUE = "#2563eb"
CYAN = "#0891b2"
GREEN = "#15803d"
AMBER = "#b45309"
RED = "#b91c1c"
INK = "#172033"
MUTED = "#64748b"
PALE = "#f1f5f9"
LINE = "#cbd5e1"


def load_json(name: str):
    with (RESULTS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


def setup_fonts():
    font_path = "/System/Library/Fonts/Supplemental/AppleGothic.ttf"
    pdfmetrics.registerFont(TTFont("Korean", font_path))
    font_manager.fontManager.addfont(font_path)
    return font_manager.FontProperties(fname=font_path)


KFONT = setup_fonts()
plt.rcParams.update(
    {
        "font.family": KFONT.get_name(),
        "axes.unicode_minus": False,
        "figure.dpi": 180,
        "savefig.dpi": 220,
    }
)


def savefig(fig, name: str) -> Path:
    path = TMP / name
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def architecture_figure() -> Path:
    fig, ax = plt.subplots(figsize=(12.5, 5.2))
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 5.2)
    ax.axis("off")
    columns = [
        (0.2, 3.75, 1.75, 1.0, "공개 실험", "IF / BF / calcium\nRNA / CyTOF / TFM / PIV", NAVY),
        (2.25, 3.75, 1.9, 1.0, "Observation Compiler", "단위·protocol·provenance\nreplicate·uncertainty", BLUE),
        (4.5, 3.75, 1.9, 1.0, "Modality Encoder", "CNN / temporal CNN\nMLP / operator features", CYAN),
        (6.75, 3.75, 1.75, 1.0, "Task Head", "분류·회귀·direction\nspatial prediction", GREEN),
        (8.85, 3.75, 1.7, 1.0, "Evidence Gate", "calibration·OOD\nholdout·refusal", AMBER),
        (10.9, 3.75, 1.35, 1.0, "출력", "posterior\n한계·다음 assay", NAVY),
    ]
    for x, y, w, h, title, body, color in columns:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.08", facecolor=color, edgecolor="none"))
        ax.text(x + w / 2, y + 0.67, title, ha="center", va="center", color="white", fontsize=10.5, fontweight="bold")
        ax.text(x + w / 2, y + 0.28, body, ha="center", va="center", color="white", fontsize=7.7, linespacing=1.25)
    for i in range(len(columns) - 1):
        x1 = columns[i][0] + columns[i][2]
        x2 = columns[i + 1][0]
        ax.annotate("", xy=(x2 - 0.05, 4.25), xytext=(x1 + 0.05, 4.25), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.7))

    heads = [
        (0.5, 2.05, "영상", "위치·형태·channel 관계"),
        (3.0, 2.05, "시간 신호", "peak·duration·waveform"),
        (5.5, 2.05, "분자 상태", "gene/pathway·cell type"),
        (8.0, 2.05, "힘·흐름", "traction·velocity·tension"),
    ]
    for x, y, title, body in heads:
        ax.add_patch(FancyBboxPatch((x, y), 2.1, 0.78, boxstyle="round,pad=0.03", facecolor="#e2e8f0", edgecolor=LINE))
        ax.text(x + 1.05, y + 0.5, title, ha="center", fontsize=9.2, color=INK, fontweight="bold")
        ax.text(x + 1.05, y + 0.19, body, ha="center", fontsize=7.6, color=MUTED)
    ax.text(6.25, 1.35, "현재는 공유 latent 하나가 아니라 task-scoped late fusion", ha="center", fontsize=10.5, color=RED, fontweight="bold")
    ax.text(6.25, 0.86, "이유: modality 간 동일 샘플 정렬이 부족해 universal latent를 만들면 leakage와 가짜 결합 위험이 큼", ha="center", fontsize=8.7, color=MUTED)
    ax.add_patch(FancyBboxPatch((3.25, 0.12), 6.0, 0.45, boxstyle="round,pad=0.03", facecolor="#fee2e2", edgecolor="#fecaca"))
    ax.text(6.25, 0.34, "권한 방화벽: 외부 근거 → 보고 가능 | 숫자 Aleph 파라미터·physics mutation → 현재 차단", ha="center", fontsize=9.3, color=RED, fontweight="bold")
    return savefig(fig, "figure_architecture.png")


def corpus_figure(audit: dict) -> Path:
    per = audit["per_modality"]
    groups = {
        "시간 신호": ["live_cell_calcium_imaging_derived", "live_cell_FRET_imaging_derived"],
        "영상/형태": ["brightfield_single_cell_patch_derived", "imaging_flow_cytometry_multichannel", "immunofluorescence_derived", "multiscale_live_cell_imaging_derived"],
        "분자 상태": ["single_cell_mass_cytometry_CyTOF_nontransformed", "bulk_RNA_expression_derived", "bulk_RNA_seq_author_processed", "bulk_RNA_seq_derived", "RT_qPCR_derived"],
        "힘/흐름": ["traction_force_microscopy_derived", "particle_image_velocimetry_derived", "particle_image_velocimetry_and_live_actin_derived", "traction_force_and_monolayer_stress_microscopy_derived", "FLIM_FRET_vinculin_tension_sensor_derived"],
        "이동/조직": ["live_cell_tracking_derived", "author_supervised_migration_mode_derived", "3D_CDM_cell_migration_trajectory_derived"],
    }
    values = [sum(per.get(k, {}).get("observations", 0) for k in keys) for keys in groups.values()]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.1), gridspec_kw={"width_ratios": [1.2, 1]})
    order = list(groups)
    colors_ = [BLUE, CYAN, GREEN, AMBER, NAVY]
    ax.barh(order[::-1], values[::-1], color=colors_[::-1])
    ax.set_xlabel("대표 modality 묶음의 canonical observations (log scale)")
    ax.set_xscale("log")
    ax.grid(axis="x", alpha=0.2)
    for idx, value in enumerate(values[::-1]):
        ax.text(value * 1.06, idx, f"{value:,}", va="center", fontsize=8.5, color=INK)
    ax.set_title("데이터는 하나의 이미지셋이 아니라 이질적 관측의 집합", loc="left", color=INK, fontweight="bold")

    cards = [
        ("38", "datasets"),
        ("37", "conservative lab groups"),
        ("193,832", "canonical samples"),
        ("329,303", "canonical observations"),
        ("30", "task heads"),
        ("9", "trained encoder heads"),
    ]
    ax2.axis("off")
    for i, (number, label) in enumerate(cards):
        row, col = divmod(i, 2)
        x, y = 0.03 + col * 0.49, 0.72 - row * 0.3
        ax2.add_patch(FancyBboxPatch((x, y), 0.44, 0.23, transform=ax2.transAxes, boxstyle="round,pad=0.02", facecolor=PALE, edgecolor=LINE))
        ax2.text(x + 0.22, y + 0.145, number, transform=ax2.transAxes, ha="center", fontsize=19, color=BLUE, fontweight="bold")
        ax2.text(x + 0.22, y + 0.055, label, transform=ax2.transAxes, ha="center", fontsize=8.2, color=MUTED)
    ax2.text(0.03, 0.98, "현재 자산", transform=ax2.transAxes, va="top", fontsize=12, color=INK, fontweight="bold")
    ax2.text(0.03, 0.02, "HPA 81,007 embedding images는 물리 관측과 분리 집계", transform=ax2.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    return savefig(fig, "figure_corpus.png")


def status_figure(multitask: dict) -> Path:
    ready = multitask["readiness"]
    labels = ["외부 근거 통과", "외부 classifier 차단", "training-only", "기타 제한 head"]
    values = [len(ready["passed_external_evidence_heads"]), len(ready["blocked_external_classifier_heads"]), len(ready["training_only_heads"]), 2]
    colors_ = [GREEN, RED, AMBER, MUTED]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [0.9, 1.5]})
    ax.bar(labels, values, color=colors_)
    ax.set_ylim(0, 22)
    ax.set_ylabel("head 수")
    ax.tick_params(axis="x", rotation=22)
    ax.grid(axis="y", alpha=0.2)
    for i, v in enumerate(values):
        ax.text(i, v + 0.6, str(v), ha="center", fontweight="bold", color=INK)
    ax.set_title("30개 head의 현재 권한", loc="left", fontweight="bold", color=INK)
    ax2.axis("off")
    rows = [
        ("PASS", "Piezo1 calcium mechanism", "외부 lab 방향성"),
        ("PASS", "paired migration state", "retrospective 외부 성능"),
        ("PASS", "fibrotic-state mechanism", "조건 수준 방향성"),
        ("PASS", "senescence transcriptomic direction", "6/6 방향, absolute state 거부"),
        ("BLOCK", "calcium cell state", "macro-F1 0.550"),
        ("BLOCK", "PBMC cell type/IFN", "class coverage 실패"),
        ("BLOCK", "raw IF localization transfer", "macro-F1 0.306"),
        ("BLOCK", "label-free organelle prediction", "Pearson/OOD 실패"),
    ]
    for i, (state, name, why) in enumerate(rows):
        y = 0.92 - i * 0.115
        color = GREEN if state == "PASS" else RED
        ax2.text(0.02, y, state, transform=ax2.transAxes, color=color, fontweight="bold", fontsize=9)
        ax2.text(0.16, y, name, transform=ax2.transAxes, color=INK, fontsize=9, fontweight="bold")
        ax2.text(0.63, y, why, transform=ax2.transAxes, color=MUTED, fontsize=8.2)
    ax2.text(0.02, 0.02, "PASS는 해당 task의 제한된 외부 근거이며 universal cell-state model의 통과가 아님", transform=ax2.transAxes, color=RED, fontsize=8.3, fontweight="bold")
    fig.tight_layout()
    return savefig(fig, "figure_head_status.png")


def validation_figure() -> Path:
    fig, ax = plt.subplots(figsize=(12.5, 5.4))
    ax.axis("off")
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 5.4)
    stages = [
        (0.2, "1. 동결", "endpoint·threshold\nsplit·no-replacement"),
        (2.55, "2. 분리", "dataset / lab / donor\nreplicate / compound"),
        (4.9, "3. 예측", "checkpoint 고정\n외부 label 미사용"),
        (7.25, "4. 불확실성", "conformal coverage\ncalibration·class coverage"),
        (9.6, "5. 거부", "OOD·integrity·scope\nauthority firewall"),
    ]
    for i, (x, title, body) in enumerate(stages):
        c = [NAVY, BLUE, CYAN, AMBER, RED][i]
        ax.add_patch(FancyBboxPatch((x, 3.45), 2.05, 1.18, boxstyle="round,pad=0.04", facecolor=c, edgecolor="none"))
        ax.text(x + 1.025, 4.25, title, ha="center", color="white", fontsize=10.5, fontweight="bold")
        ax.text(x + 1.025, 3.78, body, ha="center", color="white", fontsize=8.2, linespacing=1.3)
        if i < len(stages) - 1:
            ax.annotate("", xy=(x + 2.32, 4.04), xytext=(x + 2.08, 4.04), arrowprops=dict(arrowstyle="->", lw=1.6, color=MUTED))
    examples = [
        (0.5, 1.65, 3.55, "PBMC", "AUROC/F1가 보여도 class-conditional coverage 실패 → 차단", RED),
        (4.48, 1.65, 3.55, "Light My Cells", "SSIM·MAE·coverage 일부 통과, Pearson/OOD 실패 → 차단", RED),
        (8.47, 1.65, 3.55, "IDR0168", "공식 archive CRC 오류 → tensor 0, prediction 0, 영구 소진", RED),
    ]
    for x, y, w, title, body, c in examples:
        ax.add_patch(FancyBboxPatch((x, y), w, 0.95, boxstyle="round,pad=0.04", facecolor="#fff7ed", edgecolor="#fdba74"))
        ax.text(x + 0.18, y + 0.64, title, ha="left", color=c, fontsize=10, fontweight="bold")
        ax.text(x + 0.18, y + 0.27, body, ha="left", color=INK, fontsize=8.1)
    ax.text(6.25, 0.72, "중요: 높은 AUROC 하나, 좋은 예시 이미지 하나, nominal coverage 하나만으로는 통과하지 않는다.", ha="center", fontsize=10.5, color=INK, fontweight="bold")
    ax.text(6.25, 0.30, "실패는 삭제하지 않고 consumed diagnostic 또는 source-integrity refusal로 보존한다.", ha="center", fontsize=9, color=MUTED)
    return savefig(fig, "figure_validation.png")


def axes_figure(audit: dict) -> Path:
    rows = [
        ("Mechanosensitive calcium", "부분 외부 근거", GREEN),
        ("Migration / ROCK", "context-dependent", AMBER),
        ("Fibrotic / myofibroblast", "2-lab direction", GREEN),
        ("Senescence", "direction only", GREEN),
        ("EMT", "direction holdout", AMBER),
        ("Cell cycle", "single-provider", AMBER),
        ("Mito stress", "single-provider", AMBER),
        ("Apoptosis / death", "single-provider", AMBER),
        ("Immune / inflammation", "single annotated replicate", AMBER),
        ("Hypoxia", "MISSING", RED),
        ("Differentiation / stemness", "MISSING", RED),
        ("DNA damage / genotoxic", "MISSING", RED),
        ("Infection response", "MISSING", RED),
    ]
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    ax.axis("off")
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, len(rows) + 1)
    ax.text(0.1, len(rows) + 0.55, "Cell-state axis coverage: 넓지만 아직 완전하지 않음", fontsize=14, color=INK, fontweight="bold")
    for i, (name, status, c) in enumerate(rows):
        y = len(rows) - i - 0.1
        ax.add_patch(FancyBboxPatch((0.15, y - 0.55), 6.1, 0.62, boxstyle="round,pad=0.02", facecolor=PALE, edgecolor="none"))
        ax.text(0.35, y - 0.24, name, va="center", fontsize=9, color=INK)
        ax.add_patch(FancyBboxPatch((6.55, y - 0.55), 3.65, 0.62, boxstyle="round,pad=0.02", facecolor=c, edgecolor="none", alpha=0.92))
        ax.text(8.38, y - 0.24, status, va="center", ha="center", fontsize=8.5, color="white", fontweight="bold")
    return savefig(fig, "figure_axes.png")


class PaperDocTemplate(BaseDocTemplate):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, pagesize=A4, rightMargin=16 * mm, leftMargin=16 * mm, topMargin=17 * mm, bottomMargin=17 * mm, **kwargs)
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="normal")
        self.addPageTemplates(PageTemplate(id="paper", frames=frame, onPage=self._page))

    @staticmethod
    def _page(canvas, doc):
        canvas.saveState()
        canvas.setFont("Korean", 7.5)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.drawString(16 * mm, 9 * mm, "ALEPH OUTER LIBRARY · VISUAL PAPER · 2026-08-07")
        canvas.drawRightString(A4[0] - 16 * mm, 9 * mm, str(doc.page))
        canvas.setStrokeColor(colors.HexColor(LINE))
        canvas.line(16 * mm, 13 * mm, A4[0] - 16 * mm, 13 * mm)
        canvas.restoreState()


_ACTIVE_STYLES = None


def styles():
    global _ACTIVE_STYLES
    ss = getSampleStyleSheet()
    base = ParagraphStyle("BaseKR", parent=ss["BodyText"], fontName="Korean", fontSize=9.2, leading=14.3, textColor=colors.HexColor(INK), spaceAfter=6)
    result = {
        "body": base,
        "small": ParagraphStyle("SmallKR", parent=base, fontSize=7.7, leading=11.4, textColor=colors.HexColor(MUTED)),
        "caption": ParagraphStyle("CaptionKR", parent=base, fontSize=7.8, leading=11.2, textColor=colors.HexColor(MUTED), spaceBefore=3, spaceAfter=8),
        "h1": ParagraphStyle("H1KR", parent=base, fontSize=18, leading=23, textColor=colors.HexColor(NAVY), spaceBefore=8, spaceAfter=10),
        "h2": ParagraphStyle("H2KR", parent=base, fontSize=13, leading=17, textColor=colors.HexColor(BLUE), spaceBefore=8, spaceAfter=6),
        "h3": ParagraphStyle("H3KR", parent=base, fontSize=10.5, leading=14, textColor=colors.HexColor(INK), spaceBefore=5, spaceAfter=4),
        "title": ParagraphStyle("TitleKR", parent=base, fontSize=24, leading=32, textColor=colors.white, alignment=TA_LEFT, spaceAfter=10),
        "subtitle": ParagraphStyle("SubtitleKR", parent=base, fontSize=11, leading=17, textColor=colors.HexColor("#dbeafe")),
        "callout": ParagraphStyle("CalloutKR", parent=base, fontSize=10.2, leading=15.5, textColor=colors.HexColor(NAVY), leftIndent=8, rightIndent=8, borderColor=colors.HexColor("#93c5fd"), borderWidth=0.8, borderPadding=8, backColor=colors.HexColor("#eff6ff"), spaceBefore=5, spaceAfter=9),
        "warning": ParagraphStyle("WarningKR", parent=base, fontSize=9.4, leading=14.5, textColor=colors.HexColor(RED), leftIndent=8, rightIndent=8, borderColor=colors.HexColor("#fecaca"), borderWidth=0.8, borderPadding=8, backColor=colors.HexColor("#fef2f2"), spaceBefore=5, spaceAfter=9),
        "toc": ParagraphStyle("TOCKR", parent=base, fontSize=10, leading=18, leftIndent=10),
    }
    _ACTIVE_STYLES = result
    return result


def P(text: str, st: dict | None = None, name: str = "body") -> Paragraph:
    if st is None:
        if _ACTIVE_STYLES is None:
            raise RuntimeError("styles() must be initialized before creating paragraphs")
        st = _ACTIVE_STYLES
    return Paragraph(text, st[name])


def bullet(text: str, st: dict) -> Paragraph:
    return Paragraph(f"• {text}", st["body"])


def table(data, widths, header=True, font_size=7.7):
    body_style = ParagraphStyle(
        "TableBodyKR",
        fontName="Korean",
        fontSize=font_size,
        leading=font_size * 1.35,
        textColor=colors.HexColor(INK),
        wordWrap="CJK",
    )
    header_style = ParagraphStyle(
        "TableHeaderKR",
        parent=body_style,
        textColor=colors.white,
        leading=font_size * 1.3,
    )
    wrapped = []
    for row_idx, row in enumerate(data):
        cell_style = header_style if header and row_idx == 0 else body_style
        wrapped.append(
            [cell if isinstance(cell, Paragraph) else Paragraph(str(cell), cell_style) for cell in row]
        )
    t = Table(wrapped, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    commands = [
        ("FONT", (0, 0), (-1, -1), "Korean", font_size),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor(LINE)),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1), [colors.white, colors.HexColor("#f8fafc")]),
    ]
    if header:
        commands.extend([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(NAVY)), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white)])
    t.setStyle(TableStyle(commands))
    return t


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    if not REAL_IMAGE_SOURCE.exists():
        raise FileNotFoundError(f"Real-image panel missing: {REAL_IMAGE_SOURCE}")
    shutil.copy2(REAL_IMAGE_SOURCE, IMAGE_COPY)

    audit = load_json("adversarial_coverage_audit.json")
    multitask = load_json("multitask_outer_model_report.json")
    sweep = load_json("sweep_gate_report.json")
    idr = load_json("idr0168_confirmation_decode_refusal.json")
    real_meta = json.loads(REAL_IMAGE_META.read_text(encoding="utf-8"))

    figs = {
        "architecture": architecture_figure(),
        "corpus": corpus_figure(audit),
        "status": status_figure(multitask),
        "validation": validation_figure(),
        "axes": axes_figure(audit),
    }
    st = styles()
    story = []

    # Cover
    cover = Table(
        [[P("ALEPH OUTER LIBRARY", st, "subtitle")], [P("공개 멀티모달 세포 실험으로부터<br/>검증 가능한 생물물리 상태를 추론하는 외부 신경망", st, "title")], [P("현재 구현, 실제 영상 판독, 외부 검증, 불확실성, OOD 거부 및 Aleph 연결 경계", st, "subtitle")], [Spacer(1, 15 * mm)], [P("연구 보고서 / 논문 초안<br/>2026년 8월 7일 · 발표 목표 2026년 8월 10일", st, "subtitle")]],
        colWidths=[178 * mm], rowHeights=[16 * mm, 55 * mm, 28 * mm, 22 * mm, 28 * mm],
    )
    cover.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(NAVY)), ("LEFTPADDING", (0, 0), (-1, -1), 13 * mm), ("RIGHTPADDING", (0, 0), (-1, -1), 13 * mm), ("TOPPADDING", (0, 0), (-1, -1), 7 * mm), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
    story += [cover, Spacer(1, 12 * mm), P("핵심 결론", st, "h2"), P("Outer Library는 현재 <b>범용 세포종/세포상태 모델이나 Aleph 파라미터 추정기</b>가 아니다. 그러나 공개 실험을 provenance-aware canonical observation으로 바꾸고, modality별 신경망과 theory/operator head를 학습하며, 외부 holdout·calibration·OOD·source integrity를 함께 판정하는 <b>task-scoped evidence runtime</b>은 작동한다. 외부 근거로 통과한 head 4개와 차단된 classifier head 4개가 동시에 존재한다. 이 비대칭이 현재 시스템의 가장 중요한 과학적 결과다.", st, "callout"), P("Aleph의 전면 파라미터 스윕이 아직 계산적으로 불가능하여 독립 외부 추론망을 우선 구축했다. 목적은 포기된 것이 아니라 순서와 권한이 분리되었다. 검증된 observation operator와 실행 가능한 sweep이 생기기 전까지 숫자 Aleph 파라미터를 출력하지 않는다.", st, "warning"), PageBreak()]

    # Abstract + contents
    story += [P("초록", st, "h1"), P("세포 역학과 상태 추론에는 IF, transmitted-light, calcium/FRET time series, RNA/CyTOF, TFM/PIV처럼 단위와 오류 구조가 다른 실험이 동시에 필요하다. 본 연구는 38개 공개 데이터셋과 37개 보수적 lab/provider group에서 193,832개 canonical sample과 329,303개 canonical observation을 구성하고, 30개 task-specific head를 갖는 late-fusion 외부 추론 시스템을 구현했다. 모델은 영상의 subcellular localization과 organelle spatial prediction, 시간 신호의 기전 방향, migration/fibrotic/senescence state, force-flow observation을 서로 다른 encoder와 split semantics로 처리한다. 4개 head가 제한된 외부 evidence endpoint를 통과했으나 calcium per-cell state, PBMC uncertainty, raw IF cross-provider transfer, label-free organelle transfer는 차단되었다. IDR0168 single-use IF confirmation은 공식 archive의 CRC/Zarr integrity 실패로 prediction 0인 상태에서 영구 소진되었다. 따라서 현재 결과는 task-scoped 외부 생물물리 evidence runtime의 proof-of-system이며, universal cell-state model, 임상 도구, 또는 Aleph numeric sweep authority의 증거는 아니다.", st), P("키워드: multimodal cell state, immunofluorescence, label-free imaging, OOD detection, conformal prediction, mechanobiology, provenance, task-scoped neural network", st, "small"), P("목차", st, "h2"), P("0. 처음 보는 사람을 위한 개념 · 1. 연구 동기와 목표 · 2. 데이터 언어와 말뭉치 · 3. 신경망 구조 · 4. 실제 영상 판독 · 5. modality별 판단 범위 · 6. 검증 설계 · 7. 현재 결과 · 8. 실패와 거부 · 9. cell-state coverage · 10. Aleph 연결 · 11. 논문 가능성 · 12. 다음 단계 · 재현성·참고문헌", st, "toc"), PageBreak()]

    # Primer for a first-time reader
    story += [P("0. 처음 보는 사람을 위한 핵심 개념", st, "h1"), P("0.1 Aleph란 무엇인가", st, "h2"), P("Aleph는 세포막, cortex, cytosol, adhesion, 외부 유체와 같은 구성요소를 물리 법칙으로 계산해 세포의 위치·힘·흐름·형태가 시간에 따라 어떻게 변하는지 예측하려는 whole-cell simulation framework다. 입력은 ‘세포가 얼마나 수축하는가’, ‘막과 cortex가 얼마나 붙어 있는가’, ‘유체의 점도는 얼마인가’ 같은 숫자이며, 출력은 node 위치, force, velocity, shape 같은 simulation observable이다."), P("0.2 파라미터란 무엇인가", st, "h2"), P("<b>파라미터(parameter)</b>는 시뮬레이션의 물질·기전 특성을 정하는 숫자다. 예를 들어 cortical tension, membrane tension, viscosity, adhesion strength, actin turnover rate가 해당한다. <b>State</b>는 특정 순간의 node 위치나 calcium level처럼 시간에 따라 변하는 값이고, <b>observable</b>은 실험과 직접 비교할 수 있는 traction, velocity, fluorescence, cell area 같은 값이다."), table([["용어", "쉬운 정의", "예시"], ["Parameter", "시뮬레이터의 성질을 정하는 손잡이", "cortical tension, viscosity, adhesion"], ["State", "그 손잡이 아래에서 시간에 따라 변하는 내부 상태", "node position, actin activity"], ["Observable", "실험과 같은 단위로 비교할 수 있는 출력", "TFM traction, PIV velocity, IF intensity"], ["Mechanism", "관측을 만든 가능한 원인", "actomyosin 증가, adhesion 감소"], ["Uncertainty", "하나의 답 대신 가능한 범위와 확률", "95% interval, conformal set"]], [30*mm, 86*mm, 62*mm]), P("0.3 Parameter sweep은 왜 필요한가", st, "h2"), P("실험에서 보이는 결과가 어느 파라미터 조합에서 나오는지 알기 위해 여러 값을 반복 실행한다. 예를 들어 cortical tension을 10개 값, adhesion을 10개 값, viscosity를 10개 값으로 시험하면 이미 1,000개 조합이다. 축이 177개이면 모든 조합을 전부 도는 Cartesian sweep은 사실상 불가능하다. 더 큰 문제는 서로 다른 파라미터 조합이 같은 cell shape를 만들 수 있다는 <b>비식별성(non-identifiability)</b>이다."), P("0.4 외부 신경망이란 무엇인가", st, "h2"), P("외부 신경망은 Aleph simulation 내부에서 node를 업데이트하는 모델이 아니다. 공개 논문의 실제 이미지·곡선·힘·유전자 데이터를 읽어 ‘어떤 상태와 기전이 가능한지’를 추론하는 별도 시스템이다. 지금은 Aleph 파라미터 숫자를 출력하지 않고, 관측된 state/mechanism의 확률, 예측 가능한 실험량, 불확실성, OOD 여부, 다음에 필요한 실험을 출력한다."), P("실험 데이터 → 외부 신경망 → 가능한 상태/기전 posterior → (미래) Aleph 후보 sweep", st, "callout"), PageBreak()]

    story += [P("0.5 이미지 한 장은 어떻게 판단으로 바뀌는가", st, "h2"), table([["순서", "질문", "처리", "출력"], ["1. 입력 무결성", "파일과 channel이 맞는가?", "hash, CRC, metadata, channel audit", "accept / refusal"], ["2. 관측 변환", "픽셀에서 무엇을 측정하는가?", "normalization, texture, localization, morphology, vector field", "canonical observation"], ["3. Encoder", "어떤 반복 패턴이 있는가?", "CNN/U-Net/temporal CNN/feature encoder", "task-specific representation"], ["4. Head", "어떤 label/quantity를 예측하는가?", "classification, regression, contrast, spatial decoder", "probability/map/value"], ["5. 외부 검증", "새 lab에서도 맞는가?", "held-out F1/AP/Pearson/coverage/OOD", "pass / blocked"], ["6. 권한", "결과로 무엇을 말해도 되는가?", "evidence and Aleph firewall", "scoped claim / refusal"]], [22*mm, 39*mm, 68*mm, 49*mm]), P("예를 들어 IF에서는 DAPI가 핵의 위치를 알려주고 protein channel이 핵 내부·세포질·막·filament·puncta 중 어디에 분포하는지를 보여준다. CNN은 이 공간 패턴을 확률로 바꾸지만, 새 현미경/provider에서 같은 성능이 유지되는지는 별도의 외부 검증 문제다. TFM/PIV처럼 이미 힘·속도 field인 자료는 pixel localization CNN을 거치지 않고 물리 observation operator로 직접 정규화한다."), P("0.6 이 보고서의 색상 규칙", st, "h2"), table([["표기", "뜻"], ["PASS / 녹색", "명시된 task와 외부 endpoint 범위에서만 근거로 사용 가능"], ["TRAINING-ONLY / 주황", "representation 또는 same-provider 학습 자산; 외부 일반화 주장 불가"], ["BLOCK / 빨강", "frozen endpoint, uncertainty, OOD 또는 integrity 실패"], ["ALEPH AUTHORITY: NONE", "numeric parameter range, 축 제거, physics mutation 모두 금지"]], [52*mm, 126*mm]), P("이후 장의 모든 성능 숫자는 이 권한 규칙과 함께 읽어야 한다. 정확도가 높다는 사실만으로 parameter가 식별되거나 Aleph sweep이 줄어드는 것은 아니다.", st, "warning"), PageBreak()]

    # Motivation
    story += [P("1. 연구 동기와 목표", st, "h1"), P("원래 목표는 공개 실험과 ffn_cellsim reference/tag 정보를 이용해 Aleph의 177개 축을 좁히고 parameter sweep 비용을 줄이는 것이었다. 하지만 현재 Aleph는 고해상도 whole-cell sweep을 수행할 계산 비용과 memory scaling 문제가 해결 중이므로, 외부 신경망이 곧바로 Aleph 파라미터를 학습했다는 주장은 검증할 수 없다."), P("따라서 현재의 실행 목표는 두 층으로 분리된다."), bullet("<b>현재 층:</b> 외부 실험만으로 관측 가능한 상태, 기전 방향, 불확실성, 적용 범위, 추가 실험을 추론한다.", st), bullet("<b>조건부 미래 층:</b> Aleph forward simulation이 가능해지고 동일 observation operator가 검증되면 posterior/ambiguity set을 pre-sweep constraint로 시험한다.", st), P("이 분리는 소극적 축소가 아니다. 오히려 이미지, 힘, 흐름, 분자 상태, 시간 신호를 모두 다루되, 각 결과가 실제로 말할 수 있는 범위를 좁게 표기하는 공격적인 확장 전략이다.", st, "callout"), Image(str(figs["architecture"]), width=178 * mm, height=74 * mm), P("그림 1. 전체 시스템. 동일 환자/세포의 완전 정렬이 부족하므로 공유 latent 하나 대신 task-scoped late fusion을 사용한다. evidence gate와 Aleph authority firewall은 모델 출력의 일부다.", st, "caption"), PageBreak()]

    # Corpus
    story += [P("2. 우리의 데이터 언어와 말뭉치", st, "h1"), P("‘우리 언어’는 자연어 tag 목록이 아니라, 어떤 실험도 공통 최소 단위로 기록하는 observation schema다. 각 row에는 source dataset, lab/provider group, cell line/type, treatment, dose, time, geometry, modality, measured quantity, unit, biological/technical replicate, author label, transformation, uncertainty, licence와 content hash가 들어간다. 이 구조가 없으면 서로 다른 논문의 숫자를 같은 물리량처럼 섞거나, field 수를 biological replicate 수로 오인하게 된다."), Image(str(figs["corpus"]), width=178 * mm, height=73 * mm), P("그림 2. 현재 corpus의 규모와 대표 modality 분포. log scale은 calcium/CyTOF와 작은 force dataset을 한 그림에 함께 표시하기 위한 것이다. HPA embedding은 물리 관측과 분리 집계한다.", st, "caption"), P("데이터 계층", st, "h2"), table([["계층", "입력", "표현", "사용 가능 범위"], ["시각 입력", "IF, BF, PC, DIC, multichannel live-cell", "pixel tensor, segmentation-free patch, morphology feature", "localization, spatial prediction, morphology state"], ["시간 입력", "calcium/FRET, migration tracks", "normalized waveform, event feature, trajectory/operator", "mechanism direction, temporal state, persistence"], ["힘·흐름", "TFM, MSM, PIV, tether/AFM", "vector/scalar field, exact unit, geometry", "load transfer, contractility direction, apparent tension product"], ["분자 입력", "RNA-seq, qPCR, CyTOF, IF/WB 후보", "gene/pathway vector, marker panel", "state direction, pathway context, cross-assay consistency"]], [25*mm, 42*mm, 52*mm, 59*mm]), PageBreak()]

    # Neural architecture
    story += [P("3. 신경망은 어떻게 구성되어 있는가", st, "h1"), P("현재 신경망은 단순한 한 줄 MLP가 아니다. 9개 trainable encoder head와 rule/theory/operator head가 서로 다른 입력 구조를 처리하며, 공통 결론은 evidence gate에서만 합쳐진다. architecture report의 명시적 값은 <b>shared_latent=false</b>, <b>fusion=task_scoped_late_fusion_only</b>이다."), table([["입력", "Encoder/representation", "Head", "검증 단위"], ["Raw IF / multichannel image", "CNN, channel-aware image embedding", "multilabel localization / cell-cycle / stress", "provider, plate, image, cell"], ["BF/PC/DIC", "conditional residual U-Net", "organelle fluorescence spatial map", "study, experimenter group, acquisition"], ["Calcium curve", "temporal convolution + pooled features", "cell state / mechanism contrast", "cell, condition, independent lab"], ["RNA/CyTOF", "pathway panel MLP / ridge / contrast", "perturbation or state direction", "compound, donor, dataset/lab"], ["Morphology", "curated feature encoder / ridge", "senescence/death endpoints", "image, well, biological repeat"], ["TFM/PIV/FRET", "physics-aware observation operator", "force/flow direction and identifiability", "field, experiment, lab"]], [31*mm, 48*mm, 47*mm, 52*mm]), P("왜 universal latent를 아직 만들지 않았는가", st, "h2"), bullet("대부분 modality가 동일 세포·동일 시간점에서 정렬되어 있지 않다.", st), bullet("연구실과 protocol 정보가 modality와 강하게 confounded되어 있다.", st), bullet("무리한 shared embedding은 biology가 아니라 provider signature를 학습할 가능성이 높다.", st), bullet("task별 split unit, label meaning, calibration target이 달라 하나의 accuracy로 합칠 수 없다.", st), P("장기적으로는 modality encoder → shared biological latent → theory experts → decoder 구조가 목표지만, paired cross-modal samples, batch/lab adversarial validation, missing-modality training이 확보될 때만 승격한다.", st, "callout"), Image(str(figs["status"]), width=178 * mm, height=68 * mm), P("그림 3. 현재 30개 head의 권한 분포. 통과와 차단을 함께 제시해야 시스템의 실제 상태가 보인다.", st, "caption"), PageBreak()]

    # Real images
    story += [P("4. 실제로 어떤 이미지를 보고 어떻게 판단하는가", st, "h1"), P("아래 패널은 개념도가 아니라 저장된 실제 입력과 실제 frozen output이다. 왼쪽 HPA 예시는 개발용 localization representation을, 오른쪽 Light My Cells 예시는 외부 single-use spatial confirmation의 실패를 보여준다."), Image(str(IMAGE_COPY), width=178 * mm, height=147 * mm), P("그림 4. 실제 영상 evidence panel. HPA: DAPI와 protein channel의 공간 패턴에서 nuclear interior 확률 0.9995를 냈고 image-quality OOD는 통과했지만 같은-provider development 예시이므로 외부 권한은 없다. Light My Cells: transmitted-light에서 mitochondria fluorescence를 예측하고 author target과 error map을 비교한다. 보이는 구조가 일부 맞아도 aggregate Pearson 0.329와 embedding OOD AUROC 0.385가 frozen endpoint를 실패하여 외부 spatial transfer는 차단된다.", st, "caption"), PageBreak()]

    story += [P("4.1 IF 영상 판단 단계", st, "h2"), table([["단계", "모델이 실제로 보는 것", "계산", "허용되는 판단"], ["Channel 확인", "DAPI/nucleus, protein target, 보조 channel", "channel order·dynamic range·saturation·missingness", "이미지 품질/입력 적합성"], ["공간 표현", "핵 내부, 막 주변, filament, puncta, diffuse texture", "CNN embedding 또는 morphology/statistical feature", "localization class probability"], ["저자 label 비교", "author localization / condition", "multilabel F1, AP, calibration", "동일 task의 predictive evidence"], ["외부 이동", "새 provider/lab/optics", "OOD AUROC, conformal coverage, frozen endpoint", "통과 시에만 외부 transfer"], ["거부", "corrupt source, channel mismatch, semantic shift", "integrity and authority gate", "prediction 0 또는 blocked claim"]], [22*mm, 49*mm, 52*mm, 55*mm]), P("4.2 Label-free spatial prediction 단계", st, "h2"), P("BF/PC/DIC 입력을 target fluorescence 이미지로 decoder가 복원한다. 판단은 예쁜 결과 한 장이 아니라 (i) study-macro Pearson, (ii) SSIM, (iii) fit-mean baseline 대비 MAE 개선, (iv) pixel interval coverage, (v) embedding OOD를 동시에 본다. 현재 Light My Cells는 5개 endpoint 중 Pearson과 OOD를 실패했다."), P("4.3 이미지로 직접 말할 수 없는 것", st, "h2"), P("한 장의 IF/BF 이미지만으로 cortical tension, membrane-cortex adhesion, contractility, cell type을 고유하게 역추정할 수 없다. 동일 형태를 만드는 여러 기전이 존재하기 때문이다. 따라서 영상 결과는 force/flow/RNA/time-series expert와 결합하거나 ‘현재 관측으로 구별 불가능’한 ambiguity set을 반환해야 한다.", st, "warning"), PageBreak()]

    # Modality matrix
    story += [P("5. Modality별로 어디까지 판단할 수 있는가", st, "h1"), table([["Modality", "보는 데이터", "현재 가능한 출력", "추가되어야 강해지는 것"], ["IF / HPA / raw multichannel", "핵·표적 단백질·세포골격 위치", "localization representation, translocation, morphology state", "새 독립 provider localization holdout"], ["BF / PC / DIC", "세포 경계·texture·phase contrast", "organelle spatial prediction 후보", "새 external study, semantic OOD"], ["Calcium / FRET", "개별 세포 waveform과 조건", "Piezo1 기전 방향, temporal state 후보", "replicate ID와 새 lab per-cell endpoint"], ["Migration tracking / PIV", "track, velocity vector field, coherence", "migration regime, persistence, flow direction", "aligned biological replicates와 force pairing"], ["TFM / MSM / tether", "traction/stress/tether force와 geometry", "load-transfer/contractility direction, identifiable product", "geometry-matched Aleph operator"], ["RNA / qPCR / CyTOF", "gene/pathway/marker response", "state direction, treatment context", "donor/lab holdout, protein orthogonal validation"], ["WB / proteomics", "band/abundance/phosphorylation", "현재 registry target, 직접 head는 미완성", "정량 원자료·replicate·antibody metadata"], ["EM / metabolomics", "ultrastructure / metabolic profile", "현재 coverage gap", "공개 raw data + harmonized operator"]], [25*mm, 43*mm, 55*mm, 55*mm], font_size=7.2), P("시각 데이터 해석의 두 레이어", st, "h2"), bullet("<b>Observation layer:</b> 픽셀을 segmentation, localization, morphology, vector field, intensity/texture로 바꾸며 단위와 uncertainty를 보존한다.", st), bullet("<b>Biophysical inference layer:</b> observation을 state/mechanism posterior와 competing theory expert의 지지도로 바꾼다.", st), P("TFM처럼 이미 힘 단위인 데이터는 image encoder를 거치지 않고 observation operator로 직접 들어간다. 반대로 IF/WB/PCR은 먼저 시각·밴드·Ct 값을 정량 observation으로 compile한 뒤 기전 head로 보낸다. 두 경로를 억지로 동일 CNN에 넣지 않는다.", st, "callout"), PageBreak()]

    # Validation
    story += [P("6. 외부 검증, 불확실성, OOD와 거부", st, "h1"), Image(str(figs["validation"]), width=178 * mm, height=77 * mm), P("그림 5. Single-use 외부 검증 절차. 결과를 본 뒤 threshold나 replacement를 고치지 않는 것이 핵심이다.", st, "caption"), P("평가 원칙", st, "h2"), bullet("cell/field 수를 biological replicate 수로 바꾸지 않는다.", st), bullet("dataset, lab, donor, compound, well 등 task에 맞는 최상위 split unit을 사용한다.", st), bullet("selection, calibration, final confirmation을 분리하고 final set은 한 번만 사용한다.", st), bullet("macro-F1/AUROC와 calibration/conformal coverage/OOD를 별도 endpoint로 둔다.", st), bullet("source integrity가 깨지면 예측 성능 이전 단계에서 prediction 0으로 거부한다.", st), bullet("동일 provider의 test는 training asset일 수 있지만 independent evidence로 승격하지 않는다.", st), P("Conformal coverage가 90%에 가깝다고 자동 통과하지 않는다. prediction set이 너무 크거나 특정 class coverage가 무너지면 uncertainty가 유용하지 않다. PBMC에서 이것이 실제로 관찰되었다.", st, "warning"), PageBreak()]

    # Results
    story += [P("7. 현재까지의 외부 결과", st, "h1"), P("7.1 통과한 task-scoped evidence", st, "h2"), table([["Head", "외부 결과", "허용되는 주장", "금지되는 확대"], ["Piezo1 calcium mechanism", "equal-stratum mean diff -0.1670; 95% CI [-0.2609, -0.0874]; 2/2 direction", "활성 strata에서 기전 방향 재현", "범용 calcium cell-state classifier"], ["Paired migration state", "AUROC 0.9849; macro-F1 0.7221; coverage 1.0", "retrospective external-lab task 성능", "pristine prospective generalization"], ["Fibrotic-state mechanism", "TGFβ-control +0.8720; celastrol reversal +1.0629; direction p=1.0", "조건 수준 기전 방향", "absolute fibrotic state classifier"], ["Senescence transcriptomic direction", "두 contrast 모두 3/3 positive; combined 6/6; sign-test p=0.015625", "외부 transcriptomic direction", "absolute state; OOD에서 거부"]], [35*mm, 54*mm, 45*mm, 44*mm], font_size=7.2), P("7.2 차단된 classifier/transfer", st, "h2"), table([["Head", "좋아 보이는 부분", "실패 endpoint", "판정"], ["Calcium cell state", "AUROC 0.7617; coverage 0.9228", "macro-F1 0.5499", "blocked"], ["PBMC cell type/IFN", "Wilk macro-F1 0.7458; OOD AUROC 0.7250", "marginal coverage 0.6937; monocyte 0.2654 등 class coverage", "blocked uncertainty"], ["Raw IF localization", "macro AUROC 0.7157", "macro-F1 0.3057; AP 0.3169; ECE 0.2370; 5/5 frozen endpoints fail", "blocked domain shift"], ["Label-free organelle", "SSIM 0.2746; MAE improvement 0.1192; coverage 0.9231", "Pearson 0.3292; OOD AUROC 0.3854", "blocked spatial transfer"]], [35*mm, 51*mm, 64*mm, 28*mm], font_size=7.2), P("이 결과는 ‘학습이 끝났다/안 끝났다’의 이분법보다 더 구체적이다. 일부 task는 외부 기전 방향을 지지하지만, 범용 classifier와 cross-provider vision transfer는 아직 논문 수준의 일반화 증거가 부족하다.", st, "callout"), PageBreak()]

    # Integrity failure
    story += [P("8. 실패를 결과로 보존하는 시스템", st, "h1"), P("IDR0168 single-use confirmation은 45개 field와 9개 gene × 5개 cell line으로 사전 동결되었다. 45/45 object, 총 34,187,550,441 byte를 취득했지만 A375/RBM23 공식 929,941,454-byte tgz에서 gzip CRC mismatch가 발생했고 native ./0/0/.zarray가 없었다. 서로 다른 위치의 1 MiB HTTP range 5개가 local bytes와 일치했으므로 local download 손상으로 임의 판단하지 않았다."), table([["Integrity endpoint", "결과"], ["Frozen protocol", idr["protocol_sha256"][:20] + "…"], ["Failed object SHA-256", idr["failed_archive"]["sha256"][:20] + "…"], ["Gzip integrity", "FAIL: " + idr["source_integrity"]["gzip_error"]], ["Native Zarr metadata", "absent: ./0/0/.zarray"], ["Remote range checks", "5/5 local object와 일치"], ["Replacement field", "사용하지 않음"], ["Tensor / predictions", f"tensor={idr['tensor_written']}; predictions={idr['model_predictions_made']}"], ["Final status", idr["status"]]], [55*mm, 123*mm]), P("이 확인셋은 <b>failed_source_integrity_and_permanently_consumed</b>이며 pristine set으로 재사용할 수 없다. 실패 전에 memory에서 decode된 7개 field도 모델 evidence로 사용하지 않았다. 좋은 결과가 나올 때까지 field를 교체하는 것을 막은 실제 사례다.", st, "warning"), P("다른 실패", st, "h2"), bullet("IDR0006: phenotype schema 불일치로 pixel 전에 거부; prediction 0.", st), bullet("IDR0072: 2,593 image 진단에서 모든 frozen transfer endpoint 실패; consumed diagnostic.", st), bullet("Light My Cells: 일부 pixel metric 통과에도 cross-study correlation/OOD 실패; confirmation consumed.", st), bullet("Whole-file test에는 unrelated numerical replay drift 3건이 남아 있으며 이를 숨기거나 baseline을 재작성하지 않았다.", st), PageBreak()]

    # Axes
    story += [P("9. Cortical tension을 넘는 cell-state library", st, "h1"), Image(str(figs["axes"]), width=160 * mm, height=101 * mm), P("그림 6. 현재 cell-state axis coverage. 녹색도 universal validation을 뜻하지 않으며, amber는 training-only 또는 context-dependent 상태다.", st, "caption"), P("현재 library는 cortical tension 하나로 좁혀지지 않았다. calcium mechanosensing, migration/ROCK, fibrosis, senescence, EMT, cell cycle, mitochondrial stress, apoptosis, immune morphology, supracellular contractility와 focal-adhesion load transfer를 포함한다. 그러나 hypoxia, differentiation/stemness, DNA damage/genotoxic stress, infection response는 명시적으로 비어 있다."), P("다음 infection source로는 IDR0128 influenza/A549 공개 영상을 우선 후보로 두되, pixel을 열기 전에 author outcome, biological replicate, split, endpoint, licence, no-replacement 규칙을 preregister해야 한다. WB/proteomics는 orthogonal protein confirmation을 위해 별도 축으로 필요하다.", st, "callout"), PageBreak()]

    # Aleph connection
    story += [P("10. Aleph와의 관계 및 parameter sweep", st, "h1"), P("현재 Outer Library의 Aleph authority는 <b>none</b>이고 sweep gate status는 <b>refused</b>다. 이 시스템이 곧바로 177개 parameter의 숫자 범위를 내놓거나 physics module을 바꾸면 안 된다."), table([["단계", "필요한 증거", "현재"], ["외부 state/mechanism inference", "independent task endpoint + uncertainty/OOD", "4개 제한 task 통과"], ["Observation operator", "실험량과 simulation output의 동일 단위·geometry·protocol", "부분 구현, 미검증"], ["Identifiability", "같은 관측을 만드는 competing parameter 조합 분리", "미충족"], ["Aleph forward sensitivity", "각 축 변화가 관측에 미치는 방향/크기", "sweep 계산비로 미충족"], ["Pre-sweep reduction", "위 세 단계 통과 + physics review", "금지"], ["Numeric sweep", "엔진 cost gate + bitwise-safe physics", "현재 불가"]], [38*mm, 94*mm, 46*mm]), P("Cortical tension 예시", st, "h2"), P("TFM/PIV/contour/AFM/tether 데이터를 통해 apparent tension 또는 contractility direction을 좁힐 수는 있다. 하지만 막 장력, membrane-cortex adhesion, actomyosin contractility가 동일 관측에 섞여 있다. 따라서 현재 출력은 ‘가능한 기전 집합 + 상대 지지도 + 구별할 다음 assay’여야 한다. 향후 Aleph가 빨라지면 이 posterior가 full Cartesian sweep 대신 active learning/Bayesian optimization의 proposal distribution이 될 수 있다."), P("파동함수/양자역학 비유는 여러 가능한 상태의 분포를 표현하는 데는 유용하지만, 각 노드 위치를 물리적 wavefunction으로 바꾼다고 고전적 세포역학 sweep 비용이 자동 감소하지 않는다. 필요한 것은 surrogate, multi-fidelity simulation, amortized inference, active subspace, Bayesian experimental design이며, 불확실성은 probability distribution으로 명시하면 충분하다.", st, "warning"), PageBreak()]

    # Paper readiness
    story += [P("11. 현재 결과로 논문이 가능한가", st, "h1"), P("가능한 논문의 가장 강한 형태는 <b>multimodal evidence compiler + task-scoped neural heads + single-use external validation + refusal/authority framework</b>의 methods/resource paper다. 실제 영상 panel, 4개 외부 evidence task, 4개 차단 task, source-integrity refusal은 단순 계획이 아니라 실행 결과다."), P("현재 주장 가능한 것", st, "h2"), bullet("서로 다른 공개 실험을 provenance/unit/split-aware observation으로 통합하는 구현.", st), bullet("모든 modality를 하나의 latent로 억지 통합하지 않는 leakage-aware architecture.", st), bullet("좋은 ranking과 실패한 calibration을 분리하여 보고하는 uncertainty/OOD evaluation.", st), bullet("사전 동결된 single-use confirmation이 손상되면 prediction 0으로 멈추는 재현 가능한 refusal.", st), bullet("제한된 4개 task에서 외부 mechanism/state direction evidence.", st), P("아직 주장할 수 없는 것", st, "h2"), bullet("모든 세포종·세포상태를 일반화하는 foundation model.", st), bullet("새 연구실 raw image에서 안정적인 subcellular localization/label-free spatial transfer.", st), bullet("Aleph parameter recovery, sweep reduction 배수, physics prediction improvement.", st), bullet("임상 진단 또는 intervention efficacy prediction.", st), P("따라서 지금 제출한다면 system/resource/method 논문 초안으로는 성립 가능하지만, broad biological model 논문을 위해서는 새 pristine provider, missing axes, prospective assay와 Aleph matched-operator 실험이 추가되어야 한다.", st, "callout"), PageBreak()]

    # Roadmap
    story += [P("12. 다음 개발 순서", st, "h1"), table([["우선순위", "작업", "완료 기준", "얻는 권한"], ["P0", "새 independent-provider raw IF confirmation", "preregistered common localization task; macro-F1/AP/ECE/conformal/OOD 모두 통과", "cross-provider IF evidence"], ["P0", "PBMC domain-conditional uncertainty", "새 calibration + untouched confirmation lab; marginal/class coverage", "cell-type context uncertainty"], ["P0", "Label-free external replacement", "새 provider; Pearson·SSIM·MAE·coverage·OOD 5/5", "spatial prediction evidence"], ["P1", "IDR0128 infection axis", "replicate-aware frozen outcome and clean source", "infection-state task"], ["P1", "Hypoxia / stemness / DNA damage", "각각 ≥2 provider/lab 또는 명시적 training-only", "coverage 확대"], ["P1", "Aligned TFM/PIV/IF", "같은 perturbation, exact unit, biological replicate, independent lab", "mechanical state identifiability"], ["P2", "WB/proteomics orthogonal layer", "raw quantitative replicate data and antibody/protein provenance", "cross-assay mechanism confirmation"], ["P2", "Aleph adapter", "geometry-matched operator + forward sensitivity + executable sweep", "conditional pre-sweep proposal"]], [17*mm, 48*mm, 78*mm, 35*mm], font_size=7.1), P("발표에서 보여줄 최소 live sequence", st, "h2"), bullet("실제 HPA/Light My Cells image panel에서 pixel → prediction → error/OOD를 설명한다.", st), bullet("30개 head status에서 PASS 4와 BLOCK 4를 나란히 보여준다.", st), bullet("IDR0168 frozen protocol → acquisition → CRC failure → prediction 0 refusal을 재생한다.", st), bullet("마지막에 sweep_gate_report.json의 refused와 Aleph authority none을 보여준다.", st), P("이 순서가 ‘무엇을 만들었는가’뿐 아니라 ‘어떤 상황에서 스스로 멈추는가’를 시각적으로 입증한다.", st, "callout"), PageBreak()]

    # Reproducibility
    story += [P("재현성, artifact와 테스트 상태", st, "h1"), table([["Artifact", "역할"], ["results/adversarial_coverage_audit.json", "dataset/lab/sample/observation counts와 missing axes"], ["results/multitask_outer_model_report.json", "30 head architecture, readiness, per-head evidence"], ["results/sweep_gate_report.json", "Aleph numeric authority refusal"], ["results/idr0168_confirmation_decode_refusal.json", "single-use source-integrity refusal"], ["REPORT.md", "상세 구현·dataset·metric narrative"], ["MONDAY_PRESENTATION_PLAN.md", "2026-08-10 발표 claim discipline"], ["output/pdf/outer_library_real_image_evidence.png", "실제 HPA/Light My Cells 판독 패널"]], [70*mm, 108*mm]), P("현재 test_outer_library.py는 129개 test를 수집한다. IDR0168 refusal, summary hash, multitask/sweep authority의 targeted test 4개는 통과했다. 전체 파일 실행에서는 HPA IF conformal radius 약 3.5e-13, BBBC054 optimizer result, calcium temporal CNN optimizer result의 세 numerical replay drift가 남았다. 의도적 report update로 생긴 summary hash mismatch는 재생성되어 통과했다. 따라서 ‘whole suite green’이라고 주장하지 않는다.", st, "warning"), P("보고서 생성 규칙", st, "h2"), P("본 PDF의 corpus/head/status 수치는 위 JSON에서 읽어 생성했다. 실제 image panel은 content-hashed tensor와 frozen output에서 deterministic하게 생성된 기존 panel을 그대로 삽입했다. 본 생성기는 Aleph physics를 실행하거나 수정하지 않는다."), PageBreak()]

    # References
    story += [P("선정 데이터·방법 참고문헌", st, "h1"), P("아래는 현재 Outer Library에 직접 연결된 대표 공개 source다. 전체 38개 source의 content hash, licence, acquisition path는 acquisition_receipts.jsonl 및 dataset manifest에 있다."), table([["분야", "공개 source / primary record"], ["Mechanotransduction calcium", "Zenodo 18495219; Zenodo 19890985; Zenodo 8215150"], ["Migration", "MULTIMOT 2D, DOI 10.17044/scilifelab.21407402; eLife 71032, DOI 10.7554/eLife.71032"], ["Mechano-osmotic / tether", "eLife 72381, DOI 10.7554/eLife.72381"], ["Fibrotic state", "NCBI GEO GSE226374; 2026 tendon mechanoculture source"], ["ROCK / fibronectin", "Dryad 9jh6m, DOI 10.5061/dryad.9jh6m; eLife 11384"], ["IF benchmarks", "Broad Bioimage Benchmark Collection BBBC013, BBBC014, BBBC048, BBBC053, BBBC054"], ["Subcellular localization", "Human Protein Atlas v25.1 subcellular image embeddings; IDR0072 diagnostic; IDR0168/PUPS confirmation attempt"], ["Label-free imaging", "Light My Cells public training and single-use confirmation partitions"], ["Perturbation state", "sci-Plex3; Reactome human pathway registry"], ["Senescence", "GEO GSE250041; GSE301164 confirmation; SenSCOUT Dryad"], ["Programmed cell death", "BioImage Archive S-BIAD2515; Figshare 28202864 v2"]], [48*mm, 130*mm], font_size=7.4), P("방법론적으로 본 시스템은 split conformal prediction, dataset/lab-held-out validation, OOD scoring, multitask/late-fusion representation, physics-informed observation operator의 원칙을 사용한다. 각 head의 exact algorithm, threshold와 split contract는 개별 preregistration JSON 및 result report에 기록되어 있다."), P("결론", st, "h2"), P("Outer Library는 아직 ‘세포 상태를 모두 아는 신경망’이 아니다. 대신 어떤 공개 실험이 실제 evidence가 되는지, 어떤 영상에서 무엇을 읽고 무엇을 읽지 못하는지, domain shift와 source failure에서 언제 멈추는지를 실행 가능한 형태로 만들었다. 다음 연구의 중심은 head 수를 무작정 늘리는 것이 아니라, 비어 있는 state axis와 실패한 외부 transfer를 새 pristine lab/provider에서 다시 검증하고, force/flow/image를 Aleph와 동일한 observation operator로 연결하는 것이다.", st, "callout")]

    doc = PaperDocTemplate(str(PDF_PATH), title="Aleph Outer Library Visual Paper 2026-08-07", author="Project Aleph")
    doc.build(story)
    print(PDF_PATH)


if __name__ == "__main__":
    build()
