"""Build the evidence-bounded External Network architecture report.

The report is descriptive documentation. It reads committed Outer Library
artifacts and never runs Aleph physics, scores a validation gate, or promotes a
literature quantity to parameter authority.

Example:
    /path/to/codex-runtime/python3 build_architecture_report.py \
      --output aleph/outputs/outer/mechanics_protocol/\
external_neural_network_architecture_report_2026-08-10_v3.pdf
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    Image,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[5]
OUTER = ROOT / "aleph/outer/experiment_factory/outer_library"
MECHANICS_OUTPUT = ROOT / "aleph/outputs/outer/mechanics_protocol"

FONT_PATH = Path(
    "/System/Library/AssetsV2/com_apple_MobileAsset_Font8/"
    "7a0b5c0f3c1d41c4c52a33343496c9c65ad52c50.asset/"
    "AssetData/NanumGothic.ttc"
)

NAVY = HexColor("#0B1220")
INK = HexColor("#152238")
SLATE = HexColor("#475569")
MUTED = HexColor("#64748B")
LIGHT = HexColor("#F4F7FB")
LINE = HexColor("#D8E2EF")
CYAN = HexColor("#06B6D4")
TEAL = HexColor("#0F9F8F")
GREEN = HexColor("#16A34A")
AMBER = HexColor("#D97706")
RED = HexColor("#DC2626")
BLUE = HexColor("#2563EB")
PURPLE = HexColor("#7C3AED")
WHITE = colors.white


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _register_fonts() -> None:
    if not FONT_PATH.exists():
        raise FileNotFoundError(f"Korean font is unavailable: {FONT_PATH}")
    pdfmetrics.registerFont(TTFont("Nanum", str(FONT_PATH), subfontIndex=0))
    pdfmetrics.registerFont(TTFont("NanumBold", str(FONT_PATH), subfontIndex=1))


def _styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "ReportTitle",
            parent=base["Title"],
            fontName="NanumBold",
            fontSize=24,
            leading=31,
            textColor=NAVY,
            alignment=TA_LEFT,
            spaceAfter=8,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=base["Heading1"],
            fontName="NanumBold",
            fontSize=17,
            leading=22,
            textColor=NAVY,
            spaceBefore=0,
            spaceAfter=8,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=base["Heading2"],
            fontName="NanumBold",
            fontSize=11.5,
            leading=15,
            textColor=INK,
            spaceBefore=5,
            spaceAfter=4,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="Nanum",
            fontSize=8.4,
            leading=12.7,
            textColor=INK,
            spaceAfter=4,
            wordWrap="CJK",
        ),
        "small": ParagraphStyle(
            "Small",
            parent=base["BodyText"],
            fontName="Nanum",
            fontSize=7.1,
            leading=10.3,
            textColor=SLATE,
            wordWrap="CJK",
        ),
        "tiny": ParagraphStyle(
            "Tiny",
            parent=base["BodyText"],
            fontName="Nanum",
            fontSize=6.1,
            leading=8.2,
            textColor=MUTED,
            wordWrap="CJK",
        ),
        "metric": ParagraphStyle(
            "Metric",
            parent=base["BodyText"],
            fontName="NanumBold",
            fontSize=15,
            leading=18,
            textColor=NAVY,
            alignment=TA_CENTER,
        ),
        "metric_label": ParagraphStyle(
            "MetricLabel",
            parent=base["BodyText"],
            fontName="Nanum",
            fontSize=6.8,
            leading=9.5,
            textColor=SLATE,
            alignment=TA_CENTER,
            wordWrap="CJK",
        ),
        "code": ParagraphStyle(
            "Code",
            parent=base["Code"],
            fontName="Nanum",
            fontSize=6.6,
            leading=9.4,
            textColor=INK,
            backColor=HexColor("#EEF3F8"),
            borderPadding=6,
            borderRadius=3,
            wordWrap="CJK",
        ),
    }


class CoverFlowable(Flowable):
    """Full-page report cover."""

    def __init__(self, width: float, height: float) -> None:
        super().__init__()
        self.width = width
        self.height = height

    def draw(self) -> None:
        c = self.canv
        w, h = self.width, self.height
        c.setFillColor(NAVY)
        c.rect(0, 0, w, h, stroke=0, fill=1)
        c.setFillColor(HexColor("#102A43"))
        c.circle(w * 0.86, h * 0.80, 145, stroke=0, fill=1)
        c.setFillColor(HexColor("#0C4A6E"))
        c.circle(w * 0.93, h * 0.89, 74, stroke=0, fill=1)
        c.setStrokeColor(CYAN)
        c.setLineWidth(1.2)
        for radius in (38, 62, 92):
            c.circle(w * 0.79, h * 0.76, radius, stroke=1, fill=0)
        for angle in (0, 60, 120, 180, 240, 300):
            rad = math.radians(angle)
            x1 = w * 0.79 + 38 * math.cos(rad)
            y1 = h * 0.76 + 38 * math.sin(rad)
            x2 = w * 0.79 + 92 * math.cos(rad)
            y2 = h * 0.76 + 92 * math.sin(rad)
            c.line(x1, y1, x2, y2)
        c.setFillColor(CYAN)
        c.roundRect(38, h - 91, 124, 23, 11, stroke=0, fill=1)
        c.setFillColor(NAVY)
        c.setFont("NanumBold", 9)
        c.drawCentredString(100, h - 83, "PROJECT ALEPH / OUTER")

        c.setFillColor(WHITE)
        c.setFont("NanumBold", 26)
        c.drawString(38, h - 185, "외부 신경망 및 Cell-State 추론")
        c.drawString(38, h - 222, "구조 보고서 v3")
        c.setFont("Nanum", 11)
        c.setFillColor(HexColor("#BFD8EA"))
        c.drawString(40, h - 253, "실제 구현 · 증거 경계 · mechanics protocol · C-1 판별 실험")
        c.setFillColor(HexColor("#89D8E7"))
        c.rect(40, h - 273, 104, 3, stroke=0, fill=1)

        cards = [
            ("30 + 1", "기존 multitask heads + mechanics optional head", CYAN),
            ("9 / 19 / 196", "primary papers / cohorts / raw wires", TEAL),
            ("25× vs 1×", "a = 1→5 μm, poroelastic vs a-independent", AMBER),
        ]
        x = 39
        for value, label, color in cards:
            c.setFillColor(HexColor("#132D46"))
            c.roundRect(x, h - 400, 164, 93, 9, stroke=0, fill=1)
            c.setFillColor(color)
            c.rect(x, h - 400, 5, 93, stroke=0, fill=1)
            c.setFillColor(WHITE)
            c.setFont("NanumBold", 16)
            c.drawString(x + 15, h - 337, value)
            c.setFillColor(HexColor("#C8D7E5"))
            c.setFont("Nanum", 7.1)
            lines = label.split(" / ") if len(label) > 31 else [label]
            for index, line in enumerate(lines[:2]):
                c.drawString(x + 15, h - 358 - index * 11, line)
            x += 174

        c.setFillColor(WHITE)
        c.setFont("NanumBold", 11)
        c.drawString(40, 176, "현재 결론")
        c.setFillColor(HexColor("#CBD5E1"))
        c.setFont("Nanum", 8.5)
        lines = [
            "공통 cell-state 본체는 목표 아키텍처이며 아직 완성되지 않았다.",
            "protocol-aware mechanics representation과 exact-context cortical-tension head는 구현됐다.",
            "operator reconstruction은 성공했지만 cell-state transfer는 실패했으며, 런타임은 이를 거부한다.",
        ]
        y = 154
        for line in lines:
            c.setFillColor(CYAN)
            c.circle(44, y + 2, 2.5, stroke=0, fill=1)
            c.setFillColor(HexColor("#DCE7F1"))
            c.drawString(53, y, line)
            y -= 22
        c.setFillColor(HexColor("#8AA4B8"))
        c.setFont("Nanum", 7)
        c.drawString(40, 42, "기준일 2026-08-10  |  문서 성격: 구현·실증·로드맵 통합 보고서  |  Aleph parameter authority: none")


class ArchitectureFlowable(Flowable):
    """Architecture overview with implemented and target layers."""

    def __init__(self, width: float, height: float = 335) -> None:
        super().__init__()
        self.width = width
        self.height = height

    def _box(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        lines: list[str],
        color: colors.Color,
        dashed: bool = False,
    ) -> None:
        c = self.canv
        c.setStrokeColor(color)
        c.setLineWidth(1.1)
        c.setDash(4, 3) if dashed else c.setDash()
        c.setFillColor(colors.Color(color.red, color.green, color.blue, alpha=0.07))
        c.roundRect(x, y, w, h, 6, stroke=1, fill=1)
        c.setDash()
        c.setFillColor(color)
        c.setFont("NanumBold", 7.7)
        c.drawString(x + 8, y + h - 15, title)
        c.setFillColor(SLATE)
        c.setFont("Nanum", 5.8)
        for index, line in enumerate(lines):
            c.drawString(x + 8, y + h - 29 - index * 9, line)

    def draw(self) -> None:
        c = self.canv
        w = self.width
        c.setFillColor(LIGHT)
        c.roundRect(0, 0, w, self.height, 9, stroke=0, fill=1)
        c.setFillColor(MUTED)
        c.setFont("NanumBold", 6.4)
        c.drawString(14, self.height - 18, "INPUT LANES")
        c.drawString(150, self.height - 18, "EVIDENCE COMPILERS")
        c.drawString(312, self.height - 18, "SHARED STATE")
        c.drawString(442, self.height - 18, "OUTPUT / AUTHORITY")

        inputs = [
            ("Pixels", ["BF/PC/DIC/IF", "OME metadata"], BLUE),
            ("Dynamics", ["calcium / PIV", "curves and fields"], PURPLE),
            ("Mechanics", ["probe · geometry", "unit · time · state"], TEAL),
            ("Context", ["cell · lab", "substrate · treatment"], AMBER),
        ]
        ys = [246, 179, 112, 45]
        for (title, lines, color), y in zip(inputs, ys, strict=True):
            self._box(12, y, 105, 52, title, lines, color)
            c.setStrokeColor(LINE)
            c.line(117, y + 26, 143, y + 26)
            c.setFillColor(LINE)
            c.wedge(139, y + 22, 147, y + 30, 270, 180, stroke=0, fill=1)

        compilers = [
            ("Visual object compiler", ["segment → object tokens", "pixel provenance"], BLUE),
            ("Temporal/field encoders", ["direction · dynamics", "paired evidence only"], PURPLE),
            ("Mechanics ontology", ["method/scope route", "exact match or abstain"], TEAL),
            ("Evidence ledger", ["source audit · split", "OOD · refusal"], AMBER),
        ]
        for (title, lines, color), y in zip(compilers, ys, strict=True):
            self._box(145, y, 135, 52, title, lines, color, dashed=title == "Visual object compiler")

        c.setStrokeColor(LINE)
        for y in ys:
            c.line(280, y + 26, 305, 175)
        self._box(
            306,
            112,
            112,
            126,
            "Probabilistic cell-state latent",
            [
                "morphology · cycle",
                "polarity · migration",
                "cytoskeleton · adhesion",
                "contractility · membrane",
                "viscoelastic response",
                "uncertainty + ambiguity",
            ],
            MUTED,
            dashed=True,
        )
        c.setStrokeColor(LINE)
        c.line(418, 175, 438, 175)
        self._box(
            440,
            196,
            103,
            76,
            "Evidence posterior",
            ["supported axes", "intervals · conflicts", "next assay"],
            GREEN,
        )
        self._box(
            440,
            103,
            103,
            76,
            "Theory-gated experts",
            ["operator match", "sensitivity required", "no direct parameter"],
            AMBER,
            dashed=True,
        )
        self._box(
            440,
            16,
            103,
            70,
            "Refusal",
            ["OOD · unit mismatch", "unseen scope · blocked", "insufficient evidence"],
            RED,
        )
        c.setFillColor(GREEN)
        c.setFont("NanumBold", 6.1)
        c.drawString(16, 17, "SOLID = implemented path")
        c.setFillColor(MUTED)
        c.drawString(143, 17, "DASHED = target / operator-gated")


class CurrentNeuralRuntimeFlowable(Flowable):
    """The implemented tensor-to-head runtime, without implying a shared latent."""

    def __init__(self, width: float, height: float = 350) -> None:
        super().__init__()
        self.width = width
        self.height = height

    def _box(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        lines: list[str],
        color: colors.Color,
    ) -> None:
        c = self.canv
        c.setFillColor(colors.Color(color.red, color.green, color.blue, alpha=0.07))
        c.setStrokeColor(color)
        c.setLineWidth(1)
        c.roundRect(x, y, w, h, 5, stroke=1, fill=1)
        c.setFillColor(color)
        c.setFont("NanumBold", 6.8)
        c.drawString(x + 6, y + h - 13, title)
        c.setFillColor(SLATE)
        c.setFont("Nanum", 5.35)
        for index, line in enumerate(lines):
            c.drawString(x + 6, y + h - 25 - index * 8, line)

    def draw(self) -> None:
        c = self.canv
        c.setFillColor(LIGHT)
        c.roundRect(0, 0, self.width, self.height, 9, stroke=0, fill=1)
        headers = [(12, "INPUT TENSOR"), (110, "TRAINED ENCODER"), (257, "TASK HEAD BANK"), (373, "EVIDENCE GATE")]
        c.setFillColor(MUTED)
        c.setFont("NanumBold", 6.2)
        for x, label in headers:
            c.drawString(x, self.height - 17, label)

        lanes = [
            ("IF / BF image", ["x: B×C×H×W", "pixel + provider"], "CNN / residual U-Net", ["Conv2d or encoder-decoder", "spatial or pooled features"], BLUE),
            ("calcium curve", ["x: B×1×128", "normalized time series"], "temporal CNN", ["Conv1d 1→8, k=9", "global mean + max"], PURPLE),
            ("RNA / CyTOF", ["x: B×G", "fit-only transform"], "MLP / logistic", ["768→96 or marker panel", "calibrated logits"], TEAL),
            ("morphology", ["x: B×87", "signed log + z-score"], "linear multitask", ["β-gal logistic", "4 biomarker ridge heads"], GREEN),
            ("mechanics / field", ["x: B×F + context", "probe·unit·geometry"], "operator features", ["wire MLP / exact routing", "no universal tensor yet"], AMBER),
        ]
        ys = [273, 220, 167, 114, 61]
        for (input_title, input_lines, encoder_title, encoder_lines, color), y in zip(lanes, ys, strict=True):
            self._box(10, y, 86, 40, input_title, input_lines, color)
            self._box(108, y, 126, 40, encoder_title, encoder_lines, color)
            c.setStrokeColor(LINE)
            c.line(96, y + 20, 106, y + 20)
            c.line(234, y + 20, 251, 187)

        self._box(
            253,
            111,
            102,
            152,
            "30 task heads",
            [
                "classification / regression",
                "direction / spatial map",
                "force-flow constraints",
                "task-specific calibration",
                "task-specific split unit",
                "mechanics = optional +1",
                "shared latent = FALSE",
            ],
            CYAN,
        )
        c.setStrokeColor(RED)
        c.setDash(3, 2)
        c.line(361, 55, 361, 310)
        c.setDash()
        c.setFillColor(RED)
        c.setFont("NanumBold", 5.1)
        c.drawCentredString(361, 318, "NO SHARED LATENT")
        self._box(
            378,
            111,
            104,
            152,
            "per-head evidence gate",
            [
                "lab / provider holdout",
                "calibration + conformal",
                "semantic OOD",
                "operator compatibility",
                "PASS / training-only",
                "BLOCKED / REFUSE",
                "Aleph authority = none",
            ],
            RED,
        )
        c.setStrokeColor(LINE)
        c.line(355, 187, 376, 187)
        c.setFillColor(WHITE)
        c.setStrokeColor(NAVY)
        c.roundRect(108, 12, 374, 31, 5, stroke=1, fill=1)
        c.setFillColor(NAVY)
        c.setFont("NanumBold", 7)
        c.drawString(120, 29, "OUTPUT")
        c.setFillColor(SLATE)
        c.setFont("Nanum", 5.7)
        c.drawString(166, 29, "task posterior · calibrated interval/set · OOD score · refusal reason · next assay")
        c.drawString(166, 18, "numeric Aleph parameter range와 physics mutation은 항상 차단")
        c.setStrokeColor(LINE)
        c.line(430, 111, 430, 45)


class WireNetworkFlowable(Flowable):
    """Small shared-latent mechanics network schematic."""

    def __init__(self, width: float, height: float = 210) -> None:
        super().__init__()
        self.width = width
        self.height = height

    def draw(self) -> None:
        c = self.canv
        c.setFillColor(LIGHT)
        c.roundRect(0, 0, self.width, self.height, 8, stroke=0, fill=1)
        inputs = ["wire length", "wire diameter", "critical frequency", "elastic modulus"]
        x0 = 12
        for index, label in enumerate(inputs):
            y = 157 - index * 36
            c.setFillColor(WHITE)
            c.setStrokeColor(BLUE)
            c.roundRect(x0, y, 110, 25, 4, stroke=1, fill=1)
            c.setFillColor(INK)
            c.setFont("Nanum", 6.6)
            c.drawCentredString(x0 + 55, y + 9, label)
            c.setStrokeColor(LINE)
            c.line(x0 + 110, y + 12, 169, 104)
        c.setFillColor(HexColor("#DFF7F4"))
        c.setStrokeColor(TEAL)
        c.roundRect(169, 65, 115, 78, 7, stroke=1, fill=1)
        c.setFillColor(TEAL)
        c.setFont("NanumBold", 10)
        c.drawCentredString(226.5, 112, "shared latent")
        c.setFont("Nanum", 7)
        c.setFillColor(SLATE)
        c.drawCentredString(226.5, 92, "16 tanh units")
        c.setStrokeColor(LINE)
        c.line(284, 104, 341, 145)
        c.line(284, 104, 341, 61)
        c.setFillColor(HexColor("#E8F0FE"))
        c.setStrokeColor(BLUE)
        c.roundRect(341, 119, 174, 58, 7, stroke=1, fill=1)
        c.setFillColor(BLUE)
        c.setFont("NanumBold", 8)
        c.drawString(352, 156, "log η regression")
        c.setFillColor(INK)
        c.setFont("Nanum", 6.4)
        c.drawString(352, 140, "date-disjoint R² = 0.980")
        c.drawString(352, 128, "median absolute % error = 12.1%")
        c.setFillColor(HexColor("#FDECEC"))
        c.setStrokeColor(RED)
        c.roundRect(341, 32, 174, 58, 7, stroke=1, fill=1)
        c.setFillColor(RED)
        c.setFont("NanumBold", 8)
        c.drawString(352, 69, "cell-line classification")
        c.setFillColor(INK)
        c.setFont("Nanum", 6.4)
        c.drawString(352, 53, "balanced accuracy = 0.385")
        c.drawString(352, 41, "macro-F1 = 0.227 → runtime refusal")
        c.setFillColor(MUTED)
        c.setFont("Nanum", 6.3)
        c.drawString(14, 13, "Split unit: experimental date. One lab / one individual-level source; no claim of new-cell or new-lab transfer.")


class CorticalPipelineFlowable(Flowable):
    """Cortical-tension reference-to-engine readiness diagram."""

    def __init__(self, width: float, height: float = 190) -> None:
        super().__init__()
        self.width = width
        self.height = height

    def draw(self) -> None:
        c = self.canv
        steps = [
            ("External AFM", ["MCF-7 · suspended", "interphase · 37 C · 1 Hz"], GREEN),
            ("Reference head", ["median 270 pN/μm", "IQR 180-400 · n=27"], TEAL),
            ("Operator adapter", ["method-of-planes candidate", "geometry match required"], AMBER),
            ("Engine records", ["full native · accepted", "dt-independent · ≥3 seeds"], RED),
        ]
        gap = 15
        box_w = (self.width - gap * 3) / 4
        for index, (title, lines, color) in enumerate(steps):
            x = index * (box_w + gap)
            c.setFillColor(colors.Color(color.red, color.green, color.blue, alpha=0.08))
            c.setStrokeColor(color)
            c.setDash(4, 3) if index >= 2 else c.setDash()
            c.roundRect(x, 65, box_w, 88, 7, stroke=1, fill=1)
            c.setDash()
            c.setFillColor(color)
            c.setFont("NanumBold", 8)
            c.drawString(x + 9, 131, title)
            c.setFillColor(INK)
            c.setFont("Nanum", 6.2)
            for line_index, line in enumerate(lines):
                c.drawString(x + 9, 109 - line_index * 14, line)
            if index < 3:
                x_arrow = x + box_w
                c.setStrokeColor(LINE)
                c.line(x_arrow + 2, 109, x_arrow + gap - 3, 109)
                c.setFillColor(LINE)
                c.wedge(x_arrow + gap - 7, 105, x_arrow + gap + 1, 113, 270, 180, stroke=0, fill=1)
        c.setFillColor(HexColor("#FEF2F2"))
        c.setStrokeColor(RED)
        c.roundRect(0, 8, self.width, 35, 6, stroke=1, fill=1)
        c.setFillColor(RED)
        c.setFont("NanumBold", 7.3)
        c.drawString(12, 28, "BLOCKED")
        c.setFillColor(INK)
        c.setFont("Nanum", 6.2)
        c.drawString(
            72,
            28,
            "C-2 transaction은 PASS. 남은 blocker는 geometry-matched operator, physical inputs, active-state dt convergence, final tension seed ensemble.",
        )
        c.drawString(72, 16, "따라서 기존 engine γ magnitude와 ratio는 보고서에서도 redacted이며 draft gate를 채점하지 않는다.")


def _section_header(number: str, title: str, subtitle: str, styles: dict[str, ParagraphStyle]) -> list[Flowable]:
    number_style = ParagraphStyle(
        f"section_number_{number}",
        parent=styles["h1"],
        textColor=CYAN,
        alignment=TA_LEFT,
    )
    heading = Table(
        [[Paragraph(number, number_style), Paragraph(title, styles["h1"])]],
        colWidths=[16 * mm, 158 * mm],
        style=TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        ),
    )
    return [
        heading,
        Paragraph(subtitle, styles["small"]),
        Spacer(1, 7),
    ]


def _metric_cards(
    items: list[tuple[str, str, str]],
    styles: dict[str, ParagraphStyle],
) -> Table:
    cells = []
    for value, label, color in items:
        cells.append(
            Table(
                [[Paragraph(value, styles["metric"])], [Paragraph(label, styles["metric_label"])]],
                colWidths=[53 * mm],
                rowHeights=[13 * mm, 14 * mm],
                style=TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, -1), HexColor("#F7FAFC")),
                        ("BOX", (0, 0), (-1, -1), 0.7, HexColor(color)),
                        ("LINEABOVE", (0, 0), (-1, 0), 3, HexColor(color)),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                        ("LEFTPADDING", (0, 0), (-1, -1), 5),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                    ]
                ),
            )
        )
    return Table(
        [cells],
        colWidths=[56 * mm] * len(cells),
        style=TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 3)]),
    )


def _callout(
    title: str,
    text: str,
    styles: dict[str, ParagraphStyle],
    color: colors.Color = CYAN,
) -> Table:
    return Table(
        [[Paragraph(f"<b>{title}</b><br/>{text}", styles["body"])]],
        colWidths=[174 * mm],
        style=TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.Color(color.red, color.green, color.blue, alpha=0.07)),
                ("LINEBEFORE", (0, 0), (0, -1), 4, color),
                ("BOX", (0, 0), (-1, -1), 0.5, colors.Color(color.red, color.green, color.blue, alpha=0.25)),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        ),
    )


def _bullet(text: str, styles: dict[str, ParagraphStyle], color: str = "#06B6D4") -> Paragraph:
    return Paragraph(f'<font color="{color}">●</font> {text}', styles["body"])


def _table(
    rows: list[list[Any]],
    widths: list[float],
    styles: dict[str, ParagraphStyle],
    header: bool = True,
    font_size: float = 6.5,
) -> Table:
    formatted: list[list[Any]] = []
    for r_index, row in enumerate(rows):
        style = ParagraphStyle(
            f"table_{r_index}_{font_size}",
            parent=styles["small"],
            fontName="NanumBold" if header and r_index == 0 else "Nanum",
            fontSize=font_size,
            leading=font_size * 1.45,
            textColor=WHITE if header and r_index == 0 else INK,
            wordWrap="CJK",
        )
        formatted.append([cell if isinstance(cell, Flowable) else Paragraph(str(cell), style) for cell in row])
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.4, LINE),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]
    if header:
        commands.append(("BACKGROUND", (0, 0), (-1, 0), NAVY))
        if len(rows) > 1:
            commands.append(("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, HexColor("#F8FAFC")]))
    return Table(formatted, colWidths=widths, repeatRows=1 if header else 0, style=TableStyle(commands))


def _prepare_figures(
    scratch: Path,
    mechanics: dict[str, Any],
    vision: dict[str, Any],
) -> dict[str, Path]:
    scratch.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
        }
    )

    figures: dict[str, Path] = {}

    # Existing Allen diagnostic is intentionally cropped to three rows so it is legible on A4.
    allen_source = OUTER / "results/allencell_label_free_v2_v3_comparison.png"
    with PILImage.open(allen_source) as image:
        crop = image.crop((0, 0, image.width, int(image.height * 0.39)))
        allen_crop = scratch / "allen_v2_v3_three_rows.png"
        crop.save(allen_crop, optimize=True)
    figures["allen_crop"] = allen_crop

    # Cross-provider vision metrics.
    diagnostic = vision["post_refusal_non_authoritative_diagnostic"]
    v2 = diagnostic["v2"]
    v3 = diagnostic["v3"]
    v2_per = v2["metrics"]["per_target"]["Mitochondria"]
    v3_per = v3["metrics"]["per_target"]["Mitochondria"]
    vision_values = {
        "Pearson": [v2_per["macro_study_Pearson"], v3_per["macro_study_Pearson"]],
        "global SSIM": [v2_per["macro_study_global_SSIM"], v3_per["macro_study_global_SSIM"]],
        "MAE improvement": [
            v2["baseline"]["Mitochondria"]["relative_MAE_improvement"],
            v3["baseline"]["Mitochondria"]["relative_MAE_improvement"],
        ],
    }
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.35))
    for axis, (name, values) in zip(axes, vision_values.items(), strict=True):
        bars = axis.bar(["v2", "v3"], values, color=["#94A3B8", "#06B6D4"], width=0.62)
        axis.axhline(0, color="#CBD5E1", linewidth=0.8)
        axis.set_title(name, fontweight="bold")
        lo = min(-0.7, min(values) * 1.2) if name == "MAE improvement" else 0
        hi = max(0.55, max(values) * 1.35)
        axis.set_ylim(lo, hi)
        axis.spines[["top", "right"]].set_visible(False)
        for bar, value in zip(bars, values, strict=True):
            label = f"{100 * value:.1f}%" if name == "MAE improvement" else f"{value:.3f}"
            y = value + (0.035 if value >= 0 else -0.075)
            axis.text(bar.get_x() + bar.get_width() / 2, y, label, ha="center", va="center", fontsize=7, fontweight="bold")
    fig.suptitle("Allen hiPSC post-refusal diagnostic", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    path = scratch / "vision_metrics.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figures["vision_metrics"] = path

    # Cross-method viscosity scale. Values are deliberately not pooled.
    fig, axis = plt.subplots(figsize=(8.8, 3.0))
    rows = [
        ("37 C water (engine solvent input)", 1.0e-3, 1.0e-3, "#64748B"),
        ("MCF-7 optical Brownian probe*", 0.002, 0.016, "#D97706"),
        ("MCF-7 magnetic rotational spectroscopy", 65.9, 65.9, "#0F9F8F"),
        ("J774 surface magnetic bead", 2000.0, 2000.0, "#2563EB"),
    ]
    for index, (_label, low, high, color) in enumerate(rows):
        if low == high:
            axis.scatter([low], [index], s=72, color=color, zorder=3)
        else:
            axis.plot([low, high], [index, index], color=color, linewidth=7, solid_capstyle="round")
        axis.text(high * 1.13, index, f"{low:g}" if low == high else f"{low:g}-{high:g}", va="center", fontsize=7, color=color, fontweight="bold")
    axis.set_xscale("log")
    axis.set_xlim(6e-4, 1.2e4)
    axis.set_yticks(range(len(rows)), [row[0] for row in rows])
    axis.set_xlabel("reported / input viscosity scale (Pa s, log axis)")
    axis.grid(axis="x", which="both", color="#E2E8F0", linewidth=0.7)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    axis.set_title("One word, different measurement operators", fontweight="bold")
    fig.tight_layout()
    path = scratch / "viscosity_scale.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figures["viscosity_scale"] = path

    # Probe-size exponents from every wire and the paper's binned summaries.
    scaling = mechanics["probe_size_dependence"]["external_mrs_wire_length_scaling"]
    lines = ["MCF-10A", "MCF-7", "MDA-MB-231"]
    raw = np.array([scaling[line]["raw_wire_level_log_log_exponent"] for line in lines])
    intervals = np.array([scaling[line]["experiment_date_cluster_bootstrap_95_interval"] for line in lines])
    published = np.array([scaling[line]["published_binned_exponent"] for line in lines])
    published_se = np.array([scaling[line]["published_standard_error"] for line in lines])
    x = np.arange(len(lines))
    fig, axis = plt.subplots(figsize=(8.8, 3.2))
    axis.errorbar(
        x - 0.08,
        raw,
        yerr=np.vstack([raw - intervals[:, 0], intervals[:, 1] - raw]),
        fmt="o",
        color="#0F9F8F",
        ecolor="#0F9F8F",
        capsize=5,
        label="all-wire raw fit + date-cluster 95% interval",
    )
    axis.errorbar(
        x + 0.08,
        published,
        yerr=published_se,
        fmt="s",
        color="#2563EB",
        ecolor="#2563EB",
        capsize=5,
        label="published binned estimate ± SE",
    )
    axis.axhline(2.0, color="#D97706", linestyle="--", linewidth=1.2, label="quadratic reference")
    axis.set_xticks(x, lines)
    axis.set_ylim(0, 3.1)
    axis.set_ylabel("log-log exponent of apparent η vs wire length")
    axis.grid(axis="y", color="#E2E8F0", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper left", frameon=False, fontsize=7)
    axis.set_title("Probe-scale dependence is present; estimator choice changes the exponent", fontweight="bold")
    fig.tight_layout()
    path = scratch / "probe_exponents.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figures["probe_exponents"] = path

    # Date-disjoint confusion matrix.
    matrix = np.asarray(
        mechanics["wire_head"]["experimental_date_disjoint"]["metrics"]["test"]["cell_line"][
            "confusion_matrix_true_by_predicted"
        ],
        dtype=float,
    )
    fig, axis = plt.subplots(figsize=(5.2, 4.1))
    image = axis.imshow(matrix, cmap="Blues", vmin=0, vmax=matrix.max())
    for row in range(3):
        for col in range(3):
            axis.text(col, row, f"{int(matrix[row, col])}", ha="center", va="center", fontsize=12, color="white" if matrix[row, col] > matrix.max() * 0.52 else "#0F172A", fontweight="bold")
    labels = ["MCF-10A", "MCF-7", "MDA-MB-231"]
    axis.set_xticks(range(3), labels, rotation=20)
    axis.set_yticks(range(3), labels)
    axis.set_xlabel("predicted")
    axis.set_ylabel("true")
    axis.set_title("Experimental-date holdout", fontweight="bold")
    fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    fig.tight_layout()
    path = scratch / "cell_line_confusion.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figures["cell_line_confusion"] = path

    # Decisive C-1 radius sweep.
    radii = np.linspace(1, 5, 200)
    fig, axis = plt.subplots(figsize=(8.8, 3.25))
    axis.plot(radii, radii**2, color="#0F9F8F", linewidth=2.8, label="poroelastic: τ/τ₁ = (a/1 μm)²")
    axis.plot(radii, np.ones_like(radii), color="#D97706", linewidth=2.8, label="turnover / drainage / speed-set: τ/τ₁ = 1")
    axis.scatter([1, 5], [1, 25], color="#0F9F8F", s=45)
    axis.scatter([1, 5], [1, 1], color="#D97706", s=45)
    axis.annotate("25×", (5, 25), xytext=(4.45, 20.2), color="#0F9F8F", fontsize=10, fontweight="bold", arrowprops={"arrowstyle": "->", "color": "#0F9F8F"})
    axis.annotate("1×", (5, 1), xytext=(4.45, 4.3), color="#D97706", fontsize=10, fontweight="bold", arrowprops={"arrowstyle": "->", "color": "#D97706"})
    axis.set_xlim(0.9, 5.1)
    axis.set_ylim(0, 27)
    axis.set_xlabel("contact radius a (μm)")
    axis.set_ylabel("normalized relaxation time τ(a) / τ(1 μm)")
    axis.grid(color="#E2E8F0", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper left", frameon=False, fontsize=7)
    axis.set_title("C-1 discriminator: exponent survives magnitude uncertainty", fontweight="bold")
    fig.tight_layout()
    path = scratch / "c1_radius_sweep.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    figures["c1_radius_sweep"] = path

    return figures


def _page_frame(canvas: Canvas, doc: BaseDocTemplate) -> None:
    if doc.page == 1:
        return
    width, height = A4
    canvas.saveState()
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.5)
    canvas.line(18 * mm, height - 14 * mm, width - 18 * mm, height - 14 * mm)
    canvas.setFillColor(MUTED)
    canvas.setFont("Nanum", 6.5)
    canvas.drawString(18 * mm, height - 10.5 * mm, "PROJECT ALEPH · EXTERNAL NETWORK ARCHITECTURE v3")
    canvas.drawRightString(width - 18 * mm, 10 * mm, str(doc.page))
    canvas.drawString(18 * mm, 10 * mm, "2026-08-10 · evidence-bounded · no Aleph parameter authority")
    canvas.restoreState()


def _build_story(
    styles: dict[str, ParagraphStyle],
    mechanics: dict[str, Any],
    cortical: dict[str, Any],
    vision_v3: dict[str, Any],
    vision_external: dict[str, Any],
    multitask: dict[str, Any],
    afm_c2: dict[str, Any],
    encoder_reports: dict[str, dict[str, Any]],
    figures: dict[str, Path],
) -> list[Flowable]:
    story: list[Flowable] = [CoverFlowable(A4[0], A4[1]), NextPageTemplate("regular"), PageBreak()]

    # 1. Executive update.
    story += _section_header(
        "01",
        "이번 개정에서 달라진 것",
        "원본의 목표 아키텍처를 유지하되, 2026-08-10에 실제로 랜딩한 mechanics 및 cortical-tension 경로를 구현 증거로 추가했다.",
        styles,
    )
    story.append(
        _metric_cards(
            [
                ("30 + 1", "현재 multitask artifact 30 heads + mechanics optional head", "#06B6D4"),
                ("9 · 19 · 196", "mechanics primary sources · cohorts · raw wires", "#0F9F8F"),
                ("0.980 / 0.385", "η reconstruction R² / cell-line balanced accuracy", "#D97706"),
            ],
            styles,
        )
    )
    story.append(Spacer(1, 10))
    changes = [
        ["층", "현재 상태", "이번 개정의 표현"],
        ["공통 image → state 본체", "목표 / 미구현", "점선으로 표시하고 task head와 분리"],
        ["mechanics protocol network", "구현 / training-only", "데이터·split·성공·실패·refusal을 모두 수치화"],
        ["MCF-7 cortical tension", "exact-context reference head 구현", "270 [180, 400] pN/μm를 조건과 함께 제시; engine 비교는 차단"],
        ["C-1 probe-radius experiment", "preregistered discriminator", "a = 1→5 μm에서 25× vs 1×를 시각화; wire length와 AFM radius는 분리"],
    ]
    story.append(_table(changes, [33 * mm, 49 * mm, 92 * mm], styles, font_size=6.7))
    story.append(Spacer(1, 9))
    story.append(
        _callout(
            "핵심 판정",
            "이 자료는 첫 단계 외부 mechanics network로는 충분하다. 그러나 새로운 세포·연구실·방법을 아우르는 일반 cell-state network로 부르기에는 불충분하다. 실패 지표를 숨기지 않고 런타임 refusal로 연결한 것이 현재의 가장 중요한 성과다.",
            styles,
            AMBER,
        )
    )
    story.append(Spacer(1, 9))
    story.append(Paragraph("읽는 법", styles["h2"]))
    story.append(_bullet("초록 실선: 외부 데이터와 코드로 현재 실행 가능한 경로", styles, "#16A34A"))
    story.append(_bullet("청록 실선: 표현 학습은 가능하지만 외부 상태 권한이 없는 경로", styles, "#0F9F8F"))
    story.append(_bullet("회색 점선: 목표 아키텍처 또는 observation-operator 검증 대기", styles, "#64748B"))
    story.append(_bullet("빨강: OOD, split, protocol 또는 engine evidence가 부족해 값 공개를 거부하는 경로", styles, "#DC2626"))
    story.append(PageBreak())

    # 2. Architecture.
    story += _section_header(
        "02",
        "목표 아키텍처와 현재 구현의 경계",
        "외부 신경망의 제품은 label 하나가 아니라, 관측에 의해 제한된 posterior와 '이 관측으로는 결정할 수 없음'의 명시적 집합이다.",
        styles,
    )
    story.append(ArchitectureFlowable(174 * mm, 118 * mm))
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "현재의 정확한 이름",
            "하나의 통합 신경망이 아니라 provenance-aware task modules의 집합이다. multitask artifact에는 30개 head가 있고, mechanics report를 명시적으로 공급하면 evaluator가 31번째 fail-closed head를 추가한다. 공통 visual compiler와 shared cell-state latent는 아직 구현 대상이다.",
            styles,
            CYAN,
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        _table(
            [
                ["출력 층", "현재 허용", "현재 금지"],
                ["evidence posterior", "task-scoped direction, exact-context summary, ambiguity", "unseen scope를 가장 가까운 값으로 보간"],
                ["cell-state posterior", "검증된 축만 개별적으로 보고", "모든 modality를 하나의 disease/state label로 압축"],
                ["Aleph adapter", "future bounded proposal after operator+sensitivity", "문헌 값을 곧바로 엔진 파라미터로 대입"],
            ],
            [42 * mm, 64 * mm, 68 * mm],
            styles,
            font_size=6.6,
        )
    )
    story.append(PageBreak())

    # 2A. Implemented neural runtime.
    story += _section_header(
        "02A",
        "실제 신경망 실행 구조",
        "현재 구현은 입력 modality별 encoder와 task별 head를 독립적으로 학습하고, head 출력만 evidence gate에서 합치는 late-fusion runtime이다.",
        styles,
    )
    story.append(CurrentNeuralRuntimeFlowable(174 * mm, 111 * mm))
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "구조의 핵심",
            f"artifact가 선언한 값은 shared_latent={str(multitask['architecture']['shared_latent']).lower()}, fusion={multitask['architecture']['fusion']}, trained encoder heads={len(multitask['architecture']['trained_encoder_heads'])}, task heads={multitask['architecture']['head_count']}다. 따라서 한 modality의 embedding을 다른 modality head에 무단으로 넣지 않는다.",
            styles,
            CYAN,
        )
    )
    story.append(Spacer(1, 6))
    story.append(
        _table(
            [
                ["lane", "입력 tensor contract", "encoder output", "fusion rule"],
                ["image", "B×C×H×W + provider/context", "spatial map 또는 pooled vector", "해당 image head에만 전달"],
                ["temporal", "B×1×128 normalized curve", "8 filters × {mean,max} = 16 features", "calcium head에만 전달"],
                ["molecular", "B×G fit-only transformed panel", "hidden state 또는 calibrated logits", "dataset-specific state head"],
                ["morphology", "B×87 curated features", "β-gal + four biomarker outputs", "same-provider repeat gate"],
                ["mechanics", "B×F + method/probe/unit/context", "operator reconstruction + route", "optional 31st fail-closed head"],
            ],
            [27 * mm, 55 * mm, 50 * mm, 42 * mm],
            styles,
            font_size=5.9,
        )
    )
    story.append(PageBreak())

    # 2B. Concrete encoder inventory.
    calcium_arch = encoder_reports["calcium"]["architecture"]
    cell_cycle = encoder_reports["cell_cycle"]
    hpa_raw = encoder_reports["hpa_raw"]
    label_free = encoder_reports["label_free"]["architecture"]
    sciplex = encoder_reports["sciplex"]["model"]
    senscout = encoder_reports["senscout"]
    cell_death = encoder_reports["cell_death"]["model"]
    senescence = encoder_reports["senescence"]["model"]
    pbmc = encoder_reports["pbmc"]["implementation"]
    story += _section_header(
        "02B",
        "학습된 9개 encoder: 레이어와 차원",
        "아래 표는 계획도가 아니라 committed training report에서 직접 읽은 현재 구조다. 서로 다른 입력·loss·split을 하나의 parameter tensor로 공유하지 않는다.",
        styles,
    )
    encoder_rows = [
        ["head", "input → layers", "output / calibration", "권한"],
        [
            "calcium",
            f"1×128 → Conv1d(1,{calcium_arch['temporal_convolution']['filter_count']},k={calcium_arch['temporal_convolution']['kernel_size']}) → ReLU → global mean+max",
            f"16 → 2-class linear softmax; {calcium_arch['parameter_count']} params",
            "external classification blocked",
        ],
        [
            "cell cycle",
            f"3×8×8 = 192 → MLP 192→32→7; {cell_cycle['optimizer']}",
            f"7 phases; selected epoch {cell_cycle['best_selection_epoch']}",
            "same-provider training only",
        ],
        [
            "HPA raw IF",
            "2ch → Conv 16→32→64→96; GroupNorm+GELU+pool → GAP",
            f"{hpa_raw['architecture'][-1]} localization logits",
            "external domain-shift blocked",
        ],
        [
            "label-free",
            f"BF/PC/DIC → {label_free['family']}",
            f"organelle spatial map; {label_free['parameter_count']:,} params; study adversary={label_free['fit_study_adversary_classes']}",
            "external spatial transfer blocked",
        ],
        [
            "sci-Plex",
            f"{sciplex['selected_gene_count']} genes → MLP {sciplex['selected_gene_count']}→96",
            f"pathway {sciplex['pathway_class_count']} + cell type {sciplex['cell_type_class_count']} + time {sciplex['time_class_count']} + nuisance {sciplex['nuisance_group_count']}",
            "unseen-compound same-provider",
        ],
        [
            "SenSCOUT",
            f"{senscout['input']['feature_count']} curated morphology features → signed-log/z-score",
            "β-gal logistic + ridge(HMGB1, LMNB1, P16, P21)",
            "same-provider repeat only",
        ],
        [
            "cell death",
            "DeepProfiler well profile → train-compound PCA 32",
            f"6-class softmax; L2={cell_death['selected_l2']}",
            "unseen-compound rejected",
        ],
        [
            "senescence RNA",
            "35-marker logCPM fit-standardized panel",
            f"binary logistic; L2={senescence['selected_l2']}; calibrated direction",
            "external direction only",
        ],
        [
            "PBMC",
            f"expression → {pbmc['classifier']}",
            "cell type 6 + IFN 2; scalar temperature + conformal",
            "external uncertainty blocked",
        ],
    ]
    story.append(_table(encoder_rows, [27 * mm, 67 * mm, 52 * mm, 28 * mm], styles, font_size=5.35))
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "읽을 때 주의",
            "trained는 production-ready와 같은 말이 아니다. 9개 encoder는 실제 최적화된 모델이지만, 현재 production_eligible=false이며 external holdout·calibration·OOD 조건이 head마다 다르다.",
            styles,
            AMBER,
        )
    )
    story.append(PageBreak())

    # 2C. Head registry and authority flow.
    passed_ids = set(multitask["readiness"]["passed_external_evidence_heads"])
    blocked_ids = set(multitask["readiness"]["blocked_external_classifier_heads"])
    training_ids = set(multitask["readiness"]["training_only_heads"])
    all_ids = [head["head_id"] for head in multitask["heads"]]
    other_ids = [head_id for head_id in all_ids if head_id not in passed_ids | blocked_ids | training_ids]
    story += _section_header(
        "02C",
        "30개 head bank와 출력 권한",
        "head 수는 하나의 정확도 점수가 아니다. 각 head는 자기 modality, split, calibration, external evidence status를 소유한다.",
        styles,
    )
    head_rows = [["group", "head IDs", "count"]]
    for label, ids, color_label in (
        ("evidence-eligible", [item for item in all_ids if item in passed_ids], "PASS"),
        ("external blocked", [item for item in all_ids if item in blocked_ids], "BLOCKED"),
        ("training-only", [item for item in all_ids if item in training_ids], "TRAIN"),
        ("operator / constrained", other_ids, "LIMITED"),
    ):
        head_rows.append([f"{color_label}<br/>{label}", "<br/>".join(ids), str(len(ids))])
    story.append(_table(head_rows, [34 * mm, 128 * mm, 12 * mm], styles, font_size=5.0))
    story.append(Spacer(1, 7))
    story.append(
        _table(
            [
                ["adapter", "30-head artifact와의 관계", "현재 출력"],
                ["mechanics protocol", "report 입력 시 optional 31st fail-closed head", "η operator reconstruction, protocol routing, refusal"],
                ["cortical tension", "exact-context reference adapter; 아직 30-head bank 밖", "270 [180,400] pN/μm reference; engine magnitude redacted"],
                ["Aleph sweep gate", "모든 evidence head 다음의 mandatory firewall", "numeric parameter range=false; physics mutation=false"],
            ],
            [35 * mm, 78 * mm, 61 * mm],
            styles,
            font_size=5.8,
        )
    )
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "최종 데이터 흐름",
            "raw input → modality-specific preprocessing → trained encoder 또는 cited operator → task head → calibration/OOD/split gate → posterior·interval·refusal. mechanics와 cortical tension은 이 마지막 두 층에 들어오며, 공통 latent를 완성한 것처럼 표시하지 않는다.",
            styles,
            TEAL,
        )
    )
    story.append(PageBreak())

    # 3. Evidence map.
    story += _section_header(
        "03",
        "현재 외부 증거 지도",
        "입력 형식, split 단위, 외부 holdout, 불확실성 권한과 Aleph 권한을 분리해야 head 수가 성숙도로 오해되지 않는다.",
        styles,
    )
    evidence_rows = [
        ["lane", "대표 입력", "현재 증거", "상태", "Aleph 권한"],
        ["label-free vision", "BF → organelle map", "Allen 32 cells / 15 plates post-refusal diagnostic", "external spatial transfer BLOCKED", "none"],
        ["calcium", "single-cell time series", "독립 lab에서 GsMTx4 direction 2/2", "direction evidence only", "none"],
        ["migration / EMT", "PIV, IF, RNA, condition contrast", "task별 독립 방향 또는 representation", "mixed: some evidence, some blocked", "none"],
        ["mechanics protocol", "method · scope · probe · unit · time", "9 papers / 19 cohorts; exact route or abstain", "routing representation only", "none"],
        ["MRS wire network", "196 individual wires", "date-disjoint η reconstruction / state transfer", "operator learned; state transfer failed", "none"],
        ["cortical tension", "dynamic AFM group summary", "MCF-7 suspended interphase exact context", "reference ready; engine comparison blocked", "none"],
    ]
    story.append(_table(evidence_rows, [26 * mm, 36 * mm, 55 * mm, 41 * mm, 16 * mm], styles, font_size=5.8))
    story.append(Spacer(1, 10))
    story.append(Paragraph("세 가지 불확실성은 같은 confidence가 아니다", styles["h2"]))
    uncertainty = [
        ["aleatoric", "blur, noise, weak signal", "interval width / measurement error"],
        ["epistemic", "unseen cell, lab, microscope, method", "dataset/lab holdout, embedding OOD, refusal"],
        ["identifiability", "여러 기전이 같은 곡선·형태를 만듦", "ambiguity set + 추가 assay + in silico clock ablation"],
    ]
    cards = []
    for title, cause, action in uncertainty:
        cards.append(
            Table(
                [[Paragraph(title, styles["h2"])], [Paragraph(cause, styles["small"])], [Paragraph(action, styles["tiny"])]],
                colWidths=[55 * mm],
                style=TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, -1), LIGHT),
                        ("BOX", (0, 0), (-1, -1), 0.6, LINE),
                        ("TOPPADDING", (0, 0), (-1, -1), 6),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ]
                ),
            )
        )
    story.append(Table([cards], colWidths=[58 * mm] * 3, style=TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")])))
    story.append(Spacer(1, 10))
    story.append(
        _callout(
            "권한 원칙",
            "외부 head는 엔진 파라미터를 선택하거나 physics를 수정하지 않는다. 같은 단위라도 measurement operator가 다르면 숫자를 비교하지 않으며, 단일 논문 aggregate를 가짜 개별 세포로 복제하지 않는다.",
            styles,
            RED,
        )
    )
    story.append(PageBreak())

    # 4. Vision.
    story += _section_header(
        "04",
        "실제 픽셀 경로: 개선됐지만 외부 확인은 거부",
        "v3 내부 selection은 Pearson 0.5066, SSIM 0.4727, MAE +31.23%였고 mito/nucleus coverage는 93.65/93.16%였다. 그러나 fine TOMM20 network와 OOD 분리를 해결하지 못했다. 아래 값은 geometry refusal 이후의 비권위 진단이다.",
        styles,
    )
    story.append(Image(str(figures["allen_crop"]), width=174 * mm, height=50 * mm))
    story.append(Spacer(1, 5))
    story.append(Image(str(figures["vision_metrics"]), width=174 * mm, height=45 * mm))
    story.append(Spacer(1, 6))
    vision_cards = _metric_cards(
        [
            ("0.3679", "Allen v3 Pearson (v2: 0.1306)", "#06B6D4"),
            ("0.3590", "Allen v3 global SSIM (v2: 0.1022)", "#0F9F8F"),
            ("0.4292", "v3 external OOD AUROC: inadequate", "#DC2626"),
        ],
        styles,
    )
    story.append(vision_cards)
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "왜 PASS가 아닌가",
            "공식 crop DNA의 exact reconstruction이 실패해 preregistered confirmation status는 refused_geometry_exact_match다. v3의 +2.30% MAE improvement는 10% 기준에 못 미치고, 95.84% interval coverage는 frozen upper bound를 넘는다. 영상이 더 그럴듯해졌다는 사실은 external authority를 만들지 않는다.",
            styles,
            RED,
        )
    )
    story.append(PageBreak())

    # 5. Mechanics problem.
    story += _section_header(
        "05",
        "'세포질 점성'은 하나의 표적이 아니다",
        "프로브가 무엇을 밀고, 어느 길이·시간 스케일을 평균하며, cell surface가 어떤 상태인지가 관측량을 정의한다. 서로 다른 값을 하나의 η label로 pooling하면 모델이 방법을 생물학으로 오인한다.",
        styles,
    )
    story.append(Image(str(figures["viscosity_scale"]), width=174 * mm, height=59 * mm))
    story.append(Spacer(1, 7))
    story.append(
        _table(
            [
                ["measurement", "scope", "핵심 해석", "audit role"],
                ["water / solvent input", "solvent viscosity", "ETA_SOLVENT = 1e-3 Pa s; 엔진 입력", "input, not target"],
                ["optical Brownian probe", "solvent microviscosity", "작은 probe가 느끼는 local solvent-like response", "NOT_REGISTERED; routing context"],
                ["magnetic wire MRS", "effective cytoplasm viscoelasticity", "섬유-용매 coupling을 포함한 composite response", "OK; individual rows"],
                ["surface magnetic bead", "local apparent viscosity", "adhesion/cortex/cytoplasm이 결합된 local operator", "OK; order-of-magnitude"],
            ],
            [37 * mm, 42 * mm, 66 * mm, 29 * mm],
            styles,
            font_size=6.1,
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        _callout(
            "MCF-7에서 보존된 충돌",
            "optical Brownian probe의 reported range 0.002-0.016 Pa s와 MRS all-wire median 65.9 Pa s는 약 네 자릿수 차이다. 전자는 repo source audit에 아직 등록되지 않았으므로 route/range 맥락으로만 표시한다. 이 둘을 평균내지 않는 것이 mechanics ontology의 첫 번째 기능이다.",
            styles,
            AMBER,
        )
    )
    story.append(PageBreak())

    # 6. Corpus.
    story += _section_header(
        "06",
        "mechanics corpus: 넓은 방법, 얕은 개별 데이터",
        "catalog는 방법 간 의미 차이를 표현할 만큼 넓지만, generalization model을 학습할 만큼 각 방법의 독립 lab과 raw curve가 충분하지 않다.",
        styles,
    )
    story.append(
        _metric_cards(
            [
                ("9", "primary papers", "#06B6D4"),
                ("19", "real cohort records; no synthetic duplication", "#0F9F8F"),
                ("196", "individual magnetic wires from one paper", "#D97706"),
            ],
            styles,
        )
    )
    story.append(Spacer(1, 9))
    source_rows = [
        ["family", "independent papers", "individual-level raw", "main scope"],
        ["magnetic rheology", "Dessard; Bausch", "Dessard wires only", "effective / local apparent viscosity"],
        ["AFM relaxation / indentation", "Moeendarbary; Flormann", "not normalized here", "poroelastic D / local stiffness"],
        ["whole-cell AFM confinement", "Hosseini; Fischer-Friedrich", "aggregate summaries", "effective surface tension"],
        ["suction / electrodeformation", "Wang; Moazzeni", "aggregate summaries", "whole-cell viscosity / short-time prestress"],
        ["optical Brownian probe", "Vaippully", "reported range only", "solvent microviscosity"],
    ]
    story.append(_table(source_rows, [45 * mm, 45 * mm, 42 * mm, 42 * mm], styles, font_size=6.2))
    story.append(Spacer(1, 9))
    story.append(Paragraph("source audit status", styles["h2"]))
    audit_table = Table(
        [[
            Paragraph("OK<br/><b>5</b>", styles["metric_label"]),
            Paragraph("CHECK<br/><b>2</b>", styles["metric_label"]),
            Paragraph("NOT_REGISTERED<br/><b>2</b>", styles["metric_label"]),
        ]],
        colWidths=[58 * mm] * 3,
        rowHeights=[18 * mm],
        style=TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, 0), HexColor("#EAF7EE")),
                ("BACKGROUND", (1, 0), (1, 0), HexColor("#FEF4E8")),
                ("BACKGROUND", (2, 0), (2, 0), HexColor("#FDECEC")),
                ("BOX", (0, 0), (-1, -1), 0.5, LINE),
                ("INNERGRID", (0, 0), (-1, -1), 0.5, LINE),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ]
        ),
    )
    story.append(audit_table)
    story.append(Spacer(1, 8))
    story.append(_bullet("OK 5개 출처만 numeric evidence의 중심에 둔다.", styles, "#16A34A"))
    story.append(_bullet("CHECK 2개는 state comparison 또는 order-of-magnitude routing에만 쓴다.", styles, "#D97706"))
    story.append(_bullet("NOT_REGISTERED 2개는 catalog gap을 드러내며 point estimate authority가 없다.", styles, "#DC2626"))
    story.append(PageBreak())

    # 7. Wire model.
    story += _section_header(
        "07",
        "wire multi-head: operator는 배우고 state는 못 배웠다",
        "동일 shared latent에 회귀와 분류를 붙였기 때문에 성공과 실패의 원인이 더 선명해졌다. η는 입력과 결정적으로 연결되지만 cell line은 acquisition-date shift 아래 분리되지 않는다.",
        styles,
    )
    story.append(WireNetworkFlowable(174 * mm, 66 * mm))
    story.append(Spacer(1, 7))
    story.append(
        _table(
            [
                ["endpoint", "sealed date-disjoint result", "정당한 해석", "runtime"],
                ["log η regression", "R² 0.979995; log MAE 0.1354; MdAPE 12.06%", "measurement-operator reconstruction", "known wire scope only"],
                ["cell line", "accuracy 0.28; BA 0.385; macro-F1 0.227", "state transfer failure", "refuse new cell/lab"],
                ["protocol router", "training resubstitution 100%", "wiring diagnostic; not generalization", "exact ontology or abstain"],
            ],
            [31 * mm, 50 * mm, 60 * mm, 33 * mm],
            styles,
            font_size=6.0,
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        _callout(
            "왜 R² 0.98이 biomarker validation이 아닌가",
            "wire geometry와 critical frequency는 viscosity 계산식의 입력이다. 따라서 회귀 head의 성공은 dynamic-response → effective-viscosity operator representation이 학습됐다는 뜻이다. 독립적인 생물학적 marker로 viscosity 또는 cell state를 맞혔다는 뜻이 아니다.",
            styles,
            CYAN,
        )
    )
    story.append(PageBreak())

    # 8. Confusion.
    story += _section_header(
        "08",
        "실패를 산출물로 만든 confusion matrix",
        "실험 날짜 전체를 holdout한 50개 wire에서 분류기는 MDA-MB-231 쪽으로 무너진다. 이 실패가 바로 new-cell inference를 차단하는 근거다.",
        styles,
    )
    image = Image(str(figures["cell_line_confusion"]), width=92 * mm, height=73 * mm)
    notes = [
        Paragraph("<b>정량 결과</b>", styles["h2"]),
        _bullet("MCF-10A recall 0.105", styles, "#DC2626"),
        _bullet("MCF-7 recall 0.050", styles, "#DC2626"),
        _bullet("MDA-MB-231 recall 1.000", styles, "#D97706"),
        _bullet("random wire holdout BA도 0.464에 불과", styles, "#D97706"),
        Spacer(1, 6),
        Paragraph("<b>의미</b>", styles["h2"]),
        Paragraph("같은 논문의 196개 technical observation은 세 cell line을 외부에서 판별하는 독립 표본이 아니다. acquisition date와 wire geometry의 변동을 state 신호로 착각할 위험이 크다.", styles["body"]),
        Paragraph("<b>정책:</b> model output이 class label보다 먼저 scope와 refusal reason을 반환한다.", styles["body"]),
    ]
    story.append(Table([[image, notes]], colWidths=[96 * mm, 78 * mm], style=TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")])))
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            '{"status":"refused","reason":"cell-line transfer failed on experimental-date holdout",'
            '"permitted":"known-protocol operator reconstruction only"}',
            styles["code"],
        )
    )
    story.append(Spacer(1, 9))
    story.append(
        _callout(
            "일반화를 위해 필요한 최소 데이터",
            "각 방법별 최소 2개 독립 lab, paper-level holdout, biological replicate identifier, raw curves/fields, protocol geometry와 physical scale가 필요하다. aggregate row를 synthetic cells로 늘리는 것은 이 문제를 해결하지 않는다.",
            styles,
            AMBER,
        )
    )
    story.append(PageBreak())

    # 9. Probe dependence.
    story += _section_header(
        "09",
        "probe-size dependence: 결과와 경고를 함께 보존",
        "모든 raw wire의 log-log fit은 η가 wire length에 의존함을 보인다. 그러나 raw fit과 paper의 binned exponent가 다르므로 숫자 하나로 접지 않는다.",
        styles,
    )
    story.append(Image(str(figures["probe_exponents"]), width=174 * mm, height=63 * mm))
    story.append(Spacer(1, 6))
    scaling = mechanics["probe_size_dependence"]["external_mrs_wire_length_scaling"]
    probe_rows = [["cell line", "raw exponent", "date-cluster 95% interval", "paper binned", "η fold at 1→5 μm"]]
    for cell_line in ("MCF-10A", "MCF-7", "MDA-MB-231"):
        row = scaling[cell_line]
        interval = row["experiment_date_cluster_bootstrap_95_interval"]
        probe_rows.append(
            [
                cell_line,
                f'{row["raw_wire_level_log_log_exponent"]:.3f}',
                f"[{interval[0]:.3f}, {interval[1]:.3f}]",
                f'{row["published_binned_exponent"]:.1f} ± {row["published_standard_error"]:.1f}',
                f'{row["predicted_eta_fold_1_to_5_um_from_raw_fit"]:.1f}×',
            ]
        )
    story.append(_table(probe_rows, [34 * mm, 31 * mm, 48 * mm, 31 * mm, 30 * mm], styles, font_size=6.3))
    story.append(Spacer(1, 8))
    story.append(
        _callout(
            "중요한 비동일성",
            "Dessard의 wire length는 AFM contact radius a가 아니다. 이 분석은 apparent cytoplasmic viscosity가 probe scale에 의존한다는 외부 근거다. C-1의 τ(a) ∝ a²를 검증하지는 않는다.",
            styles,
            RED,
        )
    )
    story.append(PageBreak())

    # 10. C-1.
    story += _section_header(
        "10",
        "C-1 결정적 실험: 네 시계를 exponent로 분리",
        "단일 probe의 relaxation curve는 poroelasticity, crosslink turnover, osmotic drainage, speed-set viscous transient를 모두 '점탄성'처럼 보이게 한다. contact radius sweep만이 하나의 시계를 움직인다.",
        styles,
    )
    story.append(Image(str(figures["c1_radius_sweep"]), width=174 * mm, height=64 * mm))
    story.append(Spacer(1, 6))
    clock_rows = [
        ["clock", "현재 크기", "a dependence", "in silico ablation"],
        ["poroelastic τp", "~1 s (cross-cell D context)", "a² / D", "Biot/FSI path off"],
        ["crosslink turnover", "~15 s", "a⁰", "xl_koff_per_s=None"],
        ["membrane osmotic drainage", "~30-200 s", "a⁰", "Lp=None"],
        ["press-speed viscous transient", "speed-set", "a⁰", "v_press_um_s=None"],
    ]
    story.append(_table(clock_rows, [51 * mm, 43 * mm, 34 * mm, 46 * mm], styles, font_size=6.2))
    story.append(Spacer(1, 7))
    story.append(
        _callout(
            "사전 선언할 scalar",
            "각 a에서 동일 투영으로 stress-relaxation τ를 fit하고, log τ 대 log a의 slope를 채점한다. poroelastic hypothesis는 slope 2, probe-independent hypotheses는 slope 0이다. magnitude가 흔들려도 exponent 판정은 남는다.",
            styles,
            TEAL,
        )
    )
    story.append(PageBreak())

    # 11. Execution order.
    story += _section_header(
        "11",
        "C-2는 닫혔다: physical inputs → active state → measurement",
        "병렬 AFM branch의 membrane subdivision 8 full-native 3-seed run이 accepted physical step을 commit했다. 이제 blocker는 수치 acceptance가 아니라 apparatus/contact/exterior provenance와 active-cortex handoff다.",
        styles,
    )
    c2 = afm_c2["c2"]
    story.append(
        _callout(
            "C-2 PASS - numerical/transaction scope",
            f'{c2["total_nodes"]:,} nodes, all compartments, RTX 3090 Slurm, seeds 0/1/2가 각각 1 physical step을 accept하고 0.01 s를 commit했다. 이 결과는 C-2 readiness이며 material F-δ 또는 cortical-tension magnitude가 아니다.',
            styles,
            GREEN,
        )
    )
    story.append(Spacer(1, 7))
    sequence = [
        ["step", "must establish", "artifact / refusal condition", "status"],
        ["C-2", "full-native accepted-step transaction", "subdiv 8 · seeds 0/1/2 · apparatus CUDA seam", "PASS"],
        ["physical inputs", "apparatus, Π₀, contact, exterior viscosity", "source or explicit mechanism-demo-only scope", "open"],
        ["active-state handoff", "post-induction, equilibrated, dt-converged cortex", "STATE (c) 3/17 remain binding", "open"],
        ["observation operator", "AFM geometry/contact radius mapping", "same protocol, unit, surface state, time window", "open"],
        ["C-1 run", "a = 1, ..., 5 μm stress relaxation", "predeclared τ fit and slope", "not run"],
        ["external score", "τ ∝ a² vs a⁰", "uncertainty + seed scatter + ablations", "not permitted yet"],
    ]
    story.append(_table(sequence, [20 * mm, 48 * mm, 79 * mm, 27 * mm], styles, font_size=6.4))
    story.append(Spacer(1, 12))
    story.append(Paragraph("Gate contract - 실행 전에 고정", styles["h2"]))
    contract_items = [
        "cell line, surface state, temperature, indenter geometry, contact radius와 deformation amplitude",
        "τ estimator, fit window, sampling rate, projection, excluded transient와 missing-data rule",
        "radius별 독립 seed와 between-seed uncertainty; within-run SEM으로 대체 금지",
        "Lp=None, xl_koff_per_s=None, v_press_um_s=None의 off-is-backward-compatible control",
        "slope와 uncertainty를 보고하되 D_um2_s 또는 M_biot_Pa를 맞추기 위해 gate를 조정하지 않음",
    ]
    for item in contract_items:
        story.append(_bullet(item, styles, "#0F9F8F"))
    story.append(Spacer(1, 9))
    story.append(
        _callout(
            "현재 외부 신경망이 할 수 있는 추가 일",
            "geometry-matched curve가 들어오면 protocol route, τ(a) extraction request, hypothesis label, uncertainty와 refusal reason을 구조화할 수 있다. curve를 만들어내거나 unsourced D 값을 선택할 수는 없다.",
            styles,
            CYAN,
        )
    )
    story.append(PageBreak())

    # 12. Cortical tension.
    story += _section_header(
        "12",
        "cortical tension: reference head는 준비, engine comparison은 차단",
        "문헌 band와 엔진 값을 같은 단위로 놓는 것만으로는 비교가 아니다. surface state, cell-cycle state, temperature, frequency와 observation operator가 모두 맞아야 한다.",
        styles,
    )
    story.append(CorticalPipelineFlowable(174 * mm, 59 * mm))
    story.append(Spacer(1, 7))
    story.append(
        _metric_cards(
            [
                ("270", "median pN/μm", "#0F9F8F"),
                ("180-400", "empirical IQR pN/μm", "#06B6D4"),
                ("n = 27", "MCF-7 suspended rounded interphase", "#D97706"),
            ],
            styles,
        )
    )
    story.append(Spacer(1, 7))
    readiness = cortical["blockers"]
    story.append(Paragraph("comparison adapter가 요구하지만 아직 충족되지 않은 핵심 항목", styles["h2"]))
    for item in [
        "native-validated, geometry-matched total effective Laplace operator",
        "spherical apparatus source 또는 cited flat-cantilever geometry",
        "target-cell Π₀, contact rest/stiffness, exterior-medium viscosity provenance",
        "post-induction equilibrated and time-step-independent active-cortex state",
        "최소 3개 final stationary-tension record의 between-seed scatter",
        "exact context: MCF-7 / suspended rounded / interphase / 37 C / 1 Hz",
    ]:
        story.append(_bullet(item, styles, "#DC2626"))
    story.append(Paragraph(f"기존 external adapter artifact는 C-2 이전 상태라 {len(readiness)} blocker entries를 기록한다. 새 C-2 snapshot은 numerical transaction blocker만 닫으며, final tension record와 magnitude는 계속 redacted다.", styles["small"]))
    story.append(PageBreak())

    # 13. Runtime interface.
    story += _section_header(
        "13",
        "외부 신경망이 지금 반환할 수 있는 것",
        "현재 runtime의 강점은 무엇이든 맞히는 것이 아니라, 관측의 정의를 보존하고 권한 밖이면 정해진 이유로 멈추는 것이다.",
        styles,
    )
    story.append(
        _table(
            [
                ["request", "returned", "refused"],
                ["known mechanics protocol", "method family, mechanical scope, surface state, ontology identity", "unit/scope가 바뀐 nearest-neighbor magnitude"],
                ["MRS wire response", "effective η operator reconstruction + in-scope uncertainty note", "new-cell cell-line label"],
                ["future AFM radius sweep", "τ(a) hypothesis route and predeclared scoring request", "wire exponent를 AFM slope로 대체"],
                ["cortical-tension engine record", "missing-field checklist", "blocked γ magnitude or ratio"],
                ["image observation", "task-scoped spatial map / post-refusal diagnostic", "general cell-state or production authority"],
            ],
            [42 * mm, 69 * mm, 63 * mm],
            styles,
            font_size=6.4,
        )
    )
    story.append(Spacer(1, 10))
    story.append(Paragraph("예시: exact ontology route", styles["h2"]))
    story.append(
        Paragraph(
            '{"cell_line":"MCF-7","method":"magnetic_rotational_spectroscopy",'
            '"scope":"effective_cytoplasm_viscoelasticity","observable":"static_viscosity",'
            '"unit":"Pa s","status":"exact_catalog_route"}',
            styles["code"],
        )
    )
    story.append(Spacer(1, 8))
    story.append(Paragraph("예시: unseen scope refusal", styles["h2"]))
    story.append(
        Paragraph(
            '{"status":"refused","reason":"unseen mechanical_scope or unit",'
            '"policy":"no cross-scope interpolation; no Aleph parameter selection"}',
            styles["code"],
        )
    )
    story.append(Spacer(1, 10))
    story.append(
        _callout(
            "최종 output contract",
            "state posterior, evidence provenance, protocol context, uncertainty class, ambiguity set, required next assay, OOD/refusal reason을 함께 반환한다. point estimate 하나만 반환하는 API는 이 프로젝트의 목표와 맞지 않는다.",
            styles,
            TEAL,
        )
    )
    story.append(PageBreak())

    # 14. Sufficiency matrix.
    story += _section_header(
        "14",
        "무엇이 더 있어야 '충분하다'고 말할 수 있는가",
        "다음 단계의 병목은 모델 크기가 아니라 독립 raw data와 geometry-matched observation operator다.",
        styles,
    )
    gaps = [
        ["requirement", "current", "minimum next evidence", "why it matters"],
        ["cross-lab mechanics", "MRS raw: 1 paper / 1 lab", "≥2 labs per method", "method token과 lab signature 분리"],
        ["raw relaxation curves", "aggregate for most methods", "AFM curves + contact radius + metadata", "τ estimator와 uncertainty 재현"],
        ["paper-level holdout", "date-disjoint within one source", "sealed independent paper", "cell-state transfer claim"],
        ["C-1 radius sweep", "hypothesis only", "same cell/state, a = 1→5 μm", "a² vs a⁰ clock discrimination"],
        ["cortical-tension adapter", "C-2 PASS + exact-context external band", "physical inputs + active state + operator + final seeds", "engine magnitude comparison"],
        ["shared state latent", "task modules", "paired multimodal samples", "sample-level fusion without pseudo-pairing"],
        ["uncertainty/OOD", "head-specific, several failures", "external calibration + explicit abstention", "posterior reliability"],
    ]
    story.append(_table(gaps, [41 * mm, 40 * mm, 51 * mm, 42 * mm], styles, font_size=5.8))
    story.append(Spacer(1, 10))
    story.append(Paragraph("3090 사용 기준", styles["h2"]))
    story.append(_bullet("현재 196-row NumPy mechanics head와 PDF 생성은 CPU가 더 빠르고 재현 가능하다.", styles, "#16A34A"))
    story.append(_bullet("raw AFM curve 대규모 학습, 고해상도 image encoder 또는 multimodal transformer가 열릴 때만 GPU가 이득이다.", styles, "#D97706"))
    story.append(_bullet("그 경우에도 shared workstation의 Slurm allocation에서만, PI가 카드와 시간을 명시하고 3090이 비어 있을 때 사용한다.", styles, "#DC2626"))
    story.append(Spacer(1, 9))
    story.append(
        _callout(
            "다음 우선순위",
            "논문 aggregate row를 더 모으는 것보다, contact radius가 다른 AFM stress-relaxation raw curve와 독립 lab의 individual-level mechanics data를 확보하는 것이 정보량이 크다.",
            styles,
            AMBER,
        )
    )
    story.append(PageBreak())

    # 15. Roadmap.
    story += _section_header(
        "15",
        "실행 로드맵과 완료 정의",
        "outer network의 다음 성공 지표는 head 수가 아니라, 같은 runtime이 여러 acquisition 조건에서 object/state evidence와 정당한 refusal을 반환하는가다.",
        styles,
    )
    roadmap = [
        ["phase", "deliverable", "exit criterion"],
        ["A. protocol and visual ingestion", "OME/curve/field metadata audit + missing-modality mask", "동일 API, pixel/object/protocol provenance"],
        ["B. initial state axes", "cycle, morphology, mitochondria", "axis posterior + evidence overlay + OOD"],
        ["C. mechanobiology", "mechanics router, AFM/TFM/PIV/FRET encoders", "method/lab holdout and non-poolable scope"],
        ["D. shared multimodal latent", "paired sample fusion + unpaired direction constraints", "no pseudo-pairing; calibrated external posterior"],
        ["E. Aleph adapter", "operator and sensitivity-gated bounded proposal", "native measurement gate before any parameter sweep"],
    ]
    story.append(_table(roadmap, [37 * mm, 65 * mm, 72 * mm], styles, font_size=6.3))
    story.append(Spacer(1, 10))
    story.append(Paragraph("C-1 mechanics closeout checklist", styles["h2"]))
    checklist = [
        ["OPEN", "외부 protocol contract와 radius schedule를 run 전에 동결"],
        ["PASS", "C-2 full-native 3-seed accepted-step transaction"],
        ["OPEN", "physical inputs + active-state handoff + geometry-matched operator"],
        ["OPEN", "τ(a) slope와 uncertainty, clock-by-clock ablation"],
        ["OPEN", "result를 parameter가 아니라 evidence-bounded mechanism support로 반환"],
    ]
    story.append(_table(checklist, [13 * mm, 161 * mm], styles, header=False, font_size=6.7))
    story.append(Spacer(1, 10))
    story.append(
        _callout(
            "완료의 문장",
            "'모델이 값을 냈다'가 아니라 '관측과 protocol이 허용한 mechanism posterior를 냈고, 구별되지 않는 clock과 필요한 다음 실험을 함께 냈다'가 완료다.",
            styles,
            TEAL,
        )
    )
    story.append(PageBreak())

    # 16. References and reproducibility.
    story += _section_header(
        "16",
        "주요 근거와 재현 가능한 산출물",
        "숫자는 committed JSON에서 읽었고, source audit 상태와 authority를 보고서 안에서 함께 표시했다.",
        styles,
    )
    refs = [
        ["key", "source / role", "identifier"],
        ["[D1]", "Dessard et al. 2024 - cytoplasmic viscosity and probe-size dependence", '<link href="https://doi.org/10.1039/D4NA00003J">10.1039/D4NA00003J</link> · audit OK'],
        ["[D2]", "Moeendarbary et al. - cytoplasm as a poroelastic material", '<link href="https://doi.org/10.1038/nmat3517">10.1038/nmat3517</link> · audit OK'],
        ["[D3]", "Hosseini et al. 2020 - MCF-7 dynamic AFM effective tension", '<link href="https://doi.org/10.1002/advs.202001276">10.1002/advs.202001276</link> · audit OK'],
        ["[D4]", "Fischer-Friedrich et al. 2014 - surface tension and internal pressure", '<link href="https://doi.org/10.1038/srep06213">10.1038/srep06213</link> · audit OK'],
        ["[D5]", "Bausch et al. 1998 - local magnetic-bead microrheometry", '<link href="https://doi.org/10.1016/S0006-3495(98)77646-5">10.1016/S0006-3495(98)77646-5</link> · audit OK'],
        ["[C1]", "Vaippully et al. 2020 - optical Brownian-probe range", '<link href="https://doi.org/10.1088/1361-648X/ab76ac">10.1088/1361-648X/ab76ac</link> · NOT_REGISTERED; routing context'],
        ["[C2]", "Moazzeni et al. 2021 - electrodeformation tension", '<link href="https://doi.org/10.1103/PhysRevE.103.032409">10.1103/PhysRevE.103.032409</link> · CHECK; OOM route'],
        ["[C3]", "Flormann et al. 2024 - cortex structure and mechanics", '<link href="https://doi.org/10.1073/pnas.2320372121">10.1073/pnas.2320372121</link> · CHECK; state comparison'],
    ]
    story.append(_table(refs, [13 * mm, 91 * mm, 70 * mm], styles, font_size=5.7))
    story.append(Spacer(1, 9))
    artifact_rows = [
        ["artifact", "role"],
        ["mechanics_protocol_report.json", "split, metrics, exponents, capability and refusal source"],
        ["source_receipt.json", "publisher workbook identity, MD5/SHA-256, 196-row normalization receipt"],
        ["cortical_tension_engine_readiness.json", "exact-context reference and engine comparison blocker ledger"],
        ["afm_c2_readiness_snapshot.json", "sibling AFM branch C-2 accepted evidence, source commit and hashes"],
        ["allencell_label_free_confirmation_evaluation.json", "geometry refusal and post-refusal v2/v3 diagnostic"],
        ["multitask_outer_model_report.json", "30-head current artifact; mechanics head is optional evaluator input"],
    ]
    story.append(_table(artifact_rows, [70 * mm, 104 * mm], styles, font_size=6.0))
    story.append(Spacer(1, 9))
    story.append(
        _callout(
            "최종 상태",
            "protocol-aware mechanics representation, cortical-tension exact-context reference, C-2 full-native accepted transaction, source receipts, leakage-safe date split, explicit failure와 refusal은 준비됐다. 공통 visual compiler, shared probabilistic cell-state latent, independent mechanics lab holdout와 quantitative AFM/engine comparison은 아직 준비되지 않았다.",
            styles,
            CYAN,
        )
    )

    return story


def build_report(output: Path, scratch: Path) -> dict[str, Any]:
    """Create the revised PDF and return its reproducibility manifest."""
    _register_fonts()
    styles = _styles()

    input_paths = {
        "mechanics_report": MECHANICS_OUTPUT / "model/mechanics_protocol_report.json",
        "mechanics_receipt": MECHANICS_OUTPUT / "source_receipt.json",
        "cortical_readiness": OUTER / "results/cortical_tension_engine_readiness.json",
        "vision_v3": OUTER / "results/lightmycells_model_training_report_v3.json",
        "vision_external": OUTER / "results/allencell_label_free_confirmation_evaluation.json",
        "multitask": OUTER / "results/multitask_outer_model_report.json",
        "allen_figure": OUTER / "results/allencell_label_free_v2_v3_comparison.png",
        "afm_c2_snapshot": MECHANICS_OUTPUT / "afm_c2_readiness_snapshot.json",
        "encoder_calcium": OUTER / "results/calcium_temporal_cnn_report.json",
        "encoder_cell_cycle": OUTER / "results/bbbc048_mlp_report.json",
        "encoder_hpa_raw": OUTER / "results/hpa_raw_if_cnn_training_report.json",
        "encoder_sciplex": OUTER / "results/sciplex3_state_model_report.json",
        "encoder_senscout": OUTER / "results/senscout_morphology_model_report.json",
        "encoder_cell_death": OUTER / "results/figshare_cell_death_model_report.json",
        "encoder_senescence": OUTER / "results/external_senescence_model_report.json",
        "encoder_pbmc": OUTER / "results/pbmc_cross_lab_model_report.json",
    }
    missing = [str(path) for path in input_paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Required committed report inputs are missing: {missing}")

    mechanics = _load_json(input_paths["mechanics_report"])
    cortical = _load_json(input_paths["cortical_readiness"])
    vision_v3 = _load_json(input_paths["vision_v3"])
    vision_external = _load_json(input_paths["vision_external"])
    multitask = _load_json(input_paths["multitask"])
    afm_c2 = _load_json(input_paths["afm_c2_snapshot"])
    encoder_reports = {
        "calcium": _load_json(input_paths["encoder_calcium"]),
        "cell_cycle": _load_json(input_paths["encoder_cell_cycle"]),
        "hpa_raw": _load_json(input_paths["encoder_hpa_raw"]),
        "label_free": vision_v3,
        "sciplex": _load_json(input_paths["encoder_sciplex"]),
        "senscout": _load_json(input_paths["encoder_senscout"]),
        "cell_death": _load_json(input_paths["encoder_cell_death"]),
        "senescence": _load_json(input_paths["encoder_senescence"]),
        "pbmc": _load_json(input_paths["encoder_pbmc"]),
    }
    figures = _prepare_figures(scratch, mechanics, vision_external)

    output.parent.mkdir(parents=True, exist_ok=True)
    width, height = A4
    document = BaseDocTemplate(
        str(output),
        pagesize=A4,
        title="Project Aleph 외부 신경망 및 Cell-State 추론 구조 보고서 v3",
        author="Project Aleph",
        subject="External neural network architecture, mechanics protocol representation, cortical tension and C-1 probe-radius discriminator",
        creator="Project Aleph reproducible ReportLab builder",
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=18 * mm,
        bottomMargin=16 * mm,
    )
    cover_frame = Frame(0, 0, width, height, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0, id="cover")
    regular_frame = Frame(
        18 * mm,
        16 * mm,
        width - 36 * mm,
        height - 34 * mm,
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
        id="regular",
    )
    document.addPageTemplates(
        [
            PageTemplate(id="cover", frames=[cover_frame], onPage=_page_frame),
            PageTemplate(id="regular", frames=[regular_frame], onPage=_page_frame),
        ]
    )
    document.build(
        _build_story(
            styles,
            mechanics,
            cortical,
            vision_v3,
            vision_external,
            multitask,
            afm_c2,
            encoder_reports,
            figures,
        )
    )

    manifest = {
        "schema": "aleph.outer.external_network_architecture_report.v3",
        "report": {
            "path": str(output.relative_to(ROOT) if output.is_relative_to(ROOT) else output),
            "sha256": _sha256(output),
            "size_bytes": output.stat().st_size,
        },
        "inputs": {
            name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": _sha256(path),
            }
            for name, path in input_paths.items()
        },
        "claims": {
            "aleph_parameter_authority": "none",
            "gpu_used": False,
            "physics_run": False,
            "validation_gate_scored": False,
            "mechanics_state_transfer": "failed",
            "afm_c2": "accepted_3_seeds_numerical_transaction_scope",
            "cortical_engine_comparison": "blocked",
            "shared_cell_state_network": "not_implemented",
            "implemented_neural_runtime": "9_trained_encoders_30_task_heads_task_scoped_late_fusion",
        },
    }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=MECHANICS_OUTPUT / "external_neural_network_architecture_report_2026-08-10_v3.pdf",
    )
    parser.add_argument(
        "--scratch",
        type=Path,
        default=ROOT / "tmp/pdfs/external_nn_report_v3",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=MECHANICS_OUTPUT / "architecture_report_v3_manifest.json",
    )
    parser.add_argument(
        "--delivery-copy",
        type=Path,
        default=None,
        help="Optional exact copy, for example output/pdf/<name>.pdf.",
    )
    args = parser.parse_args()
    manifest = build_report(args.output.resolve(), args.scratch.resolve())
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.delivery_copy is not None:
        args.delivery_copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.output, args.delivery_copy)
    print(json.dumps(manifest["report"], sort_keys=True))


if __name__ == "__main__":
    main()
