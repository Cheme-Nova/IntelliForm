"""
modules/pdf_proposal.py
Branded PDF proposal generator for IntelliForm.

Full-bleed cover page + branded content pages using ReportLab.
Brand: ChemeNova navy #0A1628 · teal #0D9488 · amber #D97706.

Pages:
  1 — Full-bleed cover (canvas-drawn, no flowables)
  2 — Executive summary: blend composition + key metrics
  3 — EcoMetrics™ radar chart + axis table (if available)
  4 — Regulatory intelligence table + certification pathways (if available)
  5 — Next steps, financial summary, contact
"""
import io
import math
from datetime import datetime
from typing import Dict, Optional

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, PageBreak, KeepTogether,
)
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.spider import SpiderChart

# ── Brand palette ─────────────────────────────────────────────────────────────
NAVY    = colors.HexColor("#0A1628")
NAVY2   = colors.HexColor("#1E3A5F")   # lighter navy for cards
TEAL    = colors.HexColor("#0D9488")
AMBER   = colors.HexColor("#D97706")
WHITE   = colors.white
LIGHT   = colors.HexColor("#F0F9FF")
SLATE   = colors.HexColor("#64748B")
MUTED   = colors.HexColor("#94A3B8")
GREEN   = colors.HexColor("#00C853")
RED_C   = colors.HexColor("#EF4444")
RULE    = colors.HexColor("#E2E8F0")
DARK_TXT = colors.HexColor("#1E293B")

PW, PH = letter   # 8.5 × 11 in


# ── Styles ────────────────────────────────────────────────────────────────────

def _styles() -> dict:
    base = getSampleStyleSheet()
    extra = {
        "section_head": ParagraphStyle(
            "section_head", fontSize=13, leading=16,
            textColor=NAVY, fontName="Helvetica-Bold",
            spaceBefore=14, spaceAfter=6,
        ),
        "body": ParagraphStyle(
            "body", fontSize=9, leading=13,
            textColor=DARK_TXT, fontName="Helvetica", spaceAfter=4,
        ),
        "small": ParagraphStyle(
            "small", fontSize=7.5, leading=10,
            textColor=SLATE, fontName="Helvetica",
        ),
        "metric_label": ParagraphStyle(
            "metric_label", fontSize=8, leading=10,
            textColor=SLATE, fontName="Helvetica", alignment=TA_CENTER,
        ),
        "metric_value": ParagraphStyle(
            "metric_value", fontSize=20, leading=24,
            textColor=NAVY, fontName="Helvetica-Bold", alignment=TA_CENTER,
        ),
        "metric_unit": ParagraphStyle(
            "metric_unit", fontSize=7.5, leading=9,
            textColor=MUTED, fontName="Helvetica", alignment=TA_CENTER,
        ),
        "warn": ParagraphStyle(
            "warn", fontSize=9, leading=12,
            textColor=AMBER, fontName="Helvetica-Bold",
        ),
    }
    return {**{k: base[k] for k in base.byName}, **extra}


# ── Molecule icon ─────────────────────────────────────────────────────────────

def _draw_molecule(canvas, cx: float, cy: float, r: float):
    """Draw a simple hexagonal molecule motif centred at (cx, cy) with radius r."""
    canvas.setStrokeColor(TEAL)
    canvas.setFillColor(TEAL)
    canvas.setLineWidth(1.8)

    # Hexagon outline
    pts = [
        (cx + r * math.cos(math.pi / 2 + i * math.pi / 3),
         cy + r * math.sin(math.pi / 2 + i * math.pi / 3))
        for i in range(6)
    ]
    p = canvas.beginPath()
    p.moveTo(*pts[0])
    for x, y in pts[1:]:
        p.lineTo(x, y)
    p.close()
    canvas.drawPath(p, fill=0, stroke=1)

    # Alternate-vertex spokes
    for i in range(0, 6, 2):
        canvas.line(cx, cy, pts[i][0], pts[i][1])

    # Centre dot
    canvas.circle(cx, cy, r * 0.14, fill=1, stroke=0)


# ── Cover page ────────────────────────────────────────────────────────────────

def _make_cover(project: Dict, customer_name: str, customer_company: str,
                logo_bytes: Optional[bytes], date_str: str):
    """Return an onFirstPage callback that paints the full-bleed cover."""

    def _draw(canvas, doc):
        canvas.saveState()

        # ── Full navy background ──────────────────────────────────────────────
        canvas.setFillColor(NAVY)
        canvas.rect(0, 0, PW, PH, fill=1, stroke=0)

        # Amber bar — very top
        canvas.setFillColor(AMBER)
        canvas.rect(0, PH - 0.13 * inch, PW, 0.13 * inch, fill=1, stroke=0)

        # Teal left accent strip
        canvas.setFillColor(TEAL)
        canvas.rect(0, 0, 0.18 * inch, PH, fill=1, stroke=0)

        # ── IntelliForm logo block ────────────────────────────────────────────
        lx = 0.65 * inch   # left margin for logo area
        ly = PH - 2.1 * inch

        _draw_molecule(canvas, lx, ly + 0.3 * inch, 0.42 * inch)

        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 38)
        canvas.drawString(lx + 0.62 * inch, ly + 0.28 * inch, "IntelliForm")

        canvas.setFont("Helvetica", 9)
        canvas.setFillColor(TEAL)
        canvas.drawString(lx + 0.64 * inch, ly + 0.02 * inch,
                          "AGENTIC GREEN CHEMISTRY PLATFORM")

        # ── Divider ───────────────────────────────────────────────────────────
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(1.5)
        canvas.line(lx, ly - 0.26 * inch, PW - 0.5 * inch, ly - 0.26 * inch)

        # ── Proposal title ────────────────────────────────────────────────────
        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 28)
        canvas.drawString(lx, PH - 3.35 * inch, "Green Chemistry")
        canvas.drawString(lx, PH - 3.85 * inch, "Formulation Proposal")

        app_label = project.get("application", "").replace("_", " ").title()
        canvas.setFont("Helvetica", 14)
        canvas.setFillColor(AMBER)
        canvas.drawString(lx, PH - 4.35 * inch, app_label)

        # ── "Prepared for" card ───────────────────────────────────────────────
        card_y = PH / 2 - 0.8 * inch
        card_h = 1.45 * inch
        card_x = lx - 0.1 * inch
        card_w = PW - card_x - 0.5 * inch

        canvas.setFillColor(NAVY2)
        canvas.roundRect(card_x, card_y, card_w, card_h, 6, fill=1, stroke=0)

        # Teal left pip on card
        canvas.setFillColor(TEAL)
        canvas.roundRect(card_x, card_y, 0.06 * inch, card_h, 3, fill=1, stroke=0)

        tx = card_x + 0.2 * inch
        canvas.setFillColor(TEAL)
        canvas.setFont("Helvetica-Bold", 7.5)
        canvas.drawString(tx, card_y + card_h - 0.28 * inch, "PREPARED FOR")

        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 17)
        canvas.drawString(tx, card_y + card_h - 0.6 * inch,
                          customer_name or "Valued Customer")

        if customer_company:
            canvas.setFont("Helvetica", 11)
            canvas.setFillColor(MUTED)
            canvas.drawString(tx, card_y + card_h - 0.88 * inch, customer_company)

        # Hairline separator inside card
        canvas.setStrokeColor(colors.HexColor("#2D5A8E"))
        canvas.setLineWidth(0.5)
        sep_y = card_y + 0.48 * inch
        canvas.line(tx, sep_y, card_x + card_w - 0.2 * inch, sep_y)

        canvas.setFont("Helvetica", 8.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(tx, card_y + 0.22 * inch,
                          "Prepared by: ChemeNova LLC × ChemRich Global")

        # ── Customer logo (white-label) ───────────────────────────────────────
        if logo_bytes:
            try:
                reader = ImageReader(io.BytesIO(logo_bytes))
                iw, ih = reader.getSize()
                max_w, max_h = 1.6 * inch, 0.55 * inch
                scale = min(max_w / iw, max_h / ih)
                dw, dh = iw * scale, ih * scale
                canvas.drawImage(reader,
                                 PW - 0.5 * inch - dw,
                                 card_y + (card_h - dh) / 2,
                                 width=dw, height=dh, mask="auto")
            except Exception:
                pass

        # ── Bottom strip ──────────────────────────────────────────────────────
        canvas.setFillColor(colors.HexColor("#0A1E35"))
        canvas.rect(0, 0, PW, 0.9 * inch, fill=1, stroke=0)

        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 8)
        canvas.drawString(lx, 0.58 * inch, date_str)

        canvas.setFillColor(MUTED)
        canvas.setFont("Helvetica", 7)
        canvas.drawString(lx, 0.38 * inch,
                          "CONFIDENTIAL — proprietary formulation data, not for distribution")

        canvas.setFillColor(AMBER)
        canvas.setFont("Helvetica-Bold", 8)
        canvas.drawRightString(PW - 0.5 * inch, 0.58 * inch, "intelliform.chemenova.com")

        canvas.restoreState()

    return _draw


# ── Running header / footer ───────────────────────────────────────────────────

def _make_header_footer(date_str: str):
    def _draw(canvas, doc):
        canvas.saveState()

        # Header bar
        canvas.setFillColor(NAVY)
        canvas.rect(0, PH - 0.52 * inch, PW, 0.52 * inch, fill=1, stroke=0)

        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 9)
        canvas.drawString(0.5 * inch, PH - 0.33 * inch,
                          "Green Chemistry Formulation Proposal")
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(TEAL)
        canvas.drawRightString(PW - 0.5 * inch, PH - 0.33 * inch,
                               f"ChemeNova LLC × ChemRich Global  |  {date_str}")

        # Teal rule below header
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(1.5)
        canvas.line(0.5 * inch, PH - 0.55 * inch, PW - 0.5 * inch, PH - 0.55 * inch)

        # Footer rule
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(0.5 * inch, 0.5 * inch, PW - 0.5 * inch, 0.5 * inch)

        canvas.setFillColor(SLATE)
        canvas.setFont("Helvetica", 7)
        canvas.drawString(0.5 * inch, 0.33 * inch,
                          "CONFIDENTIAL — ChemeNova LLC × ChemRich Global")
        canvas.drawRightString(PW - 0.5 * inch, 0.33 * inch,
                               f"Page {doc.page}")
        canvas.drawCentredString(PW / 2, 0.33 * inch, "shehan@chemenova.com")

        canvas.restoreState()

    return _draw


# ── Metric boxes ──────────────────────────────────────────────────────────────

def _metric_table(metrics: list) -> Table:
    """metrics = [(label, value, unit), ...]"""
    s = _styles()
    cells = [[
        Paragraph(str(v), s["metric_value"]),
        Paragraph(u, s["metric_unit"]),
        Paragraph(lbl, s["metric_label"]),
    ] for lbl, v, u in metrics]

    col_w = (PW - inch) / len(metrics)
    t = Table([cells], colWidths=[col_w] * len(metrics))
    t.setStyle(TableStyle([
        ("BOX",           (0, 0), (-1, -1), 0.5, RULE),
        ("LINEAFTER",     (0, 0), (-2, -1), 0.5, RULE),
        ("BACKGROUND",    (0, 0), (-1, -1), LIGHT),
        ("TOPPADDING",    (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
        ("ALIGN",         (0, 0), (-1, -1), "CENTER"),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ]))
    return t


# ── Ingredient row ────────────────────────────────────────────────────────────

def _ingredient_row(name: str, pct: float, db_row: dict) -> Table:
    s = _styles()
    func = str(db_row.get("Function", "—")) if db_row else "—"
    cost = str(db_row.get("Cost_USD_kg", "—")) if db_row else "—"
    bio  = str(db_row.get("Bio_based_pct", "—")) if db_row else "—"

    data = [[
        Paragraph(f"<b>{name}</b>", s["body"]),
        Paragraph(f"{pct:.1f}%", s["body"]),
        Paragraph(func, s["small"]),
        Paragraph(f"${cost}/kg", s["small"]),
        Paragraph(f"{bio}% bio", s["small"]),
    ]]
    col_w = [2.5 * inch, 0.6 * inch, 1.4 * inch, 0.8 * inch, 0.7 * inch]
    t = Table(data, colWidths=col_w)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, -1), LIGHT),
        ("BOX",           (0, 0), (-1, -1), 0.5, RULE),
        ("TOPPADDING",    (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING",   (0, 0), (-1, -1), 6),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
        ("BACKGROUND",    (1, 0), (1, 0),   TEAL),
        ("TEXTCOLOR",     (1, 0), (1, 0),   WHITE),
        ("FONTNAME",      (1, 0), (1, 0),   "Helvetica-Bold"),
        ("ALIGN",         (1, 0), (1, 0),   "CENTER"),
    ]))
    return t


# ── Radar chart ───────────────────────────────────────────────────────────────

def _radar_drawing(eco_result) -> Optional[Drawing]:
    if eco_result is None:
        return None
    try:
        d = Drawing(280, 210)
        sp = SpiderChart()
        sp.x, sp.y, sp.width, sp.height = 10, 10, 260, 190
        sp.labels = ["Biodegradability", "Carbon\nFootprint",
                     "Ecotoxicity", "Renewability", "Regulatory"]
        sp.data = [
            [eco_result.biodegradability, eco_result.carbon_footprint,
             eco_result.ecotoxicity, eco_result.renewability, eco_result.regulatory],
            [52, 38, 41, 25, 60],
        ]
        sp.lines[0].strokeColor = GREEN
        sp.lines[0].fillColor   = colors.Color(0, 0.78, 0.32, 0.18)
        sp.lines[1].strokeColor = RED_C
        sp.lines[1].fillColor   = colors.Color(0.94, 0.27, 0.27, 0.10)
        sp.strands[0].strokeColor = GREEN
        sp.strands[1].strokeColor = RED_C
        sp.labelRadius = 1.2
        sp.spokes.strokeColor = colors.HexColor("#CBD5E1")
        d.add(sp)
        return d
    except Exception:
        return None


# ── Regulatory table ──────────────────────────────────────────────────────────

def _reg_table(reg_report) -> Table:
    s = _styles()
    flag_colors = {
        "Green": colors.HexColor("#D1FAE5"),
        "Amber": colors.HexColor("#FEF3C7"),
        "Red":   colors.HexColor("#FEE2E2"),
    }
    rows = [["Ingredient", "REACH", "EPA SC", "EU Ecolabel", "COSMOS", "Restrictions"]]
    row_styles = []
    for i, (name, p) in enumerate(reg_report.profiles.items(), start=1):
        restr = (p.restrictions[0][:55] + "…") if p.restrictions else "None"
        rows.append([
            name, p.reach_status, p.epa_safer_choice,
            "Yes" if p.eu_ecolabel else "No",
            "Yes" if p.cosmos_approved else "No",
            restr,
        ])
        row_styles.append(("BACKGROUND", (0, i), (-1, i),
                           flag_colors.get(p.reach_flag, WHITE)))

    col_w = [1.8 * inch, 0.9 * inch, 0.6 * inch, 0.85 * inch, 0.7 * inch, 2.15 * inch]
    t = Table(rows, colWidths=col_w, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR",     (0, 0), (-1, 0), WHITE),
        ("FONTNAME",      (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",      (0, 0), (-1, -1), 7.5),
        ("GRID",          (0, 0), (-1, -1), 0.3, RULE),
        ("TOPPADDING",    (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING",   (0, 0), (-1, -1), 4),
        ("ALIGN",         (2, 0), (-1, -1), "CENTER"),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ] + row_styles))
    return t


# ── Main entry point ──────────────────────────────────────────────────────────

def generate_proposal_pdf(
    project: Dict,
    eco_result=None,
    reg_report=None,
    db=None,
    customer_name: Optional[str] = None,
    customer_company: Optional[str] = None,
    logo_bytes: Optional[bytes] = None,
) -> bytes:
    """
    Generate a branded ChemeNova PDF proposal.

    Args:
        project:          project dict from session_state.projects
        eco_result:       EcoMetricsResult (optional)
        reg_report:       BlendRegulatoryReport (optional)
        db:               ingredients DataFrame (optional)
        customer_name:    recipient name (for cover page)
        customer_company: recipient company (for cover page)
        logo_bytes:       PNG/JPG bytes of customer logo for white-labelling

    Returns:
        PDF bytes (ready for st.download_button)
    """
    buf      = io.BytesIO()
    date_str = datetime.now().strftime("%B %d, %Y")
    s        = _styles()

    doc = SimpleDocTemplate(
        buf,
        pagesize=letter,
        leftMargin=0.5 * inch,
        rightMargin=0.5 * inch,
        topMargin=0.72 * inch,
        bottomMargin=0.65 * inch,
    )

    story = [PageBreak()]   # page 1 is the cover (canvas-drawn); content starts p.2

    # ── Executive summary ─────────────────────────────────────────────────────
    eco_score_val = f"{eco_result.eco_score:.0f}" if eco_result else "—"
    eco_grade_val = eco_result.grade if eco_result else "—"

    story.append(Paragraph("Executive Summary", s["section_head"]))
    story.append(_metric_table([
        ("Cost / kg",   f"${project.get('cost', 0):.2f}", "USD"),
        ("Bio-based",   f"{project.get('bio', 0):.1f}",   "%"),
        ("Performance", f"{project.get('perf', 0):.0f}",  "/ 100"),
        ("EcoScore™",   eco_score_val,                    f"Grade {eco_grade_val}"),
    ]))
    story.append(Spacer(1, 0.12 * inch))

    # Application + optimizer meta
    app_label = project.get("application", "unknown").replace("_", " ").title()
    story.append(Paragraph(
        f"Application: <b>{app_label}</b>  ·  "
        f"Optimizer: <b>{project.get('optimizer','pulp').upper()}</b>  ·  "
        f"Parser: <b>{project.get('parser','regex').upper()}</b>",
        s["body"],
    ))
    story.append(Spacer(1, 0.06 * inch))

    # ── Blend composition ─────────────────────────────────────────────────────
    story.append(HRFlowable(width="100%", thickness=1.5, color=TEAL, spaceAfter=6))
    story.append(Paragraph("Optimized Blend Composition", s["section_head"]))

    blend = project.get("blend", {})
    if db is not None:
        try:
            import pandas as pd
            idx = db.set_index("Ingredient") if isinstance(db, pd.DataFrame) else None
        except ImportError:
            idx = None
        for ing, pct in blend.items():
            row = (idx.loc[ing].to_dict()
                   if idx is not None and ing in idx.index else {})
            story.append(_ingredient_row(ing, pct, row))
            story.append(Spacer(1, 3))
    else:
        for ing, pct in blend.items():
            story.append(Paragraph(f"- <b>{ing}</b>  {pct:.1f}%", s["body"]))

    if project.get("relaxed"):
        story.append(Spacer(1, 6))
        story.append(Paragraph(
            "NOTE: Constraints were auto-relaxed to find this blend. "
            "Verify metrics against specification before filing.",
            s["warn"],
        ))

    # ── EcoMetrics page ───────────────────────────────────────────────────────
    if eco_result:
        story.append(PageBreak())
        story.append(Paragraph("EcoMetrics™ Sustainability Profile", s["section_head"]))
        story.append(Paragraph(
            "Five-axis scoring vs. a representative petrochemical surfactant baseline. "
            "Methods: OECD 301B (biodegradability) · kgCO₂eq/kg inverted (carbon footprint) · "
            "ECHA aquatic rating inverted (ecotoxicity) · ASTM D6866 composite (renewability) · "
            "REACH/EPA/EU Ecolabel composite (regulatory). "
            "Published in IntelliForm JCIM Supporting Information (2026).",
            s["small"],
        ))
        story.append(Spacer(1, 8))

        radar = _radar_drawing(eco_result)
        if radar:
            story.append(radar)
            story.append(Spacer(1, 8))

        axis_data = [
            ["Axis", "IntelliForm", "Baseline", "Δ Improvement", "Weight"],
            ["Biodegradability",
             f"{eco_result.biodegradability:.1f}", "52.0",
             f"+{eco_result.vs_baseline.get('Biodegradability', 0):.1f}", "25%"],
            ["Carbon Footprint",
             f"{eco_result.carbon_footprint:.1f}", "38.0",
             f"+{eco_result.vs_baseline.get('Carbon Footprint', 0):.1f}", "20%"],
            ["Ecotoxicity",
             f"{eco_result.ecotoxicity:.1f}", "41.0",
             f"+{eco_result.vs_baseline.get('Ecotoxicity', 0):.1f}", "20%"],
            ["Renewability",
             f"{eco_result.renewability:.1f}", "25.0",
             f"+{eco_result.vs_baseline.get('Renewability', 0):.1f}", "20%"],
            ["Regulatory",
             f"{eco_result.regulatory:.1f}", "60.0",
             f"+{eco_result.vs_baseline.get('Regulatory', 0):.1f}", "15%"],
            ["EcoScore™ (composite)",
             f"{eco_result.eco_score:.1f}", "43.2",
             f"+{eco_result.eco_score - 43.2:.1f}", "—"],
        ]
        eco_t = Table(
            axis_data,
            colWidths=[1.85 * inch, 1.1 * inch, 1.0 * inch, 1.3 * inch, 0.75 * inch],
            repeatRows=1,
        )
        eco_t.setStyle(TableStyle([
            ("BACKGROUND",    (0, 0), (-1, 0),  NAVY),
            ("TEXTCOLOR",     (0, 0), (-1, 0),  WHITE),
            ("FONTNAME",      (0, 0), (-1, 0),  "Helvetica-Bold"),
            ("FONTSIZE",      (0, 0), (-1, -1), 8),
            ("GRID",          (0, 0), (-1, -1), 0.3, RULE),
            ("TOPPADDING",    (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("ALIGN",         (1, 0), (-1, -1), "CENTER"),
            ("FONTNAME",      (0, 6), (-1, 6),  "Helvetica-Bold"),
            ("BACKGROUND",    (0, 6), (-1, 6),  colors.HexColor("#DCFCE7")),
        ]))
        story.append(eco_t)

    # ── Regulatory page ───────────────────────────────────────────────────────
    if reg_report:
        story.append(PageBreak())
        story.append(Paragraph("Regulatory Intelligence", s["section_head"]))
        story.append(Paragraph(
            f"Overall status: <b>{reg_report.overall_status}</b>  ·  "
            f"REACH: {'All Green' if reg_report.all_green else 'Review Required'}  ·  "
            f"EU Ecolabel eligible: {'Yes' if reg_report.eu_ecolabel_eligible else 'No'}  ·  "
            f"COSMOS eligible: {'Yes' if reg_report.cosmos_eligible else 'No'}  ·  "
            f"EPA Safer Choice: {'Yes' if reg_report.epa_safer_choice_eligible else 'No'}",
            s["body"],
        ))
        story.append(Spacer(1, 8))
        story.append(_reg_table(reg_report))

        if reg_report.certification_pathways:
            story.append(Spacer(1, 10))
            story.append(Paragraph("Eligible Certification Pathways", s["section_head"]))
            for p in reg_report.certification_pathways:
                story.append(Paragraph(f"· {p}", s["body"]))

        if reg_report.amber_flags:
            story.append(Spacer(1, 8))
            story.append(Paragraph("Items Requiring Review", s["section_head"]))
            for flag in reg_report.amber_flags:
                story.append(Paragraph(f"⚠ {flag}", s["body"]))

    # ── Next steps + financial ────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("Next Steps", s["section_head"]))

    next_data = [
        ["Step", "Action", "Timeline"],
        ["1", "Review this proposal with your formulation team", "This week"],
        ["2", "Confirm ingredient stock availability with ChemRich NJ", "2–3 days"],
        ["3", "Book 500 kg pilot batch", "5 business days turnaround"],
        ["4", "Receive pilot sample + QC report", "Day 10"],
        ["5", "Scale to commercial production", "30–60 days"],
    ]
    ns_t = Table(next_data, colWidths=[0.5 * inch, 4.1 * inch, 1.4 * inch], repeatRows=1)
    ns_t.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, 0),  TEAL),
        ("TEXTCOLOR",     (0, 0), (-1, 0),  WHITE),
        ("FONTNAME",      (0, 0), (-1, 0),  "Helvetica-Bold"),
        ("FONTSIZE",      (0, 0), (-1, -1), 8.5),
        ("GRID",          (0, 0), (-1, -1), 0.3, RULE),
        ("TOPPADDING",    (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("ALIGN",         (0, 0), (0, -1),  "CENTER"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [LIGHT, WHITE]),
    ]))
    story.append(ns_t)
    story.append(Spacer(1, 14))

    savings = project.get("savings", 0)
    co2     = project.get("co2_kg", 0)
    quote   = round(project.get("cost", 0) * 500 * 1.12, 0)

    fin_data = [
        ["Financial Summary", ""],
        ["Pilot quote (500 kg + 12% ChemRich fee)", f"${quote:,.0f}"],
        ["Projected savings vs. market",            f"${savings:,.0f} / 500 kg batch"],
        ["CO₂ avoided",                             f"{co2:,.0f} kg / batch"],
        ["Estimated payback period",                "< 1 production run"],
    ]
    fin_t = Table(fin_data, colWidths=[3.5 * inch, 2.5 * inch])
    fin_t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0),  AMBER),
        ("TEXTCOLOR",  (0, 0), (-1, 0),  WHITE),
        ("FONTNAME",   (0, 0), (-1, 0),  "Helvetica-Bold"),
        ("SPAN",       (0, 0), (-1, 0)),
        ("FONTSIZE",   (0, 0), (-1, -1), 9),
        ("GRID",       (0, 0), (-1, -1), 0.3, RULE),
        ("TOPPADDING",    (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING",   (0, 0), (-1, -1), 8),
        ("FONTNAME",   (0, 1), (0, -1),  "Helvetica"),
        ("FONTNAME",   (1, 1), (-1, -1), "Helvetica-Bold"),
        ("TEXTCOLOR",  (1, 1), (-1, -1), NAVY),
        ("BACKGROUND", (0, 1), (-1, -1), LIGHT),
    ]))
    story.append(fin_t)
    story.append(Spacer(1, 18))

    # ── Contact + closing ─────────────────────────────────────────────────────
    story.append(HRFlowable(width="100%", thickness=1, color=TEAL, spaceAfter=8))
    story.append(Paragraph(
        "<b>Contact:</b>  Shehan Makani, ChemeNova LLC  ·  "
        "shehan@chemenova.com  ·  intelliform.chemenova.com  ·  Pearl River, NJ",
        s["body"],
    ))
    story.append(Spacer(1, 5))
    story.append(Paragraph(
        "<i>Generated by IntelliForm AI Platform v0.9 (ChemeNova LLC). "
        "Formulation metrics computed using validated QSAR models (JCIM, 2026) "
        "and real-time ChemRich inventory data. Regulatory data sourced from ECHA, "
        "EPA Safer Choice, and COSMOS-standard databases as of January 2026. "
        "This document is confidential and prepared exclusively for the named recipient.</i>",
        s["small"],
    ))

    # ── Build ─────────────────────────────────────────────────────────────────
    cover_fn   = _make_cover(project, customer_name or "Valued Customer",
                              customer_company or "", logo_bytes, date_str)
    content_fn = _make_header_footer(date_str)

    doc.build(story, onFirstPage=cover_fn, onLaterPages=content_fn)
    return buf.getvalue()
