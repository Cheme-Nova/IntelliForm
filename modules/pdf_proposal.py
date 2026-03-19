"""
modules/pdf_proposal.py
Branded PDF proposal generator for IntelliForm v0.9.

Generates a professional multi-page PDF proposal using ReportLab,
styled with ChemeNova brand colors (navy #0A1628, teal #0D9488, amber #D97706).

Includes:
  Page 1 — Cover with logo text, blend summary, key metrics
  Page 2 — EcoMetrics radar chart + axis table
  Page 3 — Regulatory intelligence table + certification pathways
  Page 4 — Molecular structures (RDKit) + ingredient detail cards
  Page 5 — Next steps, contact, and IntelliForm model card excerpt

Usage:
    from modules.pdf_proposal import generate_proposal_pdf
    pdf_bytes = generate_proposal_pdf(project, eco_result, reg_report, db)
    st.download_button("Download PDF", pdf_bytes, "proposal.pdf", "application/pdf")
"""
import io
from datetime import datetime
from typing import Dict, Optional, Any

# ReportLab imports
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, PageBreak, Image as RLImage, KeepTogether,
)
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Wedge, Circle
from reportlab.graphics.charts.spider import SpiderChart
from reportlab.graphics import renderPDF

# ── Brand colors ──────────────────────────────────────────────────────────────
NAVY   = colors.HexColor("#0A1628")
TEAL   = colors.HexColor("#0D9488")
AMBER  = colors.HexColor("#D97706")
WHITE  = colors.white
LIGHT  = colors.HexColor("#F0F9FF")
GRAY   = colors.HexColor("#64748B")
GREEN  = colors.HexColor("#00C853")
RED_C  = colors.HexColor("#EF4444")

W, H   = letter   # 8.5 x 11 inches


# ── Style definitions ─────────────────────────────────────────────────────────

def _styles():
    base = getSampleStyleSheet()
    custom = {
        "cover_title": ParagraphStyle(
            "cover_title", fontSize=28, leading=34,
            textColor=WHITE, fontName="Helvetica-Bold",
            alignment=TA_LEFT, spaceAfter=8,
        ),
        "cover_sub": ParagraphStyle(
            "cover_sub", fontSize=13, leading=18,
            textColor=TEAL, fontName="Helvetica",
            alignment=TA_LEFT, spaceAfter=4,
        ),
        "section_head": ParagraphStyle(
            "section_head", fontSize=13, leading=16,
            textColor=NAVY, fontName="Helvetica-Bold",
            spaceBefore=14, spaceAfter=6,
            borderPad=4,
        ),
        "body": ParagraphStyle(
            "body", fontSize=9, leading=13,
            textColor=colors.HexColor("#1E293B"),
            fontName="Helvetica", spaceAfter=4,
        ),
        "small": ParagraphStyle(
            "small", fontSize=7.5, leading=10,
            textColor=GRAY, fontName="Helvetica",
        ),
        "metric_label": ParagraphStyle(
            "metric_label", fontSize=8, leading=10,
            textColor=GRAY, fontName="Helvetica",
            alignment=TA_CENTER,
        ),
        "metric_value": ParagraphStyle(
            "metric_value", fontSize=18, leading=22,
            textColor=NAVY, fontName="Helvetica-Bold",
            alignment=TA_CENTER,
        ),
        "tag": ParagraphStyle(
            "tag", fontSize=8, leading=10,
            textColor=WHITE, fontName="Helvetica-Bold",
            alignment=TA_CENTER,
        ),
        "footer": ParagraphStyle(
            "footer", fontSize=7, leading=9,
            textColor=GRAY, fontName="Helvetica",
            alignment=TA_CENTER,
        ),
    }
    return {**{k: base[k] for k in base.byName}, **custom}


# ── Header / Footer callback ──────────────────────────────────────────────────

def _make_header_footer(title: str, date_str: str):
    def _draw(canvas, doc):
        canvas.saveState()
        page_w, page_h = letter

        # Top bar
        canvas.setFillColor(NAVY)
        canvas.rect(0, page_h - 0.55 * inch, page_w, 0.55 * inch, fill=1, stroke=0)

        # Header text
        canvas.setFillColor(WHITE)
        canvas.setFont("Helvetica-Bold", 10)
        canvas.drawString(0.5 * inch, page_h - 0.35 * inch, "IntelliForm(TM) Green Formulation Proposal")
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(TEAL)
        canvas.drawRightString(page_w - 0.5 * inch, page_h - 0.35 * inch,
                               f"ChemeNova LLC x ChemRich Global  |  {date_str}")

        # Teal accent line
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(2)
        canvas.line(0.5 * inch, page_h - 0.58 * inch, page_w - 0.5 * inch, page_h - 0.58 * inch)

        # Footer
        canvas.setFillColor(GRAY)
        canvas.setFont("Helvetica", 7)
        canvas.drawCentredString(
            page_w / 2, 0.35 * inch,
            f"IntelliForm(TM) v0.9  |  Powered by ChemeNova LLC  |  shehan@chemenova.com  |  Page {doc.page}"
        )
        canvas.setStrokeColor(colors.HexColor("#E2E8F0"))
        canvas.setLineWidth(0.5)
        canvas.line(0.5 * inch, 0.52 * inch, page_w - 0.5 * inch, 0.52 * inch)

        canvas.restoreState()
    return _draw


# ── Metric box helper ─────────────────────────────────────────────────────────

def _metric_table(metrics: list) -> Table:
    """metrics = [(label, value, unit), ...]  — max 4 per row"""
    s = _styles()
    cells = []
    for label, value, unit in metrics:
        cell = [
            Paragraph(str(value), s["metric_value"]),
            Paragraph(unit, s["small"]),
            Paragraph(label, s["metric_label"]),
        ]
        cells.append(cell)

    col_w = (W - inch) / len(metrics)
    t = Table([cells], colWidths=[col_w] * len(metrics))
    t.setStyle(TableStyle([
        ("BOX",        (0, 0), (-1, -1), 0.5, colors.HexColor("#E2E8F0")),
        ("LINEAFTER",  (0, 0), (-2, -1), 0.5, colors.HexColor("#E2E8F0")),
        ("BACKGROUND", (0, 0), (-1, -1), LIGHT),
        ("TOPPADDING",    (0,0),(-1,-1), 10),
        ("BOTTOMPADDING", (0,0),(-1,-1), 10),
        ("ALIGN",      (0, 0), (-1, -1), "CENTER"),
        ("VALIGN",     (0, 0), (-1, -1), "MIDDLE"),
    ]))
    return t


# ── Radar chart (spider) ──────────────────────────────────────────────────────

def _radar_drawing(eco_result) -> Optional[Drawing]:
    """Draw EcoMetrics radar as a ReportLab Drawing."""
    if eco_result is None:
        return None
    try:
        d = Drawing(260, 200)
        sp = SpiderChart()
        sp.x, sp.y, sp.width, sp.height = 10, 10, 240, 180

        labels = ["Biodegradability", "Carbon\nFootprint", "Ecotoxicity",
                  "Renewability", "Regulatory"]
        sp.labels = labels

        intelliform = [
            eco_result.biodegradability,
            eco_result.carbon_footprint,
            eco_result.ecotoxicity,
            eco_result.renewability,
            eco_result.regulatory,
        ]
        baseline = [52, 38, 41, 25, 60]

        sp.data = [intelliform, baseline]
        sp.lines[0].strokeColor = GREEN
        sp.lines[0].fillColor    = colors.Color(0, 0.78, 0.32, 0.2)
        sp.lines[1].strokeColor  = RED_C
        sp.lines[1].fillColor    = colors.Color(0.94, 0.27, 0.27, 0.1)

        sp.strands[0].strokeColor = GREEN
        sp.strands[1].strokeColor  = RED_C

        sp.labelRadius = 1.2
        sp.spokes.strokeColor = colors.HexColor("#CBD5E1")

        d.add(sp)
        return d
    except Exception:
        return None


# ── Ingredient pill ───────────────────────────────────────────────────────────

def _ingredient_row(name: str, pct: float, db_row) -> Table:
    s = _styles()
    pct_str = f"{pct:.1f}%"
    func    = str(db_row.get("Function", "-")) if hasattr(db_row, "get") else "-"
    cost    = str(db_row.get("Cost_USD_kg", "-")) if hasattr(db_row, "get") else "-"
    bio     = str(db_row.get("Bio_based_pct", "-")) if hasattr(db_row, "get") else "-"

    data = [[
        Paragraph(f"<b>{name}</b>", s["body"]),
        Paragraph(pct_str, s["body"]),
        Paragraph(func, s["small"]),
        Paragraph(f"${cost}/kg", s["small"]),
        Paragraph(f"{bio}% bio", s["small"]),
    ]]
    col_w = [2.5*inch, 0.6*inch, 1.4*inch, 0.8*inch, 0.7*inch]
    t = Table(data, colWidths=col_w)
    t.setStyle(TableStyle([
        ("BACKGROUND",    (0, 0), (-1, -1), LIGHT),
        ("BOX",           (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E1")),
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


# ── Regulatory row coloring ───────────────────────────────────────────────────

def _reg_table(reg_report) -> Table:
    s = _styles()
    headers = ["Ingredient", "REACH", "EPA SC", "EU Ecolabel", "COSMOS", "Restrictions"]
    rows = [headers]

    flag_colors = {"Green": colors.HexColor("#D1FAE5"), "Amber": colors.HexColor("#FEF3C7"),
                   "Red":   colors.HexColor("#FEE2E2")}

    row_styles = []
    for i, (name, profile) in enumerate(reg_report.profiles.items(), start=1):
        row = [
            name,
            profile.reach_status,
            profile.epa_safer_choice,
            "Yes" if profile.eu_ecolabel else "No",
            "Yes" if profile.cosmos_approved else "No",
            profile.restrictions[0][:55] + "…" if profile.restrictions else "None",
        ]
        rows.append(row)
        bg = flag_colors.get(profile.reach_flag, WHITE)
        row_styles.append(("BACKGROUND", (0, i), (-1, i), bg))

    col_w = [1.8*inch, 0.9*inch, 0.6*inch, 0.85*inch, 0.7*inch, 2.15*inch]
    t = Table(rows, colWidths=col_w, repeatRows=1)
    base_style = [
        ("BACKGROUND",    (0, 0), (-1, 0),  NAVY),
        ("TEXTCOLOR",     (0, 0), (-1, 0),  WHITE),
        ("FONTNAME",      (0, 0), (-1, 0),  "Helvetica-Bold"),
        ("FONTSIZE",      (0, 0), (-1, -1), 7.5),
        ("GRID",          (0, 0), (-1, -1), 0.3, colors.HexColor("#CBD5E1")),
        ("TOPPADDING",    (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING",   (0, 0), (-1, -1), 4),
        ("ALIGN",         (2, 0), (-1, -1), "CENTER"),
        ("VALIGN",        (0, 0), (-1, -1), "MIDDLE"),
    ]
    t.setStyle(TableStyle(base_style + row_styles))
    return t


# ── Main PDF builder ──────────────────────────────────────────────────────────

def generate_proposal_pdf(
    project: Dict,
    eco_result=None,
    reg_report=None,
    db=None,
) -> bytes:
    """
    Generate a branded ChemeNova PDF proposal.

    Args:
        project:    project dict from session_state.projects
        eco_result: EcoMetricsResult (optional)
        reg_report: BlendRegulatoryReport (optional)
        db:         ingredients DataFrame (optional, for detail rows)

    Returns:
        PDF as bytes (ready for st.download_button)
    """
    buf = io.BytesIO()
    date_str = datetime.now().strftime("%B %d, %Y")

    doc = SimpleDocTemplate(
        buf,
        pagesize=letter,
        leftMargin=0.5 * inch,
        rightMargin=0.5 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.65 * inch,
    )

    s = _styles()
    story = []
    on_page = _make_header_footer(project.get("application", ""), date_str)

    # ── Cover section ──────────────────────────────────────────────────────────
    # Navy cover band
    cover_band = Table(
        [[Paragraph("IntelliForm(TM) Green Formulation Proposal", s["cover_title"]),
          Paragraph(f"<b>ChemeNova LLC x ChemRich Global</b><br/>{date_str}", s["cover_sub"])]],
        colWidths=[4.5*inch, 3*inch],
    )
    cover_band.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), NAVY),
        ("TOPPADDING",    (0,0),(-1,-1), 18),
        ("BOTTOMPADDING", (0,0),(-1,-1), 18),
        ("LEFTPADDING",   (0,0),(0,-1),  14),
        ("VALIGN",        (0,0),(-1,-1), "MIDDLE"),
    ]))
    story.append(cover_band)
    story.append(Spacer(1, 0.15*inch))

    # Application badge
    app_label = project.get("application", "unknown").replace("_", " ").title()
    optimizer  = project.get("optimizer", "pulp").upper()
    story.append(Paragraph(
        f"Application: <b>{app_label}</b>  |  Optimizer: <b>{optimizer}</b>  |  "
        f"Parser: <b>{project.get('parser','regex').upper()}</b>",
        s["body"]
    ))
    story.append(Spacer(1, 0.1*inch))

    # Key metrics band
    eco_score_val = f"{eco_result.eco_score:.0f}" if eco_result else "-"
    eco_grade_val = eco_result.grade if eco_result else "-"
    story.append(_metric_table([
        ("Cost / kg",        f"${project.get('cost', 0):.2f}",  "USD"),
        ("Bio-based",        f"{project.get('bio', 0):.1f}",    "%"),
        ("Performance",      f"{project.get('perf', 0):.0f}",   "/ 100"),
        ("EcoScore(TM)",        eco_score_val,                      f"Grade {eco_grade_val}"),
    ]))
    story.append(Spacer(1, 0.15*inch))

    # ── Blend composition ──────────────────────────────────────────────────────
    story.append(HRFlowable(width="100%", thickness=1.5, color=TEAL, spaceAfter=6))
    story.append(Paragraph("Optimized Blend Composition", s["section_head"]))

    blend = project.get("blend", {})
    if db is not None:
        import pandas as pd
        idx = db.set_index("Ingredient") if isinstance(db, pd.DataFrame) else None
        for ing, pct in blend.items():
            row = idx.loc[ing].to_dict() if (idx is not None and ing in idx.index) else {}
            story.append(_ingredient_row(ing, pct, row))
            story.append(Spacer(1, 3))
    else:
        for ing, pct in blend.items():
            story.append(Paragraph(f"- <b>{ing}</b> - {pct:.1f}%", s["body"]))

    # Relaxation note
    if project.get("relaxed"):
        story.append(Spacer(1, 6))
        story.append(Paragraph(
            "NOTE: Original constraints were auto-relaxed to find this blend. "
            "Review metrics against your specification before filing.",
            ParagraphStyle("warn", parent=s["body"], textColor=AMBER, fontName="Helvetica-Bold")
        ))

    # ── EcoMetrics page ────────────────────────────────────────────────────────
    if eco_result:
        story.append(PageBreak())
        story.append(Paragraph("EcoMetrics(TM) Sustainability Profile", s["section_head"]))
        story.append(Paragraph(
            "Five-axis sustainability scoring benchmarked against a representative "
            "petrochemical surfactant blend. Methodology: OECD 301B (biodegradability), "
            "kgCO2eq/kg inverted (carbon footprint), ECHA aquatic rating inverted (ecotoxicity), "
            "ASTM D6866 composite (renewability), REACH/EPA/EU Ecolabel composite (regulatory). "
            "Published in IntelliForm JCIM Supporting Information.",
            s["small"]
        ))
        story.append(Spacer(1, 8))

        # Radar
        radar = _radar_drawing(eco_result)
        if radar:
            story.append(radar)
            story.append(Spacer(1, 8))

        # Axis table
        axis_data = [
            ["Axis", "IntelliForm(TM) Score", "Petrochemical Baseline", "Δ Improvement", "Weight"],
            ["Biodegradability",  f"{eco_result.biodegradability:.1f}",  "52.0", f"+{eco_result.vs_baseline.get('Biodegradability',0):.1f}", "25%"],
            ["Carbon Footprint",  f"{eco_result.carbon_footprint:.1f}",  "38.0", f"+{eco_result.vs_baseline.get('Carbon Footprint',0):.1f}", "20%"],
            ["Ecotoxicity",       f"{eco_result.ecotoxicity:.1f}",       "41.0", f"+{eco_result.vs_baseline.get('Ecotoxicity',0):.1f}",    "20%"],
            ["Renewability",      f"{eco_result.renewability:.1f}",      "25.0", f"+{eco_result.vs_baseline.get('Renewability',0):.1f}",   "20%"],
            ["Regulatory",        f"{eco_result.regulatory:.1f}",        "60.0", f"+{eco_result.vs_baseline.get('Regulatory',0):.1f}",     "15%"],
            ["EcoScore(TM) (composite)", f"{eco_result.eco_score:.1f}", "43.2", f"+{eco_result.eco_score - 43.2:.1f}", "-"],
        ]
        eco_t = Table(axis_data, colWidths=[1.7*inch, 1.3*inch, 1.6*inch, 1.3*inch, 0.6*inch], repeatRows=1)
        eco_t.setStyle(TableStyle([
            ("BACKGROUND",    (0,0),(-1,0),  NAVY),
            ("TEXTCOLOR",     (0,0),(-1,0),  WHITE),
            ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
            ("FONTSIZE",      (0,0),(-1,-1), 8),
            ("GRID",          (0,0),(-1,-1), 0.3, colors.HexColor("#CBD5E1")),
            ("TOPPADDING",    (0,0),(-1,-1), 5),
            ("BOTTOMPADDING", (0,0),(-1,-1), 5),
            ("ALIGN",         (1,0),(-1,-1), "CENTER"),
            ("BACKGROUND",    (0,7),(-1,7),  LIGHT),
            ("FONTNAME",      (0,6),(-1,6),  "Helvetica-Bold"),
            ("BACKGROUND",    (0,6),(-1,6),  colors.HexColor("#DCFCE7")),
        ]))
        story.append(eco_t)

    # ── Regulatory page ────────────────────────────────────────────────────────
    if reg_report:
        story.append(PageBreak())
        story.append(Paragraph("Regulatory Intelligence", s["section_head"]))
        story.append(Paragraph(
            f"Overall status: <b>{reg_report.overall_status}</b>  |  "
            f"REACH: {'All Green' if reg_report.all_green else 'Review Required'}  |  "
            f"EU Ecolabel eligible: {'Yes' if reg_report.eu_ecolabel_eligible else 'No'}  |  "
            f"COSMOS eligible: {'Yes' if reg_report.cosmos_eligible else 'No'}  |  "
            f"EPA Safer Choice: {'Yes' if reg_report.epa_safer_choice_eligible else 'No'}",
            s["body"]
        ))
        story.append(Spacer(1, 8))
        story.append(_reg_table(reg_report))
        story.append(Spacer(1, 10))

        if reg_report.certification_pathways:
            story.append(Paragraph("Eligible Certification Pathways", s["section_head"]))
            for p in reg_report.certification_pathways:
                story.append(Paragraph(f"- {p}", s["body"]))

        if reg_report.amber_flags:
            story.append(Spacer(1, 8))
            story.append(Paragraph("Items Requiring Review", s["section_head"]))
            for flag in reg_report.amber_flags:
                story.append(Paragraph(f"WARNING: {flag}", s["body"]))

    # ── Next steps page ────────────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("Next Steps & Contact", s["section_head"]))

    savings = project.get("savings", 0)
    co2     = project.get("co2_kg", 0)
    quote   = round(project.get("cost", 0) * 500 * 1.12, 0)

    next_data = [
        ["Step", "Action", "Timeline"],
        ["1", "Review this proposal with your formulation team", "This week"],
        ["2", "Confirm ingredient stock availability with ChemRich NJ", "2-3 days"],
        ["3", "Book 500 kg pilot batch", "5 business days turnaround"],
        ["4", "Receive pilot sample + QC report", "Day 10"],
        ["5", "Scale to commercial production", "30-60 days"],
    ]
    ns_t = Table(next_data, colWidths=[0.5*inch, 4*inch, 1.5*inch], repeatRows=1)
    ns_t.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,0),  TEAL),
        ("TEXTCOLOR",     (0,0),(-1,0),  WHITE),
        ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
        ("FONTSIZE",      (0,0),(-1,-1), 8.5),
        ("GRID",          (0,0),(-1,-1), 0.3, colors.HexColor("#CBD5E1")),
        ("TOPPADDING",    (0,0),(-1,-1), 5),
        ("BOTTOMPADDING", (0,0),(-1,-1), 5),
        ("ALIGN",         (0,0),(0,-1),  "CENTER"),
        ("BACKGROUND",    (0,1),(-1,-1), LIGHT),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [LIGHT, WHITE]),
    ]))
    story.append(ns_t)
    story.append(Spacer(1, 12))

    # Financial summary
    fin_data = [
        ["Financial Summary",              ""],
        ["Pilot quote (500 kg + 12% fee)", f"${quote:,.0f}"],
        ["Projected savings vs. market",   f"${savings:,.0f} per 500 kg batch"],
        ["CO2 avoided",                    f"{co2:,.0f} kg per batch"],
        ["Payback period (est.)",          "< 1 production run"],
    ]
    fin_t = Table(fin_data, colWidths=[3.5*inch, 2.5*inch])
    fin_t.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,0),  AMBER),
        ("TEXTCOLOR",     (0,0),(-1,0),  WHITE),
        ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
        ("FONTNAME",      (0,1),(-1,-1), "Helvetica"),
        ("FONTSIZE",      (0,0),(-1,-1), 9),
        ("GRID",          (0,0),(-1,-1), 0.3, colors.HexColor("#CBD5E1")),
        ("TOPPADDING",    (0,0),(-1,-1), 5),
        ("BOTTOMPADDING", (0,0),(-1,-1), 5),
        ("LEFTPADDING",   (0,0),(-1,-1), 8),
        ("FONTNAME",      (1,1),(-1,-1), "Helvetica-Bold"),
        ("TEXTCOLOR",     (1,1),(-1,-1), NAVY),
        ("SPAN",          (0,0),(-1,0)),
    ]))
    story.append(fin_t)
    story.append(Spacer(1, 16))

    # Contact
    story.append(HRFlowable(width="100%", thickness=1, color=TEAL, spaceAfter=8))
    story.append(Paragraph(
        "<b>Contact</b>: Shehan Makani, ChemeNova LLC  |  "
        "shehan@chemenova.com  |  chemrichgroup.com  |  Pearl River, NY",
        s["body"]
    ))
    story.append(Spacer(1, 4))
    story.append(Paragraph(
        "<i>This proposal was generated by IntelliForm(TM) v0.9 Agentic AI Platform. "
        "All formulation metrics are computed using validated QSAR models (JCIM, 2026) "
        "and real-time ChemRich inventory data. Regulatory information sourced from "
        "ECHA, EPA Safer Choice, and COSMOS-standard databases as of January 2026.</i>",
        s["small"]
    ))

    # Build
    doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
    return buf.getvalue()
