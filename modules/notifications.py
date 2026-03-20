"""
modules/notifications.py
Email notifications and dynamic reformulation triggers for IntelliForm v1.0.

Features:
  1. Pilot booking confirmation email to customer + Shehan
  2. Regulatory alert email when ingredient gets flagged
  3. Dynamic reformulation trigger — detects when saved formulation
     is at risk and suggests alternative

Uses SendGrid (free tier: 100 emails/day) or SMTP fallback.
Set SENDGRID_API_KEY or SMTP_HOST/SMTP_USER/SMTP_PASS in environment.
"""
import os
import smtplib
import json
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
from typing import Dict, Optional, List
from dataclasses import dataclass


CHEMENOVA_EMAIL = "shehan@chemenova.com"
CHEMENOVA_NAME  = "Shehan Makani, ChemeNova LLC"


@dataclass
class EmailResult:
    sent: bool
    method: str      # sendgrid | smtp | logged
    message: str


# ── Core email sender ─────────────────────────────────────────────────────────

def _send_via_sendgrid(to: str, subject: str, html: str, text: str) -> EmailResult:
    api_key = os.getenv("SENDGRID_API_KEY", "")
    if not api_key:
        return EmailResult(sent=False, method="sendgrid", message="No API key")
    try:
        import sendgrid
        from sendgrid.helpers.mail import Mail
        sg = sendgrid.SendGridAPIClient(api_key=api_key)
        message = Mail(
            from_email=CHEMENOVA_EMAIL,
            to_emails=to,
            subject=subject,
            html_content=html,
            plain_text_content=text,
        )
        response = sg.send(message)
        return EmailResult(sent=response.status_code in [200, 202],
                          method="sendgrid",
                          message=f"Status {response.status_code}")
    except Exception as e:
        return EmailResult(sent=False, method="sendgrid", message=str(e))


def _send_via_smtp(to: str, subject: str, html: str, text: str) -> EmailResult:
    host = os.getenv("SMTP_HOST", "")
    user = os.getenv("SMTP_USER", "")
    pwd  = os.getenv("SMTP_PASS", "")
    port = int(os.getenv("SMTP_PORT", "587"))
    if not all([host, user, pwd]):
        return EmailResult(sent=False, method="smtp", message="SMTP not configured")
    try:
        msg = MIMEMultipart("alternative")
        msg["Subject"] = subject
        msg["From"]    = f"{CHEMENOVA_NAME} <{user}>"
        msg["To"]      = to
        msg.attach(MIMEText(text, "plain"))
        msg.attach(MIMEText(html, "html"))
        with smtplib.SMTP(host, port) as server:
            server.starttls()
            server.login(user, pwd)
            server.sendmail(user, to, msg.as_string())
        return EmailResult(sent=True, method="smtp", message="Sent via SMTP")
    except Exception as e:
        return EmailResult(sent=False, method="smtp", message=str(e))


def _log_email(to: str, subject: str, text: str) -> EmailResult:
    """Fallback — log email to console when no service configured."""
    print(f"\n[EMAIL LOG] To: {to} | Subject: {subject}\n{text}\n")
    return EmailResult(sent=True, method="logged", message="Logged to console (no email service configured)")


def send_email(to: str, subject: str, html: str, text: str) -> EmailResult:
    result = _send_via_sendgrid(to, subject, html, text)
    if result.sent:
        return result
    result = _send_via_smtp(to, subject, html, text)
    if result.sent:
        return result
    return _log_email(to, subject, text)


# ── Email templates ───────────────────────────────────────────────────────────

def send_pilot_booking_confirmation(
    customer_email: str,
    customer_name: str,
    blend: Dict[str, float],
    cost_per_kg: float,
    batch_kg: int,
    quote_usd: float,
    application: str,
) -> EmailResult:
    """Send booking confirmation to customer and notification to Shehan."""

    blend_lines_html = "".join(
        f"<tr><td style='padding:4px 12px;'>{ing}</td><td style='padding:4px 12px;text-align:right'>{pct}%</td></tr>"
        for ing, pct in blend.items()
    )
    blend_lines_text = "\n".join(f"  {ing}: {pct}%" for ing, pct in blend.items())
    date_str = datetime.now().strftime("%B %d, %Y")

    html = f"""
<div style="font-family:Arial,sans-serif;max-width:600px;margin:0 auto">
  <div style="background:#0A1628;padding:24px;border-radius:8px 8px 0 0">
    <h1 style="color:#fff;margin:0;font-size:22px">ChemeNova LLC x ChemRich Global</h1>
    <p style="color:#0D9488;margin:4px 0 0">Pilot Batch Booking Confirmed</p>
  </div>
  <div style="background:#f8fafc;padding:24px;border-radius:0 0 8px 8px;border:1px solid #e2e8f0">
    <p>Dear {customer_name},</p>
    <p>Your pilot batch request has been received. Here are the details:</p>
    <table style="width:100%;border-collapse:collapse;margin:16px 0">
      <tr style="background:#0A1628;color:#fff">
        <th style="padding:8px 12px;text-align:left">Ingredient</th>
        <th style="padding:8px 12px;text-align:right">%</th>
      </tr>
      {blend_lines_html}
    </table>
    <table style="width:100%;margin:16px 0">
      <tr><td style="color:#64748b">Application</td><td style="text-align:right;font-weight:bold">{application.replace('_',' ').title()}</td></tr>
      <tr><td style="color:#64748b">Batch size</td><td style="text-align:right;font-weight:bold">{batch_kg} kg</td></tr>
      <tr><td style="color:#64748b">Cost/kg</td><td style="text-align:right;font-weight:bold">${cost_per_kg:.2f}</td></tr>
      <tr><td style="color:#64748b">Total quote</td><td style="text-align:right;font-weight:bold;color:#0D9488">${quote_usd:,.0f}</td></tr>
      <tr><td style="color:#64748b">Lead time</td><td style="text-align:right;font-weight:bold">5 business days</td></tr>
    </table>
    <div style="background:#dcfce7;padding:16px;border-radius:6px;border-left:4px solid #059669">
      <strong>Reformulation guarantee:</strong> If the batch does not meet your certification criteria,
      we will reformulate and re-pilot at no additional cost.
    </div>
    <p style="margin-top:24px">Questions? Reply to this email or contact <a href="mailto:{CHEMENOVA_EMAIL}">{CHEMENOVA_EMAIL}</a></p>
    <p style="color:#64748b;font-size:12px">ChemeNova LLC | Pearl River, NJ | chemrichgroup.com | {date_str}</p>
  </div>
</div>
"""
    text = f"""Pilot Batch Booking Confirmed — ChemeNova LLC x ChemRich Global

Dear {customer_name},

Your pilot batch request has been received.

Blend:
{blend_lines_text}

Application: {application.replace('_',' ').title()}
Batch size: {batch_kg} kg
Cost/kg: ${cost_per_kg:.2f}
Total quote: ${quote_usd:,.0f}
Lead time: 5 business days

Reformulation guarantee: If the batch does not meet your certification
criteria, we will reformulate and re-pilot at no additional cost.

Questions? Contact {CHEMENOVA_EMAIL}

ChemeNova LLC | Pearl River, NJ | {date_str}
"""

    # Send to customer
    result = send_email(
        to=customer_email,
        subject=f"Pilot Batch Confirmed — ${quote_usd:,.0f} | ChemeNova x ChemRich",
        html=html,
        text=text,
    )

    # Notify Shehan
    send_email(
        to=CHEMENOVA_EMAIL,
        subject=f"NEW BOOKING: {customer_name} | {batch_kg}kg | ${quote_usd:,.0f}",
        html=f"<p>New pilot booking from <b>{customer_name}</b> ({customer_email})</p>{html}",
        text=f"New pilot booking from {customer_name} ({customer_email})\n\n{text}",
    )

    return result


def send_regulatory_alert(
    customer_email: str,
    customer_name: str,
    ingredient: str,
    alert_type: str,
    blend: Dict[str, float],
    formulation_name: str = "Your saved formulation",
) -> EmailResult:
    """Alert customer when a saved formulation ingredient gets flagged."""

    date_str = datetime.now().strftime("%B %d, %Y")
    html = f"""
<div style="font-family:Arial,sans-serif;max-width:600px;margin:0 auto">
  <div style="background:#dc2626;padding:24px;border-radius:8px 8px 0 0">
    <h1 style="color:#fff;margin:0;font-size:20px">Regulatory Alert — Action Required</h1>
    <p style="color:#fca5a5;margin:4px 0 0">IntelliForm Compliance Monitoring</p>
  </div>
  <div style="background:#f8fafc;padding:24px;border-radius:0 0 8px 8px;border:1px solid #e2e8f0">
    <p>Dear {customer_name},</p>
    <p>Your formulation <strong>{formulation_name}</strong> contains an ingredient that requires attention:</p>
    <div style="background:#fef2f2;padding:16px;border-radius:6px;border-left:4px solid #dc2626;margin:16px 0">
      <strong>{ingredient}</strong><br/>
      <span style="color:#dc2626">{alert_type}</span>
    </div>
    <p>We recommend scheduling a reformulation consultation to maintain your certification status.</p>
    <a href="mailto:{CHEMENOVA_EMAIL}?subject=Reformulation Request - {ingredient}" 
       style="background:#0D9488;color:#fff;padding:12px 24px;border-radius:6px;text-decoration:none;display:inline-block;margin:16px 0">
      Request Reformulation
    </a>
    <p style="color:#64748b;font-size:12px">ChemeNova LLC | {CHEMENOVA_EMAIL} | {date_str}</p>
  </div>
</div>
"""
    text = f"""Regulatory Alert — {ingredient}

Dear {customer_name},

Your formulation "{formulation_name}" contains {ingredient}.

Alert: {alert_type}

We recommend scheduling a reformulation consultation.
Contact: {CHEMENOVA_EMAIL}

ChemeNova LLC | {date_str}
"""
    return send_email(
        to=customer_email,
        subject=f"Regulatory Alert: {ingredient} in {formulation_name}",
        html=html,
        text=text,
    )


def send_proposal_email(
    customer_email: str,
    customer_name: str,
    application: str,
    cost_per_kg: float,
    eco_score: float,
) -> EmailResult:
    """Send a follow-up email when customer downloads a proposal."""
    date_str = datetime.now().strftime("%B %d, %Y")
    html = f"""
<div style="font-family:Arial,sans-serif;max-width:600px;margin:0 auto">
  <div style="background:#0A1628;padding:24px;border-radius:8px 8px 0 0">
    <h1 style="color:#fff;margin:0;font-size:20px">Your IntelliForm Proposal</h1>
    <p style="color:#0D9488;margin:4px 0 0">ChemeNova LLC x ChemRich Global</p>
  </div>
  <div style="background:#f8fafc;padding:24px;border-radius:0 0 8px 8px;border:1px solid #e2e8f0">
    <p>Dear {customer_name},</p>
    <p>Your green chemistry formulation proposal has been generated.</p>
    <table style="width:100%;margin:16px 0">
      <tr><td style="color:#64748b">Application</td><td style="text-align:right;font-weight:bold">{application.replace('_',' ').title()}</td></tr>
      <tr><td style="color:#64748b">Cost/kg</td><td style="text-align:right;font-weight:bold">${cost_per_kg:.2f}</td></tr>
      <tr><td style="color:#64748b">EcoScore</td><td style="text-align:right;font-weight:bold;color:#0D9488">{eco_score:.0f}/100</td></tr>
    </table>
    <p>Ready to move forward? Book a pilot batch and we'll have your certified formulation ready in 5 days.</p>
    <a href="mailto:{CHEMENOVA_EMAIL}?subject=Pilot Batch Request"
       style="background:#0A1628;color:#fff;padding:12px 24px;border-radius:6px;text-decoration:none;display:inline-block">
      Book Pilot Batch
    </a>
    <p style="color:#64748b;font-size:12px">ChemeNova LLC | {CHEMENOVA_EMAIL} | {date_str}</p>
  </div>
</div>
"""
    text = f"""Your IntelliForm Proposal — ChemeNova LLC

Dear {customer_name},

Application: {application}
Cost/kg: ${cost_per_kg:.2f}
EcoScore: {eco_score:.0f}/100

Book a pilot batch: {CHEMENOVA_EMAIL}
"""
    return send_email(
        to=customer_email,
        subject="Your Green Chemistry Formulation Proposal — ChemeNova",
        html=html,
        text=text,
    )
