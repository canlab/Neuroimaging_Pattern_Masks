# -*- coding: utf-8 -*-
# Build docs/multivariate_signature_taxonomy.svg  (Figure 2: signature taxonomy).
# Run build_thumbnails.py first, then:  python3 build_taxonomy.py
#
# To re-classify a signature, move / add / edit an entry in DOMAINS below.
# A signature may appear in more than one branch (e.g. Coll 2022 is under both
# Pain > Modulation & value and Appetitive > Reward & craving). The chip's
# thumbnail is the 3rd tuple field (a key in build_thumbnails.SIGNATURES), or
# None to draw a plain vector brain glyph in the domain colour.
import os
import common as C
from PIL import ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
STH = os.path.join(HERE, "_thumbs", "signatures")
OUT = os.path.join(C.REPO_ROOT, "docs", "multivariate_signature_taxonomy.svg")
TARGET_WIDTH_PX = 1180.0    # ~12.3 in @ 96 DPI  (see common.bake_scale / README)

# Arial gives predictable metrics that match PowerPoint; fall back to PIL's
# default font if Arial is unavailable (wrapping stays close enough).
_CANDIDATES = {
    True:  ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "Arial Bold.ttf", "arialbd.ttf"],
    False: ["/System/Library/Fonts/Supplemental/Arial.ttf", "Arial.ttf", "arial.ttf"],
}
_fb = {}


def _font(size, bold):
    key = (size, bold)
    if key not in _fb:
        for path in _CANDIDATES[bold]:
            try:
                _fb[key] = ImageFont.truetype(path, size)
                break
            except OSError:
                continue
        else:
            _fb[key] = ImageFont.load_default()
    return _fb[key]


def measure(text, size, bold=True):
    return _font(size, bold).getlength(text)


def wrap(text, size, maxw, bold=True, maxlines=2):
    lines, cur = [], ""
    for w in text.split():
        t = (cur + " " + w).strip()
        if measure(t, size, bold) <= maxw or not cur:
            cur = t
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    if len(lines) > maxlines:
        lines = lines[:maxlines - 1] + [" ".join(lines[maxlines - 1:])]
    return lines


def durl(key):
    return C.datauri(os.path.join(STH, key + ".png"))


esc = C.esc

# ---------------- taxonomy ----------------
# leaf: (display name, "Author · year", thumbnail key or None, tag or None)
DOMAINS = [
 dict(key="pain", title="Pain", color="#E05C4E", dark="#B23F33", tint="#FBEAE7",
   subs=[("Nociceptive systems", [
            ("NPS — Neurologic Pain Signature", "Wager · 2013", "nps", "on request"),
            ("SIIPS1 — cerebral pain", "Woo · 2017", "siips", None),
            ("Facial expression of pain (FEPS)", "Picard · 2024", "feps", None),
            ("Somatovisceral pain", "Van Oudenhove · 2020", None, None),
            ("Pain mediation (cPDM)", "Geuter · 2020", "cpdm", None)]),
         ("Vicarious / observed", [
            ("Vicarious pain (VPS)", "Krishnan · 2016", "vps", None),
            ("General vicarious pain", "Zhou · 2020", "genvic", None)]),
         ("Modulation & value", [
            ("Placebo analgesia", "Wager · 2011", "placebo", None),
            ("Pain × cognitive control", "Silvestrini · 2020", "paincog", None),
            ("Pain decision value", "Coll · 2022", "painval", None)])]),
 dict(key="clinical", title="Clinical", color="#5AA46E", dark="#3B7A50", tint="#E6F2EA",
   subs=[("Patient & chronic-pain markers", [
            ("Fibromyalgia — pain & multisensory", "López-Solà · 2017", "fibro", None),
            ("Chronic back pain (S1)", "Lee · 2019", "backpain", None),
            ("Fear of neck movement (PiFoneM)", "Murillo · 2026", "pifonem", None)])]),
 dict(key="aversive", title="Aversive / Negative affect", color="#8A6FC0", dark="#5F4A94", tint="#EEEAF7",
   subs=[("Threat", [
            ("Threat conditioning (ImEx)", "Reddan · 2018", "threat", None),
            ("Subjective fear (VIFS)", "Zhou · 2021", "vifs", None),
            ("Multiaversive (MPA2)", "Čeko · 2021", "mpa2", None)]),
         ("Negative emotion", [
            ("Negative affect (PINES)", "Chang · 2015", "pines", None),
            ("Negative emotions (BPLS)", "Kragel · 2015", "bpls_neg", None),
            ("Negative emotion schemas", "Kragel · 2019", "schema_neg", None)]),
         ("Social cognition", [
            ("Empathic distress", "Ashar · 2017", "distress", None),
            ("Interpersonal guilt", "Yu · 2019", "guilt", None),
            ("Romantic rejection (dpSP)", "Woo · 2015", "rejection", None)])]),
 dict(key="physio", title="Physiology", color="#3FA0A8", dark="#256F76", tint="#E2F1F2",
   subs=[("Autonomic responses", [
            ("Skin conductance (SCL)", "Eisenbarth · 2016", "gsr", None),
            ("Heart rate", "Eisenbarth · 2016", "hr", None),
            ("Skin-conductance + pain", "Matthewson · 2019", None, None)])]),
 dict(key="appetitive", title="Appetitive / Reward", color="#E3A33C", dark="#B57516", tint="#FBF2E0",
   subs=[("Reward & craving", [
            ("Brain Reward Signature (BRS)", "Speer · 2023", "brs", None),
            ("Craving — drug & food (NCS)", "Koban · 2022", "ncs", None),
            ("Sexual-image classifier (BASIC)", "van 't Hof · 2021", "basic", None),
            ("Reward decision value", "Coll · 2022", "rewval", None),
            ("Empathic care", "Ashar · 2017", "care", None)]),
         ("Positive emotion", [
            ("Positive emotions (BPLS)", "Kragel · 2015", "bpls_pos", None),
            ("Positive emotion schemas", "Kragel · 2019", "schema_pos", None)])]),
 dict(key="cognitive", title="Cognitive & Social", color="#3E86C0", dark="#28608F", tint="#E5EFF8",
   subs=[("Attention & control", [
            ("Sustained attention (saCPM)", "Rosenberg · 2017", None, None),
            ("Cognitive control (Stroop)", "Silvestrini · 2020", "stroop", None),
            ("MFC generalizability", "Kragel · 2018", "mfc", None)]),
         ("Social cognition", [
            ("Mentalizing — self / other", "Açıl · 2026", "mentalizing", None)])]),
]
# masonry packing: two stacked domains per column balances the heights
COLUMNS = [["pain", "physio"], ["aversive", "clinical"], ["appetitive", "cognitive"]]
DMAP = {d["key"]: d for d in DOMAINS}

# ---------------- layout constants ----------------
MARGIN = 44
COLGAP = 30
NCOL = 3
PANEL_W = 530
W = MARGIN * 2 + NCOL * PANEL_W + (NCOL - 1) * COLGAP
HUBY = 132
HDR_TOP = 248
HDR_H = 64
AFTER_HDR = 22
SUBHEAD_H = 30
CHIP_GAP = 9
SECTION_GAP = 16
NAME_SZ = 21
AUTH_SZ = 15
LINEH = 24
THUMB_R = 26
TEXT_L = THUMB_R * 2 + 22
CHIP_INSET = 8
CHIP_W = PANEL_W - 2 * CHIP_INSET
WRAP_W = 330    # wrap threshold; leaves slack before the chip edge


def colx(i):
    return MARGIN + i * (PANEL_W + COLGAP)


for d in DOMAINS:
    for sub in d["subs"]:
        for i, lf in enumerate(sub[1]):
            sub[1][i] = lf + (wrap(lf[0], NAME_SZ, WRAP_W),)


def chip_h(lines):
    return 12 + len(lines) * LINEH + 18 + 8


col_bottom = [0] * NCOL
layout = {}
for ci, dks in enumerate(COLUMNS):
    x = colx(ci)
    y = HDR_TOP
    for dk in dks:
        d = DMAP[dk]
        block = {"x": x, "hdr_y": y}
        y2 = y + HDR_H + AFTER_HDR
        secs = []
        for (subname, leaves) in d["subs"]:
            sec = {"sub_y": y2, "subname": subname, "chips": []}
            y2 += SUBHEAD_H + 2
            for lf in leaves:
                h = chip_h(lf[4])
                sec["chips"].append({"y": y2, "h": h, "lf": lf})
                y2 += h + CHIP_GAP
            y2 += SECTION_GAP
            secs.append(sec)
        block["secs"] = secs
        block["bottom"] = y2 - SECTION_GAP
        layout[dk] = block
        y = block["bottom"] + 34
    col_bottom[ci] = y
H = max(col_bottom) + 50

# ---------------- render ----------------
S = []
S.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{int(H)}" viewBox="0 0 {W} {int(H)}" font-family="Arial, Helvetica, sans-serif">')
S.append('<defs>')
S.append('<linearGradient id="bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#fbfcfe"/><stop offset="1" stop-color="#eef1f7"/></linearGradient>')
S.append('<radialGradient id="rootfill" cx="0.4" cy="0.3" r="0.9"><stop offset="0" stop-color="#2C4A63"/><stop offset="0.6" stop-color="#1E3346"/><stop offset="1" stop-color="#152535"/></radialGradient>')
S.append('<filter id="chipsh" x="-8%" y="-15%" width="116%" height="140%"><feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="#243049" flood-opacity="0.13"/></filter>')
S.append('<filter id="hdrsh" x="-15%" y="-25%" width="130%" height="160%"><feDropShadow dx="0" dy="4" stdDeviation="6" flood-color="#243049" flood-opacity="0.20"/></filter>')
S.append('<filter id="rootsh" x="-40%" y="-40%" width="180%" height="180%"><feDropShadow dx="0" dy="6" stdDeviation="10" flood-color="#10202f" flood-opacity="0.28"/></filter>')
for d in DOMAINS:
    S.append(f'<clipPath id="cp_{d["key"]}"><circle cx="0" cy="0" r="{THUMB_R}"/></clipPath>')
S.append('</defs>')
S.append(f'<rect width="{W}" height="{int(H)}" fill="url(#bg)"/>')

S.append(f'<text x="{MARGIN+6}" y="60" font-size="34" font-weight="700" fill="#22314a">Multivariate Signature Patterns</text>')
S.append(f'<text x="{MARGIN+8}" y="90" font-size="16" fill="#6b7089">29 studies · 36 pre-trained neuromarkers, grouped by what they predict · first author &amp; year below each</text>')

RX = W / 2
S.append(f'<circle cx="{RX}" cy="{HUBY}" r="52" fill="url(#rootfill)" filter="url(#rootsh)"/>')
S.append(C.brain_glyph(RX, HUBY - 10, 24, "#7FD3CE", "#2FA4A0", sulci="#eafffd"))
S.append(f'<text x="{RX}" y="{HUBY+26}" font-size="15" font-weight="800" fill="#ffffff" text-anchor="middle" letter-spacing="1">SIGNATURES</text>')


def cubicV(x1, y1, x2, y2):
    my = (y1 + y2) / 2
    return f'M {x1:.1f} {y1:.1f} C {x1:.1f} {my:.1f}, {x2:.1f} {my:.1f}, {x2:.1f} {y2:.1f}'


for ci, dks in enumerate(COLUMNS):
    d = DMAP[dks[0]]
    b = layout[dks[0]]
    cx = b["x"] + PANEL_W / 2
    S.append(f'<path d="{cubicV(RX, HUBY+52, cx, HDR_TOP)}" fill="none" stroke="{d["color"]}" stroke-width="3" opacity="0.5"/>')

for dk, b in layout.items():
    d = DMAP[dk]
    x = b["x"]
    color = d["color"]
    dark = d["dark"]
    nsig = sum(len(s[1]) for s in d["subs"])
    S.append(f'<rect x="{x-4}" y="{b["hdr_y"]-4}" width="{PANEL_W+8}" height="{b["bottom"]-b["hdr_y"]+12:.0f}" rx="18" fill="{d["tint"]}" opacity="0.5"/>')
    S.append(f'<g filter="url(#hdrsh)"><rect x="{x}" y="{b["hdr_y"]}" width="{PANEL_W}" height="{HDR_H}" rx="15" fill="{color}"/></g>')
    S.append(C.brain_glyph(x + 32, b["hdr_y"] + HDR_H / 2, 19, "#ffffff", dark, sulci=color, opacity=0.95))
    htl = wrap(d["title"], 17, PANEL_W - 96, bold=True, maxlines=2)
    if len(htl) == 1:
        S.append(f'<text x="{x+58}" y="{b["hdr_y"]+30}" font-size="17" font-weight="700" fill="#ffffff">{esc(htl[0])}</text>')
        S.append(f'<text x="{x+58}" y="{b["hdr_y"]+50}" font-size="12" fill="#ffffff" opacity="0.85">{nsig} signatures</text>')
    else:
        S.append(f'<text x="{x+58}" y="{b["hdr_y"]+25}" font-size="15" font-weight="700" fill="#ffffff">{esc(htl[0])}</text>')
        S.append(f'<text x="{x+58}" y="{b["hdr_y"]+43}" font-size="15" font-weight="700" fill="#ffffff">{esc(htl[1])}</text>')
        S.append(f'<text x="{x+58}" y="{b["hdr_y"]+58}" font-size="11" fill="#ffffff" opacity="0.85">{nsig} signatures</text>')
    for sec in b["secs"]:
        sy = sec["sub_y"]
        S.append(f'<text x="{x+10}" y="{sy+20}" font-size="12.5" font-weight="700" fill="{dark}" letter-spacing="0.4">{esc(sec["subname"].upper())}</text>')
        lblw = measure(sec["subname"].upper(), 13) + 4
        S.append(f'<line x1="{x+16+lblw}" y1="{sy+15}" x2="{x+PANEL_W-8}" y2="{sy+15}" stroke="{color}" stroke-width="1.5" opacity="0.35"/>')
        for chip in sec["chips"]:
            y = chip["y"]
            h = chip["h"]
            name, cite, thumb, tag, lines = chip["lf"]
            cxl = x + CHIP_INSET
            S.append(f'<g filter="url(#chipsh)"><rect x="{cxl}" y="{y}" width="{CHIP_W}" height="{h}" rx="12" fill="#ffffff"/></g>')
            S.append(f'<rect x="{cxl}" y="{y}" width="5" height="{h}" rx="2.5" fill="{color}"/>')
            gcx = cxl + 18 + THUMB_R
            gcy = y + h / 2
            S.append(f'<circle cx="{gcx}" cy="{gcy}" r="{THUMB_R+3}" fill="{d["tint"]}"/>')
            if thumb:
                S.append(f'<g transform="translate({gcx},{gcy})"><image href="{durl(thumb)}" x="{-THUMB_R}" y="{-THUMB_R}" width="{2*THUMB_R}" height="{2*THUMB_R}" clip-path="url(#cp_{d["key"]})" preserveAspectRatio="xMidYMid slice"/></g>')
            else:
                S.append(C.brain_glyph(gcx, gcy, THUMB_R * 0.8, color, dark, sulci="#ffffff"))
            S.append(f'<circle cx="{gcx}" cy="{gcy}" r="{THUMB_R+3}" fill="none" stroke="{color}" stroke-width="2" opacity="0.8"/>')
            tx = cxl + TEXT_L
            ty = y + 12 + (0 if len(lines) == 2 else 6)
            for li, ln in enumerate(lines):
                S.append(f'<text x="{tx}" y="{ty+18+li*LINEH}" font-size="{NAME_SZ}" font-weight="700" fill="#22314a">{esc(ln)}</text>')
            ay = ty + 18 + len(lines) * LINEH - 2
            S.append(f'<text x="{tx}" y="{ay}" font-size="{AUTH_SZ}" font-style="italic" fill="#7c8598">{esc(cite)}</text>')
            if tag:
                aw = measure(cite, AUTH_SZ, bold=False)
                tgw = measure(tag, 11, bold=True) + 16
                S.append(f'<rect x="{tx+aw+12}" y="{ay-13}" width="{tgw:.0f}" height="17" rx="8.5" fill="{color}" opacity="0.16"/>')
                S.append(f'<text x="{tx+aw+12+tgw/2:.0f}" y="{ay-1}" font-size="10" font-weight="700" fill="{dark}" text-anchor="middle">{esc(tag)}</text>')

S.append(f'<text x="{MARGIN+6}" y="{int(H)-22}" font-size="13" fill="#8a93a6">Some papers contribute more than one signature and appear in multiple branches. · CANlab · github.com/canlab/Neuroimaging_Pattern_Masks</text>')
S.append('</svg>')

out = C.bake_scale("\n".join(S), TARGET_WIDTH_PX)
with open(OUT, "w") as fh:
    fh.write(out)
print(f"wrote {OUT}  native {W}x{int(H)} -> width {TARGET_WIDTH_PX/96:.1f} in")
