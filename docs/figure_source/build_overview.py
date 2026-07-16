# Build docs/npm_repository_overview.svg  (Figure 1: repository mind-map).
# Run build_thumbnails.py first, then:  python3 build_overview.py
import os
import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
OVDIR = os.path.join(HERE, "_thumbs", "overview")
OUT = os.path.join(C.REPO_ROOT, "docs", "npm_repository_overview.svg")
TARGET_WIDTH_PX = 1180.0    # ~12.3 in @ 96 DPI  (see common.bake_scale / README)

W, H = 1720, 1180
HCX, HCY, HR = 860, 610, 94

cards = [
 dict(key="sig",   x=430,  y=360, color="#E36A5C", dark="#B8463B",
      title="Multivariate Signatures", count="29 predictive brain patterns",
      eg="NPS · SIIPS1 · PINES · VPS · FEPS", thumb="sig.png"),
 dict(key="atlas", x=860,  y=250, color="#2FA4A0", dark="#1F6E6B",
      title="Atlases &amp; Parcellations", count="33 labeled atlases &amp; networks",
      eg="CANlab · Glasser · Schaefer · Tian", thumb="atlas.png"),
 dict(key="meta",  x=1290, y=360, color="#E0A44E", dark="#B57C25",
      title="Meta-Analysis Maps", count="23 consensus map sets",
      eg="Emotion · Working memory · Placebo", thumb="meta_kober.png"),
 dict(key="ns",    x=255,  y=650, color="#7E6BC4", dark="#584A91",
      title="Neurosynth Maps", count="100 topics + MKDA term maps",
      eg="Automated large-scale meta-analysis", thumb="neurosynth_iso.png"),
 dict(key="study", x=1465, y=650, color="#4C86C6", dark="#2F5E93",
      title="Individual Study Maps", count="Per-study effect maps",
      eg="Bayes-factor evidence maps", thumb="study.png"),
 dict(key="basis", x=540,  y=965, color="#5AAE7A", dark="#3B7E54",
      title="Spatial Basis Functions", count="5 gradient / ICA map sets",
      eg="Cortical gradients · ICAs · transcriptomic", thumb=None),
 dict(key="tmpl",  x=1180, y=965, color="#8792A6", dark="#5B6070",
      title="Templates", count="MNI reference brains",
      eg="T1 templates · cerebellum · transforms", thumb=None),
]

CW, CH = 350, 150   # card size
TR = 52             # thumbnail radius


def edge_point_card(c):
    cx, cy = c["x"], c["y"]
    dx, dy = HCX - cx, HCY - cy
    if abs(dx) < 1e-6:
        dx = 1e-6
    scale = min((CW / 2) / abs(dx), (CH / 2) / abs(dy)) if dy != 0 else (CW / 2) / abs(dx)
    return cx + dx * scale, cy + dy * scale


def hub_edge(px, py):
    import math
    a = math.atan2(py - HCY, px - HCX)
    return HCX + HR * math.cos(a), HCY + HR * math.sin(a)


def thumb_uri(name):
    return C.datauri(os.path.join(OVDIR, name))


svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">')
svg.append('<defs>')
svg.append('<linearGradient id="bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#fbfcfe"/><stop offset="1" stop-color="#eef1f7"/></linearGradient>')
svg.append('<radialGradient id="hubglow" cx="0.5" cy="0.5" r="0.5"><stop offset="0" stop-color="#2FA4A0" stop-opacity="0.20"/><stop offset="1" stop-color="#2FA4A0" stop-opacity="0"/></radialGradient>')
svg.append('<radialGradient id="hubfill" cx="0.38" cy="0.32" r="0.85"><stop offset="0" stop-color="#2C4A63"/><stop offset="0.55" stop-color="#1E3346"/><stop offset="1" stop-color="#152535"/></radialGradient>')
for c in cards:
    svg.append(f'<clipPath id="clip_{c["key"]}"><circle cx="{c["x"]-CW/2+TR+18}" cy="{c["y"]}" r="{TR}"/></clipPath>')
    svg.append(f'<linearGradient id="conn_{c["key"]}" x1="0" y1="0" x2="1" y2="0"><stop offset="0" stop-color="{c["color"]}" stop-opacity="0.30"/><stop offset="1" stop-color="{c["color"]}" stop-opacity="0.95"/></linearGradient>')
svg.append('<filter id="cardshadow" x="-20%" y="-20%" width="140%" height="150%"><feDropShadow dx="0" dy="6" stdDeviation="9" flood-color="#243049" flood-opacity="0.16"/></filter>')
svg.append('<filter id="hubshadow" x="-40%" y="-40%" width="180%" height="180%"><feDropShadow dx="0" dy="8" stdDeviation="14" flood-color="#10202f" flood-opacity="0.30"/></filter>')
svg.append('</defs>')

svg.append(f'<rect width="{W}" height="{H}" fill="url(#bg)"/>')
svg.append(f'<circle cx="{HCX}" cy="{HCY}" r="430" fill="url(#hubglow)"/>')

# connectors (behind cards)
for c in cards:
    px, py = edge_point_card(c)
    hx, hy = hub_edge(px, py)
    d = C.cubic(hx, hy, px, py, curv=0.42)
    svg.append(f'<path d="{d}" fill="none" stroke="url(#conn_{c["key"]})" stroke-width="5" stroke-linecap="round"/>')
    svg.append(f'<circle cx="{hx:.1f}" cy="{hy:.1f}" r="5.5" fill="{c["color"]}"/>')
    svg.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6" fill="#ffffff" stroke="{c["color"]}" stroke-width="3"/>')

# title
svg.append(f'<text x="70" y="86" font-size="40" font-weight="700" fill="#22314a">Neuroimaging Pattern Masks</text>')
svg.append(f'<text x="72" y="120" font-size="19" fill="#6b7089">A CANlab open library of brain signatures, atlases &amp; meta-analytic maps</text>')


def card(c):
    x0, y0 = c["x"] - CW / 2, c["y"] - CH / 2
    g = ['<g filter="url(#cardshadow)">']
    g.append(f'<rect x="{x0}" y="{y0}" width="{CW}" height="{CH}" rx="22" fill="#ffffff"/>')
    g.append(f'<rect x="{x0}" y="{y0}" width="10" height="{CH}" rx="5" fill="{c["color"]}"/>')
    g.append('</g>')
    tcx = c["x"] - CW / 2 + TR + 18
    g.append(f'<circle cx="{tcx}" cy="{c["y"]}" r="{TR+5}" fill="{c["color"]}" opacity="0.14"/>')
    if c["thumb"]:
        g.append(f'<image href="{thumb_uri(c["thumb"])}" x="{tcx-TR}" y="{c["y"]-TR}" width="{2*TR}" height="{2*TR}" clip-path="url(#clip_{c["key"]})" preserveAspectRatio="xMidYMid meet"/>')
    else:
        g.append(C.brain_glyph(tcx, c["y"], TR * 0.74, c["color"], c["dark"], sulci="#ffffff"))
    g.append(f'<circle cx="{tcx}" cy="{c["y"]}" r="{TR+5}" fill="none" stroke="{c["color"]}" stroke-width="3" opacity="0.85"/>')
    txl = tcx + TR + 26
    g.append(f'<text x="{txl}" y="{c["y"]-20}" font-size="21" font-weight="700" fill="#22314a">{c["title"]}</text>')
    g.append(f'<text x="{txl}" y="{c["y"]+7}" font-size="14.5" font-weight="600" fill="{c["dark"]}">{c["count"]}</text>')
    g.append(f'<text x="{txl}" y="{c["y"]+30}" font-size="13" fill="#8a93a6">{c["eg"]}</text>')
    return "".join(g)


for c in cards:
    svg.append(card(c))

# hub
svg.append(f'<circle cx="{HCX}" cy="{HCY}" r="{HR}" fill="url(#hubfill)" filter="url(#hubshadow)"/>')
svg.append(f'<circle cx="{HCX}" cy="{HCY}" r="{HR}" fill="none" stroke="#ffffff" stroke-width="2" opacity="0.18"/>')
svg.append(C.brain_glyph(HCX, HCY - 14, 40, "#7FD3CE", "#2FA4A0", sulci="#eafffd", opacity=0.95))
svg.append(f'<text x="{HCX}" y="{HCY+46}" font-size="27" font-weight="800" fill="#ffffff" text-anchor="middle" letter-spacing="2">NPM</text>')
svg.append(f'<text x="{HCX}" y="{HCY+68}" font-size="12" fill="#9fb4c4" text-anchor="middle" letter-spacing="1.5">OPEN  LIBRARY</text>')

# footer
svg.append(f'<text x="70" y="{H-40}" font-size="14" fill="#8a93a6">Cognitive &amp; Affective Neuroscience Lab (Tor Wager, PI)  ·  github.com/canlab/Neuroimaging_Pattern_Masks</text>')
svg.append(f'<text x="{W-70}" y="{H-40}" font-size="14" fill="#aab2c2" text-anchor="end">MATLAB · NIfTI · CanlabCore</text>')
svg.append('</svg>')

out = C.bake_scale("\n".join(svg), TARGET_WIDTH_PX)
with open(OUT, "w") as fh:
    fh.write(out)
print(f"wrote {OUT}  (native width {TARGET_WIDTH_PX/96:.1f} in)")
