# Shared helpers for the NPM overview figures.
#   - brain_glyph : stylised vector brain icon (recolourable)
#   - cubic       : smooth S-curve connector between two points
#   - datauri     : base64-embed a PNG so the SVG stays self-contained
#   - esc         : XML-escape text (&, <, >)
#   - bake_scale  : uniformly shrink a finished SVG so PowerPoint "Convert to
#                   Shape" keeps text inside its box (see README.md for why)
import base64, math, os, re
import xml.dom.minidom as MD

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def datauri(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ---- brain glyph: organic wobbly-blob cortex + midline + sulci ----
def brain_path(cx, cy, R):
    pts = []
    N = 72
    for i in range(N):
        t = 2 * math.pi * i / N
        rr = R * (1 + 0.06 * math.sin(3 * t + 0.6) + 0.045 * math.sin(5 * t) + 0.03 * math.sin(2 * t))
        x = cx + rr * 1.06 * math.cos(t)
        y = cy - rr * 0.92 * math.sin(t)
        pts.append((x, y))
    d = "M {:.2f} {:.2f} ".format(*pts[0])
    for i in range(1, N):
        d += "L {:.2f} {:.2f} ".format(*pts[i])
    d += "Z"
    return d


def brain_glyph(cx, cy, R, fill, stroke, sulci="#ffffff", sw=None, opacity=1.0):
    if sw is None:
        sw = max(1.0, R * 0.09)
    d = brain_path(cx, cy, R)
    mid = 'M {:.2f} {:.2f} C {:.2f} {:.2f}, {:.2f} {:.2f}, {:.2f} {:.2f}'.format(
        cx, cy - R * 0.86, cx + R * 0.10, cy - R * 0.3, cx - R * 0.10, cy + R * 0.3, cx, cy + R * 0.82)

    def wave(x0, x1, y, amp):
        return 'M {:.2f} {:.2f} Q {:.2f} {:.2f}, {:.2f} {:.2f} T {:.2f} {:.2f}'.format(
            x0, y, (x0 + x1) / 2, y - amp, (x0 + x1) / 2, y, x1, y + amp * 0.4)

    s = [
        wave(cx - R * 0.78, cx - R * 0.18, cy - R * 0.32, R * 0.16),
        wave(cx - R * 0.72, cx - R * 0.2, cy + R * 0.22, R * 0.14),
        wave(cx + R * 0.18, cx + R * 0.78, cy - R * 0.30, R * 0.16),
        wave(cx + R * 0.2, cx + R * 0.72, cy + R * 0.24, R * 0.14),
    ]
    out = '<g opacity="{:.2f}">'.format(opacity)
    out += '<path d="{}" fill="{}" stroke="{}" stroke-width="{:.2f}" stroke-linejoin="round"/>'.format(d, fill, stroke, sw)
    out += '<path d="{}" fill="none" stroke="{}" stroke-width="{:.2f}" stroke-linecap="round" opacity="0.85"/>'.format(mid, sulci, sw * 0.85)
    for w in s:
        out += '<path d="{}" fill="none" stroke="{}" stroke-width="{:.2f}" stroke-linecap="round" opacity="0.7"/>'.format(w, sulci, sw * 0.75)
    out += '</g>'
    return out


def cubic(x1, y1, x2, y2, curv=0.5):
    dx, dy = x2 - x1, y2 - y1
    c1x, c1y = x1 + dx * curv, y1 + dy * 0.08
    c2x, c2y = x2 - dx * curv, y2 - dy * 0.08
    return 'M {:.1f} {:.1f} C {:.1f} {:.1f}, {:.1f} {:.1f}, {:.1f} {:.1f}'.format(x1, y1, c1x, c1y, c2x, c2y, x2, y2)


# ---- scale-bake so the SVG converts cleanly to shapes in PowerPoint ----
# PowerPoint reads font-size (px) as points at 96 DPI (a FIXED value) but sizes
# each text box from the figure's on-slide display size. If the picture is
# shrunk before "Convert to Shape", the boxes shrink while the font does not, so
# long names overflow. Authoring the SVG at ~slide-width native size (~1180 px =
# 12.3 in @ 96 DPI) means it is placed ~1:1 and text stays inside its box.
_SKIP_TAGS = {"lineargradient", "radialgradient", "stop"}
_LEN_ATTRS = {"x", "y", "cx", "cy", "r", "rx", "ry", "x1", "y1", "x2", "y2",
              "width", "height", "font-size", "stroke-width", "letter-spacing",
              "stdDeviation", "dx", "dy"}


def _scale_num(v, S):
    try:
        return f"{float(v) * S:.3f}".rstrip("0").rstrip(".")
    except ValueError:
        return v


def _scale_numbers_in(s, S):
    return re.sub(r"-?\d*\.?\d+", lambda m: _scale_num(m.group(0), S), s)


def _walk(el, S):
    if el.nodeType == 1 and el.tagName.lower() not in _SKIP_TAGS:
        for a in list(el.attributes.keys()):
            if a in _LEN_ATTRS:
                el.setAttribute(a, _scale_num(el.getAttribute(a), S))
            elif a == "d":
                el.setAttribute("d", _scale_numbers_in(el.getAttribute(a), S))
            elif a == "transform":
                tv = re.sub(r"translate\(([^)]*)\)",
                            lambda m: "translate(" + _scale_numbers_in(m.group(1), S) + ")",
                            el.getAttribute("transform"))
                el.setAttribute("transform", tv)
    for c in list(el.childNodes):
        if c.nodeType == 1:
            _walk(c, S)


def bake_scale(svg_string, target_width_px=1180.0):
    """Return svg_string uniformly scaled so its width == target_width_px.
    Gradient elements (objectBoundingBox coords, 0..1) are left untouched."""
    doc = MD.parseString(svg_string)
    svg = doc.documentElement
    w = float(svg.getAttribute("width"))
    h = float(svg.getAttribute("height"))
    S = target_width_px / w
    _walk(svg, S)
    svg.setAttribute("width", f"{w * S:.1f}")
    svg.setAttribute("height", f"{h * S:.1f}")
    svg.setAttribute("viewBox", f"0 0 {w * S:.1f} {h * S:.1f}")
    return doc.toxml()
