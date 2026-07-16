# Overview-figure source

Scripts that generate the two editable SVG figures used in the repo READMEs:

| Output (in `docs/`) | Built by | What it shows |
| --- | --- | --- |
| `npm_repository_overview.svg` | `build_overview.py` | Radial mind-map of the seven repository collections |
| `multivariate_signature_taxonomy.svg` | `build_taxonomy.py` | The 36 signatures grouped into six domains + sub-branches |

Both embed small circular **brain renders taken from the repo's own
`png_images/`** folders (base64-embedded, so the SVGs are self-contained) plus a
stylised vector brain glyph, curved connectors, and a modern seaborn-ish
palette. They are plain `<text>`/`<rect>`/`<path>`/`<image>` — fully editable in
PowerPoint or Illustrator.

## Regenerate

```bash
cd docs/figure_source
python3 build_thumbnails.py     # crops png_images renders -> _thumbs/  (run first)
python3 build_overview.py       # -> ../npm_repository_overview.svg
python3 build_taxonomy.py       # -> ../multivariate_signature_taxonomy.svg
```

Dependencies: Python 3 + **Pillow** (`pip install pillow`). Text is measured
with **Arial** for wrap widths that match PowerPoint (falls back to a default
font if Arial is missing — install it for best results). `_thumbs/` is
regenerable and git-ignored.

Preview on macOS without opening PowerPoint:
`qlmanage -t -s 1400 -o . ../multivariate_signature_taxonomy.svg`
(qlmanage centre-crops to a square; pad the short axis in a scratch copy to see
the whole thing.)

## Editing guide

- **Re-classify / add / rename a signature** — edit the `DOMAINS` list in
  `build_taxonomy.py`. Each leaf is
  `("Display name", "First-author · year", "<thumb-key>" | None, "<tag>" | None)`.
  A signature may appear in more than one branch (e.g. Coll 2022 and Silvestrini
  2020 each appear twice). `<thumb-key>` must exist in `build_thumbnails.SIGNATURES`;
  use `None` to draw a plain vector brain in the domain colour.
- **Add a thumbnail** — add `"<thumb-key>": "<repo-relative *_surface.png>"` to
  `SIGNATURES` in `build_thumbnails.py`, rerun it, then rebuild the figure.
- **Domain colours / sub-branches** — the `color`/`dark`/`tint` fields and the
  `subs` structure in `DOMAINS`. Columns are packed two domains high in
  `COLUMNS` to balance height.
- **Overview cards / counts** — the `cards` list in `build_overview.py`.

## ⚠️ Why the figures are ~12.3 in wide (PowerPoint "Convert to Shape")

PowerPoint converts an imported SVG's `font-size` to points at **96 DPI as a
fixed value** (e.g. `font-size="14.26"` → ~10.7 pt) but sizes each text *box*
from the figure's **on-slide display size**. If you shrink the picture before
converting, the boxes shrink while the font does not, so long names spill out of
their box.

The fix, applied by `common.bake_scale()`, is to author each figure at a native
size of ~**1180 px = 12.3 in @ 96 DPI** (`TARGET_WIDTH_PX`), i.e. roughly slide
width. Then it imports ~1:1 and text stays inside its box. **Keep this scale**
if you edit the builders — regenerating without `bake_scale` reintroduces the
overflow.

Reliable PowerPoint workflow:

1. Drop the SVG on the slide (it lands ~12 in wide — about slide width).
2. **Convert to Shape *first*, at that size** — do not shrink the picture beforehand.
3. Only then resize, by selecting the resulting **group** and dragging a corner.
   Resizing a group scales the font and the box together, so text cannot spill.

The one thing to avoid is shrinking the *picture* before converting.
