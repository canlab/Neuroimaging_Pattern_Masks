# CANlab Neuromarker Gallery

A static, map-based gallery for GitHub Pages, alongside the existing repository
README and documentation. Search by domain, sensory modality, target, name,
acronym, or citation; explore slices with CANlab2024 labels and HCP cortical
surfaces. Separate positive/negative percentile or magnitude controls support
optional synchronization. Map downloads preserve source values; share and embed
URLs encode the selected map and display settings.

## Build and preview

Requires Python 3.9+ and Node 20+ (Node is used for tests, not production).
From the repository root:

```sh
python3 -m venv /tmp/neuromarker-venv
/tmp/neuromarker-venv/bin/pip install -r site/requirements.txt
/tmp/neuromarker-venv/bin/python site/scripts/build.py
/tmp/neuromarker-venv/bin/python site/scripts/validate.py
npm ci --prefix site
npm test --prefix site
python3 -m http.server 8765 --directory site/dist
```

Open http://localhost:8765. For browser integration tests, with this server running:
`node site/tests/browser.mjs` uses installed Chrome locally; CI uses Playwright
Chromium (`npx playwright install --with-deps chromium`).

## Catalog maintenance

Edit `data/catalog.json`; it is the curated source of truth. Each study has
publication metadata, domain/modality/target tags, and explicit maps with stable
IDs, source paths, map roles, descriptions and map-specific target tags. Source
paths are relative to the repository root. Do not infer map role from a file's
extension. The one-time `seed_catalog.py` documents initial extraction but is not
part of the build: running it would overwrite subsequent curation.

The build validates dimensions and affine, expands 4D maps into components,
precomputes independent sign quantiles and trilinear surface samples, and creates
axial thumbnails and static detail pages for social metadata. Existing 3D NIfTI
files are retained byte-for-byte inside gzip; Analyze conversions and 4D splits
preserve physical values and affine. `validate.py` checks all exported maps
against sources, including every frame, and validates surface payloads.

No connectivity-edge models are included. Lee's perfusion map is included; its
connectivity predictor is excluded. Matthewson has no local maps; Van Oudenhove
provides network coefficients. Fibromyalgia is not packaged because the top-level
README describes request-only availability. Supporting masks/p-value maps are
excluded. Supplied statistical-support maps are explicitly labeled; an empty
thresholded component is retained with an explanatory note.

Exact MNI template variants remain undocumented for many source maps. The HCP
midthickness projection is an approximate world-coordinate visualization, not
nonlinear registration. It is not a replacement predictive model. The
`methods.html` page exposes these limitations to users.

## Deployment

`.github/workflows/neuromarker-gallery.yml` builds, validates and publishes
`site/dist` using GitHub Actions. Configure repository Settings → Pages → Source
as **GitHub Actions**. Publishing is triggered by relevant changes on `master`
or manual dispatch; pull requests build/test without publishing. The expected
URL is https://canlab.github.io/Neuroimaging_Pattern_Masks/.

Only generated public web assets enter the Pages artifact. PDFs, videos, local
notes, node_modules and repository documentation are not copied. Shared assets
and map assets have content-derived URLs. Build output is ignored by git.
The runtime has no npm dependencies, CDN imports, backend, MATLAB dependency,
analytics, or uploads. Local overlays remain in the browser and cannot be shared
as reproducible gallery links. NiiVue uses WebGL2; a static preview/download
remains available if initialization fails.

Rollback: revert the relevant source commit and rerun the publishing workflow.
For custom hosting set `GALLERY_ORIGIN` during the build to the canonical base URL.

## Validation coverage

- Node tests: independent percentiles, zero/full tails, ties, nonfinite values,
  absolute thresholds, sign controls, URL parsing and filter semantics.
- Python: all map values and affines, all percentile metadata, unique IDs,
  publisher links, static detail pages, surface sizes, finite values and sampling.
- Browser: search, CANlab atlas controls/readout, volumes, true cortical meshes,
  independent controls, saved URL state, absolute/sign controls, hemispheres,
  NIfTI and PNG downloads, mobile width, embed and empty-results states.

The optional feature-detected WebMCP search tool mirrors the visible search.
Browsers without a WebMCP context ignore it; the standard browser does not verify
native WebMCP registration.

The gallery defaults to tiles with a domain-tree alternative, independently adjustable 35% sign thresholds, no atlas outline, and matte cortical lighting. The color scale is drawn outside the slice canvas and included in PNG exports. Thumbnail overlays are masked with the repository MNI152NLin2009cAsym brain mask; source downloads remain unchanged.

`data/shared/atlas-descriptions.json` and `.csv` provide all 518 short/full parcel names and row-level provenance. `scripts/prepare_labels.py` regenerates these from repository atlas dictionaries; this preparation step is separate from the network-free site build. Cortex names follow Glasser’s neuroanatomical supplement, including “Area posterior 24.” Coarse subcortical components are preserved in the table. The public lookup is `atlas.html`.

Branding uses the unmodified `ncs_thumbnail.png` supplied by Tor Wager, served as `brand/neuromarkers-ncs.png` for the hero, logo, favicon, and social previews. CSS presents it in grayscale at rest and reveals its original colors on hover, focus, or tap. Previous NPS renderings and their offline preparation script are retained as backups. Map-specific catalog tile previews remain tied to their respective maps.
