# Thiebaut de Schotten white-matter / functional component atlas (2020)

## Overview

This folder ships the **BCBlab functional-component white-matter atlas**
described by Thiebaut de Schotten and colleagues (2020). Rather than
parcellating fibre tracts per se, the atlas decomposes white-matter
anatomy into functional **components** that account for the
functional segregation of the cortex. Maps are distributed by the
BCBlab and through NeuroVault collections; this folder is mostly a
**code wrapper** (Python `bcblib`) used to download / preprocess
those collections — it does not redistribute the NIfTIs themselves.

> See [`README.rst`](./README.rst) for canonical download links and
> the NeuroVault collection IDs (per-term-letter splits A–C, D–H,
> I–N, O–R, S–U, V–Z, plus the main 2020 Nat Comms components and
> the BCBlab Atlas of Human Brain Connections).

There is no local `atlas` object to load; treat this folder as a
pointer to the upstream maps.

## Primary reference

- Thiebaut de Schotten, M., Foulon, C., & Nachev, P. (2020).
  *Brain disconnections link structural connectivity with function
  and behaviour.* **Nature Communications, 11**, 5094.
  [doi:10.1038/s41467-020-18920-9](https://doi.org/10.1038/s41467-020-18920-9)

No local PDF is checked in. See the DOI link above.

## Key images

No PNGs are bundled here. Pre-rendered map visualisations are
available directly on the upstream NeuroVault collections (see
[`README.rst`](./README.rst) for IDs). Once you download a component
NIfTI locally, [`visualize_contents.m`](./visualize_contents.m) will
write a montage + isosurface into `png_images/` if you point it at
the downloaded file.

## How to load

There is no `load_atlas` keyword. Download a component NIfTI from the
NeuroVault collections listed in [`README.rst`](./README.rst) and:

```matlab
obj = fmri_data('path/to/component.nii.gz');
```

The BCBtoolkit (<http://toolkit.bcblab.com>) is the preferred way to
work with the BCBlab Atlas of Human Brain Connections.

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `README.rst` | reStructuredText | **Authoritative download / source links.** |
| `LICENSE.txt` | text | Licence. |
| `setup.py`, `__init__.py` | Python | Package metadata for the `bcblib` wrapper. |
| `bcblib/` | dir | BCBlab Python helper library. |
| `bash/` | dir | Bash scripts used by `bcblib`. |
| `visualize_contents.m` | MATLAB | Optional renderer once components are downloaded locally. |

## Citations

- Thiebaut de Schotten M, Foulon C, Nachev P. (2020). Brain
  disconnections link structural connectivity with function and
  behaviour. *Nat Commun* 11:5094.
  [doi:10.1038/s41467-020-18920-9](https://doi.org/10.1038/s41467-020-18920-9)
- Rojkova K, Volle E, Urbanski M, Humbert F, Dell'Acqua F,
  Thiebaut de Schotten M. (2016). Atlasing the frontal lobe
  connections and their variability due to age and education.
  *Brain Struct Funct* 221:1751–1766.
  [doi:10.1007/s00429-015-1001-3](https://doi.org/10.1007/s00429-015-1001-3)
