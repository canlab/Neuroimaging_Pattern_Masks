# Peripheral-inflammation neuroimaging meta-analysis (Kraynak et al. 2018)

## Overview

Coordinate database for a meta-analysis of human fMRI / PET studies of
**peripheral inflammatory physiology** — immunomodulatory challenges
(e.g. endotoxin, vaccination) and brain–immune correlations. The folder
ships the original study foci (`Kraynak_2018_Meta_DB_foci.mat`) and the
MKDA `SETUP.mat`; no voxelwise NIfTI consensus maps are included. To
regenerate consensus maps, re-run the MKDA pipeline against the bundled
SETUP / foci.

See [`Readme.md`](./Readme.md) for the authoritative author write-up.

## Primary reference

Kraynak, T. E., Marsland, A. L., Wager, T. D., & Gianaros, P. J. (2018).
Functional neuroanatomy of peripheral inflammatory physiology: a
meta-analysis of human neuroimaging studies. *Neuroscience &
Biobehavioral Reviews*, 94, 76–92.
[doi:10.1016/j.neubiorev.2018.07.013](https://doi.org/10.1016/j.neubiorev.2018.07.013)

(No paper PDF is bundled in this folder.)

## Key images

No NIfTI result maps are bundled — `visualize_contents.m` is a stub
that warns the user. To produce maps, run the CanlabCore MKDA
pipeline (`Meta_Activation_FWE` etc.) against the database in
`Kraynak_2018_Meta_DB_foci.mat`.

## How to load

Not registered in `load_image_set`. Load the database / SETUP directly:

```matlab
DB    = load(which('Kraynak_2018_Meta_DB_foci.mat'));
SETUP = load(which('SETUP.mat'));
```

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `Kraynak_2018_Meta_DB_foci.mat` | MAT | Coordinate database for the inflammation meta-analysis (study foci + study metadata). |
| `SETUP.mat` | MAT | MKDA analysis setup structure. |
| `Readme.md` | text | Author readme (citation + scope). |
| `visualize_contents.m` | MATLAB | Stub renderer; emits a warning since no NIfTIs are bundled. |

## Citations

- Kraynak TE, Marsland AL, Wager TD, Gianaros PJ (2018). Functional
  neuroanatomy of peripheral inflammatory physiology: a meta-analysis
  of human neuroimaging studies. *Neurosci Biobehav Rev* 94:76–92.
  [doi:10.1016/j.neubiorev.2018.07.013](https://doi.org/10.1016/j.neubiorev.2018.07.013)
- Eisenberger NI, Moieni M, Inagaki TK, Muscatell KA, Irwin MR (2017).
  In sickness and in health: the co-regulation of inflammation and
  social behavior. *Neuropsychopharmacology* 42:242–253.
  [doi:10.1038/npp.2016.141](https://doi.org/10.1038/npp.2016.141)
