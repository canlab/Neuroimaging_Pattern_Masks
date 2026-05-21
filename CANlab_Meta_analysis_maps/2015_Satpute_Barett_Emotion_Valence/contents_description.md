# Emotion valence meta-analysis database (Satpute, Barrett & Wager 2015)

## Overview

Coordinate database and MKDA pipeline setup used by the Satpute / Barrett
/ Wager line of work on the **neural representation of emotional
valence**. The folder ships only the analysis setup (`SETUP.mat`) and
the underlying emotion meta-analysis coordinate database
(`Emotion_Meta_DB_all_3_30_13.mat`, 30 March 2013 snapshot) — no
voxelwise NIfTI results are stored here. Re-run the MKDA pipeline (see
the CanlabCore `meta_analysis_mask_tool` workflow) against `SETUP.mat`
to regenerate consensus maps.

## Primary reference

Lindquist, K. A., Satpute, A. B., Wager, T. D., Weber, J., & Barrett,
L. F. (2016). The brain basis of positive and negative affect: evidence
from a meta-analysis of the human neuroimaging literature. *Cerebral
Cortex*, 26(5), 1910–1922.
[doi:10.1093/cercor/bhv001](https://doi.org/10.1093/cercor/bhv001)

(No PDF is included in this folder.)

## Key images

No NIfTI results are bundled here — `visualize_contents.m` will print a
warning that no images were found. To produce maps, run the MKDA
pipeline (`Meta_Activation_FWE` etc., see CanlabCore) on the database
in `Emotion_Meta_DB_all_3_30_13.mat`.

## How to load

Not registered in `load_image_set` (the **related** 5-emotion meta-set
*is* registered as `'emometa'` — see
[`../2015_Wager_Kang_etal_Emotion_Meta_BSPP`](../2015_Wager_Kang_etal_Emotion_Meta_BSPP)).
Load the database / SETUP directly:

```matlab
DB    = load(which('Emotion_Meta_DB_all_3_30_13.mat'));
SETUP = load(which('SETUP.mat'));
```

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `Emotion_Meta_DB_all_3_30_13.mat` | MAT | Emotion meta-analysis coordinate database (30 March 2013 snapshot). |
| `SETUP.mat` | MAT | MKDA analysis setup structure for the valence meta-analysis. |
| `visualize_contents.m` | MATLAB | Stub renderer; emits a warning since no NIfTIs are bundled. |

## Citations

- Lindquist KA, Satpute AB, Wager TD, Weber J, Barrett LF (2016). The
  brain basis of positive and negative affect: evidence from a
  meta-analysis of the human neuroimaging literature. *Cereb Cortex*
  26:1910–1922.
  [doi:10.1093/cercor/bhv001](https://doi.org/10.1093/cercor/bhv001)
- Wager TD, Kang J, Johnson TD, Nichols TE, Satpute AB, Barrett LF
  (2015). A Bayesian model of category-specific emotional brain
  responses. *PLoS Comput Biol* 11:e1004066.
  [doi:10.1371/journal.pcbi.1004066](https://doi.org/10.1371/journal.pcbi.1004066)
