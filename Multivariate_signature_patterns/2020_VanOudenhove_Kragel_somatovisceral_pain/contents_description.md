# Somatovisceral pain classifier (Van Oudenhove, Kragel et al. 2020)

## Overview

A **support-vector-classifier** trained to distinguish **somatic
(cutaneous heat) pain** from **visceral (rectal-distension) pain** based
on fMRI activity within the Yeo 7-network parcellation. Provided as
Yeo-network beta coefficients in a single `.mat`, plus an example
application script.

**Primary reference (open access).** Van Oudenhove, L., Kragel, P. A.,
Dupont, P., Ly, H. G., Pazmany, E., Enzlin, P., Rubio, A., Delon-Martin,
C., Bonaz, B., Aziz, Q., Tack, J., Fukudo, S., Kano, M., & Wager, T. D.
(2020). *Common and distinct neural representations of aversive
somatic and visceral stimulation in healthy individuals.* **Nature
Communications, 11**, 5939.
[doi:10.1038/s41467-020-19688-8](https://doi.org/10.1038/s41467-020-19688-8)
· [local PDF](./VanOudenhove_2020_NatComms_somatovisceral_pain.pdf)

## Key images

![Van Oudenhove 2020 inventory / Yeo-network betas](./png_images/VanOudenhove2020_inventory.png)

This pattern uses **network-level coefficients** rather than
voxelwise weights, so the standard surface/montage/isosurface trio
does not apply. [`visualize_contents.m`](./visualize_contents.m)
renders this inventory / per-Yeo-network coefficient panel into
`png_images/`.

## How to load

Not registered in `load_image_set`. Load the `.mat`:

```matlab
S = load(which('Visceral_vs_Somatic_betas_Yeo_Networks.mat'));
% S contains the per-Yeo-network betas of the SVC.
```

Apply with the worked example:

```matlab
edit classify_somatovisceral_pain.m
```

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `Visceral_vs_Somatic_betas_Yeo_Networks.mat` | MAT | Yeo 7-network coefficients of the visceral-vs-somatic SVC. |
| `classify_somatovisceral_pain.m` | MATLAB | Worked example applying the classifier to a new image. |
| `VanOudenhove_2020_NatComms_somatovisceral_pain.pdf` | PDF | Primary reference (OA). |
| `visualize_contents.m` | MATLAB | Renders Yeo-network bar chart into `png_images/`. |

## Citations

- Van Oudenhove L, Kragel PA, Dupont P, et al. (2020). Common and
  distinct neural representations of aversive somatic and visceral
  stimulation in healthy individuals. *Nat Commun* 11:5939.
  [doi:10.1038/s41467-020-19688-8](https://doi.org/10.1038/s41467-020-19688-8)
