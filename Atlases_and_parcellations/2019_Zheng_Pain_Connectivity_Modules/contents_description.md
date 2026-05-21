# Pain-evoked functional connectivity modules (Zheng et al. 2019)

## Overview

Functional connectivity **modules** derived from pain-evoked fMRI in a
large multi-study CANlab dataset (Zheng et al. 2019 *Cerebral Cortex*).
Cross-subject networks were constructed at the Brainnetome (273-region)
parcellation level for noxious / innocuous and painful / non-painful
stimulus conditions; community detection then produced module
assignments that show how brain networks reorganise during pain.
This folder contains the **module label NIfTIs**, the **network
adjacency matrices**, the underlying Brainnetome template, and the
multi-study metadata MAT used to drive the original analyses.

> See [`Readme.txt`](./Readme.txt) for the authoritative
> methods walkthrough (preparation, mean-data extraction, network
> construction). Note: this folder ships data and code; there is no
> packaged `atlas` object — the modules live as NIfTIs in
> [`module/`](./module).

## Primary reference

- Zheng, W., Woo, C.-W., Yao, Z., Goldstein, P., Atlas, L. Y.,
  Roy, M., Schmidt, L., Krishnan, A., Jepma, M., Hu, B., & Wager, T. D.
  (2019). *Pain-evoked reorganization in functional brain networks.*
  **Cerebral Cortex, bhz260.**
  [doi:10.1093/cercor/bhz260](https://doi.org/10.1093/cercor/bhz260)

Local preprint: [`2019_Zheng_Draft_Preprint.pdf`](./2019_Zheng_Draft_Preprint.pdf).

## Key images

No pre-rendered PNGs are bundled with this folder.
[`visualize_contents.m`](./visualize_contents.m) renders the
Brainnetome template (used as the parcellation backbone in the
paper) plus the four module label volumes (Noxious, Innocuous,
Painful, Non-painful) into `png_images/`.

## How to load

This is a code + data release, not a packaged `atlas` object:

```matlab
% Brainnetome parcellation used as the backbone
atl = load_atlas('brainnetome');

% Module label volumes per stimulus condition (NIfTIs)
mod_noxious = fmri_data(fullfile('module', 'Noxious', '<file>.nii'));
% (see contents of module/Innocuous|Noxious|Painful|Non-painful)

% Multi-study metadata + group-level data used to build the networks
S = load('canlab_datasets_and_metadata.mat');
```

To run the original pipeline, follow the four steps in
[`Readme.txt`](./Readme.txt) and the scripts under [`code/`](./code).

## File inventory

| File / Folder | Type | What it is |
| --- | --- | --- |
| `Readme.txt` | text | **Authoritative methods + usage walkthrough.** |
| `2019_Zheng_Draft_Preprint.pdf` | PDF | Local draft preprint. |
| `canlab_datasets_and_metadata.mat` | MAT | Multi-study CANlab metadata + extracted regional activity. |
| `code/` | dir | MATLAB scripts (preparation, network construction, etc.). |
| `module/` | dir | Module label NIfTIs grouped by condition (Noxious / Innocuous / Painful / Non-painful). |
| `Network/` | dir | Cross-subject network matrices (`Net_noxious_innocoous.mat`, `Net_pain_nopain.mat`). |
| `Template/` | dir | Brainnetome 273-region atlas NIfTI + label files used as the backbone. |
| `visualize_contents.m` | MATLAB | Renders Brainnetome backbone + module label volumes into `png_images/`. |

## Citations

- Zheng W, Woo C-W, Yao Z, et al. (2019). Pain-evoked reorganization
  in functional brain networks. *Cereb Cortex* bhz260.
  [doi:10.1093/cercor/bhz260](https://doi.org/10.1093/cercor/bhz260)
- Fan L, Li H, Zhuo J, et al. (2016). The human Brainnetome atlas:
  a new brain atlas based on connectional architecture.
  *Cereb Cortex* 26:3508–3526.
  [doi:10.1093/cercor/bhw157](https://doi.org/10.1093/cercor/bhw157)
- Rubinov M, Sporns O. (2010). Complex network measures of brain
  connectivity: uses and interpretations. *NeuroImage* 52:1059–1069.
  [doi:10.1016/j.neuroimage.2009.10.003](https://doi.org/10.1016/j.neuroimage.2009.10.003)
