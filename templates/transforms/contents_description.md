# Cross-template warps (MNI152NLin6Asym ↔ MNI152NLin2009cAsym ↔ Colin27)

## Overview

**Inter-template warp fields** for moving image data between the three
MNI spaces commonly used in CANlab work — **MNI152NLin6Asym** (FSL
default, also the canlab2024 grayordinate space),
**MNI152NLin2009cAsym** (fmriprep default), and **Colin27v1998** — in
three software-native formats:

- `ants/` — ANTs `.h5` composite transforms (use with
  `antsApplyTransforms`).
- `fsl/` — Equivalent affine + displacement-field pairs in FSL format
  (use with `applywarp`).
- `spm/` — Equivalent affine CSVs + SPM y-field deformations (use with
  `spm_get_space` / `Normalise: Write`, or the helper
  [`code/apply_spm_warp.m`](./code/apply_spm_warp.m)).

The original ANTs transforms were produced by running each template
through **fmriprep 20.2.3** with surface reconstruction enabled; a
**subcortically-weighted** set (`*_subctx_*`) was generated with more
aggressive SyN parameters and a cifti-volumetric mask metric for
improved subcortical alignment, since the default fmriprep warps prioritise
cortical correspondence. See [`code/README.md`](./code/README.md) for
the full methods narrative and example invocation lines for each
toolchain.

A `leads_dbs/` subfolder ships the Lead-DBS community's
MNI152NLin6Asym <-> MNI152NLin2009bAsym warps (mass-univariate
DBS-targeting use case) — see its
[`README.txt`](./ants/leads_dbs/README.txt) and `COPYING` for
provenance and licence.

[`download_warpfield.m`](./download_warpfield.m) is a one-line helper
that fetches the large SPM `y_*` deformation files from figshare on
demand, so the repository checkout itself can stay small.

## Primary reference

Derived from publicly available MNI / Colin27 templates using
[fmriprep](https://fmriprep.org/) (ANTs SyN registration) as part of
the CANlab build pipeline. For the underlying registration method
see:

- Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008).
  *Symmetric diffeomorphic image registration with cross-correlation.*
  **Medical Image Analysis, 12**(1), 26–41.
  [doi:10.1016/j.media.2007.06.004](https://doi.org/10.1016/j.media.2007.06.004)
- Esteban O, Markiewicz CJ, Blair RW, et al. (2019). *fMRIPrep: A
  robust preprocessing pipeline for functional MRI.* **Nature Methods,
  16**(1), 111–116.
  [doi:10.1038/s41592-018-0235-9](https://doi.org/10.1038/s41592-018-0235-9)

## How to use

ANTs (recommended — single resampling):

```bash
antsApplyTransforms -d 3 \
    -i input_MNI152NLin6Asym.nii.gz \
    -r MNI152NLin2009cAsym_T1_1mm.nii.gz \
    -t ants/MNI152NLin6Asym_to_MNI152NLin2009cAsym.h5 \
    -o output_in_fmriprep_space.nii.gz
```

FSL (apply both transforms at once — see `code/README.md` for premat
ordering):

```bash
applywarp -i MNI152NLin6Asym_T1_1mm.nii.gz \
          -r MNI152NLin2009cAsym_T1_1mm.nii.gz \
          -o MNI152NLin6Asym_to_MNI152NLin2009cAsym_fsl.nii.gz \
          --premat=fsl/00_fsl_to_fmriprep_AffineTransform_mod.mat \
          -w     fsl/01_fsl_to_fmriprep_DisplacementFieldTransform_mod.nii.gz
```

SPM (MATLAB, via the helper `code/apply_spm_warp.m`):

```matlab
apply_transforms('MNI152NLin2009cAsym_T1_1mm.nii', ...
                 'MNI152NLin6Asym_T1_1mm.nii', ...
                 [], ...
                 'spm/01_fmriprep_to_fsl_AffineTransform.csv', ...
                 'spm/y_00_fmriprep_to_fsl_DisplacementFieldTransform.nii', ...
                 'MNI152NLin2009cAsym_to_MNI152NLin6Asym_spm.nii', 1);
```

If the SPM `y_*` deformation files are not present locally:

```matlab
download_warpfield('MNI152NLin6Asym', 'MNI152NLin2009cAsym', 'spm');
```

## File inventory

| File / dir | What it is |
| --- | --- |
| `ants/MNI152NLin6Asym_to_MNI152NLin2009cAsym.h5` | ANTs composite transform, FSL to fmriprep space. |
| `ants/MNI152NLin2009cAsym_to_MNI152NLin6Asym.h5` | Inverse direction. |
| `ants/MNI152NLin6Asym_to_MNI152NLin2009cAsym_subctx.h5` | Subcortically-weighted version (and inverse pair). |
| `ants/MNI152NLin2009cAsym_to_MNI152NLin6Asym_subctx.h5` | Inverse subcortically-weighted. |
| `ants/MNI152NLin6Asym_to_MNI152NLin2009bAsym.h5` | FSL to MNI152NLin2009bAsym (and inverse). |
| `ants/MNI152NLin2009bAsym_to_MNI152NLin2009cAsym.h5` | Between 2009b and 2009c (and inverse). |
| `ants/MNI152NLin6Asym_to_MNIColin27v1998.h5` | FSL to Colin27 (and inverse). |
| `ants/MNI152NLin2009cAsym_to_MNIColin27v1998.h5` | fmriprep to Colin27 (and inverse). |
| `ants/leads_dbs/` | Lead-DBS MNI152NLin6Asym <-> MNI152NLin2009bAsym warps (see `README.txt`, `COPYING`). |
| `fsl/00_*` / `01_*` `_AffineTransform_mod.mat` / `_DisplacementFieldTransform_mod.nii.gz` | FSL-format affine + warp pairs (cortical and `_subctx`). |
| `spm/00_*` / `01_*` `_AffineTransform.csv` | Affines as CSV for SPM. |
| `spm/y_01_fsl_to_fmriprep_DisplacementFieldTransform.nii` | SPM y-field, FSL to fmriprep (cortical). |
| `spm/y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii` | SPM y-field, subcortically weighted. |
| `spm/README` | Example invocation. |
| `code/README.md` | Full methods write-up + per-toolchain examples. |
| `code/apply_spm_warp.m` | MATLAB helper to apply affine + y-field in one shot. |
| `code/ea_antsmat2mat.m` | Convert ANTs affine to plain `.mat`. |
| `code/ants_to_spm_warps.py`, `code/2_ants_to_spm_warps.py` | Convert ANTs `.h5` to SPM y-field. |
| `code/fmriprep_to_fsl.sh`, `code/fsl_to_fmriprep.sh` | Populate the `ants/` and `fsl/` folders from a single `.h5`. |
| `code/subctx_alignment.sh` | Builds the subcortically-weighted SyN. |
| `code/0_populate_ants_and_convert_to_fsl.sh`, `code/1_ants_to_spm_mat.sh` | Driver scripts. |
| `download_warpfield.m` | Fetches large SPM y-fields from figshare on demand. |

## Citations

- Avants BB, Epstein CL, Grossman M, Gee JC (2008). Symmetric
  diffeomorphic image registration with cross-correlation. *Med Image
  Anal* 12:26–41.
  [doi:10.1016/j.media.2007.06.004](https://doi.org/10.1016/j.media.2007.06.004)
- Esteban O, Markiewicz CJ, Blair RW, et al. (2019). fMRIPrep: A
  robust preprocessing pipeline for functional MRI. *Nat Methods*
  16:111–116.
  [doi:10.1038/s41592-018-0235-9](https://doi.org/10.1038/s41592-018-0235-9)
- Fonov VS, Evans AC, Botteron K, Almli CR, McKinstry RC, Collins DL,
  Brain Development Cooperative Group (2011). Unbiased average
  age-appropriate atlases for pediatric studies. *NeuroImage*
  54:313–327.
  [doi:10.1016/j.neuroimage.2010.07.033](https://doi.org/10.1016/j.neuroimage.2010.07.033)
- Holmes CJ, Hoge R, Collins L, Woods R, Toga AW, Evans AC (1998).
  Enhancement of MR images using registration for signal averaging.
  *J Comput Assist Tomogr* 22:324–333.
  [doi:10.1097/00004728-199803000-00032](https://doi.org/10.1097/00004728-199803000-00032)
- Lead-DBS warps: Horn A, Li N, Dembek TA, et al. (2019). Lead-DBS
  v2: Towards a comprehensive pipeline for deep brain stimulation
  imaging. *NeuroImage* 184:293–316.
  [doi:10.1016/j.neuroimage.2018.08.068](https://doi.org/10.1016/j.neuroimage.2018.08.068)
