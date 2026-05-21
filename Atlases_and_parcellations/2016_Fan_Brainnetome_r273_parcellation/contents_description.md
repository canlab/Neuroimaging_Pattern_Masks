# Brainnetome 273-region atlas (Fan et al. 2016)

## Overview

The **Brainnetome atlas** is a connectivity-based parcellation of the
human brain into **246 cortical + subcortical regions** (123 per
hemisphere) derived from in-vivo HCP diffusion and resting-state
data, with each parcel anchored to its anatomical connectivity
fingerprint. This folder ships the Brainnetome cortical/subcortical
parcellation extended with the **CANlab "Wani" brainstem labels** to
reach **273–280 regions**, packaged as a CANlab `atlas` object in MNI
space.

The primary PDF is included as
[`Fan_2016_CerebralCortex.pdf`](./Fan_2016_CerebralCortex.pdf).
Build helpers include the constructor
[`Brainnetome_create_atlas_object.m`](./Brainnetome_create_atlas_object.m)
and the region-lookup utility
[`get_brainnetome_regions_by_name.m`](./get_brainnetome_regions_by_name.m).

## Primary reference

Fan, L., Li, H., Zhuo, J., Zhang, Y., Wang, J., Chen, L., Yang, Z., et
al. (2016). *The Human Brainnetome Atlas: a new brain atlas based on
connectional architecture.* **Cerebral Cortex, 26**(8), 3508–3526.
[doi:10.1093/cercor/bhw157](https://doi.org/10.1093/cercor/bhw157)

Local PDF: [`Fan_2016_CerebralCortex.pdf`](./Fan_2016_CerebralCortex.pdf).

## Key images

[`visualize_contents.m`](./visualize_contents.m) writes axial+sagittal
montage and 3-D isosurface PNGs of the Brainnetome atlas into
`png_images/`.

## How to load

Use the CANlab Core
[`load_atlas`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Data_extraction/load_atlas.m)
keyword:

```matlab
atl = load_atlas('brainnetome');   % Brainnetome_atlas_object.mat
```

Or load the `.mat` directly:

```matlab
S = load('Brainnetome_atlas_object.mat');
atl = S.atlas_obj;
```

Or read one of the raw NIfTIs:

```matlab
obj = fmri_data('BN_Atlas_274_plus_brainstem_r280_wani.nii.gz');
```

Look up a region by its abbreviated name:

```matlab
r = get_brainnetome_regions_by_name('A4hf');  % e.g. area 4 head/face
```

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `Brainnetome_atlas_object.mat` | MAT (`atlas`) | CANlab atlas object. `load_atlas('brainnetome')`. |
| `Brainnetome_atlas_regions.mat` | MAT (`region`) | Per-region `region` array. |
| `Brainnetome_create_atlas_object.m` | MATLAB | Constructor script. |
| `BN_Atlas_274_noCb_uint16.nii.gz` | NIfTI | Original Brainnetome 246-cortical + 28-subcortical (no cerebellum). |
| `BN_Atlas_274_plus_brainstem_r280_wani.nii.gz` | NIfTI | Brainnetome extended with CANlab "Wani" brainstem labels (~280 regions). |
| `BN_Atlas_274_T1_HCP40_MNI_1.25mm.nii.gz` | NIfTI | 1.25 mm T1 template (HCP40) the labels were drawn on. |
| `cluster_names.mat` | MAT | Region-name lookup table. |
| `get_brainnetome_regions_by_name.m` | MATLAB | Helper to retrieve regions by abbreviated label. |
| `get_brainnetome_regions_by_name_copy.m` | MATLAB | Backup copy of the helper. |
| `Fan_2016_CerebralCortex.pdf` | PDF | Primary reference (open access). |
| `visualize_contents.m` | MATLAB | Writes `png_images/`. |

## Citations

- Fan L, Li H, Zhuo J, et al. (2016). The Human Brainnetome Atlas: a
  new brain atlas based on connectional architecture. *Cereb Cortex*
  26:3508–3526.
  [doi:10.1093/cercor/bhw157](https://doi.org/10.1093/cercor/bhw157)
- Eickhoff SB, Yeo BTT, Genon S (2018). Imaging-based parcellations of
  the human brain. *Nat Rev Neurosci* 19:672–686.
  [doi:10.1038/s41583-018-0071-7](https://doi.org/10.1038/s41583-018-0071-7)
