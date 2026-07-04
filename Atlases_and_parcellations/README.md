# Atlases and parcellations

Brain atlases — labelled regions, networks, and parcellations —
distributed as NIfTI files and CANlab
[`atlas`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40atlas/atlas.m)
objects. Each sub-folder hosts one published atlas (or a CANlab-built
volumetric projection of a surface atlas) with the data, a
[`contents_description.md`](./), a `visualize_contents.m`, and
(where redistributable) the primary-reference PDF.

For the general loading API see
[`CanlabCore/Data_extraction/load_atlas.m`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Data_extraction/load_atlas.m).

```matlab
atl = load_atlas('canlab2024');             % most recent CANlab combined atlas
atl = load_atlas('canlab2023');             % previous combined atlas
atl = load_atlas('glasser_fmriprep20');     % 360 cortical parcels (HCP MMP1)
atl = load_atlas('schaefer400');            % 400 cortical parcels (Schaefer)
atl = load_atlas('cit168');                 % subcortex (CIT168)
atl = load_atlas('julichbrain');            % cytoarchitectonic (Julich-Brain)
```

## Sub-folders

### Combined / whole-brain CANlab atlases

| Year | Folder | Description |
| --- | --- | --- |
| 2024 | [2024_CANLab_atlas](./2024_CANLab_atlas) | Latest CANlab combined atlas (cortex + subcortex + brainstem + cerebellum). |
| 2023 | [2023_CANLab_atlas](./2023_CANLab_atlas) | Previous combined atlas. |
| 2023 | [2023_CANlab_atlas_CIFTI](./2023_CANlab_atlas_CIFTI) | CIFTI / surface-grayordinate version. |
| 2018 | [2018_Wager_combined_atlas](./2018_Wager_combined_atlas) | Original Wager-lab combined atlas (Glasser cortex + Pauli BG + Bianciardi brainstem + Morel thalamus). |

### Cortical parcellations

| Year | Folder | Description |
| --- | --- | --- |
| 2018 | [2018_Schaefer_Yeo_multires_cortical_parcellation](./2018_Schaefer_Yeo_multires_cortical_parcellation) | Multi-resolution (100–1000 parcels) cortex parcellation. |
| 2016 | [2016_Glasser_Nature_HumanConnectomeParcellation](./2016_Glasser_Nature_HumanConnectomeParcellation/contents_description.md) | **HCP MMP1** — 360 multi-modal cortical parcels (volumetric registration-fusion build). |
| 2016 | [2016_Fan_Brainnetome_r273_parcellation](./2016_Fan_Brainnetome_r273_parcellation) | Brainnetome 273-region parcellation. |
| 2013 | [2013_Shen_Constable_NIMG_268_parcellation](./2013_Shen_Constable_NIMG_268_parcellation) | Shen 268-node functional parcellation. |
| 2012 | [2012_desikan_killiany_tourville](./2012_desikan_killiany_tourville) | DKT cortical atlas (FreeSurfer). |
| 2011 | [2011_Buckner_7networks](./2011_Buckner_7networks) | 7-network resting-state parcellation. |
| 2011 | [2011_Yeo_17networks](./2011_Yeo_17networks) | 17-network resting-state parcellation. |
| 2009 | [2009_destrieux](./2009_destrieux) | Destrieux gyral / sulcal cortical atlas. |
| 2006 | [2006_desikan_killiany](./2006_desikan_killiany) | Desikan-Killiany cortical atlas. |

### Subcortex

| Year | Folder | Description |
| --- | --- | --- |
| 2020 | [2020_Tian_subcortical_v1.1](./2020_Tian_subcortical_v1.1) | Tian multi-scale subcortical parcellation (S1–S4). |
| 2018 | [2018_CIT168_Reinf_Learn_v1.1.0](./2018_CIT168_Reinf_Learn_v1.1.0) | CIT168 reinforcement-learning subcortical atlas v1.1.0. |
| 2018 | [2018_CIT168_Reinf_Learn_v1.0.0](./2018_CIT168_Reinf_Learn_v1.0.0) | CIT168 v1.0.0 (legacy). |
| 2016 | [2016_CIT168_Amygdala_v1.0.3](./2016_CIT168_Amygdala_v1.0.3) | CIT168 amygdala subnuclei v1.0.3. |
| 2014 | [2014_Keuken_7T_subcortex](./2014_Keuken_7T_subcortex) | 7-T high-resolution subcortical atlas. |
| 2020 | [2020_iglesias_hypothalamus](./2020_iglesias_hypothalamus) | Iglesias hypothalamic subregion atlas. |
| 2018 | [2018_Iglesias_thalamic_reconstruction](./2018_Iglesias_thalamic_reconstruction) | Iglesias thalamic nuclei reconstruction. |

### Brainstem

| Year | Folder | Description |
| --- | --- | --- |
| 2023 | [2023_Bianciardi_BrainstemNavigatorV0.9](./2023_Bianciardi_BrainstemNavigatorV0.9) | Bianciardi Brainstem Navigator v0.9. |
| 2023 | [2023_harvard_aan_brainstem_atlas](./2023_harvard_aan_brainstem_atlas) | Harvard Ascending Arousal Network brainstem atlas. |
| 2023 | [2023_levinson_bari_limbic_brainstem_atlas](./2023_levinson_bari_limbic_brainstem_atlas) | Levinson / Bari limbic brainstem atlas. |
| 2019 | [2019_Kragel_PAG](./2019_Kragel_PAG) | Periaqueductal gray subregions. |

### Cytoarchitecture / mixed

| Year | Folder | Description |
| --- | --- | --- |
| 2020 | [2020_JulichBrain_v3.0.3](./2020_JulichBrain_v3.0.3) | Julich-Brain cytoarchitectonic atlas v3.0.3. |

### Targeted regions / pathways

| Year | Folder | Description |
| --- | --- | --- |
| 2019 | [2019_Cartmell_NAcCoreShell](./2019_Cartmell_NAcCoreShell) | NAcc core/shell. |
| 2019 | [2019_Wager_pain_pathways](./2019_Wager_pain_pathways) | Canonical ascending/descending pain pathways. |
| 2019 | [2019_Zheng_Pain_Connectivity_Modules](./2019_Zheng_Pain_Connectivity_Modules) | Pain-related connectivity modules. |
| 2022 | [2022_Hansen_PET_tracer_maps](./2022_Hansen_PET_tracer_maps) | Hansen PET-tracer neurotransmitter maps. |
| 2026 | [2026_Pourmajidian — Mitochondrial energetic capacity](./2026_Pourmajidian_mitochondrial%20energetic%20capacity_Map) | Mitochondrial energetic capacity map. |

### White-matter atlases

| Year | Folder | Description |
| --- | --- | --- |
| 2020 | [2020_Pandora_white_matter_atlas](./2020_Pandora_white_matter_atlas) | Pandora white-matter bundle atlas. |
| 2020 | [2020_Thiebaut_de_Schotten_white_matter_atlas](./2020_Thiebaut_de_Schotten_white_matter_atlas) | Thiebaut de Schotten white-matter atlas. |

### Utilities

- [`readme_and_summary/`](./readme_and_summary) — older descriptive
  materials / cross-atlas summary files.

## Conventions

See the [docs README](../docs/README.md). Each folder has the same
structure: `contents_description.md` (overview / loading / inventory /
citations), `visualize_contents.m` (regenerates `png_images/`), the
`.nii` / `.mat` data files, and (where redistributable) the primary-
reference PDF. Existing per-folder `README.md`s are preserved
verbatim — `contents_description.md` links to them.
