# Neuroimaging Pattern Masks — Documentation

This `docs/` directory is the documentation hub for the
[`canlab/Neuroimaging_Pattern_Masks`](https://github.com/canlab/Neuroimaging_Pattern_Masks)
repository. It indexes every study folder, links to per-folder documentation
(`contents_description.md`), and explains how to load each map with the
[CANlab Core Tools](https://github.com/canlab/CanlabCore).

> GitHub renders this page automatically when you open the `docs/` folder.
> The file is also compatible with GitHub Pages (Jekyll) — to enable a public
> site, turn on Pages with source `docs/` in repository settings.

For the project-level overview, see the [top-level README](../README.md).

---

## What is in this repository

The repository collects three kinds of brain maps, each with its own
top-level landing page:

| Type | Top-level folder | Landing page |
| --- | --- | --- |
| **Multivariate signatures** — pre-trained predictive patterns | [`Multivariate_signature_patterns/`](../Multivariate_signature_patterns) | [README](../Multivariate_signature_patterns/README.md) |
| **Atlases & parcellations** — labeled regions / networks | [`Atlases_and_parcellations/`](../Atlases_and_parcellations) | [README](../Atlases_and_parcellations/README.md) |
| **CANlab meta-analyses** — consensus maps across studies | [`CANlab_Meta_analysis_maps/`](../CANlab_Meta_analysis_maps) | [README](../CANlab_Meta_analysis_maps/README.md) |
| **Neurosynth maps** | [`Neurosynth_maps/`](../Neurosynth_maps) | [README](../Neurosynth_maps/README.md) |
| **Individual-study maps** | [`Individual_study_maps/`](../Individual_study_maps) | [README](../Individual_study_maps/README.md) |
| **Spatial basis functions** — continuous gradients / ICA components | [`spatial_basis_functions/`](../spatial_basis_functions) | [README](../spatial_basis_functions/README.md) |
| **Templates** — reference T1 templates | [`templates/`](../templates) | [README](../templates/README.md) |

A web companion with figures and references is hosted at the
[CANlab Brain Patterns site](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home).

## Loading maps

Most maps can be loaded by keyword through
[CanlabCore](https://github.com/canlab/CanlabCore):

```matlab
% Multivariate signatures (returns fmri_data + names)
[obj, networknames, imagenames] = load_image_set('siips');     % SIIPS1
[obj, networknames, imagenames] = load_image_set('npsplus');   % NPS, SIIPS, PINES, Rejection, VPS, GSR...
[obj, networknames, imagenames] = load_image_set('painsig');   % NPS + SIIPS only

% Atlases (returns atlas object)
atl = load_atlas('canlab2024');
atl = load_atlas('canlab2023');
atl = load_atlas('glasser_fmriprep20');
atl = load_atlas('schaefer400');
```

See `load_image_set.m` and `load_atlas.m` in CanlabCore for the full keyword list.

## Applying signatures to new data

The repo ships two convenience wrappers at the root of
[`Multivariate_signature_patterns/`](../Multivariate_signature_patterns):

- [`apply_all_signatures.m`](../Multivariate_signature_patterns/apply_all_signatures.m)
  — applies a configurable set of signatures (NPS, SIIPS, PINES, VPS,
  Rejection, GSR, HR, …) to a cell array of `fmri_data` objects and
  returns a results table. Supports `similarity_metric` (dot-product /
  cosine / correlation), `image_scaling` (none / center / z-score /
  l2-norm), and a named `image_set` to choose which signature bundle
  to apply. Example:

  ```matlab
  [SIG, sigtable] = apply_all_signatures(DATA_OBJ, ...
      'similarity_metric', 'cosine_similarity', ...
      'image_scaling', 'l2norm_images', ...
      'image_set', 'npsplus');
  ```

- [`apply_siips.m`](../Multivariate_signature_patterns/apply_siips.m)
  — SIIPS1-specific wrapper that accepts wildcards, filename lists, or
  `fmri_data` objects and returns whole-signature responses plus local
  responses for each of the 44 FDR-thresholded SIIPS subregions.

Both helpers internally call CanlabCore's `apply_mask` /
`canlab_pattern_similarity`, but they are easier to drop into a batch
pipeline than wiring those up by hand.

## Per-study documentation

Each study folder contains:

- A primary-reference PDF named `<lastauthor>_<year>_<keywords>.pdf` (where available).
- A `contents_description.md` that describes the folder (overview, key images, loading instructions, file inventory, citations).
- A `visualize_contents.m` MATLAB script that regenerates cortical-surface, isosurface (where relevant), and montage figures into `png_images/`.
- The image files themselves (NIfTI / `atlas` `.mat` objects / signature `.mat` objects).
- Any existing author-authored `README.md` is preserved verbatim and linked from `contents_description.md`.

Shared rendering helpers used by every per-folder script live in
[`docs/`](./):

- [`canlab_render_patterns.m`](./canlab_render_patterns.m) — surface + montage + isosurface for `fmri_data` patterns.
- [`canlab_render_atlas.m`](./canlab_render_atlas.m) — montage + isosurface for `atlas` objects.
- [`run_all_multivariate.m`](./run_all_multivariate.m) — batch-re-renders every folder in a chosen category.

---

### Multivariate signatures — all entries documented

| Year | Study | Topic |
| --- | --- | --- |
| 2011 | [Wager — Placebo prediction (J Neurosci)](../Multivariate_signature_patterns/2011_Wager_JNeuro_placebo_prediction/contents_description.md) | Placebo response |
| 2015 | [Chang — PINES (PLoS Biol)](../Multivariate_signature_patterns/2015_Chang_PLoSBiology_PINES/contents_description.md) | Picture-induced negative affect |
| 2015 | [Kragel — Emotion BPLS (SCAN)](../Multivariate_signature_patterns/2015_Kragel_emotionClassificationBPLS/contents_description.md) | 7 emotion categories |
| 2015 | [Woo — Romantic Rejection (Nat Comms)](../Multivariate_signature_patterns/2015_Woo_NatureComms_Rejection/contents_description.md) | Social rejection |
| 2016 | [Eisenbarth — Autonomic GSR/HR (J Neurosci)](../Multivariate_signature_patterns/2016_Eisenbarth_JNeuro_autonomic_patterns/contents_description.md) | Skin conductance, heart rate |
| 2016 | [Krishnan — VPS (eLife)](../Multivariate_signature_patterns/2016_Krishnan_eLife_VPS/contents_description.md) | Vicarious pain |
| 2017 | [Ashar — Care / Distress (Neuron)](../Multivariate_signature_patterns/2017_Ashar_care_distress/contents_description.md) | Empathic care vs distress |
| 2017 | [Rosenberg — saCPM (Nat Neurosci)](../Multivariate_signature_patterns/2017_Rosenberg_sustained_attention/contents_description.md) | Sustained attention (connectome) |
| 2017 | [Woo — SIIPS1 (Nat Comms)](../Multivariate_signature_patterns/2017_Woo_SIIPS1/contents_description.md) | Cerebral pain (beyond nociception) |
| 2018 | [Kragel — MFC generalizability (Nat Neurosci)](../Multivariate_signature_patterns/2018_Kragel_MFC_Generalizability/contents_description.md) | Pain × emotion × cognitive control |
| 2018 | [Reddan — Threat ImEx (Neuron)](../Multivariate_signature_patterns/2018_Reddan_Threat_Conditioning_ImEx/contents_description.md) | CS+ vs CS− threat conditioning |
| 2019 | [Kragel — Emotion Schemas (Sci Adv)](../Multivariate_signature_patterns/2019_Kragel_Emotion_Schemas/contents_description.md) | 20 emotion categories |
| 2019 | [Lee — Back pain (PAIN)](../Multivariate_signature_patterns/2019_Lee_JPain_backpain/contents_description.md) | Chronic-back-pain markers |
| 2019 | [Matthewson/Woo — SCR pain (PAIN)](../Multivariate_signature_patterns/2019_Matthewson_Woo_SCR_pain/contents_description.md) | Skin-conductance + pain |
| 2019 | [Yu/Koban — Guilt (Cereb Cortex)](../Multivariate_signature_patterns/2019_Yu_Koban_Guilt/contents_description.md) | Interpersonal guilt |
| 2020 | [Geuter — Pain PDM mediation (Cereb Cortex)](../Multivariate_signature_patterns/2020_Geuter_pain_multivariate_mediation_PDM/contents_description.md) | Pain mediation directions |
| 2020 | [Silvestrini/Rainville — aMCC pain × cognitive control (NeuroImage)](../Multivariate_signature_patterns/2020_Silvestrini_Rainville_Pain_CogControl_interaction_aMCC/contents_description.md) | dACC pain / Stroop |
| 2020 | [Van Oudenhove/Kragel — Somatovisceral (Nat Comms)](../Multivariate_signature_patterns/2020_VanOudenhove_Kragel_somatovisceral_pain/contents_description.md) | Visceral vs somatic |
| 2020 | [Zhou — General vicarious pain (eLife)](../Multivariate_signature_patterns/2020_Zhou_general_vicarious_pain/contents_description.md) | NS / FE / general vicarious pain |
| 2021 | [Čeko — MPA2 multiaversive (Nat Neurosci 2022)](../Multivariate_signature_patterns/2021_Ceko_MPA2_multiaversive/contents_description.md) | 5 aversive modalities |
| 2021 | [van 't Hof — BASIC (Cereb Cortex)](../Multivariate_signature_patterns/2021_vantHoff_BASIC_sexual_image_classifier/contents_description.md) | Sexual-image classifier |
| 2021 | [Zhou — Subjective fear VIFS (Nat Comms)](../Multivariate_signature_patterns/2021_Zhou_Subjective_Fear/contents_description.md) | Subjective fear |
| 2022 | [Coll — Pain × money decision value](../Multivariate_signature_patterns/2022_coll_pain_monetary_reward_decision_value/contents_description.md) | Pain / money / shock value |
| 2022 | [Koban — NCS Craving (Nat Neurosci 2023)](../Multivariate_signature_patterns/2022_Koban_NCS_Craving/contents_description.md) | Drug & food craving |
| 2023 | [Speer — Brain Reward Signature (BRS)](../Multivariate_signature_patterns/2023_Speer_Brain_Reward_Signature_BRS/contents_description.md) | Reward signature |
| 2024 | [FEPS — Facial Expressions of Pain (eLife)](../Multivariate_signature_patterns/2024_FEPS_Facial_Expressions_of_Pain_Signature/contents_description.md) | Facial-expression-based pain |
| 2026 | [Murillo — PiFoneM (J Pain)](../Multivariate_signature_patterns/2026_Murillo_PiFoneM/contents_description.md) | Pain-induced fear / negative-affect |

### Atlases and parcellations — all entries documented

| Year | Folder | Description |
| --- | --- | --- |
| 2006 | [Desikan-Killiany](../Atlases_and_parcellations/2006_desikan_killiany/contents_description.md) | Desikan-Killiany cortical atlas |
| 2009 | [Destrieux](../Atlases_and_parcellations/2009_destrieux/contents_description.md) | Destrieux gyral/sulcal cortical atlas |
| 2011 | [Buckner 7networks](../Atlases_and_parcellations/2011_Buckner_7networks/contents_description.md) | 7 resting-state networks |
| 2011 | [Yeo 17networks](../Atlases_and_parcellations/2011_Yeo_17networks/contents_description.md) | 17 resting-state networks |
| 2012 | [Desikan-Killiany-Tourville (DKT)](../Atlases_and_parcellations/2012_desikan_killiany_tourville/contents_description.md) | DKT cortical atlas |
| 2013 | [Shen/Constable 268](../Atlases_and_parcellations/2013_Shen_Constable_NIMG_268_parcellation/contents_description.md) | 268-node functional parcellation |
| 2014 | [Keuken 7T subcortex](../Atlases_and_parcellations/2014_Keuken_7T_subcortex/contents_description.md) | 7T high-resolution subcortex |
| 2016 | [CIT168 amygdala v1.0.3](../Atlases_and_parcellations/2016_CIT168_Amygdala_v1.0.3/contents_description.md) | Amygdala subnuclei |
| 2016 | [Fan Brainnetome 273](../Atlases_and_parcellations/2016_Fan_Brainnetome_r273_parcellation/contents_description.md) | Brainnetome parcellation |
| 2016 | [Glasser HCP-MMP1](../Atlases_and_parcellations/2016_Glasser_Nature_HumanConnectomeParcellation/contents_description.md) | 360 multimodal cortical parcels (pilot) |
| 2018 | [CIT168 reinf-learn v1.0.0](../Atlases_and_parcellations/2018_CIT168_Reinf_Learn_v1.0.0/contents_description.md) | CIT168 subcortical v1.0.0 |
| 2018 | [CIT168 reinf-learn v1.1.0](../Atlases_and_parcellations/2018_CIT168_Reinf_Learn_v1.1.0/contents_description.md) | CIT168 subcortical v1.1.0 |
| 2018 | [Iglesias thalamic](../Atlases_and_parcellations/2018_Iglesias_thalamic_reconstruction/contents_description.md) | Iglesias thalamic nuclei |
| 2018 | [Schaefer/Yeo multires](../Atlases_and_parcellations/2018_Schaefer_Yeo_multires_cortical_parcellation/contents_description.md) | 100–1000 parcel multi-resolution cortex |
| 2018 | [Wager combined atlas](../Atlases_and_parcellations/2018_Wager_combined_atlas/contents_description.md) | Original CANlab whole-brain composite |
| 2019 | [Cartmell NAc core/shell](../Atlases_and_parcellations/2019_Cartmell_NAcCoreShell/contents_description.md) | NAcc core / shell |
| 2019 | [Kragel PAG](../Atlases_and_parcellations/2019_Kragel_PAG/contents_description.md) | Periaqueductal gray subregions |
| 2019 | [Wager pain pathways](../Atlases_and_parcellations/2019_Wager_pain_pathways/contents_description.md) | Canonical ascending/descending pain |
| 2019 | [Zheng pain connectivity modules](../Atlases_and_parcellations/2019_Zheng_Pain_Connectivity_Modules/contents_description.md) | Pain connectivity modules |
| 2020 | [Iglesias hypothalamus](../Atlases_and_parcellations/2020_iglesias_hypothalamus/contents_description.md) | Hypothalamic subregions |
| 2020 | [JulichBrain v3.0.3](../Atlases_and_parcellations/2020_JulichBrain_v3.0.3/contents_description.md) | Cytoarchitectonic |
| 2020 | [Pandora WM atlas](../Atlases_and_parcellations/2020_Pandora_white_matter_atlas/contents_description.md) | White-matter bundles |
| 2020 | [Thiebaut de Schotten WM](../Atlases_and_parcellations/2020_Thiebaut_de_Schotten_white_matter_atlas/contents_description.md) | White-matter atlas |
| 2020 | [Tian subcortical v1.1](../Atlases_and_parcellations/2020_Tian_subcortical_v1.1/contents_description.md) | Multi-scale subcortex |
| 2022 | [Hansen PET tracer maps](../Atlases_and_parcellations/2022_Hansen_PET_tracer_maps/contents_description.md) | Neurotransmitter PET maps |
| 2023 | [Bianciardi Brainstem Navigator v0.9](../Atlases_and_parcellations/2023_Bianciardi_BrainstemNavigatorV0.9/contents_description.md) | Brainstem nuclei |
| 2023 | [CANlab atlas (2023)](../Atlases_and_parcellations/2023_CANLab_atlas/contents_description.md) | CANlab combined atlas (2023) |
| 2023 | [CANlab atlas CIFTI](../Atlases_and_parcellations/2023_CANlab_atlas_CIFTI/contents_description.md) | CIFTI grayordinate version |
| 2023 | [Harvard AAN brainstem](../Atlases_and_parcellations/2023_harvard_aan_brainstem_atlas/contents_description.md) | Ascending arousal network |
| 2023 | [Levinson/Bari limbic brainstem](../Atlases_and_parcellations/2023_levinson_bari_limbic_brainstem_atlas/contents_description.md) | Limbic brainstem |
| 2024 | [CANlab atlas (2024)](../Atlases_and_parcellations/2024_CANLab_atlas/contents_description.md) | Latest CANlab combined atlas |
| 2026 | [Pourmajidian mitochondrial energetic capacity](../Atlases_and_parcellations/2026_Pourmajidian_mitochondrial%20energetic%20capacity_Map/contents_description.md) | Mitochondrial energetic capacity |
| util | [readme_and_summary](../Atlases_and_parcellations/readme_and_summary/contents_description.md) | Cross-atlas summary materials |

### CANlab meta-analysis maps — all entries documented

| Year | Folder |
| --- | --- |
| 2003 | [Wager — Emotion 64 studies (NeuroImage)](../CANlab_Meta_analysis_maps/2003_Wager_Emotion_64_studies/contents_description.md) |
| 2003 | [Wager — Working memory 60 studies (CABN)](../CANlab_Meta_analysis_maps/2003_Wager_Working_memory_60_studies/contents_description.md) |
| 2004 | [Wager — Attention switching 31 studies](../CANlab_Meta_analysis_maps/2004_Wager_Attention_switching_31_studies/contents_description.md) |
| 2007 | [Etkin — Anxiety disorders (AJP)](../CANlab_Meta_analysis_maps/2007_Etkin_AJP_Anxiety_disorders/contents_description.md) |
| 2007 | [Nee — Inhibition 47 studies](../CANlab_Meta_analysis_maps/2007_Nee_Inhibition_47_studies/contents_description.md) |
| 2008 | [Kober — Emotion 163 studies (NeuroImage)](../CANlab_Meta_analysis_maps/2008_Kober_Emotion_163_studies/contents_description.md) |
| 2011 | [Agency meta-analysis](../CANlab_Meta_analysis_maps/2011_Agency_Meta_analysis/contents_description.md) |
| 2011 | [Meissner — Placebo masks (J Neurosci)](../CANlab_Meta_analysis_maps/2011_Meissner_Placebo_meta_analysis_masks/contents_description.md) |
| 2011 | [Yarkoni — Neurosynth (Nat Methods)](../CANlab_Meta_analysis_maps/2011_Yarkoni_Neurosynth_Original/contents_description.md) |
| 2012 | [Denny — SOMA self/other](../CANlab_Meta_analysis_maps/2012_Denny_SOMA_self_other_meta/contents_description.md) |
| 2014 | [Buhle/Silvers — Reappraisal meta](../CANlab_Meta_analysis_maps/2014_BuhleSilvers_Reappraisal_Meta/contents_description.md) |
| 2015 | [Satpute/Barrett — Emotion valence](../CANlab_Meta_analysis_maps/2015_Satpute_Barett_Emotion_Valence/contents_description.md) |
| 2015 | [Wager/Kang — Emotion meta BSPP](../CANlab_Meta_analysis_maps/2015_Wager_Kang_etal_Emotion_Meta_BSPP/contents_description.md) |
| 2016 | [de la Vega — Neurosynth mPFC parcellation](../CANlab_Meta_analysis_maps/2016_delaVega_JN_neurosynth-mfc_parcellation/contents_description.md) |
| 2016 | [Neurosynth/Wager — Social-affective](../CANlab_Meta_analysis_maps/2016_Neurosynth_Wager_SocAffective/contents_description.md) |
| 2016 | [Pauli — Basal-ganglia parcels (PNAS)](../CANlab_Meta_analysis_maps/2016_Pauli_Basal_Ganglia_Parcels/contents_description.md) |
| 2017 | [Ashar — Placebo review](../CANlab_Meta_analysis_maps/2017_Ashar_Placebo_Review/contents_description.md) |
| 2017 | [de la Vega — Neurosynth cortical parcellation](../CANlab_Meta_analysis_maps/2017_delaVega_Neurosynth_cortical_parcellation/contents_description.md) |
| 2018 | [Kraynak/Gianaros — Immune meta](../CANlab_Meta_analysis_maps/2018_Kraynak_Gianaros_immunemeta/contents_description.md) |
| 2018 | [Sha — Common networks of psychopathology](../CANlab_Meta_analysis_maps/2018_Sha_BiolPsych_Common_Networks_Psychopathology/contents_description.md) |
| 2021 | [Zunhammer — N=603 pain placebo](../CANlab_Meta_analysis_maps/2021_Zunhammer_n603_Pain_Placebo/contents_description.md) |
| 2024 | [Quah/Saggar — RDoC factor maps (Nat Comms 2025)](../CANlab_Meta_analysis_maps/2024_Quah_Saggar_factor_maps/contents_description.md) |

### Neurosynth maps

- [2016 Neurosynth 100 topics (v4/v5)](../Neurosynth_maps/2016_Neurosynth_100_topics/contents_description.md) — pilot
- [mkda/](../Neurosynth_maps/mkda/contents_description.md) — MKDA per-term kernel-density maps
- [scripts/](../Neurosynth_maps/scripts/contents_description.md) — Helper scripts
- [topic_terms_csv/](../Neurosynth_maps/topic_terms_csv/contents_description.md) — Per-topic term-frequency CSVs

### Individual-study maps

- [2024 Bo — Emotion regulation Bayes factor (Nat Neurosci)](../Individual_study_maps/2024_Bo_EmotionRegulation_BayesFactor/contents_description.md) — pilot
- [2026 Miao — Social/ToM Bayes factor](../Individual_study_maps/2026_Miao_Social-ToM_BayesFactor/contents_description.md)

### Spatial basis functions

- [Margulies — Principal cortical gradient (PNAS 2016)](../spatial_basis_functions/margulies/contents_description.md) — pilot
- [HCP 91k](../spatial_basis_functions/hcp_91k/contents_description.md)
- [HCP group ICAs](../spatial_basis_functions/hcp_groupICAs/contents_description.md)
- [Transcriptomic gradients](../spatial_basis_functions/transcriptomic_gradients/contents_description.md)
- [Mitochondrial profile maps (Mosharov 2025)](../spatial_basis_functions/mitochondrial_profile_maps/contents_description.md)

### Templates

- [templates landing page](../templates/README.md) — overview of MNI templates included.
- [cerebellum/](../templates/cerebellum/contents_description.md)
- [transforms/](../templates/transforms/contents_description.md)

---

## Conventions used in this documentation

Per the [contribution guidelines](../Steps_to_add_a_new_atlas.rtf), every study folder
follows the same layout. The `contents_description.md` in each folder is the
canonical entry point and contains:

1. **Overview** — what the resource is, in a few lines, with the primary reference
   linked at the publisher's site (and citations to any secondary references found
   in the folder).
2. **Key images** — 1–2 representative figures (montage, surface, or isosurface) taken
   from the CANlab Brain Patterns site or generated by `visualize_contents.m`.
   Pre-rendered figures already present in folders (e.g.,
   `images_and_regions/`, `figures/`, `png_images/`) are embedded directly.
3. **How to load** — the exact `load_atlas()` / `load_image_set()` call, or the
   filename to pass to `fmri_data()` / `atlas()` directly.
4. **Construction scripts** — links to constructor scripts in the folder or in CanlabCore.
5. **File inventory** — every NIfTI / `.mat` / `.m` in the folder with one line each.
   Individual `.png` files are *not* listed (they are embedded as key images instead).
6. **Citations** — primary and secondary papers, plus links.

Existing author-authored `README.md` / `Readme.md` / `Readme.txt` files in
sub-folders are preserved verbatim and linked from the new
`contents_description.md`, mirroring the Glasser-folder pattern.

### Regenerating figures

After modifying any folder, regenerate its PNGs:

```matlab
cd /path/to/Atlases_and_parcellations/2016_Glasser_Nature_HumanConnectomeParcellation
visualize_contents();   % writes png_images/<label>_{surface,montage,isosurface}.png
```

To re-render every signature in one MATLAB session:

```matlab
addpath('/path/to/Neuroimaging_Pattern_Masks/docs');
run_all_multivariate('Multivariate_signature_patterns');    % or any other category
```

### Adding a new study folder

1. Create the folder under the appropriate category with name `YYYY_<lastauthor>_<keywords>`.
2. Drop in the NIfTI / `.mat` files plus a copy of the primary reference PDF
   (filename: `<lastauthor>_<year>_<keywords>.pdf`).
3. Copy a `visualize_contents.m` from a similar folder and adapt the `imgs` cell.
4. Run it once locally to generate `png_images/`.
5. Write `contents_description.md` following the sections above.
6. Add the entry to the appropriate category list in the top-level
   landing page (e.g., `Atlases_and_parcellations/README.md`) and to
   the master list in this `docs/README.md`.
7. If applicable, add a keyword to the switch in
   `CanlabCore/Data_extraction/load_image_set.m` or `load_atlas.m`.

---

## Documentation status

Every study folder in the repository now has a `contents_description.md`
and (in nearly all cases) a `visualize_contents.m`. Top-level landing
pages list each subfolder and link to its docs.

**PDF coverage:** 30+ PDFs are checked in directly. For a few
paywalled papers (Reddan 2018 Neuron, Speer 2023 J Neurosci, Coll 2022,
Kragel 2015 SCAN) the DOI link is recorded but no local PDF is present.

**MATLAB / PNG coverage:** every `visualize_contents.m` writes its
output to `png_images/`. The full batch is reproducible by running
[`docs/run_all_multivariate.m`](./run_all_multivariate.m) once per
category.

If a folder you need is missing something, the underlying image files
and any pre-existing `README.md` / `create_*.m` script in that folder
remain authoritative — the new docs are an additional indexing /
discovery layer, not a replacement.
