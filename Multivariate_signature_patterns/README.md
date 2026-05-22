# Multivariate signature patterns

Pre-trained multivariate fMRI brain patterns ("signatures", "neuromarkers")
that can be applied to new individual-subject contrast maps to generate
predicted scores (e.g., pain intensity, emotion rating, craving). Each
sub-folder contains one published signature with its weights, supporting
files, and a per-folder [`contents_description.md`](./).

See the [docs README](../docs/README.md) for the documentation conventions
used throughout the repo, and the
[CanlabCore Tools](https://github.com/canlab/CanlabCore) for the
object-oriented MATLAB toolbox these signatures are designed to be used
with.

## Applying signatures to new data

Two convenience wrappers live at the root of this folder:

- [`apply_all_signatures.m`](./apply_all_signatures.m) — applies a
  configurable set of signatures (NPS, SIIPS, PINES, VPS, Rejection,
  GSR, HR, …) to a cell array of `fmri_data` objects and returns a
  results table. Supports `similarity_metric` (dot-product /
  cosine / correlation), `image_scaling` (none / center / z-score /
  l2-norm), and a named `image_set` to choose which signature bundle
  to apply.

  ```matlab
  [SIG, sigtable] = apply_all_signatures(DATA_OBJ, ...
      'similarity_metric', 'cosine_similarity', ...
      'image_scaling', 'l2norm_images', ...
      'image_set', 'npsplus');
  ```

- [`apply_siips.m`](./apply_siips.m) — SIIPS1-specific wrapper that
  accepts wildcards, filename lists, or `fmri_data` objects and
  returns whole-signature responses plus local responses for each of
  the 44 FDR-thresholded SIIPS subregions. See the SIIPS folder
  ([`2017_Woo_SIIPS1/contents_description.md`](./2017_Woo_SIIPS1/contents_description.md))
  for examples.

Both helpers internally call CanlabCore's `apply_mask` /
`canlab_pattern_similarity` but are easier to drop into a batch
pipeline than wiring those up by hand.

For a single signature, the most direct call is:

```matlab
[obj, networknames, imagenames] = load_image_set('siips');     % keyword from load_image_set.m
new_data = fmri_data('my_contrast.nii');
sig_response = apply_mask(new_data, obj, ...
                         'pattern_expression', 'ignore_missing');
```

## Sub-folders

Every sub-folder below contains `contents_description.md` (per-study
overview, references, file inventory, and loading instructions),
`visualize_contents.m` (regenerates `png_images/`), and the signature
files themselves.

| Year | Study | Topic | CanlabCore keyword |
| --- | --- | --- | --- |
| 2011 | [Wager — Placebo prediction (J Neurosci)](./2011_Wager_JNeuro_placebo_prediction/contents_description.md) | Placebo response | (no keyword — load `.img` directly) |
| 2015 | [Chang — PINES (PLoS Biol)](./2015_Chang_PLoSBiology_PINES/contents_description.md) | Picture-induced negative affect | `pines` |
| 2015 | [Kragel — Emotion BPLS (SCAN)](./2015_Kragel_emotionClassificationBPLS/contents_description.md) | 7 emotion categories | `kragelemotion` |
| 2015 | [Woo — Romantic Rejection (Nat Comms)](./2015_Woo_NatureComms_Rejection/contents_description.md) | Social rejection | `rejection` |
| 2016 | [Eisenbarth — Autonomic GSR/HR (J Neurosci)](./2016_Eisenbarth_JNeuro_autonomic_patterns/contents_description.md) | Skin conductance, heart rate | `gsr`, `hr` |
| 2016 | [Krishnan — VPS (eLife)](./2016_Krishnan_eLife_VPS/contents_description.md) | Vicarious pain | `vps` |
| 2017 | [Ashar — Care / Distress (Neuron)](./2017_Ashar_care_distress/contents_description.md) | Empathic care vs distress | (load `.nii` directly) |
| 2017 | [Rosenberg — saCPM (Nat Neurosci)](./2017_Rosenberg_sustained_attention/contents_description.md) | Sustained attention (connectome) | (load `.mat` directly) |
| 2017 | [Woo — SIIPS1 (Nat Comms)](./2017_Woo_SIIPS1/contents_description.md) | Cerebral pain (beyond nociception) | `siips` |
| 2018 | [Kragel — MFC generalizability (Nat Neurosci)](./2018_Kragel_MFC_Generalizability/contents_description.md) | Pain × emotion × cognitive control in MFC | `kragel18`, `pain_cog_emo` |
| 2018 | [Reddan — Threat ImEx (Neuron)](./2018_Reddan_Threat_Conditioning_ImEx/contents_description.md) | CS+ vs CS− threat conditioning | `csplus` |
| 2019 | [Kragel — Emotion Schemas (Sci Adv)](./2019_Kragel_Emotion_Schemas/contents_description.md) | 20 emotion categories | `kragelschemas` |
| 2019 | [Lee — Back pain (PAIN)](./2019_Lee_JPain_backpain/contents_description.md) | Chronic-back-pain markers (S1, PCASL, HFHRV) | (load directly) |
| 2019 | [Matthewson/Woo — SCR pain (PAIN)](./2019_Matthewson_Woo_SCR_pain/contents_description.md) | Skin-conductance + pain | (see `Readme.rtf`) |
| 2019 | [Yu/Koban — Guilt (Cereb Cortex)](./2019_Yu_Koban_Guilt/contents_description.md) | Interpersonal guilt | `guilt` |
| 2020 | [Geuter — Pain PDM mediation (Cereb Cortex)](./2020_Geuter_pain_multivariate_mediation_PDM/contents_description.md) | Pain mediation directions | `pdm`, `pain_pdm` |
| 2020 | [Silvestrini/Rainville — aMCC pain × cognitive control (NeuroImage)](./2020_Silvestrini_Rainville_Pain_CogControl_interaction_aMCC/contents_description.md) | dACC pain / Stroop patterns | `stroop` |
| 2020 | [Van Oudenhove/Kragel — Somatovisceral (Nat Comms)](./2020_VanOudenhove_Kragel_somatovisceral_pain/contents_description.md) | Visceral vs somatic pain | (load `.mat` directly) |
| 2020 | [Zhou — General vicarious pain (eLife)](./2020_Zhou_general_vicarious_pain/contents_description.md) | NS / FE / general vicarious pain | (helper `load_zhouvps`) |
| 2021 | [Čeko — MPA2 multiaversive (Nat Neurosci 2022)](./2021_Ceko_MPA2_multiaversive/contents_description.md) | 5 aversive modalities | `mpa2`, `multiaversive` |
| 2021 | [van 't Hof — BASIC (Cereb Cortex)](./2021_vantHoff_BASIC_sexual_image_classifier/contents_description.md) | Sexual-image classifier | (load `.nii` directly) |
| 2021 | [Zhou — Subjective fear VIFS (Nat Comms)](./2021_Zhou_Subjective_Fear/contents_description.md) | Fear from video | (helper `load_vifs`) |
| 2022 | [Coll — Pain × money decision value](./2022_coll_pain_monetary_reward_decision_value/contents_description.md) | Pain / money / shock-intensity value | (load directly) |
| 2022 | [Koban — NCS Craving (Nat Neurosci 2023)](./2022_Koban_NCS_Craving/contents_description.md) | Drug & food craving | `ncs` |
| 2023 | [Speer — Brain Reward Signature (BRS)](./2023_Speer_Brain_Reward_Signature_BRS/contents_description.md) | Reward signature | (load `.nii` directly) |
| 2024 | [FEPS — Facial Expressions of Pain (eLife)](./2024_FEPS_Facial_Expressions_of_Pain_Signature/contents_description.md) | Facial-expression-based pain | (load `.nii` directly) |
| 2026 | [Açıl — Mentalizing: Self / Other (Nat Comms)](./2026_Acil_Mentalizing_Self_Other/contents_description.md) | Mentalizing about self vs. other (MS, MS-Self, MS-Other, MS-SvO) | `selfother` |
| 2026 | [Murillo — PiFoneM (J Pain)](./2026_Murillo_PiFoneM/contents_description.md) | Pain-induced fear of neck movement | `pifonem`, `fearofneckpain` |

## Conventions

- Each sub-folder has the same shape: `contents_description.md`,
  `visualize_contents.m`, NIfTI / `.mat` data, a PDF of the primary
  reference (where redistributable), and a `png_images/` directory
  with regenerated surface / montage / isosurface figures.
- The shared rendering helpers live in
  [`../docs/canlab_render_patterns.m`](../docs/canlab_render_patterns.m)
  (and `canlab_render_atlas.m`).
- To regenerate every signature's PNGs in one MATLAB session, run
  [`../docs/run_all_multivariate.m`](../docs/run_all_multivariate.m).
