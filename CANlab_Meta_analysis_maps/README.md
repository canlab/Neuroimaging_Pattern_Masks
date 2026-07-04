# CANlab meta-analysis maps

Consensus brain-activation maps from meta-analyses of fMRI studies,
covering emotion, cognition, pain, placebo, mentalising, anxiety, and
more. Each sub-folder contains one meta-analysis with its result NIfTIs,
[`contents_description.md`](./), a `visualize_contents.m`, and (where
redistributable) the primary-reference PDF.

Most CANlab meta-analyses were generated with
[MKDA](https://github.com/canlab/Canlab_MKDA_MetaAnalysis) (Multilevel
Kernel Density Analysis); the result maps here are the resulting
activation-density / probability maps.

## Sub-folders

| Year | Folder | Topic |
| --- | --- | --- |
| 2024 | [Quah / Saggar — RDoC factor maps](./2024_Quah_Saggar_factor_maps/contents_description.md) | Data-driven latent factor maps across RDoC tasks. |
| 2021 | [Zunhammer — N=603 pain placebo](./2021_Zunhammer_n603_Pain_Placebo) | Mega-analytic placebo brain map. |
| 2018 | [Sha — Common networks of psychopathology (BiolPsych)](./2018_Sha_BiolPsych_Common_Networks_Psychopathology) | Shared abnormalities across psychiatric disorders. |
| 2018 | [Kraynak / Gianaros — Immune meta](./2018_Kraynak_Gianaros_immunemeta) | Brain-immune meta-analytic maps. |
| 2017 | [Ashar — Placebo review](./2017_Ashar_Placebo_Review) | Placebo review brain maps. |
| 2017 | [de la Vega — Neurosynth cortical parcellation (Cereb Cortex)](./2017_delaVega_Neurosynth_cortical_parcellation) | Neurosynth-derived cortical parcellation. |
| 2016 | [Neurosynth / Wager — Social-affective](./2016_Neurosynth_Wager_SocAffective) | Social/affective Neurosynth maps. |
| 2016 | [Pauli — Basal-ganglia parcels (PNAS)](./2016_Pauli_Basal_Ganglia_Parcels) | Basal-ganglia parcellation. |
| 2016 | [de la Vega — Neurosynth mPFC parcellation (J Neurosci)](./2016_delaVega_JN_neurosynth-mfc_parcellation) | mPFC parcellation. |
| 2015 | [Wager / Kang — Emotion meta (PLoS Comp Biol)](./2015_Wager_Kang_etal_Emotion_Meta_BSPP) | Emotion meta-analysis (BSPP). |
| 2015 | [Satpute / Barrett — Emotion valence](./2015_Satpute_Barett_Emotion_Valence) | Emotion-valence meta-analysis. |
| 2014 | [Buhle / Silvers — Reappraisal meta (Cereb Cortex)](./2014_BuhleSilvers_Reappraisal_Meta) | Cognitive-reappraisal meta-analysis. |
| 2012 | [Denny — SOMA self/other (J Cog Neuro)](./2012_Denny_SOMA_self_other_meta) | Self vs other meta-analysis. |
| 2011 | [Yarkoni — Neurosynth (original; Nat Methods)](./2011_Yarkoni_Neurosynth_Original) | Original Neurosynth dataset/maps. |
| 2011 | [Meissner — Placebo masks (J Neurosci)](./2011_Meissner_Placebo_meta_analysis_masks) | Placebo meta-analytic masks. |
| 2011 | [Agency meta-analysis](./2011_Agency_Meta_analysis) | Sense of agency meta-analysis. |
| 2008 | [Kober — Emotion 163 (NeuroImage)](./2008_Kober_Emotion_163_studies) | Emotion meta of 163 studies. |
| 2007 | [Nee — Inhibition 47](./2007_Nee_Inhibition_47_studies) | Response-inhibition meta. |
| 2007 | [Etkin — Anxiety disorders (AJP)](./2007_Etkin_AJP_Anxiety_disorders) | Anxiety-disorder meta. |
| 2004 | [Wager — Attention switching 31](./2004_Wager_Attention_switching_31_studies) | Attention-switching meta. |
| 2003 | [Wager — Working memory 60 (CABN)](./2003_Wager_Working_memory_60_studies) | Working-memory meta. |
| 2003 | [Wager — Emotion 64 (NeuroImage)](./2003_Wager_Emotion_64_studies) | Emotion meta of 64 studies. |

Auxiliary files at this level:

- `Description_of_images.doc` — legacy summary of the image set.
- `scripts_summary/` — analysis-script snippets and helpers.

## Conventions

See the [docs README](../docs/README.md). Each folder follows the
standard layout (`contents_description.md`, `visualize_contents.m`,
data files, PDF, `png_images/`). Existing `README.md` files in
sub-folders are preserved verbatim and linked from
`contents_description.md`.
