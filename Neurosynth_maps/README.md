# Neurosynth maps

Maps and topic-model outputs derived from the
[Neurosynth](https://neurosynth.org) automated meta-analytic
database. The headliner here is the **v4 100-topic LDA topic model**
(in the `2016_Neurosynth_100_topics/` sub-folder), but the folder also
holds several auxiliary resources: MKDA outputs, topic embeddings,
co-activation matrices, and helper scripts.

## Sub-folders

| Folder | Description |
| --- | --- |
| [2016_Neurosynth_100_topics](./2016_Neurosynth_100_topics/contents_description.md) | **Pilot** — Neurosynth v4 100-topic LDA model: 50 topics × 2 inference directions (pAgF + pFgA) at FDR q<0.01, plus combined `fmri_data` objects and GPT-5–derived topic labels. |
| [mkda](./mkda) | MKDA per-term meta-analytic kernel-density maps from the Neurosynth corpus. |
| [scripts](./scripts) | Helper scripts used to build and re-derive the above. |
| [topic_terms_csv](./topic_terms_csv) | Per-topic term-frequency CSVs. |

## Auxiliary files at this level

| File | What it is |
| --- | --- |
| `neurosynth_data_obj.mat` | Aggregate `fmri_data` object covering the Neurosynth term/topic maps. |
| `neurosynth_topics_grouped_by_clique.mat` | Topics grouped by maximal-clique structure. |
| `neurosynth_interregion_coactivation_canlab2018_2mm.mat` | Inter-region co-activation matrix at 2 mm using the 2018 CANlab atlas. |
| `canlab_combined_atlas_2018_parcel_means_neurosynth.mat` | Parcel-mean Neurosynth maps under the 2018 combined atlas. |
| `canlab_combined_atlas_2018_resorted_brainpathway_obj.mat` | Brain-pathway object derived from those parcel means. |
| `topic_embeddings_text-embedding-3-large.{csv,mat}` | Per-topic vector embeddings (OpenAI text-embedding-3-large). |
| `current_data.tar.gz` | Compressed snapshot of the raw Neurosynth data used. |
| `assign_unique_cliques_and_average.m`, `assign_unique_cliques_from_maximal_cliques.m`, `maximalCliques.m` | Clique-finding helpers used to group similar topics. |
| `neurosynth_interregion_coactivation.m`, `neurosynth_seed_coactivation_map.m` | Build / use co-activation matrices. |
| `2026_Neurosynth_topics_LLM_vector_embedding.pptx` | Slide deck describing the LLM-embedding work. |

## Loading

The topic maps are registered in
[`load_image_set.m`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Data_extraction/load_image_set.m):

```matlab
[obj, ~, ~] = load_image_set('neurosynth_topics_forwardinference');
[obj, ~, ~] = load_image_set('neurosynth_topics_reverseinference');
[obj, ~, ~] = load_image_set('neurosynth');     % original term-level
```

To annotate an arbitrary brain map with its top Neurosynth topic:

```matlab
new_data = fmri_data('my_contrast.nii');
[~, topics] = neurosynth_feature_labels(new_data, 'topics_ri');
disp(topics{1}.Term_or_Topic_highest(1:5));
```

## Conventions

See the [docs README](../docs/README.md).
