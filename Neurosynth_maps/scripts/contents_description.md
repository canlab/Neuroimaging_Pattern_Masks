# Neurosynth build / analysis scripts

## Overview

Helper scripts used to build and analyse the CANlab Neurosynth data
products in `Neurosynth_maps/` — preparing the MKDA inputs from the
Neurosynth coordinate database, computing co-activation / connectivity
networks across CANlab parcels, generating per-topic OpenAI embeddings,
and producing montage figures of provincial / connector hubs.

This folder is **scripts-only** — there are no NIfTIs or `.mat` images
in it. The outputs live in sibling folders
(`../mkda/`, `../2016_Neurosynth_100_topics/`,
`../topic_terms_csv/`, plus the `.mat` files in `Neurosynth_maps/`).

## Primary reference

Derived from the [Neurosynth](https://neurosynth.org) coordinate
database (Yarkoni et al. 2011 *Nat Methods*) as part of the CANlab
build pipeline.

## How to use

These are stand-alone scripts intended to be opened and run / adapted
in MATLAB or Python. They reference paths inside `Neurosynth_maps/`
and assume CanlabCore + the Neurosynth `current_data.tar.gz` snapshot
are available.

## File inventory

| File | Language | What it is |
| --- | --- | --- |
| `ns_matlab_prep_MKDA.m` | MATLAB | Prepares the MKDA `SETUP.mat` / kernel inputs from the Neurosynth coordinate database. Outputs land in `../mkda/`. |
| `tor_script_prepdata_2022_stub_(use_older_database).m` | MATLAB | Older data-prep stub (kept for reproducibility); use the `2022` workflow above for new builds. |
| `generate_neurosynth_atlases` | Shell / pipeline | Driver that re-derives Neurosynth atlas outputs end-to-end. |
| `neurosynth_default_mode_analysis_tor_dec_2019.m` | MATLAB | Default-mode-network-focused Neurosynth analysis (Dec 2019). |
| `neurosynth_cluster_defModeA_by_connectivity.m` | MATLAB | Clusters default-mode subregion A by connectivity profile. |
| `neurosynth_viz_overall_prob_activation.m` | MATLAB | Visualises the overall activation-probability map. |
| `montage_tables_provincial_connector_hubs.m` | MATLAB | Builds montage tables splitting parcels into provincial vs connector hubs. |
| `topic_embedding_similarity_analysis.m` | MATLAB | Loads the OpenAI text-embedding-3-large vectors and computes / visualises pairwise cosine similarity. |
| `generate_topic_embeddings.py` | Python | Calls the OpenAI API to embed each per-topic term list (see `codex_embedding_instructions.txt`). |
| `codex_embedding_instructions.txt` | text | Prompt / spec for the OpenAI embedding pipeline. |
| `__pycache__/` | dir | Python bytecode cache (auto-generated). |

## Citations

- Yarkoni T, Poldrack RA, Nichols TE, Van Essen DC, Wager TD (2011).
  Large-scale automated synthesis of human functional neuroimaging
  data. *Nat Methods* 8:665–670.
  [doi:10.1038/nmeth.1635](https://doi.org/10.1038/nmeth.1635)
- Wager TD, Lindquist M, Kaplan L (2007). Meta-analysis of functional
  neuroimaging data. *SCAN* 2:150–158.
  [doi:10.1093/scan/nsm015](https://doi.org/10.1093/scan/nsm015)
- OpenAI `text-embedding-3-large` model — see
  [platform.openai.com/docs/guides/embeddings](https://platform.openai.com/docs/guides/embeddings).
