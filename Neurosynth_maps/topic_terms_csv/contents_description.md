# Per-topic term-frequency CSVs (Neurosynth v4 100-topic LDA)

## Overview

One CSV per topic from the Neurosynth **v4 100-topic LDA model**, used
as the "document" inputs to the OpenAI `text-embedding-3-large`
embedding pipeline (see
[`../scripts/generate_topic_embeddings.py`](../scripts/generate_topic_embeddings.py)
and
[`../scripts/codex_embedding_instructions.txt`](../scripts/codex_embedding_instructions.txt)).
Each CSV is named `<topic_id>_<short_label>.csv`, where the integer
`topic_id` matches the topic index in `v4-topics-100_<id>_<terms>_*`
NIfTI filenames in
[`../2016_Neurosynth_100_topics/`](../2016_Neurosynth_100_topics/).
The file contains a single `Term` column listing the ~40 highest-
probability terms for that topic, one per row.

Only a subset of the 100 LDA topics survives the CANlab curation step
(topics dominated by methods / nuisance terms are dropped), so the
filenames here are sparse over `[0, 100)`.

This folder is **CSV-only** — no images, no MATLAB scripts.

## Primary reference

Derived from the
[Neurosynth v4 LDA-100 topic model](https://github.com/neurosynth/neurosynth-data)
as part of the CANlab build pipeline; the underlying corpus and topic
model are described in Yarkoni et al. 2011 *Nat Methods*.

## How to use

Read any one CSV with MATLAB's `readtable` or Python's `pandas`:

```matlab
T = readtable(which('60_Emotion_Processing.csv'));   % 1-column 'Term' table
top_terms = T.Term;
```

```python
import pandas as pd
terms = pd.read_csv('60_Emotion_Processing.csv')['Term'].tolist()
```

Each topic's CSV becomes a "document" that is fed to OpenAI's
`text-embedding-3-large` endpoint; the resulting vectors are saved
into
[`../topic_embeddings_text-embedding-3-large.csv`](../topic_embeddings_text-embedding-3-large.csv)
and the matching MAT-file.

## File inventory

49 CSVs of the form `<id>_<Label>.csv`. The numeric prefix matches
the topic index in the v4-100 model and in the matching NIfTI
filenames. Representative examples:

| File | Topic |
| --- | --- |
| `0_Sensory_Stimulation.csv` | Sensory / somatosensory stimulation, TMS, tDCS. |
| `13_Executive_Function.csv` | Executive control / cognitive control. |
| `60_Emotion_Processing.csv` | Emotion processing. |
| `61_Pain_Perception.csv` | Pain perception. |
| `68_Working_Memory.csv` | Working memory. |
| `97_Fear_conditioning.csv` | Fear conditioning. |
| `99_Reward_processing.csv` | Reward processing. |

Run `ls Neurosynth_maps/topic_terms_csv` for the full list.

## Citations

- Yarkoni T, Poldrack RA, Nichols TE, Van Essen DC, Wager TD (2011).
  Large-scale automated synthesis of human functional neuroimaging
  data. *Nat Methods* 8:665–670.
  [doi:10.1038/nmeth.1635](https://doi.org/10.1038/nmeth.1635)
- Poldrack RA, Mumford JA, Schonberg T, Kalar D, Barman B, Yarkoni T
  (2012). Discovering relations between mind, brain, and mental
  disorders using topic mapping. *PLoS Comput Biol* 8:e1002707.
  [doi:10.1371/journal.pcbi.1002707](https://doi.org/10.1371/journal.pcbi.1002707)
- OpenAI `text-embedding-3-large` model — see
  [platform.openai.com/docs/guides/embeddings](https://platform.openai.com/docs/guides/embeddings).
