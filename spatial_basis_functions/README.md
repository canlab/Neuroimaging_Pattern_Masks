# Spatial basis functions

> Spatial basis functions should be graded orthogonal maps, as opposed
> to discrete parcellations.
>
> A handful are provided here for ease of use with CANlab tools, but
> for a more comprehensive or rigorous analysis you should look into
> NeuroMaps:
> [https://netneurolab.github.io/neuromaps/usage.html](https://netneurolab.github.io/neuromaps/usage.html).

Each sub-folder contains one set of orthogonal maps (a "basis") with
its NIfTI(s), `contents_description.md`, `visualize_contents.m`, and a
short methods write-up.

## Sub-folders

| Folder | Description | Loader keyword |
| --- | --- | --- |
| [margulies](./margulies/contents_description.md) | **Pilot** — Margulies et al. 2016 first principal gradient of cortical organisation (unimodal → transmodal). CANlab volumetric build via registration fusion. | `marg`, `transmodal`, `principalgradient`, `margfsl` |
| [transcriptomic_gradients](./transcriptomic_gradients) | First three principal components of the Allen Human Brain Atlas gene-expression data (Burt et al. / Anderson et al.). | `transcriptomic_gradients` |
| [hcp_91k](./hcp_91k) | HCP 91k-grayordinate spatial bases (CIFTI). | — |
| [hcp_groupICAs](./hcp_groupICAs) | HCP group-ICA spatial maps at multiple model orders. | — |
| [mitochondrial_profile_maps](./mitochondrial_profile_maps) | Mitochondrial-pathway profile maps (Mosharov / Picard 2025 *Nature*). | — |

## Conventions

See the [docs README](../docs/README.md). Each folder retains any
existing `README.md` verbatim; the `contents_description.md` adds the
standardised overview / inventory / citation section and links back to
that README.
