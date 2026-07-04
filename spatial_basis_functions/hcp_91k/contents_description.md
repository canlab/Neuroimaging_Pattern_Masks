# HCP 91k-grayordinate spectral basis functions

## Overview

**Spectral (graph-Laplacian eigenmode) basis functions** in the HCP
91k-grayordinate space — 32k_fs_LR cortical surfaces plus 2 mm
subcortical voxels. The first 200 eigenmodes of each structure are
stacked into a single CIFTI dscalar file. The bases were generated for
a CANlab journal-club discussion of Pang et al. (2023) *Nature* on
geometric eigenmodes of brain dynamics; they provide a compact
orthonormal basis on which arbitrary brain maps in the HCP
grayordinate space can be projected.

> See [`README`](./README) for the authoritative methods write-up,
> including the choice of inverse-Euclidean edge weights, the
> 26-connectivity adjacency used for the subcortical graph, and
> comparison notes against Pang et al.'s finite-element approach.

**Primary reference (motivating paper).** Pang, J. C., Aquino, K. M.,
Oldehinkel, M., Robinson, P. A., Fulcher, B. D., Breakspear, M., &
Fornito, A. (2023). *Geometric constraints on human brain function.*
**Nature, 618**(7965), 566–574.
[doi:10.1038/s41586-023-06098-1](https://doi.org/10.1038/s41586-023-06098-1)

The CANlab spectral bases here are an instructional re-implementation
using `spharapy` on the HCP S1200 midthickness surfaces from the
[`hcp_utils`](https://github.com/rmldj/hcp-utils) repository, with
subcortical eigenmodes computed via a hand-built `networkx` adjacency
graph (see `README`).

## Key images

The file is a CIFTI dscalar with 200 maps; for best results visualise
in the
[HCP Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench)
loaded onto the matching 32k_fs_LR midthickness/inflated surfaces from
`hcp_utils`.
[`visualize_contents.m`](./visualize_contents.m) writes a small subset
of mode-by-mode PNGs into `png_images/` for quick inspection from
MATLAB.

## How to load

Not registered in `load_image_set`. Load with CanlabCore's CIFTI
support (which wraps `cifti-matlab`):

```matlab
cii = cifti_read(which('spectral_bases_200.dscalar.nii'));
% cii.cdata is grayordinates x 200
```

Or, if you only need the volumetric subcortical portion, use Workbench
to split the CIFTI:

```bash
wb_command -cifti-separate spectral_bases_200.dscalar.nii COLUMN \
    -volume-all subctx_spectral_bases_200.nii.gz
```

## File inventory

| File | Type | What it is |
| --- | --- | --- |
| `spectral_bases_200.dscalar.nii` | CIFTI dscalar | 200 spectral (graph-Laplacian) eigenmodes on HCP 91k grayordinates. |
| `README` | text | Methods write-up (B. Petre, 7/30/23). |
| `visualize_contents.m` | MATLAB | Renders selected modes to `png_images/` (after CIFTI to volume split). |

## Citations

- Pang JC, Aquino KM, Oldehinkel M, Robinson PA, Fulcher BD,
  Breakspear M, Fornito A (2023). Geometric constraints on human
  brain function. *Nature* 618:566–574.
  [doi:10.1038/s41586-023-06098-1](https://doi.org/10.1038/s41586-023-06098-1)
- Graff K, Tansey R, Iaria G, Bray S (2022). spharapy: Spatial
  harmonic analysis in Python.
  [docs](https://spharapy.readthedocs.io/en/latest/)
- HCP S1200 midthickness 32k_fs_LR surfaces via
  [hcp_utils](https://github.com/rmldj/hcp-utils).
