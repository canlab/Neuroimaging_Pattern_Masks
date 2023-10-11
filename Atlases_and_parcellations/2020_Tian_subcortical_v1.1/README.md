## Tian 2020

This atlas was adapted from a different repo,

https://github.com/yetianmed/subcortex/tree/e80ee787732536e9e89534c9a623b10aff7928f4

Group-Parcellation files were copied here for posterity, but the main contribution of 
this repo is the contribution of canlabCore atlas objects.

This atlas was generated in a hierarchical fashion with finer and finer parcellations
at each hierarchy until a natural stopping condition in the parcellation algorithm
was reached (see Tian et al \[2020\] Nat Neuro for details). Different parcellations
are indicated by the "scale", and there are 4. Furthermore, This atlas is provided 
registered to two different reference templates for use with FSL/HCP or fMRIprep data
as the case may be, which is likewise specified by FSL6 or fmriprep20 (the latest LTS)
in the name.

## Atlas space

The authors of this atlas are professionally invested in deep brain stimulation 
treatments and require greater precision from their segmentations than is typical
in fMRI research. Consequently they make a point of specifying their reference spaces.
I've made efforts to retain this specificity here, so multiple versions of this atlas
are provided which correspond to multiple common reference spaces.

As of this writing, most wagerlab users are relying on fMRIprep for spatial normalization.
This aligns all data to the MNI152NLin2009cAsym reference space. To use this atlas please
use the Tian_3T_S4_2mm_fmriprep20 segmentation, and use the appropriate underlay
for visualization. Please see the "templates" folder of this repo for more details.

## 

Bogdan Petre <br />
10/10/23
