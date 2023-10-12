## Julich Brain Atlas v3.0.3

The Julich Brain Atlas is an ongoing effort to develop a complete digital atlas of human brain cytoarchitecture
measured using histological methods. It is the basis of the SPM Anatomy Toolbox, but it has grown and developed
considerably since its original version. Canlab Tools have copies of earlier versions of this atlas, including 
the spm anatomy toolbox v2.2. This atlas should supercede those, and also supercedes the contemporaneous SPM
Anatomy Toolbox, which has fewer regions than this (e.g. CA2/3 are missing). I don't know why.

The most up to date version of the atlas is stored in a digital repository called EBRAINS. This version was found
here,

https://search.kg.ebrains.eu/?category=Dataset#d69b70e2-3002-4eaf-9c61-9c56f019bbc8

Despite decades of development, the atlas is not complete. Areas that still need work are filled in with "GapMaps".
These are removed from the visualizations in the png* folder, but are present in the atlas objects here. They
have no meaningful probabilities, so I manually set them to 0.2, which is the threshold used for generating the
*.dat fields of the atlas objects. If you use this atlas you should probably treat these GapMaps differently than
the other parcels.

## Reference spaces

The atlas is originally registered to colin27 space in MNI coordinates, which is a subject specific template (from
a guy named Colin), which was distributed with SPM96 but still remains useful for precision mapping applications
because subject specific templates are more detailed that group mean templates. Beginning with v3 the authors also
distribute this atlas in MNI152NLin2009cAsym space, used by fmriprep and qsiprep, and map cortical areas to a
surface model in fsaverage space. Although the subject specific atlas is highly detailed, and useful when registering
a micron scale histological stain, it is not representative. For instance, Colin has a somewhat squatter brain than
average, and MTL structures are more laterally displaced than in the group mean template, which will lead to systematic
errors in the mapping of regions if used naively. I have not even bothered including the Colin27 version of this
atlas here, but I have converted the MNI152NLin2009cAsym and provide the fsaverage5 surface files.

I have additionally created a mapping to MNI152NLin6Asym space, which is used by FSL, SPM12, and HCP data as a default
template in spatial normalization. The mapping was performed using transformations computed by Lead-DBS in conjunction
with Vladimir Fonov from the MNI. The transforms were obtained here,

https://figshare.com/articles/dataset/MNI_T1_6thGen_NLIN_to_MNI_2009b_NLIN_ANTs_transform/3502238

antsApplyTransform from the fmriprep 20.0.3 LTS singularity container was used to perform the alignment on the to 
probablistic versions of the parcels and was implemented by the script warp_to_MNI152NLin6Asym0.sh. This was in turn
invoked by warp_to_MNI152NLin6Asym.sh on the Dartmouth HPC system. The probablistic versions of the parcels were then used
to regenerate the atlas in MNI152NLin6Asym space, and the result is also provided in this directory.

## References

The histological atlas is a collaborative global effort and spans multiple teams and decades of work. Different
brain areas have different publications associated with them, but there are a number of specific publications
pertaining to the creation of this combined atlas in particular, listed in order of relevance,

* Amuts K, Mohlberg H, Bludau S, Zilles K. Julich-Brain: A 3D probablistic atlas of the human brain's cytoarchitecture. 
Science 369(6506) 988-992, 2020
* Amunts K, Zilles K. Architectonic Mapping of the Human Brain beyond Broadmann. Neuron 88(6), 1086-1107, 2015
* Eikhoff S, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K. A new SPM toolbox for combining probablistic
cytoarchitectonic maps and functional imaging data. NeuroImage 25(4), 1325-1335, 2005

For citations pertaining to specific parcellation or other details please refer to the Julich Brain Atlas website or 
this page, https://www.fz-juelich.de/en/inm/inm-7/resources/jubrain-anatomy-toolbox

## 
Bogdan Petre
10/12/23
