script_2023_wagerlab_combined_atlas_prep.sh should be run before matlab or other bash scripts

script_2023_wagerlab_combined_atlas_cifti.m will produce subcortical volumetric structures and two nifti 
atlases to accompany the subsequently generated CIFTI atlas. The parcel labels and order of the nifti atlases
will match the cifti atlas, unlike canlab2018. Some canlab2018 regions will be missing though due to
the way subcortical structures are segmented in grayordinate space. One cifti atlas will be 1mm, the other
will be 2mm.

script_2013_wagerlab_combined_atlas_cifti.sh does three things,
1) produces the cifti dlabel file using Diedrichsen's version of the Glasser atlas. The subcortical structures
   in the CIFTI dlabel file should perfectly match the subcortical structures of the 2mm volumetric atlas.
2) It also outputs a config file needed for using the canlab2023 volumetric atlas in qsiprep. 
3) It creates a modified copie of the 2023 1mm volumetric atlas that lacks sform data. This is required for
   use of this atlas as a custom atlas in qsirecon.

The Diedrichsen atlas files can be obtained here: https://github.com/DiedrichsenLab/fs_LR_32

The glasser atlas can also be found in the hcp_utils package here as the "MMP" atlas:
    https://rmldj.github.io/hcp-utils/

cifti visualization is best done using connectome workbench which you can get here:
    https://www.humanconnectome.org/software/get-connectome-workbench

Data files and surface models are stored separately for cifti files. There's no structure in my dlabel.nii atlas file for
instance. It can be visualized on any surface model you like. Sometimes it's helpful to visualize things on an inflated
surface to be able to get a better view into sulci. So in order to visualize any data in connectome workbench you also need
the surface mesh models to which your data maps. I've found the hcp_utils git repo helpful for this.

Note: 2mm and 1mm atlases have different orientation. This is because the 2mm atlas is designed for use with 
HCP grayordinate files which are oriented right-left, posterior-anterior and inferior-superior (LAS+), while the 1mm
files are designed for use with qsirecon, which requires right-left, anterior-posterior and infermor-superior (LPS+)
orientation. See here,
https://qsiprep.readthedocs.io/en/latest/reconstruction.html
Most software should handle this under the hood, but be aware all the same in case you get funny results when using these.
