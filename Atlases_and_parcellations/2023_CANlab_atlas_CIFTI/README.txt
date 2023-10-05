script_2023_wagerlab_combined_atlas_prep.sh should be run before matlab or other bash scripts

script_2023_wagerlab_combined_atlas_cifti.m will produce subcortical volumetric structures and a nifti 
atlas to accompany the subsequently generated CIFTI atlas. The parcel labels and order of the nifti atlas
will match the cifti atlas, unlike canlab2018. Some canlab2018 regions will be missing though due to
the way subcortical structures are segmented in grayordinate space.

script_2013_wagerlab_combined_atlas_cifti.sh will produce the cifti dlabel file using Diedrichsen's
version of the Glasser atlas. It also outputs a config file needed for using the canlab2023 volumetric
atlas in qsiprep

The Diedrichsen atlas files can be obtained here: https://github.com/DiedrichsenLab/fs_LR_32

The glasser atlas can also be found in the hcp_utils package here as the "MMP" atlas:
    https://rmldj.github.io/hcp-utils/

cifti visualization is best done using connectome workbench which you can get here:
    https://www.humanconnectome.org/software/get-connectome-workbench

Data files and surface models are stored separately for cifti files. There's no structure in my dlabel.nii atlas file for
instance. It can be visualized on any surface model you like. Sometimes it's helpful to visualize things on an inflated
surface to be able to get a better view into sulci. So in order to visualize any data in connectome workbench you also need
the surface mesh models to which your data maps. I've found the hcp_utils git repo helpful for this.

