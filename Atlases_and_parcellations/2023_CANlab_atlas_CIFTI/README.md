## Overview

This directory contains an updated version of the canlab2018 atlas that is designed to be more general purpose.
Three equivalent atlases are provided, which differ in format but identify identical regions ordered in the
same way and indexed by the same values using the same labels.

1) 2mm volumetric
* canlab_2023_2mm.mat    - a drop in replacement for the workhorse canlab2018 atlas. Contains a canlabCore/@atlas object
* canlab_2023_2mm.nii.gz - same as the above, except in nifti format (no metadata)

2) 1mm volumetric
* canlab_2023_1mm.mat    - a higher resolution version of the above. The former were generated from this by downsampling.
                           This also forms the basis of the QSIprep compatable atlas, and therefore is in LPS+ orientation,
                           while the 2mm versions are in the LAS+ orientation which seems more common in this repo.
* canlab_2023_MNI152NLin2009cAsym_1mm_lps.nii.gz 
                         - same as the above, but nifti format (no metadata), and stripped of sform data too. This version
                           is compatable with qsirecon
* atlas_config.json      - qsirecon atlas config file.

3) 2mm grayordinate
* canlab_2023_2mm.dlable.nii
                         - CIFTI version of the above atlases.

## Canlab 2023 vs. Canlab 2018

* 2023 has more lateralized structures. Except for some ambiguous ones (e.g. thalamic midline), most structures that were previously bilateral (hippocampus, amygdala, subiculum, all thalamic structures, etc.) have been split into left and right lateralized regions
* 2023 omits white matter structures. In particular some cerebellar structures that do not overlap meaningfully with grayordinate space were removed. This ensures consistency across grayordinate and volumetric atlas versions.
* 2023 was designed for use with qsirecon in mind, and includes example configuration files

## QSIPrep use

QSIprep is to DWI analysis what fMRIPrep is to fMRI analysis and some. Whereas fMRIprep only performs preprocessing, QSIPrep also performs what might be a the DWI analog of a first level analysis. In the case of tractography, the question of interest is what regions are connected to what other regions. The regions are defined by atlas files. Several default atlases are provided with qsiprep, but none are more detailed or exhaustive than the canlab atlas. Thankfully qsiprep allows for use of custom atlases if appropriately formatted and accompanied by the necessary metadata (specified in a json file). To use the canlab 2023 atlas with qsiprep you need to pass a path to a directory containing the atlas and its atlas_config.json file in during the QSIprep run. If running qsiprep in a singularity container you do this by mounting this folder at a particular location within the virtual environment by using the --bind <full_atlas_directory_path>:/atlas/qsirecon_atlases. For use without a singularity container refer to the QSIprep docs,
https://qsiprep.readthedocs.io/en/latest/reconstruction.html

For the most basic usage (e.g. you only want to use the canlab2023 atlas in your qsirecon pipelines, and no other atlases) you can simply mount this directory in your singularity container and it should work out of the box. You will only need to create your own custom atlases folder and merge atlas_config.json with the stock atlas_config.json file if you wish to combine multiple atlases into a single reconstruction workflow. In the latter case you will need the stock qsiprep files. They're linked to in the qsiprep docs, but for your convenience they can be found here,
https://upenn.box.com/shared/static/8k17yt2rfeqm3emzol5sa0j9fh3dhs0i.xz

## Developer notes

If you need to modify this atlas or regenerate it these notes may be instructive.

script_2023_wagerlab_combined_atlas_prep.sh should be run before matlab or other bash scripts. This extracts cortical labels, among other things. You will need a copy of the Glasser cortical surface atlas. One place to find it is the Diedrichsen github repo.

The Diedrichsen atlas files can be obtained here: https://github.com/DiedrichsenLab/fs_LR_32

The glasser atlas can also be found in the hcp_utils package here as the "MMP" atlas:
    https://rmldj.github.io/hcp-utils/

Once the above is run, script_2023_wagerlab_combined_atlas_cifti.m can be executed in matlab. This will produce subcortical volumetric structures and two nifti atlases to accompany the subsequently generated CIFTI atlas. The parcel labels and order of the nifti atlases will match the cifti atlas. One nifti atlas will be 1mm, the other will be 2mm.

script_2013_wagerlab_combined_atlas_cifti.sh does three things,
1) produces the cifti dlabel file using Diedrichsen's version of the Glasser atlas. The subcortical structures in the CIFTI dlabel file should perfectly match the subcortical structures of the 2mm volumetric atlas.
2) It also outputs a config file needed for using the canlab2023 volumetric atlas in qsiprep. 
3) It creates a modified copie of the 2023 1mm volumetric atlas that lacks sform data. This is required for use of this atlas as a custom atlas in qsirecon. It is outputted as canlab_2023_MNI152NLin2009cAsym_1mm_lps.nii.gz

cifti visualization is best done using connectome workbench which you can get here:
    https://www.humanconnectome.org/software/get-connectome-workbench

Data files and surface models are stored separately for cifti files. There's no structure in my dlabel.nii atlas file for
instance. It can be visualized on any surface model you like. Sometimes it's helpful to visualize things on an inflated
surface to be able to get a better view into sulci. So in order to visualize any data in connectome workbench you also need
the surface mesh models to which your data maps. I've found the hcp_utils git repo helpful for this,
https://rmldj.github.io/hcp-utils/

The differing orientations of 2mm and 1mm files are a consequence of the fact that the 2mm grayordinate files were derived from HCP data, which is in LAS+ orientation, while the 1mm files were designed for use with qsiprep which requires LPS+ orietnation. Most software should handle this under the hood, but be aware all the same in case you get funny results when using these. If you want to query the orientation you can use fslhd. LPS+ orientation should result in left-right, anterior-posterior and inferior-superior designations for the x y and z axes, respectively.

##
Bogdan Petre
10/8/2023
