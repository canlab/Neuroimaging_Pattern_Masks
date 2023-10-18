#!/bin/bash

# fmriprep was used to align MNI152NLin6Asym 1mm T1 images (non-skull stripped) to the MNI152NLin2009cAsym:res-01 space,
# which is 1mm space in templateFlow that fmriprep uses to pull its templates. This produced ants format *.h5 files which
# are composite affine/warp files.
# these files were created by converting ants *.h5 files into fsl formated files using code like this:

# this is my fmriprep output. Copies of the necessary files are also in the sister directory of this one, "ants"
srcRoot=/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas/preproc/fmriprep/sub-MNI152NLin6Asym/anat

# source data we passed for alignment to fmriprep
# this file can be obtained here
# https://templateflow.s3.amazonaws.com/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-01_T1w.nii.gz
REF_FSL=~/software/canlab/Neuroimaging_Pattern_Masks/templates/MNI152NLin6Asym_T1_1mm.nii.gz

# this file can be obtained here
# https://templateflow.s3.amazonaws.com/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz
REF_FMRIPREP=~/software/canlab/Neuroimaging_Pattern_Masks/templates/MNI152NLin2009cAsym_T1_1mm.nii.gz

# these operations break down the h5 files into constituent affine matrices and warp field files
CompositeTransformUtil --disassemble $srcRoot/sub-MNI152NLin6Asym_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 \
    fmriprep_to_fsl

mv [0-9][0-9]*fmriprep_to_fsl* ../ants/

# Both Affine matrix files and warp fields need to be modified to work with FSL. We can do this with
# two different tools, c3d_affine_tool obtained here,
# http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D
# and connectome workbench, available here:
# https://www.humanconnectome.org/software/get-connectome-workbench
 
c3d=~/software/c3d-1.1.0-Linux-x86_64/bin/c3d_affine_tool

$c3d -src $REF_FMRIPREP -ref $REF_FSL -itk ../ants/01_fmriprep_to_fsl_AffineTransform.mat -ras2fsl \
    -o ../fsl/01_fmriprep_to_fsl_AffineTransform_mod.mat
wb_command -convert-warpfield -from-itk ../ants/00_fmriprep_to_fsl_DisplacementFieldTransform.nii.gz \
    -to-fnirt ../fsl/00_fmriprep_to_fsl_DisplacementFieldTransform_mod.nii.gz $REF_FMRIPREP
