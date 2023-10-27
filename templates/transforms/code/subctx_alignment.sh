#!/bin/bash
# this version was based on fmriprep output 

ROOT=/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas
RES_DIR=$(readlink -f $ROOT/_resources)

NPM_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks

#make subcortical mask
# we doubly dilate it so that it encompases the subcortex in the refernece space in spite of whatever
# minor misalignment we may have between the two spaces
RES_DIR=$(readlink -f ../../_resources);

module load fsl 
fslmaths $RES_DIR/hcp_cifti_subctx_labels.nii.gz -bin -dilD -dilD cifti_mask.nii.gz


MIMG=$NPM_DIR/templates/MNI152NLin6Asym_T1_1mm.nii.gz
FIMG=$NPM_DIR/templates/MNI152NLin2009cAsym_T1_1mm.nii.gz

singularity exec --bind $ROOT,$RES_DIR,$NPM_DIR --cleanenv \
    /dartfs-hpc/rc/lab/C/CANlab/modules/fmriprep-20.2.3-LTS.sif \
    antsRegistration -o [$ROOT/sandbox/subctx_alignment/MNI152NLin6Asym_to_MNI152NLin2009cAsym_subctx_,1,1] \
        --dimensionality 3 \
        --collapse-output-transforms \
        --write-composite-transform \
        --interpolation LanczosWindowedSinc \
        --transform Rigid[0.05] \
        --metric Mattes[$FIMG,$MIMG,1,32,Regular,0.25] \
        --convergence [1000x500x250x100,0.00000001,10] \
        --shrink-factors 8x4x2x1 \
        --smoothing-sigmas 4x2x1x0vox \
        --use-histogram-matching 1 \
        --masks [cifti_mask.nii.gz,NULL] \
        --transform Affine[0.1] \
        --metric Mattes[$FIMG,$MIMG,1,32,Regular,0.25] \
        --convergence [1000x500x250x100,1e-8,10] \
        --shrink-factors 8x4x2x1 \
        --smoothing-sigmas 4x2x1x0vox \
        --use-histogram-matching 1 \
        --masks [cifti_mask.nii.gz,NULL] \
        --transform SyN[0.05,3,0] \
        --metric CC[$FIMG,$MIMG,1,4,None,1] \
        --convergence [100x100x50,1e-9,15] \
        --shrink-factors 4x2x1 \
        --smoothing-sigmas 2x1x0vox \
        --use-histogram-matching 1 \
        --masks [cifti_mask.nii.gz,NULL] \
        --winsorize-image-intensities [0.025,0.975] \
        --verbose 1
