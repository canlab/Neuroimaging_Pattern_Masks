#!/bin/bash

# This script should be run inside the fmriprep singularity container by warp_to_MNI152NLin6Asym.sh

set -x

TEMPLATE_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/templates
XFM_DIR=$TEMPLATE_DIR/transforms/ants

for INPUT in hcp_cifti_subctx_labels.nii.gz; do
    INPUT=$(readlink -f $INPUT)
    OUTPUT=${INPUT%%.*}_MNI152NLin2009cAsym.nii

    echo "Aligning $INPUT"

    antsApplyTransforms \
        -i $INPUT -r $TEMPLATE_DIR/MNI152NLin2009cAsym_T1_1mm.nii.gz \
        -o $OUTPUT -n NearestNeighbor \
        -t $XFM_DIR/MNI152NLin6Asym_to_MNI152NLin2009cAsym.h5
done

