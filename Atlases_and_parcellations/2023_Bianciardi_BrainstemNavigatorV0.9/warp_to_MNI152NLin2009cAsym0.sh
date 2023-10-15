#!/bin/bash

# This script should be run inside the fmriprep singularity container by warp_to_MNI152NLin6Asym.sh

set -x

TEMPLATE_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/templates
XFM_DIR=$TEMPLATE_DIR/MNI152NLin6_to_MNI152NLin2009b

for INPUT in $(find BrainstemNavigator/0.9/2[ab]*_MNI/labels_probabilistic/ -type f -name "*nii.gz"); do
    INPUT=$(readlink -f $INPUT)
    INPUT_DIR=$(dirname $INPUT)
    OUTPUT=$(echo ${INPUT} | sed 's/.nii.gz/_MNI152NLin2009cAsym.nii.gz/')

    echo "Aligning $INPUT"

    antsApplyTransforms \
        -i $INPUT -r $TEMPLATE_DIR/MNI152NLin2009cAsym_T1_1mm.nii.gz \
        -o $OUTPUT -n LanczosWindowedSinc \
        -t $XFM_DIR/MNI_6thgen_2_MNI2009b.h5
done

