#!/bin/bash

# This script should be run inside the fmriprep singularity container by warp_to_MNI152NLin6Asym.sh

TEMPLATE_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/templates
XFM_DIR=$TEMPLATE_DIR/transforms/ants

for INPUT in $(find probabilistic_maps_pmaps_157areas/ -type f -name "*nii.gz"); do
    INPUT=$(readlink -f $INPUT)
    INPUT_DIR=$(dirname $INPUT)
    OUTPUT=$(echo $INPUT | sed 's/nlin2ICBM152asym2009c/nlin2ICBM152asym6/')

    echo "Aligning $INPUT"

    antsApplyTransforms \
        -i $INPUT -r $TEMPLATE_DIR/MNI152NLin6Asym_T1_1mm.nii.gz \
        -o $OUTPUT -n LanczosWindowedSinc \
        -t $XFM_DIR/MNI152NLin2009cAsym_to_MNI152NLin6Asym.h5
done

