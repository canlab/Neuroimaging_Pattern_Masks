#!/bin/bash

# this is an example script of what was done to assemble the registration fusion results
# into a single file. The resultant file is what's shared on figshare

set -x
OUTDIR=/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas/results/glasser/subject_parcellations/

for SPACE in MNI152NLin2009cAsym MNI152NLin6Asym; do

    lh_list_paingen=($(for sid in $(cat paingen_iid_mrns.csv); do
        if [ -e $OUTDIR/lh_${SPACE}_glasser_sub-${sid}.nii.gz ] && [ -e $OUTDIR/rh_${SPACE}_glasser_sub-${sid}.nii.gz ]; then
            echo -n $OUTDIR/lh_${SPACE}_glasser_sub-${sid}.nii.gz" "
        fi
    done))

    rh_list_paingen=($(for sid in $(cat paingen_iid_mrns.csv); do
        if [ -e $OUTDIR/lh_${SPACE}_glasser_sub-${sid}.nii.gz ] && [ -e $OUTDIR/rh_${SPACE}_glasser_sub-${sid}.nii.gz ]; then
            echo -n $OUTDIR/rh_${SPACE}_glasser_sub-${sid}.nii.gz" "
        fi
    done))

    fslmerge -t $OUTDIR/lh_paingen_${SPACE}.nii.gz ${lh_list_paingen[@]} &
    fslmerge -t $OUTDIR/rh_paingen_${SPACE}.nii.gz ${rh_list_paingen[@]} &
done

wait

for SPACE in MNI152NLin2009cAsym MNI152NLin6Asym; do
    fslmaths $OUTDIR/rh_paingen_${SPACE}.nii.gz \
	-add 180 \
	-thr 181 \
	-add $OUTDIR/lh_paingen_${SPACE}.nii.gz \
	-uthr 360 \
	$OUTIDR/paingen_${SPACE}.nii.gz &
done

wait
