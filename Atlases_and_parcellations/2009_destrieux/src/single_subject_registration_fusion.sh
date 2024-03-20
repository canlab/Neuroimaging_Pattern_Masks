#!/bin/bash

#SBATCH --job-name spacetop_dk_proj
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --hint=nomultithread
#SBATCH --output dk_spacetop_logs/atlas_%a.out
#SBATCH --error dk_spacetop_logs/atlas_%a.err
#SBATCH --account dbic
#SBATCH --array=1-113
#SBATCH --dependency=4798298
#SBATCH --exclude=s11

module load freesurfer/7.4.1

ROOT=/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas/
SRCROOT=/dartfs-hpc/rc/lab/C/CANlab/labdata/data/spacetop_data/derivatives/fmriprep/
TEMPLATE_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/templates

RES_DIR=$ROOT/_resources/
SCRATCH_DIR=/dartfs-hpc/scratch/$user/$(uuidgen)
mkdir -p $SCRATCH_DIR

OUTDIR=$ROOT/results/desikan_killiany/subject_parcellations_construction

sid_list=($(for file in $(find $SRCROOT/results/fmriprep/sub-*/*/anat/ -name "*to-MNI152NLin2009c*h5"); do basename $(dirname $(dirname $(dirname $file))); done))

iter=$[${SLURM_ARRAY_TASK_ID} - 1]

set -x

export FS_LICENSE=$RES_DIR/freesurfer_license.txt

SID=${sid_list[$iter]}

for atlas in dk destrieux; do
    if [ $atlas == "dk" ]; then
        src=""
    elif [ $atlas == "destrieux" ]; then
        src=.a2009s
    else
        >&2 echo "ERROR: did not understand atlas=$atlas"
        exit
    fi

    for TEMPLATE in MNI152NLin2009cAsym MNI152NLin6Asym; do
        for h in lh rh; do
            transform=$(readlink -f $SRCROOT/results/fmriprep/$SID/ses-01/anat/${SID}*_from-T1w_to-${TEMPLATE}_mode-image_xfm.h5 | head -n 1)

            # proj annot to vol
            mri_surf2vol --so $SRCROOT/results/fmriprep/sourcedata/freesurfer/$SID/surf/${h}.white \
                    $SRCROOT/results/fmriprep/sourcedata/freesurfer/$SID/label/${h}.aparc${src}.annot \
                --subject $SID \
                --sd $SRCROOT/results/fmriprep/sourcedata/freesurfer/ \
                --o $SCRATCH_DIR/${h}.${atlas}.${SID}.nii.gz

            # now let's project to MNI152NLin2009cAsym space from T1w space
            singularity exec --bind $ROOT,$SRCROOT,$TEMPLATE_DIR,$RES_DIR,$SCRATCH_DIR,$OUTDIR --cleanenv \
                /dartfs-hpc/rc/lab/C/CANlab/modules/fmriprep-20.2.3-LTS.sif \
                antsApplyTransforms -i $SCRATCH_DIR/${h}.${atlas}.${SID}.nii.gz \
                    -o $OUTDIR/${h}_${TEMPLATE}_${atlas}_${SID}.nii.gz \
                    -r $TEMPLATE_DIR/${TEMPLATE}_T1_1mm.nii.gz \
                    -t $transform \
                    -n NearestNeighbor -v
        done
    rm -rf $SCRATCH_DIR/*
    done
done
    
rm -rf $SCRATCH_DIR
