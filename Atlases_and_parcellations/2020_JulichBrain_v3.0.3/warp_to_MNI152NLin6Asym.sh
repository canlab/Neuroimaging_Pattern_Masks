#!/bin/bash

#SBATCH --job-name antsApplyTransforms
#SBATCH --time 7-00:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --hint=nomultithread
#SBATCH --output warp_logs1/warp_%a.out
#SBATCH --error warp_logs1/warp_%a.err
#SBATCH --account dbic

# this script maps from MNI152NLin2009cAsym space to MNI152NLin6Asym space

THIS_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/

singularity exec --bind $THIS_DIR --cleanenv \
    /dartfs-hpc/rc/lab/C/CANlab/modules/fmriprep-20.2.3-LTS.sif \
    $THIS_DIR/Atlases_and_parcellations/2020_JulichBrain_v3.0.3/warp_to_MNI152NLin6Asym0.sh

