#!/bin/bash

# This is an example of our "registration fusion" implementation, modeled after
# Wu, Ngo, Greve, et al. (2018) Human Brain Mapping.
#
# this script is run within an fmriprep 20.2.3 (LTS) singularity container.
# it produces left and right hemispheric projections into volumetric space
# for a specified subject that's already had fmriprep and freesurfer recon-all
# run on its T1w volume. It is only an example script which was specifically
# designed for use with the PainGen study in the Wager lab. It may need
# slight modifications to some paths to run with other studies, but they
# would be minor.
#
# Wrap this in a loop or invoke it in job array batch script on an HPC system
# to run.
#
# The result will be two nifti files with a subjects left hemisphere glasser
# parcels in volumetric space and the right hemisphere glasser parcels in a second
# volumetric file. The alignment will be subject specific and conform to the
# particular subjetcs cortical folding patterns.

set -x

SID=$1          # e.g. sub-M80318096
SCRATCH_DIR=$2  # everything here will be erased after the run, should be unique for each invocation
OUTDIR=$3       # location to write [rl]h_MNI152NLin2009cAsym_${SID}.nii.gz outputs


# this is an fmriprep output directory which should contain transformations
# from T1w space to your target space, in this case MNI152NLin2009cAsym, and
# freesurfer outputs in preproc/freesufer
SRCROOT=/dartfs-hpc/rc/lab/C/CANlab/labdata/data/Pain_Gen/Imaging/preproc

# this should contain a reference volume, e.g. MNI152NLin2009cAsym 1mm brain
# from templateFlow (which is where fmriprep gets all its templates)
TEMPLATE_DIR=/dartfs-hpc/rc/lab/C/CANlab/modules/Neuroimaging_Pattern_Masks/templates

# RES_DIR should contain lctx_labels.txt and rctx_labels.txt, ascii files that
# contain glasser label names in the first column and indices in the second column
# (1-180, not 1-360, although it may not matter). It should also contain the
# HCP annotation files, obtained here:
# https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446
RES_DIR=/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_atlas/_resources/

export FS_LICENSE=$RES_DIR/freesurfer_license.txt

transform=$(readlink -f $SRCROOT/fmriprep/$SID/ses-*/anat/${SID}_ses-*_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 | head -n 1)

# first thing's first, we have annotation files from here:
# https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446
# but we need surfaces, so we invert the steps they list to move back on stage
# We use the same surfaces they use, which I got here:
# https://github.com/Washington-University/HCPpipelines/blob/master/global/templates/standard_mesh_atlases/fs_R/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii
# I don't know why they don't use the alternative surface in the fs_[RL] folders, 
# I'm just blindly reversing what they did trusting that it was correct

# convert annot to fsaverage6 spherical surface
mris_convert --annot $RES_DIR/lh.HCP-MMP1.annot \
    --parcstats $RES_DIR/lctx_labels.txt \
    $RES_DIR/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii \
    $SCRATCH_DIR/lh.parcstats

# the above is now a surface with vertices in the same locations as the 164k fsaverage sphere, except
# that's not where we need the vertices to project this onto a subject specific sphere. The subject
# specific sphere has different vertices (in the same coordinate system), so let's resample

# resample fsaverage6 sphere onto aligned subject sphere's vertices
mri_surf2surf --srcsubject ico \
    --sval $SCRATCH_DIR/lh.parcstats \
    --trgsubject ${SID} \
    --tval $SCRATCH_DIR/lh.gii \
    --hemi lh \
    --mapmethod nnfr \
    --sd $SRCROOT/freesurfer/

# now that we have the labels in the subject sphere's space, let's go from that to the subject's native space volume
mri_surf2vol --surfval $SCRATCH_DIR/lh.gii \
    --hemi lh \
    --template $SRCROOT/freesurfer/$SID/mri/norm.mgz \
    --identity $SID \
    --sd $SRCROOT/freesurfer \
    --fillribbon \
    --o $SCRATCH_DIR/lh.${SID}.nii.gz

# now let's project to MNI152NLin2009cAsym space from T1w space
antsApplyTransforms -i $SCRATCH_DIR/lh.${SID}.nii.gz \
    -o $OUTDIR/lh_MNI152NLin2009cAsym_glasser_${SID}.nii.gz \
    -r $TEMPLATE_DIR/MNI152NLin2009cAsym_T1_1mm.nii.gz \
    -t $transform \
    -n NearestNeighbor -v


# now repeat the above for the right hemisphere

mris_convert --annot $RES_DIR/rh.HCP-MMP1.annot \
    --parcstats $RES_DIR/rctx_labels.txt \
    $RES_DIR/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii \
    $SCRATCH_DIR/rh.parcstats

mri_surf2surf --srcsubject ico \
    --sval $SCRATCH_DIR/rh.parcstats \
    --trgsubject ${SID} \
    --tval $SCRATCH_DIR/rh.gii \
    --hemi rh \
    --mapmethod nnfr \
    --sd $SRCROOT/freesurfer/

mri_surf2vol --surfval $SCRATCH_DIR/rh.gii \
    --hemi rh \
    --template $SRCROOT/freesurfer/$SID/mri/norm.mgz \
    --identity $SID \
    --sd $SRCROOT/freesurfer \
    --fillribbon \
    --o $SCRATCH_DIR/rh.${SID}.nii.gz

antsApplyTransforms -i $SCRATCH_DIR/rh.${SID}.nii.gz \
    -o $OUTDIR/rh_MNI152NLin2009cAsym_glasser_${SID}.nii.gz \
    -r $TEMPLATE_DIR/MNI152NLin2009cAsym_T1_1mm.nii.gz \
    -t $transform -n NearestNeighbor -v
