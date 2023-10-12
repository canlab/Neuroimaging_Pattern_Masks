#!/bin/bash

# this script exports labels from the glasser atlas in the same order they're in in the atlas. We will use this to
# confirm the nifti atlas parcels are ordered correctly. This should be run before
# script_2023_wagerlab_combined_atlas_prep.sh

wb_command -label-export-table /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii \
    lctx_labels.txt

wb_command -label-export-table /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.R.label.gii \
    rctx_labels.txt
