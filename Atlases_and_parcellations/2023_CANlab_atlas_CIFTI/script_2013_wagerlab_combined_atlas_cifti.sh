#!/bin/bash

wb_command -volume-label-import subctx_atlas.nii subctx_atlas.txt subctx_atlas.label.nii
wb_command -cifti-create-label canlab_2023_2mm.dlabel.nii \
    -volume subctx_atlas.label.nii hcp_cifti_subctx_labels.nii \
    -left-label /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii \
    -right-label /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.R.label.gii
