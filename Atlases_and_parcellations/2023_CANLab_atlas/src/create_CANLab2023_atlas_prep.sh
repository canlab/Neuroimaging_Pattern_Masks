#!/bin/bash

# this script exports labels from the glasser atlas in the same order they're in in the atlas. We will use this to
# confirm the nifti atlas parcels are ordered correctly. This should be run before
# create_CANLab2023_atlas_unrestricted.m or create_CANLab2023_atlas_cifti.sh scripts

wb_command -label-export-table ../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.L.label.gii \
    lctx_labels.txt

wb_command -label-export-table ../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.R.label.gii \
    rctx_labels.txt
