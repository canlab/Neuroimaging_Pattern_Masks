#!/bin/bash

wb_command -volume-label-import subctx_atlas.nii subctx_atlas.txt subctx_atlas.label.nii
wb_command -cifti-create-label canlab_2023_2mm.dlabel.nii \
    -volume subctx_atlas.label.nii hcp_cifti_subctx_labels.nii \
    -left-label /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii \
    -right-label /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.R.label.gii


# create qsiprep atlas config file
cat <<END > atlas_config.json
{
  "canlab2023": {
    "file": "canlab_2023_2mm.nii.gz",
    "node_names": [
END

wb_command -cifti-label-export-table canlab_2023_2mm.dlabel.nii 1 canlab_2023_2mm_labels.txt

# the way we handle "," and newlines below is convoluted. It's intention is to
# have a "," after every newline except for the last in a list
n_ROI=$(echo $(cat canlab_2023_2mm_labels.txt | wc -l)/2 | bc)
for ROI in $(cat canlab_2023_2mm_labels.txt | grep -v ^[0-9]); do 
    echo ","
    echo -n "      \"$ROI\""
done | tail -n $n_ROI >> atlas_config.json

cat <<END >> atlas_config.json

    ],
    "node_ids": [
END

for i in $(seq 1 $n_ROI); do 
    echo ","
    echo -n "      $i"
done | tail -n $n_ROI >> atlas_config.json
    
cat <<END >> atlas_config.json

    ]
  }
}
END

