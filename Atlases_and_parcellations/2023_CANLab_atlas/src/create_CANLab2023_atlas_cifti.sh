#!/bin/bash

# this should run after you've evaluated this in matlab
# create_CANLab2023_CIFTI_subctx('MNI512NLin2009cAsym','fine',2,load_atlas('canlab2023_coarse_fsl6_2mm'))
# which will create a nifti file in this folder with all necessary subcortical volumes

SCALE=coarse
SPACE=MNI152NLin6Asym
res=2

wb_command -volume-label-import CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_vols.nii.gz \
    CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_vols.txt \
    subctx_atlas.label.nii

if [ $SPACE == "MNI152NLin6Asym" ]; then
    hcp_labels=hcp_cifti_subctx_labels.nii.gz
else
    hcp_labels=hcp_cifti_subctx_labels_${SPACE}.nii.gz
    if [ ! -e $hcp_labels ]; then
        echo "Could not find ${hcp_labels}"
        exit
    fi
fi
# Note to self: the label files below were copied from
# /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii
wb_command -cifti-create-label canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii \
    -volume subctx_atlas.label.nii ${hcp_labels} \
    -left-label ../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.L.label.gii \
    -right-label ../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.R.label.gii

# Remove Glasser Hippocampus
wb_command -cifti-label-export-table canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii 1 \
    canlab2023_${SCALE}_${SPACE}_${res}mm.txt

L_ind=$(cat -n canlab2023_${SCALE}_${SPACE}_${res}mm.txt | grep L_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat canlab2023_${SCALE}_${SPACE}_${res}mm.txt | wc -l)

cat canlab2023_${SCALE}_${SPACE}_${res}mm.txt | head -n $[${L_ind}-1] > canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt
cat canlab2023_${SCALE}_${SPACE}_${res}mm.txt | tail -n	$[$len - $[$L_ind+1]] >> canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt

R_ind=$(cat -n canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt | grep R_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt | wc -l)

cat canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt | head -n $[${R_ind}-1] > canlab2023_${SCALE}_${SPACE}_${res}mm.txt
cat canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt | tail -n $[$len - $[$R_ind+1]] >> canlab2023_${SCALE}_${SPACE}_${res}mm.txt

wb_command -cifti-label-import canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii \
    canlab2023_${SCALE}_${SPACE}_${res}mm.txt \
    canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii \
    -discard-others -drop-unused-labels

# Remap indices to match volumetric atlas
if [ -e CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_glasser_labels.txt ]; then
    src_ind=($(cat canlab2023_${SCALE}_${SPACE}_${res}mm.txt | awk '{if (NR%2==0) {print $1}}'))
    src_lbl=($(cat canlab2023_${SCALE}_${SPACE}_${res}mm.txt | awk '{if (NR%2==1) {print $0}}'))

    rm -f canlab2023_${SCALE}_${SPACE}_${res}mm_atlasObj_to_cifti_mapping.txt
    for i in $(seq 0 $[${#src_ind[@]}-1]); do
        this_lbl=${src_lbl[$i]}
        target_ind=$(cat CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_index_labels.txt | grep $this_lbl\$ | awk '{print $1}')
        echo "${src_ind[$i]} ${target_ind}" >> canlab2023_${SCALE}_${SPACE}_${res}mm_atlasObj_to_cifti_mapping.txt
    done

    currentDir=$(dirname $(readlink -f $0))
    wb_command -cifti-label-modify-keys canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii \
        canlab2023_${SCALE}_${SPACE}_${res}mm_atlasObj_to_cifti_mapping.txt \
        canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii

    # replace glasser labels with canlab style labels
    wb_command -cifti-label-export-table canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii 1 \
        canlab2023_${SCALE}_${SPACE}_${res}mm_label_table.txt

    glasser_ind=($(cat CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_glasser_labels.txt | awk '{print $1}'))
    glasser_lbl=($(cat CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_glasser_labels.txt | awk '{print $2}'))
    rm -f canlab2023_${SCALE}_${SPACE}_${res}mm_label_table_new.txt 
    for i in $(seq 0 $[${#glasser_lbl[@]}-1]); do
        canlab_lbl=$(cat CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_canlab_labels.txt | \
                     grep "^${glasser_ind[i]} " | awk '{print $2}')
        cat canlab2023_${SCALE}_${SPACE}_${res}mm_label_table.txt | grep -A1 ${glasser_lbl[$i]}\$ | \
            sed "s/${glasser_lbl[$i]}/$canlab_lbl/" >> canlab2023_${SCALE}_${SPACE}_${res}mm_label_table_new.txt 
    done
    wb_command -cifti-label-import canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii \
        canlab2023_${SCALE}_${SPACE}_${res}mm_label_table_new.txt \
        canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii

else
    echo "Warning: Could not find key mapping file CANLab2023_${SCALE}_${SPACE}_${res}mm_cifti_index_labels.txt. CIFTI keys won't match canlab atlas."
fi

# copy file to main folder
currentDir=$(dirname $(readlink -f $0))
mv -v canlab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii $currentDir/../CANLab2023_${SCALE}_${SPACE}_${res}mm.dlabel.nii

# garbage collection
rm -f canlab2023_${SCALE}_${SPACE}_${res}mm_abr.txt \
    canlab2023_${SCALE}_${SPACE}_${res}mm.txt \
    canlab2023_${SCALE}_${SPACE}_${res}mm_atlasObj_to_cifti_mapping.txt
#    canlab2023_${SCALE}_${SPACE}_${res}mm_label_table.txt

