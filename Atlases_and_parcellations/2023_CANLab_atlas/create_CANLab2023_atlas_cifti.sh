#!/bin/bash

# this should run after you've evaluated this in matlab
# create_CANLab2023_CIFTI_subctx('MNI152NLin6Asym','coarse',2,load_atlas('canlab2023_coarse_fsl6_2mm'))
# which will create a nifti file in this folder with all necessary subcortical volumes
#
# It assumes it's located in Atlases_and_parcellations/2023_CANLab_atlas/src/ and will not work if placed
# elsewhere. To move its location update the WD (working directory) variable so that it continues to point to
# the location of all of the dependencies used below (which are currently referenced w.r.t. the $WD variable

SCALE=coarse
SPACE=MNI152NLin6Asym
res=2

WBCMD=/home/bogdan/Downloads/workbench/bin_rh_linux64/wb_command

WD=$(cd $(dirname $(readlink -f $0)) && pwd)

$WBCMD -volume-label-import $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_vols.nii.gz \
    $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_vols.txt \
    $WD/src/subctx_atlas.label.nii

if [ $SPACE == "MNI152NLin6Asym" ]; then
    hcp_labels=$WD/src/hcp_cifti_subctx_labels.nii.gz
else
    hcp_labels=$WD/src/hcp_cifti_subctx_labels_${SPACE}.nii.gz
    if [ ! -e $hcp_labels ]; then
        echo "Could not find ${hcp_labels}"
        exit
    fi
fi
# Note to self: the label files below were copied from
# /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii
$WBCMD -cifti-create-label $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    -volume $WD/src/subctx_atlas.label.nii ${hcp_labels} \
    -left-label $WD/../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.L.label.gii \
    -right-label $WD/../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.R.label.gii

# Remove Glasser Hippocampus
$WBCMD -cifti-label-export-table $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii 1 \
    $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt

L_ind=$(cat -n $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | grep L_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | wc -l)

cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | head -n $[${L_ind}-1] > $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt
cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | tail -n	$[$len - $[$L_ind+1]] >> $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt

R_ind=$(cat -n $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt | grep R_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt | wc -l)

cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt | head -n $[${R_ind}-1] > $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt
cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt | tail -n $[$len - $[$R_ind+1]] >> $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt

$WBCMD -cifti-label-import $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt \
    $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    -discard-others -drop-unused-labels

# Remap indices to match volumetric atlas
if [ -e $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt ]; then
    src_ind=($(cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | awk '{if (NR%2==0) {print $1}}'))
    src_lbl=($(cat $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt | awk '{if (NR%2==1) {print $0}}'))

    rm -f $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt
    for i in $(seq 0 $[${#src_ind[@]}-1]); do
        this_lbl=$(echo ${src_lbl[$i]} | sed 's/\([LR]\)_\(.*\)_ROI/Ctx_\2_\1/' | sed 's/-/_/g')
        target_ind=$(cat $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt | grep \ $this_lbl\$ | awk '{print $1}')
        echo "${src_ind[$i]} ${target_ind}" >> $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt
    done

    currentDir=$(dirname $(readlink -f $0))
    $WBCMD -cifti-label-modify-keys $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
        $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt \
        $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii

    # replace glasser labels with canlab style labels
    $WBCMD -cifti-label-export-table $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii 1 \
        $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_label_table.txt

    glasser_ind=($(cat $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt | awk '{print $1}'))
    glasser_lbl=($(cat $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt | awk '{print $2}'))
    rm -f $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_label_table_new.txt 
    for i in $(seq 0 $[${#glasser_lbl[@]}-1]); do
        canlab_lbl=$(cat $WD/src/CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt | \
                     grep "^${glasser_ind[i]} " | awk '{print $2}')
        cat $WD/src/canlab2023_${SPACE}_${SCALE}_${res}mm_label_table.txt | grep -A1 ${glasser_lbl[$i]}\$ | \
            sed "s/${glasser_lbl[$i]}/$canlab_lbl/" >> $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_label_table_new.txt 
    done
    $WBCMD -cifti-label-import $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
        $WD/src/canlab2023_${SPACE}_${SCALE}_${res}mm_label_table_new.txt \
        $WD/src/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii

else
    echo "Warning: Could not find key mapping file CANLab2023_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt. CIFTI keys won't match canlab atlas."
fi

# copy file to main folder
mv -v $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii $WD/../CANLab2023_${SPACE}_${SCALE}_${res}mm.dlabel.nii

# garbage collection
rm -f $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_abr.txt \
    $WD/canlab2023_${SPACE}_${SCALE}_${res}mm.txt \
    $WD/canlab2023_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt \
    canlab2023_${SPACE}_${SCALE}_${res}mm_label_table.txt

