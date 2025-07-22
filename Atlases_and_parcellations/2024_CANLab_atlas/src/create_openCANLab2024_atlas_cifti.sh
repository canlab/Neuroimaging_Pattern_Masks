#!/bin/bash

# this should run after you've invoked this in matlab
# create_openCANLab2024_CIFTI_subctx('MNI152NLin6Asym','coarse',2,load_atlas('opencanlab2024_coarse_fsl6_2mm'))
# which will create a nifti file in this folder with all necessary subcortical volumes
#
# It assumes it's located in Atlases_and_parcellations/2024_CANLab_atlas/src/ and will not work if placed
# elsewhere. To move its location update the WD (working directory) variable so that it continues to point to
# the location of all of the dependencies used below (which are currently referenced w.r.t. the $WD variable

SCALE=coarse
SPACE=MNI152NLin6Asym
res=2

#WBCMD=/home/bogdan/Downloads/workbench/bin_rh_linux64/wb_command
WBCMD=wb_command

WD=$(cd $(dirname $(readlink -f $0)) && pwd)

$WBCMD -volume-label-import $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_vols.nii.gz \
    $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_vols.txt \
    $WD/subctx_atlas.label.nii

if [ $SPACE == "MNI152NLin6Asym" ]; then
    hcp_labels=$WD/hcp_cifti_subctx_labels.nii.gz
else
    hcp_labels=$WD/hcp_cifti_subctx_labels_${SPACE}.nii.gz
    if [ ! -e $hcp_labels ]; then
        echo "Could not find ${hcp_labels}"
        exit
    fi
fi
# Note to self: the label files below were copied from
# /dartfs-hpc/rc/home/m/f0042vm/software/diedrichsen_fs_LR_32/Glasser_2016.32k.L.label.gii
# the atlasROI files were obtained from
# https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases
$WBCMD -cifti-create-label $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    -volume $WD/subctx_atlas.label.nii ${hcp_labels} \
    -left-label $WD/../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.L.label.gii \
    -roi-left $WD/L.atlasroi.32k_fs_LR.shape.gii \
    -right-label $WD/../../2016_Glasser_Nature_HumanConnectomeParcellation/Glasser_2016.32k.R.label.gii \
    -roi-right $WD/R.atlasroi.32k_fs_LR.shape.gii

# Remove Glasser Hippocampus
$WBCMD -cifti-label-export-table $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii 1 \
    $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt

L_ind=$(cat -n $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | grep L_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | wc -l)

cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | head -n $[${L_ind}-1] > $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt
cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | tail -n	$[$len - $[$L_ind+1]] >> $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt

R_ind=$(cat -n $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt | grep R_H_ROI -B1 | tail -n 1 | awk '{print $1}')
len=$(cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt | wc -l)

cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt | head -n $[${R_ind}-1] > $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt
cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt | tail -n $[$len - $[$R_ind+1]] >> $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt

$WBCMD -cifti-label-import $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt \
    $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
    -discard-others -drop-unused-labels

# Remap indices to match volumetric atlas
if [ -e $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt ]; then
    src_ind=($(cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | awk '{if (NR%2==0) {print $1}}'))
    src_lbl=($(cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt | awk '{if (NR%2==1) {print $0}}'))

    rm -f $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt
    for i in $(seq 0 $[${#src_ind[@]}-1]); do
        this_lbl=$(echo ${src_lbl[$i]} | sed 's/\([LR]\)_\(.*\)_ROI/Ctx_\2_\1/' | sed 's/-/_/g')
        target_ind=$(cat $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt | grep \ $this_lbl\$ | awk '{print $1}')
        echo "${src_ind[$i]} ${target_ind}" >> $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt
    done

    currentDir=$(dirname $(readlink -f $0))
    $WBCMD -cifti-label-modify-keys $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
        $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt \
        $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii

    # replace glasser labels with canlab style labels
    $WBCMD -cifti-label-export-table $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii 1 \
        $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table.txt

    glasser_ind=($(cat $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt | awk '{print $1}'))
    glasser_lbl=($(cat $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_glasser_labels.txt | awk '{print $2}'))
    rm -f $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table_new.txt 
    for i in $(seq 0 $[${#glasser_lbl[@]}-1]); do
        canlab_lbl=$(cat $WD/openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt | \
                     grep "^${glasser_ind[i]} " | awk '{print $2}')
        cat $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table.txt | grep -A1 ${glasser_lbl[$i]}\$ | \
            sed "s/${glasser_lbl[$i]}/$canlab_lbl/" >> $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table_new.txt 
    done
    $WBCMD -cifti-label-import $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii \
        $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table_new.txt \
        $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii

else
    echo "Warning: Could not find key mapping file openCANLab2024_${SPACE}_${SCALE}_${res}mm_cifti_canlab_labels.txt. CIFTI keys won't match canlab atlas."
fi

# move file to main folder
mv -v $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii $WD/../openCANLab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii

# create fsaverage GIFTIs (for completion, despite being mostly redundant with glasser)
$WBCMD -cifti-separate $WD/../openCANLab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii COLUMN -label CORTEX_LEFT $WD/openCANLab2024_fs-LR-32k_hemi-L_${SCALE}.label.gii
$WBCMD -cifti-separate $WD/../openCANLab2024_${SPACE}_${SCALE}_${res}mm.dlabel.nii COLUMN -label CORTEX_RIGHT $WD/openCANLab2024_fs-LR-32k_hemi-R_${SCALE}.label.gii

# resample to surface using files derived from the fsaverage folder and a hcp 32k surface
# registered to the fsaverage surface (by default it's rotated for som ereason). The registered
# surface was obtained from 
# https://github.com/Washington-University/HCPpipelines/blob/master/global/templates/standard_mesh_atlases/fs_R/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii
# This can be improved by using adaptive barycentric interpolation but then you would need to
# generate area files for the HCP surface and that's more work than I want to do right now
# Seems the HCP people are in the same boat:
# https://www.mail-archive.com/hcp-users%40humanconnectome.org/msg02890.html
for H in L R; do
    $WBCMD -label-resample $WD/openCANLab2024_fs-LR-32k_hemi-${H}_${SCALE}.label.gii \
        $WD/S1200.${H}.sphere.32k_fs_LR.surf.gii \
        $WD/fs_${H}-to-fs_LR_fsaverage.${H}_LR.spherical_std.164k_fs_${H}.surf.gii \
        BARYCENTRIC \
        $WD/../openCANLab2024_fsaverage-164k_hemi-${H}_${SCALE}.label.gii
done

# garbage collection
rm -f $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_abr.txt \
    $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm.txt \
    $WD/opencanlab2024_${SPACE}_${SCALE}_${res}mm_atlasObj_to_cifti_mapping.txt \
    opencanlab2024_${SPACE}_${SCALE}_${res}mm_label_table.txt


