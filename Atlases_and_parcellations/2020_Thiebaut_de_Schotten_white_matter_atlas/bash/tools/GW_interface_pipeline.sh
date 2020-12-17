#! /bin/bash

/home/tolhs/fslcompiled/fsl/bin/flirt \
  -in /data/BCBToolKit/Tools/extraFiles/Priors/brainWithSkullTemplate.nii.gz \
  -ref /home/tolhs/fslcompiled/fsl/data/standard/MNI152 \
  -out /data/Chris/ants_NKI_priors_mni \
  -omat /data/Chris/ants_NKI_priors_mni.mat \
  -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 \
  -dof 12  -interp trilinear

for f in /data/BCBToolKit/Tools/extraFiles/Priors/*nii*;
do
  /home/tolhs/fslcompiled/fsl/bin/flirt -in $f \
  -ref /home/tolhs/fslcompiled/fsl/data/standard/MNI152lin_T1_1mm.nii.gz \
  -out /data/Chris/ants_NKI_priors_mni/$(basename $f) -applyxfm \
  -init /data/Chris/ants_NKI_priors_mni.mat -interp trilinear;
done;

fslmaths priors2.nii.gz -thr 0.5 gray_05
fslmaths priors3.nii.gz -thr 0.5 white_05

fslmaths gray_05 -dilM dil_gray_05
fslmaths white_05 -dilM dil_white_05
fslmaths dil_white_05.nii.gz -mas dil_gray_05.nii.gz interface

# THR 0.3
fslmaths priors2.nii.gz -thr 0.3 gray_03
fslmaths priors3.nii.gz -thr 0.3 white_03

fslmaths gray_03 -dilM dil_gray_03
fslmaths white_03 -dilM dil_white_03
fslmaths dil_white_03.nii.gz -mas dil_gray_03.nii.gz interface_03
