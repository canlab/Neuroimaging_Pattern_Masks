#! /bin/bash
[ $# -lt 2 ] && { echo 'Usage :
  $1 = fmri folder
  $2 = result folder'; exit 1; }

mkdir -p $2

for d in $1/*/;
do
  name=$(basename $d)
  imcp $d/RS_clean_plusmean_bptf_MNI_2mm.nii.gz $2/$name
done;
