#! /bin/bash
[ $# -lt 3 ] && { echo 'Usage :
  $1 = input folder (With the / after the name)
  $2 = seed mask name
  $3 = target mask name
  $4 = result folder
  $5 = other mask name (Will mask both seed and target)
  ...'; exit 1; }

###############################################################################
# Here we create a seed entirely contained in the target                      #
###############################################################################

mkdir -p $4

seed=${1}$2
target=${1}$3
if [[ $5 != "" ]];
then
  fslmaths $seed -mas ${1}$5 $4/masked_$(basename $seed)
  seed=$4/masked_$(basename $seed)
  fslmaths $target -mas ${1}$5 $4/target_$(basename $target)
  target=$4/target_$(basename $target)
fi;

fslmaths $seed -mas $target $4/seed_$(basename $2)
