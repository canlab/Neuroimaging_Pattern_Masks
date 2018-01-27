#!/bin/bash
# ANTs bimodal SyN warp CIT168 T1w and T2w head templates to MNI152 2009c nonlin asym T1w and T2w templates
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2017-10-31 JMT Adapt from mirror_warp_calc.sh

# Key directories
cit_dir=./CIT168
mni_dir=./mni_icbm152_nlin_asym_09c

# Key files
cit_t1=${cit_dir}/CIT168_T1w_head_700um.nii.gz
cit_t2=${cit_dir}/CIT168_T2w_head_700um.nii.gz
mni_t1=${mni_dir}/mni_icbm152_t1_tal_nlin_asym_09c.nii
mni_t2=${mni_dir}/mni_icbm152_t2_tal_nlin_asym_09c.nii
cit2mni_t1=CIT168toMNI152_T1w_head_700um.nii.gz
cit2mni_t2=CIT168toMNI152_T2w_head_700um.nii.gz

# ANTs output prefix
ants_prefix=cit2mni_
cit2mni_affine=${ants_prefix}Affine.txt
cit2mni_warp=${ants_prefix}Warp.nii.gz

# Bivariate ANTs SyN registration options
ants_opts="-i 30x90x20 -t SyN[0.25] -r Gauss[3,0] --use-Histogram-Matching --number-of-affine-iterations 10000x10000x1000 --MI-option 32x16000"

# Calculate CIT168 to MNI152 SyN warp
if [ ! -s ${cit2mni_warp} ]; then
    echo "Calculating CIT168 to MNI152 SyN warp"
    ANTS 3 -m CC[${mni_t1},${cit_t1},1,5] -m CC[${mni_t2},${cit_t2},1,5] ${ants_opts} -o ${ants_prefix}
fi

# Apply warp to T1w and T2w templates
if [ ! -s ${cit2mni_t1} ]; then
    echo "Warping CIT168 T1 to MNI152 space"
    WarpImageMultiTransform 3 ${cit_t1} ${cit2mni_t1} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline
fi
if [ ! -s ${cit2mni_t2} ]; then
    echo "Warping CIT168 T2 to MNI152 space"
    WarpImageMultiTransform 3 ${cit_t2} ${cit2mni_t2} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline
fi
