#!/usr/bin/env python3
# -*-coding:Utf-8 -*

import sys
import nibabel as nib
import numpy as np
import math
from scipy.special import xlogy  # xlogy returns 0 if x,
                                 # which here is also y, is=0


if len(sys.argv) < 4:
    # fmri_4D.nii.gz
    print("1 : fMRI data\n")
    # MNI152_T1_4mm_brain_mask.nii.gz
    print("2 : the brainmask, to know on which voxels to calculate H\n")
    # The result file (the parent folders should exist)
    print("3 : result file path\n")
    sys.exit()

# The bash $1 is sys.argv[1]
# We load the image with nibabel

img4D = nib.load(sys.argv[1])
# img4D = nib.load("fmri_4D.nii.gz")
# We load the data array
data4D = img4D.get_data()

# Load the mask
brainmask = nib.load(sys.argv[2])
# brainmask = nib.load("MNI152_T1_4mm_brain_mask.nii.gz")
# Load its data array
data_brainmask = brainmask.get_data()

# Check if the shape of the mask is the same as the first 3D of the fMRI
if data4D.shape[0:3] != data_brainmask.shape:
    print("ERROR : the shape of the mask should be the same as the 4D image")
    sys.exit()

data_H = np.zeros(brainmask.shape)
# where there is other value than 0 (np.where returns a tuple)
ind_mask = np.array(np.where(data_brainmask))
#  Number of voxels
nvox = ind_mask.shape[1]
# To study !!!
nb_bins=round(math.sqrt(data4D.shape[3]))
# We will calculate the entropy for EACH voxel
for i in range(0, nvox):
    # coordinates of the voxel in the mask
    vox = ind_mask[:,i]
    # We read the time course (x, y, z and the time points)
    tc_vox = data4D[vox[0], vox[1], vox[2], :]
    # Entropy is calculated on the histogram of the distribution
    counts, bin_edges = np.histogram(tc_vox, bins=nb_bins)
    # Tranform counts into probability
    p = counts / np.sum(counts, dtype=float)
    binWidth = np.diff(bin_edges)
    # xlogy is with natural logarithm but we should use log base 2 (TO STUDY !!)
    H = -np.sum(xlogy(p, p/binWidth))

    data_H[vox[0], vox[1], vox[2]] = H

image_H = nib.Nifti1Image(data_H, brainmask.affine)
nib.save(image_H, sys.argv[3])
# nib.save(image_H, "result_H.nii.gz")
