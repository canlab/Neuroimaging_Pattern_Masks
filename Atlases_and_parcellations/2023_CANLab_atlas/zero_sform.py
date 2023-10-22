#!/usr/bin/env python

# this code is copied from ChatGPT and zeros out sform 
# from a NIFTI file's header. This is required for use
# of nifti files as custom atlases with qsirecon. See
# here under "Using custom atlases",
# https://qsiprep.readthedocs.io/en/latest/reconstruction.html

import nibabel as nib
import numpy as np

# Load the NIFTI file
nii = nib.load('canlab_2023_1mm.nii.gz')

# Get the affine matrix for the sform
sform_affine = nii.get_sform()

# Print the original sform data
print("Original sform:")
print(sform_affine)

# Zero out the sform data
nii.set_sform(np.zeros((4,4)), code=0)

# Check if the sform data has been zeroed out
new_sform_affine = nii.get_sform()
print("New sform:")
print(new_sform_affine)

# Save the new NIFTI file
nib.save(nii, 'canlab_2023_MNI152NLin2009cAsym_1mm_lps.nii.gz')

