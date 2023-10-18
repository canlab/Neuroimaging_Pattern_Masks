#!/usr/bin/env python

# ironically I had trouble converting ANTs warps into SPM format in matlab because I couldn't
# easily work with 5d nifti files using spm tools, so we're doing this in python instead.

import numpy as np
import nibabel as ni

# handle fsl to fmriprep
orig_nii = ni.load('../ants/01_fsl_to_fmriprep_DisplacementFieldTransform.nii.gz')
aff = orig_nii.header.get_best_affine() # all offset by 1 for zero indexed python format

dat = orig_nii.get_fdata()

coord = np.zeros(dat.shape)
for i in range(dat.shape[0]):
    for j in range(dat.shape[1]):
        for k in range(dat.shape[2]):
            this_coord = np.array([i,j,k,1]).T
            coord[i,j,k,0,:] = this_coord[:3]

vec_dat_ijk = np.ones((np.prod(np.shape(coord)[:3]),4))
vec_dat_ijk[:,:3] = coord.reshape((-1,3))
vec_dat_mm = (aff@vec_dat_ijk.T).T
vec_dat_mm = vec_dat_mm[:,:3].reshape(coord.shape)

# Displacements are for data in LPS orientation, but SPM uses RAS, so flip x and y
dat[:,:,:,0,:2] = -1*dat[:,:,:,0,:2]
new_nii = ni.Nifti1Image(vec_dat_mm + dat, aff, header=orig_nii.header)
ni.save(new_nii,'../spm/y_01_fsl_to_fmriprep_DisplacementFieldTransform.nii')


# handle fmriprep to fsl
orig_nii = ni.load('../ants/00_fmriprep_to_fsl_DisplacementFieldTransform.nii.gz')
aff = orig_nii.header.get_best_affine() # all offset by 1 for zero indexed python format

dat = orig_nii.get_fdata()

coord = np.zeros(dat.shape)
for i in range(dat.shape[0]):
    for j in range(dat.shape[1]):
        for k in range(dat.shape[2]):
            this_coord = np.array([i,j,k,1]).T
            coord[i,j,k,0,:] = this_coord[:3]

vec_dat_ijk = np.ones((np.prod(np.shape(coord)[:3]),4))
vec_dat_ijk[:,:3] = coord.reshape((-1,3))
vec_dat_mm = (aff@vec_dat_ijk.T).T
vec_dat_mm = vec_dat_mm[:,:3].reshape(coord.shape)

# Displacements are for data in LPS orientation, but SPM uses RAS, so flip x and y
dat[:,:,:,0,:2] = -1*dat[:,:,:,0,:2]
new_nii = ni.Nifti1Image(vec_dat_mm + dat, aff, header=orig_nii.header)
ni.save(new_nii,'../spm/y_00_fmriprep_to_fsl_DisplacementFieldTransform.nii')


