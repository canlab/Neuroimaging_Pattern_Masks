#!/usr/bin/env python3
# -*-coding:Utf-8 -*
#don't forget to install Anaconda, Atom, in in atom install hydrogen, autocomplete python, python tools, python indent

#import the required library
#os for interaction with the system
#sys to script argv
#numpy to do basic maths renamed np
#nibabel to import and export nifti
#scipy to do advanced statistics
import os
import sys

import numpy as np
import nibabel as nib
from scipy.stats import spearmanr
import dask.array as da
#Defines the variable
# note that os.path.join add the path before the file

# verify if we have the right number of parameters
if len(sys.argv) < 5:
    print("usage: python funcon.py resting seed mask res_folder")
    print("resting: the path to your resting state image")
    print("seed: the path to the seed image")
    print("mask: the path to the mask image")
    print("res_folder: the path to your results folder")
    exit()

resting = sys.argv[1]
seed = sys.argv[2]
mask = sys.argv[3]
res_folder = sys.argv[4]

if os.path.exists(resting):
    rest_name = os.path.basename(resting).split(sep='.')[0]
else:
    print("ERROR: the resting state file does not exist.")
    exit()

if os.path.exists(seed):
    seed_name = os.path.basename(seed).split(sep='.')[0]
else:
    print("ERROR: the seed image file does not exist.")
    exit()

#import the nifti
i_resting = nib.load(resting)
i_seed = nib.load(seed)
i_mask = nib.load(mask)

#optional image dimension info
i_resting.shape
i_seed.shape
i_mask.shape

#check that the three shapes are the same == or different !=
i_resting.shape[:3]
if i_resting.shape[:3] != i_seed.shape:
    print("ERROR: seed and restings are not in the same space")
    exit()
if i_resting.shape[:3] != i_mask.shape:
    print("ERROR: mask and restings are not in the same space")
    exit()

#Extract the data from the images
data_resting = i_resting.get_data()
data_seed = i_seed.get_data()
data_mask = i_mask.get_data()

#Extract time course for seed in resting
#Where there is the seed in coordinates
ind_seed = np.where(data_seed)
#create an empty table for the values
seedtc = np.zeros(data_resting.shape[3])
#We read the time course (x, y, z and the time points)
for i in range(0,data_resting.shape[3]):
    #extract 3D volume out of 4D
    volume3d = data_resting[:,:,:,i]
    #average signal in seed
    seedtc[i] = np.mean(volume3d[ind_seed])

#Where is the mask in coordinates
ind_mask = np.array(np.where(data_mask))

#  Number of voxels
nvox = ind_mask.shape[1]
#create an empty table to contain spearman
table = np.zeros(i_mask.shape)

def spearman_tc(seed_tc, vox_tc):
    cor,_ = spearmanr(seed_tc, vox_tc)
    return cor

def corr(data, ind_mask, vox_ind, seed_tc):
    vox = ind_mask[:, vox_ind]
    tcvox = data_resting[vox[0], vox[1], vox[2], :]
    # test if the voxel from the mask is outside the brain in the resting state
    if all(tcvox == 0):
        # set the result to 0 (to avoid errors and nan from spearmanr function)
        cor = 0
    else:
        #calculate the Spearman
        cor = spearman_tc(seed_tc, tcvox)
        # cor,_ = spearmanr(seedtc,tcvox)
    table[vox[0], vox[1], vox[2]] = cor

for i in range(0, nvox):
    corr(data_resting, ind_mask, i, seedtc)
# delayed_task = []
# for i in range(0, nvox):
#     delayed_task.append(corr(data_resting, ind_mask, i, seedtc))
# dask.compute(delayed_task)

# # We will calculate the Spearman for EACH voxel
# for i in range(0, nvox):
#     # coordinates of the voxel in the mask
#     vox = ind_mask[:,i]
#     # We read the time course (x, y, z and the time points)
#     tcvox = data_resting[vox[0], vox[1], vox[2], :]
#     # test if tcvox is only zeros.
#     # meaning voxels outside or the resting brain image
#     if all(tcvox == 0):
#         # set the result to 0 (to avoid errors and nan from spearmanr function)
#         cor = 0
#     else:
#         #calculate the Spearman
#         cor,_ = spearmanr(seedtc,tcvox)
#     print(cor)
#     table[vox[0], vox[1], vox[2]] = cor

#write the output images
spearman = nib.Nifti1Image(table, i_mask.affine)
res_file = os.path.join(res_folder,
                        'spearman_' + rest_name + '_' + seed_name + '.nii.gz')
nib.save(spearman, res_file)
