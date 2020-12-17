# -*- coding: utf-8 -*-

import os
import numpy as np

import nibabel as nib


AFFINE_MNI_1MM = np.array([[  -1.,    0.,    0.,   90.],
       [   0.,    1.,    0., -126.],
       [   0.,    0.,    1.,  -72.],
       [   0.,    0.,    0.,    1.]])

AFFINE_MNI_2MM = np.array([[  -2.,    0.,    0.,   90.],
       [   0.,    2.,    0., -126.],
       [   0.,    0.,    2.,  -72.],
       [   0.,    0.,    0.,    1.]])

SHAPE_MNI_1MM = tuple((182, 218, 182))
SHAPE_MNI_2MM = tuple((91, 109, 91))


def empty_MNI1MM():
    data = np.zeros(SHAPE_MNI_1MM)
    return nib.Nifti1Image(data, AFFINE_MNI_1MM)

def empty_MNI2MM():
    data = np.zeros(SHAPE_MNI_2MM)
    return nib.Nifti1Image(data, AFFINE_MNI_2MM)

def get_data_folder():
    return os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(__file__))), "Data")

def get_ants_priors_folder():
    return os.path.join(get_data_folder(), "ants_priors")
