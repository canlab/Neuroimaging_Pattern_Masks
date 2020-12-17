# -*- coding: utf-8 -*-
import os
import six

import nibabel as nib
import numpy as np


def simple_stats(image, method='mean', axis=3):
    """
    apply a numpy function method to the data of image on the given axis
    Parameters
    ----------
    image: str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be averaged or
        the image already loaded through nibabel
    method: {'mean', 'median', 'std', 'average' ...}, optional
        method used to average the 4D image.
        most of the functions from numpy having only one array in mandatory
        parameter and an axis option will work.
        default is 'mean'
    axis: int, optional
        default is 3 (for the temporal resolution of fMRI for example)

    Returns
    -------
    output: nibabel.Nifti1Image

    Notes
    -----
    This function has been mostly created to use some of the functions from
    the Statistics module of numpy:
    https://docs.scipy.org/doc/numpy/reference/routines.statistics.html

    Examples
    --------
    >>> import bcblib.tools.nii_stats as st
    >>> fmri_img_path = 'sub1/func/sub1_task-rest_run-1_bold.nii.gz'
    >>> out_nii = st.simple_stats(fmri_img_path, method='median')

    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        raise TypeError("Image can be either a string or a nifti1.Nifti1Image")

    data = img.get_data()
    out_data = getattr(np, method)(data, axis)

    return nib.nifti1.Nifti1Image(out_data, img.affine)
