# -*- coding: utf-8 -*-
import os
import errno
import numpy as np

import nibabel as nib


def float_in_filename(num):
    """ Function to remove the "." of the float to avoid it in filenames.
    It convert the Numeric input into a str
    If the input is simply an int, it just converts it into str
    Parameters
    ----------
    num: Numeric
    Returns:
    str
        the string "dot" will replace the dot of the float
    """
    spl = str(num).split(".")
    if len(spl) > 1:
        return spl[0] + "dot" + spl[1]
    elif len(spl) == 1:
        return spl[0]
    else:
        print("uknown value: " + str(num))


def split_clusters(nii, res_folder, name):
    """ Split a 3D image into several nifti files. Each one containing a
    distinct value from the source image.
    Parameters
    ----------
    nii: nibabel.Nifti1Image
        the image to split
    res_folder: str
        the path to the result folder. The function will create a folder
        which will contain all the splitted images
    name: str
        suffix for the result files, it will also be used to name the result
        folder
    """
    # extract the needed informations from the source image
    data = nii.get_data()
    affine = nii.affine

    folder = os.path.join(res_folder, name)
    # Try to create the folder and ignore the error in the case it alread exists
    try:
        os.mkdir(folder)
    # note that all the other errors like permissions error will be caught
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    # TODO change the amax for a walk through an array of the unique values
    # find the maximum value of the source
    o_max = np.amax(data)
    if np.amin(data) < 0:
        print("The source image contains negative values")
        return

    while o_max > 0:
        # we create an empty (full of 0) array with the same shape as data
        mask = np.zeros(data.shape)
        # save the cluster in mask with its original value
        mask[np.where(data == o_max)] = o_max
        # we remove the cluster from data
        data[np.where(data == o_max)]  = 0
        # we save the mask
        img_ROIs = nib.Nifti1Image(mask, affine)
        path = os.path.join(
            folder, "clu_" + float_in_filename(o_max) + "_" + name + ".nii.gz")
        nib.save(img_ROIs, path)
        # The maximum of the remaining values of data
        o_max = np.amax(data)
    print("All the cluster has been splitted in " + folder)
    return
    # clu = np.array(np.where(data == i))
    # mask[clu[0,], clu[1,], clu[2,]] = i
    # tt.append(len(clu))
