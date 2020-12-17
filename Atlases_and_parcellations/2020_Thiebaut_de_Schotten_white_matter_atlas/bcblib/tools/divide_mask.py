# -*- coding: utf-8 -*-

import numpy as np

import nibabel as nib
from scipy.spatial import distance_matrix


def find_seed(coords, direc):
    """ Find (one of) the most distant voxels in a given direction and
    return its coordinates
    Parameters
    ----------
    coords: np.array (coordinates of each voxel on lines)
        coordinates of voxels in the mask
    direc: int
        0: x
        1: -x
        2: y
        3: -y
        4: z
        5: -z
    Returns
    -------
    ext: [(int)x, (int)y, (int)z]
        the coordinates of the most extreme voxels in the direction direc
    """
    side = direc % 2
    axis = direc // 2
    a_side = ["np.argmin", "np.argmax"]
    # print("The direction is: " + ["+", "-"][side] + ["x", "y", "z"][axis])
    ext = eval(a_side[side])(coords[:, axis])
    return coords[ext]


def gather_round(seed, coords, size):
    """ Find the nearest voxels from the seed and return an array of their
    coordinates
    Parameters
    ----------
    seed: np.array
        array([x,y,z]) an array of the coordinates of the seed voxel
        of the cluster
    coords: np.array (coordinates of each voxel on lines)
        coordinates of voxels in the mask
    size: int
        size of the cluster. If there isn't enough coordinates, the function
        will still return a cluster but with less voxels
    Returns
    -------
    np.array
        array with the indixes of the nearest voxels(in coords) from seed
        (the array contains seed)
    """
    dist_mat = distance_matrix(np.array(seed), coords)
    ind_sort = np.argsort(dist_mat[0], 0, 'mergesort')
    # neighbours = [coords[i] for i in ind_sort[0:size]]
    return np.array(ind_sort[0:size])


def divide_compactor(img, size):
    """ Cluster img in groups of a given number of neighbour voxels
    Parameters
    ----------
    img: Nifti1Image
        The nifti mask of non-zero voxels to cluster
    size: int
        The size of each clutser (The last cluster can have a lower number
        of voxels)
    Returns
    -------
    res_img: Nifti1Image
        An image with the same dimension than img and its voxels labelled with
        their cluster number
    """
    coords = np.asarray(np.where(img.get_data())).T
    res_data = np.zeros(img.get_data().shape)
    clu_lbl = 0
    while len(coords) > 0:
        direc = clu_lbl % 6
        # We want to start by the +x direction but with the cluster label 1
        clu_lbl = clu_lbl + 1
        seed = find_seed(coords, direc)
        tmp_clu = gather_round([seed], coords, size)
        for i in tmp_clu:
            v = coords[i]
            res_data[v[0], v[1], v[2]] = clu_lbl
        coords = np.delete(coords, tmp_clu, 0)
    res_img = nib.Nifti1Image(res_data, img.affine)
    return res_img
