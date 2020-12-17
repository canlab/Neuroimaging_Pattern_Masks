# -*- coding: utf-8 -*-
import numpy as np
import sys

import nibabel as nib
import pandas as pd
from sklearn import decomposition
import sklearn.covariance as sco
import tools.mat_transform as mt
# import mat_transform as mt

path = "/data/neurosynth_data/Wholebrain.csv"
path2 = "/data/neurosynth_data/unrotated_components.csv"
mat = pd.read_csv(filepath_or_buffer=path, header=0)

pmat = mat[mat.columns[4:]]
cov = pmat.cov()
type(cov)
g, omega = np.linalg.eig(cov)
# So g is the vector containing the eigenvalues and SPSS take their absolute val
g = np.abs(g)
# According to SPSS gamma is a diagonal matrix
gamma = np.diag(g)
# Communalities (it is a vector)
# h0 = np.sum(g*(omega[0]**2))
# len(range(len(omega)))
h = [np.sum(g*(omega[i]**2)) for i in range(len(omega))]
np.testing.assert_almost_equal(np.diag(cov), h, decimal=5)
np.isclose(np.diag(cov), h, rtol=1e-07).all()
diag_h = np.diag(h)
# The matrix of factor laodings
lam = np.dot(omega, np.sqrt(gamma))
pd.DataFrame(lam)
norm_lam = np.dot(np.sqrt(diag_h), lam)
pd.DataFrame(norm_lam)
lam[:,0]
lam.shape
diag_h.shape
norm_lam.shape
n = len(h)
m = len(norm_lam)

# svi = np.sum(n * np.sum())
def n_pow_four(mat, j, n):
    return n * np.sum(mat[:,j]**4)

def pow2_pow2(mat, j, n):
    return np.sum(mat[:,j]**2)**2
len(range(10))
len(np.array([[1,2,3],[4,5,6]]))
norm_lam[1]
for j in range(m):
    n * np.sum(norm_lam[:,j]**4) - np.sum(norm_lam[:,j]**2)**2

svi = np.sum([n * np.sum(norm_lam[:,j]**4) - np.sum(norm_lam[:,j]**2)**2  for j in range(len(norm_lam))])/(n**2)
# corr = pmat.corr()
# test_cov = np.cov(np.array(pmat).T)
# sys.getsizeof(cov)
# tmat = pmat.T
# terms = tmat.index
# tcov = tmat.cov()
# skcov = sco.empirical_covariance(np.array(pmat.T))
#
# raw = np.array([[ 10, 10, 10],
#                 [ 5, 20, 30],
#                 [ 6, 15, 7],
#                 [ 6, 21, 29],
#                 [ 31, 11, 4]])
# cov = np.cov(raw)
# np.dot(raw, raw.T)/5
# matrix = np.array([[ 0.05,  0.05,  0., 0.75, 0.5],
#                 [ 0.05,  0.95,  0.05, 0.5, 0.05],
#                 [ 0.05,  0.95,  0.05, 0.5, 0.05],
#                 [ 0.25,  0.5,  0.75, 0.5, 0.25],
#                 [ 0.25,  0.25,  0.25, 0.25, 0.25]])
# sco.empirical_covariance(matrix)
# np.cov(matrix)
# pd.DataFrame(matrix).cov()
# matrix[2, 1:]
# matrix = cov
# path_pref = "osef"
# rot='quartimax'
d, v = np.linalg.eig(cov)
# Sort values avec vectors in decreasing order
indsorted = d.argsort()[::-1]
d = d[indsorted]
v = v[:,indsorted]
l = v * np.sqrt(d)
l.shape
po = np.power(l, 2)

com_full = np.sum(po, axis=0)

eigval_thr = 1
comp_thr = eigval_thr * np.mean(d)
ind = 0
while d[ind] > comp_thr:
    ind += 1
# we have the number of values > comp_thr but the index of the last is minus 1
# ind -= 1
l_trunc = l[:,0:ind]
l_trunc.shape
d[68]
comp_thr
pd.DataFrame(l_trunc)

com = np.sum(np.power(l_trunc, 2), axis=0)

com_scaled = com/com_full[0:len(com)]
# We could run loads of PCA and then use what Slava did to correlate the values
# across the trials to display how much confident we are for each component


def parcellate_PCA(matrix, mat_type, path_pref, rot='quartimax', eigval_thr=1):
        """ Parellate a 2D similarity matrix with the PCA algorithm
        Parameters
        ----------
        similarity_matrix: 2D np.array
            square similarity_matrix, e.g. correlation matrix
            (It is assumed that the original 2D_connectivity_matrix was
            normalized across ROIs)
        mat_type: str ['covariance', 'correlation']
            the type of similarity matrix. The threshold of eigenvalues will be
            eigval_thr for correlation matrices but
            eigval_thr * mean(eigevalues) for covariance matrices
        rot: str ['quartimax', 'varimax']
            Type of factor rotation
        path_pref: str
            the path and the prefixes to add before the name of files which
            will be created by the function
        eigval_thr: int

        Returns
        -------
        labels: np.array
            Nseed labels (integers) which can be used to assign to each seed
            ROI the value associated to a certain cluster
        """
        if rot == 'quartimax':
            rotation = 0.0
        elif rot == 'varimax':
            rotation = 1.0
        else:
            raise Exception('This factor rotation type is not handled')

        # To have more than just a reference of matrix  in mat
        mat = matrix + 0

        # tests
        mat = cov
        rotation = 1.0
        mat_type = 'covariance'
        eigval_thr = 1
        # Get the eigenvalues and eigenvectors of the
        # mat = cov(2D_connectivity_matrix)
        gamma_eigval, omega_eigvec = np.linalg.eig(mat)
        # # Just in case, we do the absolute value of the eigenvalues
        # gamma_eigval = np.abs(gamma_eigval)
        # Calculate the threshold to reduce the number of components we will
        # keep.
        if mat_type == "covariance":
            comp_thr = eigval_thr * np.mean(gamma_eigval)
        elif mat_type == "correlation":
            comp_thr = eigval_thr
        else:
            raise Exception('This factor rotation type is not handled')

        unro_cmp = pd.read_csv(filepath_or_buffer=path2, header=None)
        npunrot = np.array(unro_cmp)

        pd.DataFrame(npunrot)

        # Sort the Gamma_eigval in decreasing order of magnitude, and sort
        # the order of the eigenvectors accordingly
        indsort = np.argsort(gamma_eigval)[::-1]

        # The SSQ_loadings is equal to the eigenvalues of the SM (cov(data))
        # They correspond to the values in the 'Extraction Sum of Squared
        # loadings' in SPSS
        gamma_eigval_sort = gamma_eigval[indsort]
        omega_eigvec_sort = omega_eigvec[:,indsort]

        # We keep only the components which have an eigenvalue above comp_thr
        keep = np.where(gamma_eigval_sort > comp_thr)
        ind = 0
        while gamma_eigval_sort[ind] > comp_thr:
            ind += 1
        gamma_eigval_sort = gamma_eigval_sort[:ind]
        omega_eigvec_sort = omega_eigvec_sort[:,:ind]

        SSQ_loadings = gamma_eigval_sort

        # Force all eigenvalues to be > 0, as rounding errors can yield very
        # small eigenvalues to be < 0
        # https://de.mathworks.com/matlabcentral/newsreader/view_thread/322098
        # ind_negative_eigenvalues = np.where(gamma_eigval_sort < 0)
        # gamma_eigval_sort[ind_negative_eigenvalues] = np.abs(
        #     gamma_eigval_sort[ind_negative_eigenvalues])

        pd.DataFrame(omega_eigvec_sort)
        # --------------------------------------------------------------------------
        # (2) Calculate Factor loadings Lambda = Eigvecs * sqrt(Eigvals)
        # These correspond to the values in the 'Component Matrix' of SPSS
        # --------------------------------------------------------------------------
        Lambda = omega_eigvec_sort.dot(np.diag(np.sqrt(gamma_eigval_sort)))


        lam = pd.DataFrame(Lambda)
        # --------------------------------------------------------------------------
        # (3-4) Perform Orthomax rotation of Lambda -> Lambda_rot
        # Returns the Lambda_rot and its Sum of Squared loadings
        # This correspond to the values in the 'Rotated Component Matrix' of SPSS
        # up to a constant rotation factor.
        # --------------------------------------------------------------------------
        # Rotate factor loadings using the function in do_PCA_utilities.py
        # https://en.wikipedia.org/wiki/Talk:Varimax_rotation
        lambda_rot = mt.rotate_components(Lambda, q = 1000, gamma=rotation)
        pd.DataFrame(lambda_rot)
        # Get sum of squared loadings
        SSQ_loadings_rot = np.sum(lambda_rot**2, axis=0)
        # print(np.sort(SSQ_loadings_rot)[::-1])

        # Sort the SSQ_loadings_rot in descending order to prepare for the
        # power fitting
        SSQ_loadings_rot_sorted = np.sort(SSQ_loadings_rot)
        SSQ_loadings_rot_sorted_descending = SSQ_loadings_rot_sorted[::-1]
        # plt.plot(SSQ_loadings_rot_sorted_descending); plt.show()

        # --------------------------------------------------------------------------
        # (5) Fit a power law to the sorted SSQ_Loadings_rot to Estimate
        #     the number of relevant factors Npc using the fitpower function in
        #     do_PCA_utilities.py (only the first 50 SSQ_Loadings are considered).
        # Returns the number of components to consider: Npc
        # --------------------------------------------------------------------------
        npc = mt.fit_power(SSQ_loadings_rot_sorted_descending)
        print('\n Power fitting of the eigenvalues associated with the rotated')
        print('loadings estimated the presence of ' + str(npc) + ' clusters \n')



        # --------------------------------------------------------------------------
        # (6) Rotate Lambda_Npc = Lambda[:,Npc]
        # Returns the final Factor loadings, defining the clusters
        # --------------------------------------------------------------------------
        lambda_npc = Lambda[:, 0:npc]
        pd.DataFrame(lambda_npc)
        lambda_npc_rot = mt.rotate_components(lambda_npc, gamma=rotation)

        pd.DataFrame(lambda_npc_rot)


        # --------------------------------------------------------------------------
        # (7) Sort for visualization and return cluster labels
        # --------------------------------------------------------------------------
        # Sort the clusters of ROIs and count the number of ROIs per clusters
        absV = np.abs(lambda_npc_rot)
        # plt.imshow(PC, aspect='auto');
        # plt.show()

        # Find the maximum value for each ROI across principal components.
        # In other words, find to which principal component / cluster each ROI
        # should be assigned
        labels = np.argmax(absV,axis=1)

        ROI_clu = np.zeros(lambda_npc_rot.shape)

        # Display the assignment of each ROIs to its pc/cluster
        for i in np.arange(npc):
            factor_idx = np.where(labels == i)
            ROI_clu[factor_idx, i] = 1

        # plt.imshow(ROI_clu, aspect='auto', interpolation='none');
        # plt.show()

        # 'labels' is already the vector of cluster number for each ROI
        # Here we just visualize it as in SPSS
        ROI_clu_sort = np.zeros(ROI_clu.shape)

        startrow = 0

        for i in np.arange(ROI_clu_sort.shape[1]):
            tmp = np.where(ROI_clu[:,i])
            nROIs_ith_clu = np.array(tmp).shape[1]
            stoprow = startrow + nROIs_ith_clu
            ROI_clu_sort[startrow : (nROIs_ith_clu + startrow), i] = 1
            startrow = stoprow

        # plt.imshow(ROI_clu_sort, aspect='auto', interpolation='none');
        # plt.show()
        np.save(path_pref + 'ROI_clu.npy', ROI_clu)
        np.save(path_pref + "ROI_clu_sort.npy", ROI_clu_sort)

        # Show how clusters look like in the similarity_matrix
        idxsort = np.argsort(labels)
        mat_clusters = mat[idxsort,:][:,idxsort]

        return labels
