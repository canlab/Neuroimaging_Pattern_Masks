import pandas
import os
import sys
import numpy as np
import nibabel as ni
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
import abagen
from nilearn.datasets import fetch_atlas_aal

nii = ni.load('~/.matlab/canlab/Neuroimaging_Pattern_Masks/templates/MNI152NLin6Asym_T1_1mm.nii.gz')

bigdf, coords = abagen.get_samples_in_mask(mask=None)

# load expression data
# big file -- can take quite awhile to load!
gdf = pandas.read_csv('data/gxp_correlation_wholebrain_results_NEW.csv',index_col=0)

# harmonize indices of gene expression dataframe and meta dataframe
gdf.loc[:,'ind'] = gdf.index.values
gdf.index = gdf.well_id.values
gdf = gdf.loc[bigdf.well_id.values]
gdf.index = gdf.ind.values

# test
assert all(gdf.well_id.values == bigdf.well_id.values)

# harmonize indices cont.
bigdf.index = gdf.index
bigdf.drop('well_id',axis=1,inplace=True)

from scipy.interpolate import griddata

pcamod = PCA(n_components=100, random_state=123).fit(bigdf)

# Transform gene expression data based on PCA model
pca_tfm = pandas.DataFrame(pcamod.transform(bigdf),index = gdf.index)

pls_mod = PLSRegression(n_components=3)
gdf.loc[:,'abs_mni_nlin_x'] = abs(gdf.mni_nlin_x.values) # this creates symmetry around midline
full_y = gdf[['mni_nlin_y','mni_nlin_z','abs_mni_nlin_x']]
pls_mod.fit(pca_tfm, full_y)

data = nii.get_fdata()
affine = nii.affine

brain_mask = data != 0

x_grid = np.arange(data.shape[0])  # Voxel indices along x
y_grid = np.arange(data.shape[1])  # Voxel indices along y
z_grid = np.arange(data.shape[2])  # Voxel indices along z

X_grid, Y_grid, Z_grid = np.meshgrid(x_grid, y_grid, z_grid, indexing="ij")
voxel_indices = np.vstack([X_grid.ravel(), Y_grid.ravel(), Z_grid.ravel(), np.ones(X_grid.size)]).T
world_coordinates = np.dot(affine, voxel_indices.T).T[:, :3]  # Drop homogeneous coordinate
Y_prime = world_coordinates  # Regular grid in world space

Y1 = np.vstack([full_y['abs_mni_nlin_x'], full_y['mni_nlin_y'], full_y['mni_nlin_z']]).T
Y2 = np.vstack([-1*full_y['abs_mni_nlin_x'], full_y['mni_nlin_y'], full_y['mni_nlin_z']]).T
Y = np.vstack([Y1, Y2])
X = np.vstack([pls_mod.x_scores_, pls_mod.x_scores_])

# Interpolate the irregular data (Y, X) onto the regular grid (Y')
X_prime = griddata(Y, X, Y_prime, method='linear')  # Use 'linear', 'nearest', or 'cubic'

# Reshape the interpolated values back to the NIfTI grid shape
X_prime_grid = X_prime.reshape(data.shape + (3,))

X_prime_grid[np.isnan(X_prime_grid)] = 0

interpolated_nii = ni.Nifti1Image(X_prime_grid, affine)
ni.save(interpolated_nii, "transcriptomic_gradients.nii.gz")
