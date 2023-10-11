% subcortex+cortex atlas are in CIFTI file format. So will load in
% subcortical files only. Will load in atlases in 2009cAsym space, as this
% is the standard space used by fMRIprep

% atlas file names
basedir = fileparts(which('Tian_Subcortex_S4_3T_2009cAsym.nii.gz'))

atlas_fnames = filenames(fullfile(basedir, 'Tian_Subcortex_S*_3T_2009cAsym.nii.gz'))

label_fnames = filenames(fullfile(basedir, 'Tian_Subcortex_S*_label.txt'))

reference = 'Tian et al., 2021, Nat Neuro; see https://github.com/yetianmed/subcortex';

% scales 1 to 4
opts = {'references', reference, 'space_description', 'MNI152 Nonlinear 2009cAsym'};
    
tian_subcortical_S1 = atlas(atlas_fnames(1), 'labels', importdata(label_fnames{1}), opts{:});
tian_subcortical_S2 = atlas(atlas_fnames(2), 'labels', importdata(label_fnames{2}), opts{:});
tian_subcortical_S3 = atlas(atlas_fnames(3), 'labels', importdata(label_fnames{3}), opts{:});
tian_subcortical_S4 = atlas(atlas_fnames(4), 'labels', importdata(label_fnames{4}), opts{:});
    
%% save to mat files
save(fullfile(basedir, 'Tian2020_subcortical_S1.mat'), 'tian_subcortical_S1')
save(fullfile(basedir, 'Tian2020_subcortical_S2.mat'), 'tian_subcortical_S2')
save(fullfile(basedir, 'Tian2020_subcortical_S3.mat'), 'tian_subcortical_S3')
save(fullfile(basedir, 'Tian2020_subcortical_S4.mat'), 'tian_subcortical_S4')
