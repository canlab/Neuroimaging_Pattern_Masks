% This function applies the S1-connectivity-based marker described in Lee et al 2018 J Pain to
% fMRI data, most commonly a resting state run.
% 
% input: 
%   dat: a 4D fmri_data object, most commonly a (preprocessed, denoised)
%   resting state run
%
% output:
%   marker_val: the marker value
function marker_val = apply_LeeCBP_S1_marker(dat)

    % load in seed and weight map
    s1back_roi = atlas(which('S1back_Lee_2018coords.nii'), 'labels', 'S1back', 'space_description', 'MNI152 space');
    LeeS1_SVM = fmri_data(which('nilearn_pairedSVM_W_S1conn53.nii'));
    LeeS1_intercept_fname = fullfile(fileparts(which('nilearn_pairedSVM_W_S1conn53.nii')), 'intercept.txt');
    LeeS1_intercept = textread(LeeS1_intercept_fname, '%f');

    % load subject's data and resample space
    b = brainpathway(s1back_roi);    
    dat_rs = resample_space(dat, b.region_atlas);

    % load into brainpathway object
    b.voxel_dat = dat_rs.dat;

    % estimate connectivity to S1 seed
    connectivity_map = seed_connectivity(b, {'S1back'});

    % apply marker
    marker_val = apply_mask(connectivity_map, LeeS1_SVM, 'pattern_expression') + LeeS1_intercept;

end