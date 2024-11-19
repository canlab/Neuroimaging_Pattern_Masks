pla_fixed_neg = fmri_data('/Users/torwager/Dropbox (Dartmouth College)/SHARED_DATASETS/P_Zunhammer_Bingel_Placebo_Boulder_Essen/Analysis/G_meta_analysis_whole_brain/GIV/nii_results_tamas/full/pla/g/fixed/full_pla_g_pperm_tfce_FWE05_neg.nii.gz');

create_figure('placebo decreases fixed FWE 05 TFCE', 1, 3);
surface(pla_fixed_neg, 'colormap', 'winter');

subplot(1, 3, 2)
h = addbrain('right_insula_slab');

subplot(1, 3, 3)
h = [h addbrain('left_insula_slab')];

surface(pla_fixed_neg, 'surface_handles', h, 'colormap', 'winter');

%% Correlations with individual diffs

pla_fixed_neg_rrating = fmri_data('/Users/torwager/Dropbox (Dartmouth College)/SHARED_DATASETS/P_Zunhammer_Bingel_Placebo_Boulder_Essen/Analysis/G_meta_analysis_whole_brain/GIV/nii_results_tamas/full/pla/rrating/fixed/full_pla_rrating_pperm_tfce_FWE05_neg.nii.gz');

% pla_fixed_neg_rrating = fmri_data('/Users/torwager/Dropbox (Dartmouth College)/SHARED_DATASETS/P_Zunhammer_Bingel_Placebo_Boulder_Essen/Analysis/G_meta_analysis_whole_brain/GIV/nii_results_tamas/full/pla/rrating/fixed/full_pla_rrating_pperm_tfce_FDR_neg.nii.gz');

create_figure('placebo decreases fixed rrating FWE 05 TFCE', 1, 3);
surface(pla_fixed_neg_rrating, 'colormap', 'winter');

subplot(1, 3, 2)
h = addbrain('right_insula_slab');

subplot(1, 3, 3)
h = [h addbrain('left_insula_slab')];

surface(pla_fixed_neg_rrating, 'surface_handles', h, 'colormap', 'winter');

%% Study as fixed effect, placebo vs. control effects


pla_fixed = fmri_data('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2021_Zunhammer_n603_Pain_Placebo/placebo_effect_study_as_fixed_effect/full_pla_g_pperm_tfce_FWE05.nii.gz');

create_figure('placebo effects fixed FWE 05 TFCE', 1, 3);
surface(pla_fixed, 'colormap', 'winter');

subplot(1, 3, 2)
h = addbrain('right_insula_slab');

subplot(1, 3, 3)
h = [h addbrain('left_insula_slab')];

surface(pla_fixed, 'surface_handles', h, 'colormap', 'winter');

figure; montage(pla_fixed);

[rpos rneg] = table(region(pla_fixed));

%% Study as fixed effect, placebo vs. control correlations with analgesia


pla_fixed = fmri_data('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2021_Zunhammer_n603_Pain_Placebo/placebo_rating_correlation_study_as_fixed_effect/full_pla_rrating_pperm_tfce_FWE05.nii.gz');

create_figure('placebo effects fixed r rating FWE 05 TFCE', 1, 3);
surface(pla_fixed, 'colormap', 'winter');

subplot(1, 3, 2)
h = addbrain('right_insula_slab');

subplot(1, 3, 3)
h = [h addbrain('left_insula_slab')];

surface(pla_fixed, 'surface_handles', h, 'colormap', 'winter');

figure; montage(pla_fixed);

[rpos rneg] = table(region(pla_fixed));

