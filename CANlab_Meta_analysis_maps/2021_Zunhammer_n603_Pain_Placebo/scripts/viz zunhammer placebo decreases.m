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