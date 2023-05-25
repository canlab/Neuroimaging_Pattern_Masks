% Threshold the NCS in several different ways, and save .nii files

% load NCS weight map statistic image (with P-values)
tmp = load('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/Multivariate_signature_patterns/2022_Koban_NCS_Craving/Data_code_public/NCS_weightmaps/Bootstrapping_10K/results.mat')
orthviews(stats.weight_obj)

% P < .005, k = 10+ contiguous voxels
w = apply_mask(stats.weight_obj, which('gray_matter_mask.nii'));
w = threshold(stats.weight_obj, .005, 'unc', 'k', 10);
orthviews(w)
w.fullpath = fullfile(pwd, 'NCS_gray_p005_k10.nii');
write(w);

% FDR

w = apply_mask(stats.weight_obj, which('gray_matter_mask.nii'));
w = threshold(stats.weight_obj, .05, 'fdr');
w.fullpath = fullfile(pwd, 'NCS_gray_p05fdr.nii');
write(w)

% multi-threshold

w = apply_mask(stats.weight_obj, which('gray_matter_mask.nii'));
[o2, sig, poscl, negcl] = multi_threshold(w, 'nodisplay', 'thresh', [.001 .005 .05], 'k', [5 1 1]);
orthviews(sig)
sig.dat(~sig.sig) = 0;
sig.fullpath = fullfile(pwd, 'NCS_multithr_001k5_005_05_pruned.nii');
write(sig);

