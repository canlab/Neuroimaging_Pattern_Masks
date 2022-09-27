function [pexp] = apply_ncs(dat)
% applies NCS to fmri_data object dat
% Leonie Koban, 2022

ncs_filepath = which('NCS_craving_wmapN99_boot10K_02-May-2022.img');
if isempty(ncs_filepath)
    error('NCS weight map not found')
end

dat = rescale(dat, 'l2norm_images');

pexp = apply_mask(dat, ncs_filepath, 'pattern_expression', 'ignore_missing') + 3.029; % apply pattern plus intercept (3.029)

