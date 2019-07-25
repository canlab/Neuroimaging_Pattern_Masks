load data/pain_pathways_single_trial_data_2019 ST_cleaned
load pain_pathways_region_obj_with_local_patterns.mat
load pain_pathways_atlas_obj.mat

%% Questions:
%
% dividing pattern by norm(weights) makes response invariant wrt scale of weights.
% can we express intercept in units of the norm of training images in order
% to make the intercept scale-free?
%
% To preserve intensity info:
% Instead of cosine sim, can we normalize by the pattern weights and the
% WHOLE image norm? But how to deal with missing data, which will affect
% the norm?

%%
% post-windsorizing (latest)
% VPLVPM_R connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_Ig_R Ctx_PoI2_R Ctx_MI_R
% VPLVPM_L connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Amygdala_AStr__L Ctx_Ig_L
% IntralamMidline_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Ctx_a24pr_L Ctx_a32pr_R
% Thal_MD_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_a24pr_L Ctx_a24pr_R Ctx_a32pr_R
% Hythal connected to  Hythal Bstem_PAG Amygdala_CM__L Amygdala_CM__R Amygdala_SF__R Amygdala_SF__L Amygdala_LB__L
% pbn_R connected to  Hythal pbn_R pbn_L Bstem_PAG rvm_R Ctx_PoI2_R Ctx_MI_R
% pbn_L connected to  VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R Ctx_PoI2_R
% Bstem_PAG connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Hythal pbn_R pbn_L Bstem_PAG
% rvm_R connected to  pbn_R pbn_L rvm_R Ctx_PoI2_R Ctx_MI_L Ctx_MI_R Ctx_p32pr_R

% VPLVPM to ipsiateral Ctx_Ig_L Ctx_MI_L maybe Ctx_PoI2_R

% IntralamMidline_M to Ctx_a24pr_L Ctx_a32pr_R
% Thal_MD_M to Ctx_a24pr_L Ctx_a24pr_R Ctx_a32pr_R

%% Examine associations with key brainstem regions: PDM1

% rmat = corr([ST_cleaned.pdm1 ST_cleaned.pain_rating]);
% labels = format_strings_for_legend({pain_regions_pdm1.shorttitle ['pain_rating']});

%rmat = corr(ST_cleaned.pdm1);
rmat = partialcorr(ST_cleaned.pdm1, 'type', 'Spearman');
labels = format_strings_for_legend({pain_regions_pdm1.shorttitle});

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('PDM1 inter-region single-trial correlations');
saveas(gcf, 'figures/PDM1_local_pattern_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end
