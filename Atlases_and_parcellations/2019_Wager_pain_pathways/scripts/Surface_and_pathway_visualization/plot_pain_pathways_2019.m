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
% Uses ST_cleaned

p1 = ST_cleaned.gray_white_csf;
p2 = ST_cleaned.pdm1;
p3 = [ST_cleaned.rel_temp ST_cleaned.pain_rating];
X = [p1 p2 p3];
Xlabels = [ST_cleaned.graywhitecsf_labels ST_cleaned.pdm1_labels {'Temp' 'Pain'}];
Xpartitions = [ones(size(p1, 2), 1); 2 * ones(size(p2, 2), 1); 3 * ones(size(p3, 2), 1)]; 
partitionlabels = {'Nuisance' 'PDM1 local patterns' 'Temp Pain'};

create_figure('PDM1 correlations', 1, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
    'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

title('Full pairwise Spearman Zwithin');
drawnow

subplot(1, 2, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'dopartial', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
    'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

title('Partial r Spearman Zwithin');
drawnow

saveas(gcf, 'figures/PDM1_local_pattern_ST_connectivity_matrix_Zwithin.png');

[b, dev, stats] = glmfit(X(:, 1:end-1), ST_cleaned.pain_rating);
glm_table(stats, Xlabels(1:end-1));


% for i = 1:9 % do for brainstem regions 
%     
%     fprintf('%s connected to ', labels{i});
%     wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
%     for j = 1:length(lab), fprintf(' %s', lab{j}); end
%     fprintf('\n');
%     
% end

%%
% Correlations between raw extracted ROI values and potential nuisance components
% Question is: can we identify nuisance variables that are correlated with
% regional brain activity but not pain/temp?  And thus controlling for
% these, the relationship with pain becomes stronger. 
% We should consider whether we are dealing with raw or z-scored-within
% person/cleaned data.

X = scale(ST_cleaned.pdm1, 1);
pX = pinv(X);

pdm1fit = @(y) X * pX * y;  % no intercept... ?? assuming all vars are z-scored...?
pdm1varexplained = @(y) var(pdm1fit(y)) ./ var(y);

% behavior to nuisance 
X = scale([ST_cleaned.rel_temp ST_cleaned.pain_rating], 1);
pX = pinv(X);

paintempfit = @(y) X * pX * y;
paintempvarexplained = @(y) var(paintempfit(y)) ./ var(y);

% we want to find nuisance covariates where v1/v2 is high

for i = 1:3
    
    v1 = pdm1varexplained(ST_cleaned.gray_white_csf(:, i)); % brain to nuisance
    v2 = paintempvarexplained(ST_cleaned.gray_white_csf(:, i)); % behavior to nuisance
    fprintf('%s: brain %3.2f%% behavior %3.2f%% ratio (higher is better): %3.2f\n', ST_cleaned.graywhitecsf_labels{i}, 100*v1, 100*v2, v1/v2);
    
end

% These have already had white and CSF removed in the 'cleaned' version: 
%
% Gray_mean: brain 27% behavior 0.01% ratio (higher is better): 32.00
% White_mean: brain 00% behavior 0.00% ratio (higher is better): 0.50
% CSF_mean: brain 00% behavior 0.01% ratio (higher is better): 0.12

for i = 1:15
    
    v1 = pdm1varexplained(ST_cleaned.gray_white_csf_components(:, i)); % brain to nuisance
    v2 = paintempvarexplained(ST_cleaned.gray_white_csf_components(:, i)); % behavior to nuisance
    fprintf('%s: brain %3.2f%% behavior %3.2f%% ratio (higher is better): %3.2f\n', ST_cleaned.graywhitecsfcomponents_labels{i}, 100*v1, 100*v2, v1/v2);
    
end

% Gray1: brain 7.50% behavior 0.63% ratio (higher is better): 12.00
% Gray2: brain 1.32% behavior 0.02% ratio (higher is better): 76.13
% Gray3: brain 0.84% behavior 0.02% ratio (higher is better): 44.90
% Gray4: brain 0.58% behavior 0.02% ratio (higher is better): 27.66
% Gray5: brain 0.37% behavior 0.01% ratio (higher is better): 32.22
% White1: brain 2.05% behavior 0.14% ratio (higher is better): 14.41
% White2: brain 0.68% behavior 0.08% ratio (higher is better): 8.00
% White3: brain 0.88% behavior 0.13% ratio (higher is better): 6.60
% White4: brain 0.51% behavior 0.06% ratio (higher is better): 7.93
% White5: brain 0.32% behavior 0.01% ratio (higher is better): 26.38
% CSF1: brain 1.68% behavior 0.75% ratio (higher is better): 2.23
% CSF2: brain 1.32% behavior 0.02% ratio (higher is better): 57.36
% CSF3: brain 0.69% behavior 0.01% ratio (higher is better): 101.02
% CSF4: brain 0.41% behavior 0.01% ratio (higher is better): 38.63
% CSF5: brain 0.31% behavior 0.03% ratio (higher is better): 11.35

% All components seem to be viable nuisance components - except maybe CSF1
X = scale(ST_cleaned.gray_white_csf_components(:, [1:10 12:15]), 1);
X(:, end + 1) = 1; % intercept
remove_nuisance = @(y) y - X * pinv(X) * y;

%% PDM1, removing nuisance

X = scale(ST_cleaned.pdm1, 1);
X = remove_nuisance(X);

p1 = ST_cleaned.gray_white_csf;
p2 = X;
p3 = [ST_cleaned.rel_temp ST_cleaned.pain_rating];
X = [p1 p2 p3];
Xlabels = [ST_cleaned.graywhitecsf_labels ST_cleaned.pdm1_labels {'Temp' 'Pain'}];
Xpartitions = [ones(size(p1, 2), 1); 2 * ones(size(p2, 2), 1); 3 * ones(size(p3, 2), 1)]; 
partitionlabels = {'Nuisance' 'PDM1 local patterns' 'Temp Pain'};

create_figure('PDM1 correlations', 1, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
    'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

title('Full pairwise Spearman Zwithin GWCSFcomp removed');
drawnow

subplot(1, 2, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'dopartial', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
    'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

title('Partial r Spearman Zwithin GWCSFcomp removed');
drawnow

saveas(gcf, 'figures/PDM1_local_pattern_ST_connectivity_matrix_Zwithin_GWCSFcompRem.png');

[b, dev, stats] = glmfit(X(:, 1:end-1), ST_cleaned.pain_rating);
glm_table(stats, Xlabels(1:end-1));

% Multiple correlation: PDM1 without and with components removed
X = scale(ST_cleaned.pdm1, 1);

pdm1fit = @(y) X * pinv(X) * y;  % no intercept... ?? assuming all vars are z-scored...?
pdm1varexplained = @(y) var(pdm1fit(y)) ./ var(y);

pdm_explained1 = pdm1varexplained(ST_cleaned.pain_rating);

X = remove_nuisance(X);
pdm1fit = @(y) X * pinv(X) * y;  % no intercept... ?? assuming all vars are z-scored...?
pdm1varexplained = @(y) var(pdm1fit(y)) ./ var(y);

pdm_explained2 = pdm1varexplained(ST_cleaned.pain_rating);

fprintf('Var explained in pain before = %3.2f%% and after GWCSF component removal = %3.2f%%\n', 100*pdm_explained1, 100*pdm_explained2);
% Var explained in pain before = 7.11% and after GWCSF component removal = 6.72%


%%

      c = []; c.covs_nointerest = [SETUP.X SETUP.covariates]; c.outcome = SETUP.Y; doranks = 1;
      for i = 1:length(cl), c.names{i} = cl(i).shorttitle; end
      c = nmdsfig_tools('get_data_matrix',c,cl,'timeseries',1,[],doranks);
      c = nmdsfig_tools('get_correlations',c);
      [c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);
      c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:5, 1000, []);
      c = nmdsfig_tools('apply_clusters',c);
      c = nmdsfig_tools('predict_behavior',c,'classes');
      cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'sagittal', 0);
      nmdsfig_tools('nmdsfig_plot',c, 0, 0, 'fill');
  c_complete = nmdsfig_tools('nmdsfig_plot_withcovs',c, cl, doranks)
  
  %%
  remove_nuisance = @(y) y - X * pinv(X) * y;
X = scale(ST_cleaned.pdm1, 1);
X = remove_nuisance(X);

  clustercolors = scn_standard_colors(max(Xpartitions));
  
  [stats, handles] = canlab_force_directed_graph(X, 'names', Xlabels, ...
'rset', Xpartitions, 'setcolors', clustercolors, ...
'sizescale', 'custom', 'sizes', 18 * ones(size(X, 2), 1), 'LineWidth', 2, ...
'cl', pain_regions_pdm1, 'nograph');


%% Examine associations: Fine-grained atlas regions

rmat = partialcorr(ST_cleaned.fine_regions, 'type', 'Spearman');
labels = format_strings_for_legend(pain_pathways_finegrained.labels);

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('Fine regions inter-region single-trial correlations');
saveas(gcf, 'figures/Fine_regions_local_pattern_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% VPLVPM R connected to  VPLVPM R VPLVPM L IntralamMidline M pbn R rvm R Amygdala CM  R *Ctx RI R Ctx Ig R Ctx PoI2 R Ctx 5mv R Ctx 2 R
% VPLVPM L connected to  VPLVPM R VPLVPM L IntralamMidline M Amygdala CM  L Amygdala AStr  L Ctx RI L Ctx Ig L Ctx p24pr L Ctx 2 L
% IntralamMidline M connected to  VPLVPM R VPLVPM L IntralamMidline M Thal MD M Hythal Bstem PAG Ctx RI R Ctx Ig L Ctx Ig R Ctx AVI R Ctx 1 R Ctx 3a L
% Thal MD M connected to  IntralamMidline M Thal MD M Amygdala CM  L Amygdala CM  R Amygdala SF  L Ctx OP1 R Ctx a24pr L Ctx a32pr R
% Hythal connected to  IntralamMidline M Hythal pbn R pbn L Bstem PAG Amygdala SF  R Amygdala SF  L Ctx 45 L Ctx AVI L Ctx a24pr L Ctx a32pr R
% pbn R connected to  VPLVPM R VPLVPM L Hythal pbn R pbn L Bstem PAG rvm R Amygdala AStr  R Ctx PFcm R Ctx PoI2 R Ctx MI L Ctx AVI R Ctx p32pr L
% pbn L connected to  VPLVPM L Hythal pbn R pbn L Bstem PAG rvm R Amygdala LB  L Ctx PoI2 L Ctx FOP3 R Ctx AVI L Ctx a24pr L Ctx 1 L
% Bstem PAG connected to  IntralamMidline M Hythal pbn R pbn L Bstem PAG Amygdala CM  R Amygdala AStr  L Ctx PFcm L Ctx FOP2 L Ctx 33pr R Ctx a24pr R
% rvm R connected to  VPLVPM R VPLVPM L pbn R pbn L rvm R Ctx OP1 L Ctx RI L Ctx PoI2 R Ctx FOP5 L Ctx p32pr R Ctx 5mv R

create_figure('r');
imagesc(rmat(1:9, 10:end), [-.4 .4])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:9,  'YTickLabel', labels, 'XTick', 1:length(rmat)-9, 'XTickLabel', labels(10:end), 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)

%% cluster non-brainstem/thal regions by associations with bstem/thal
X = rmat(1:9, 10:end)';
c = clusterdata(X,'linkage','ward','savememory','on','maxclust',6);
[c_sorted, indx] = sort(c, 'ascend');

create_figure('r');
imagesc(X(indx, :)', [-.4 .4])
axis tight
mylabels = labels(10:end); mylabels = mylabels(indx);

set(gca, 'YDir', 'reverse', 'YTick', 1:9,  'YTickLabel', labels, 'XTick', 1:length(rmat)-9, 'XTickLabel', mylabels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)

xbreaks = find(diff(c_sorted)) + 0.5;
for i =  1:length(xbreaks)
    plot_vertical_line(xbreaks(i));
end

title('Fine regions inter-region single-trial correlations');
saveas(gcf, 'figures/Fine_regions_ST_connectivity_matrix_brainstem_detail.png');

%%

x_cell = canlab_mat2cell(x, s)
out = cellfun(@(x) corr(x, 'Type', 'Spearman'), x_cell, 'UniformOutput', false);
r3d = cat(3, out{:});
out = ttest3d(r3d);
mean(~out.sig) = 0; % get rid of non-sig edges

mean = (mean' + mean) ./ 2;
A = graph(mean, ST_cleaned.fine_labels, 'OmitSelfLoops');

A.Nodes = table(ST_cleaned.fine_labels', 'VariableNames', {'Labels'});

plot(A)

figure; plot_correlation_matrix(mean, 'docircles', false, 'doimage', true, 'dotext', false);

