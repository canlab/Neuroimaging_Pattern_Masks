% predict Craving with/without bootstrapping (10000 samples)
% and generate ROC plots
% Leonie Koban, 2022

clear all
close all
clc

%%
load(which('AllCues_NCS_fMRIdata_byCravingLevel.mat'))

%% cross-validated LASSO-PCR

[cverr, stats, optout] = predict(Craving_5levels, 'algorithm_name', 'cv_lassopcr', 'nfolds', Craving_5levels.metadata_table.CV10Folds);
% [cverr, stats, optout] = predict(bigdat99, 'algorithm_name', 'cv_lassopcr', 'nfolds', Craving_5levels.metadata_table.CV10Folds, 'error_type', 'mse', 'bootweights', 'bootsamples', 10000);

%% visualization

ncs = stats.weight_obj;
canlab_results_fmridisplay(ncs) % unthresholded display of the weight map
% canlab_results_fmridisplay(threshold(ncs, 0.05, 'fdr'))     % thresholded weight maps (require bootstrapping)
% canlab_results_fmridisplay(threshold(ncs, 0.005, 'unc'))

%% ROC plots and stats

figure; 
scatter(stats.Y, stats.yfit, 3, 'k.'); hold on
lineplot_columns({stats.yfit(stats.Y==1), stats.yfit(stats.Y==2), stats.yfit(stats.Y==3), stats.yfit(stats.Y==4), stats.yfit(stats.Y==5)}, 'shade'); hold on
xlim([0.5 5.5]); ylim([-2.5 8])

% find lowest and highest Y for each subject
for s = 1:99
    subj_data_results.Y{s} = stats.Y(Craving_5levels.metadata_table.SubjNumber==s);
    subj_data_results.Yfit{s} = stats.yfit(Craving_5levels.metadata_table.SubjNumber==s);
    
    Ymin(s,1) = min(subj_data_results.Y{s});
    YminFit(s,1) = subj_data_results.Yfit{s}(subj_data_results.Y{s}==min(subj_data_results.Y{s}));
    Ymax(s,1) = max(subj_data_results.Y{s});
    YmaxFit(s,1) = subj_data_results.Yfit{s}(subj_data_results.Y{s}==max(subj_data_results.Y{s}));
end

figure('Name', 'ROC Craving levels discrimination', 'Color', [1 1 1])
% Single interval
LASSO_ROC_5_1 = roc_plot([stats.yfit(stats.Y==5); stats.yfit(stats.Y==1)], [stats.Y(stats.Y==5); stats.Y(stats.Y==1)]==5, 'color', 'b'); 
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	3.37	Sens:	 64% CI(53%-75%)	Spec:	 80% CI(71%-87%)	PPV:	 74% CI(64%-84%)	
% Nonparametric AUC:	0.76	Parametric d_a:	0.94	  Accuracy:	 72% +- 3.4% (SE), P = 0.000000

% forced/two-choice
LASSO_ROC_Max_Min = roc_plot([YmaxFit; YminFit], [ones(numel(Ymax),1); -ones(numel(Ymin),1)]==1, 'color', 'k', 'twochoice'); 
set(gca, 'FontSize', 18, 'LineWidth', 2); box('off')
% ROC_PLOT Output: Two-alternative forced choice, A priori threshold
% Threshold:	0.00	Sens:	 81% CI(73%-89%)	Spec:	 81% CI(73%-89%)	PPV:	 81% CI(73%-89%)	
% Nonparametric AUC:	0.91	Parametric d_a:	1.27	  Accuracy:	 81% +- 4.0% (SE), P = 0.000000
