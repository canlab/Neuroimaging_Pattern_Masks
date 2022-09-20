%% control analysis: NCS pattern expression, plots by study (Fig. 3b), and individual differences
% Leonie Koban, 2018-2022

clear all
close all
clc

cd ~/Dropbox/work/BOULDER_PROJECTS/15_Craving/Results/N99_controlanalyses/BothCues
load results_N99_10folds

dsnum99 = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'F2:F100');
user = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'C2:C100');
scanner = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'G2:G100');
sex = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'H2:H100');
age = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'I2:I100');
bmi = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'J2:J100');
edu = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'K2:K100');
pexp_D = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'P2:P100');
pexp_F = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'Q2:Q100');
pexp_DF = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'R2:R100');
motion = xlsread('CV_wmaps/Table_PatternExp_byGroupCondition_L2normed_pattern_NEW_April22_N99.xlsx', 'Subjects99_new_order', 'U2:U100');

%%

figure; 
scatter(stats.Y, stats.yfit, 3, 'k.'); hold on
lineplot_columns({stats.yfit(stats.Y==1), stats.yfit(stats.Y==2), stats.yfit(stats.Y==3), stats.yfit(stats.Y==4), stats.yfit(stats.Y==5)}, 'shade'); hold on
xlim([0.5 5.5]); ylim([-2.5 8])

% find lowest and highest Y for each subject
for s = 1:99
    subj_data_results.Y{s} = stats.Y(allsubnum99==subnums99(s));
    subj_data_results.Yfit{s} = stats.yfit(allsubnum99==subnums99(s));
    
    Ymin(s,1) = min(subj_data_results.Y{s});
    YminFit(s,1) = subj_data_results.Yfit{s}(subj_data_results.Y{s}==min(subj_data_results.Y{s}));
    Ymax(s,1) = max(subj_data_results.Y{s});
    YmaxFit(s,1) = subj_data_results.Yfit{s}(subj_data_results.Y{s}==max(subj_data_results.Y{s}));
end

%% %% GLMfit multilevel to get stats

clc

disp('=============================')
disp('All data')
disp('=============================')
stats_all = glmfit_multilevel(subj_data_results.Yfit, subj_data_results.Y, [], 'robust', 'plots'); % p < 0.0001 
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-3 9]); title('');
ylim([1 5])
xlabel('Craving rating'); ylabel('NCS response')

% Cig smokers
disp('=============================')
disp('Cig users')
disp('=============================')
stats_cigU = glmfit_multilevel(subj_data_results.Yfit(dsnum99==4), subj_data_results.Y(dsnum99==4), [], 'plots'); % p = 0.02366*	
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-2 8]); title('')

% Cig controls
disp('=============================')
disp('Cig non-users')
disp('=============================')
stats_cigN = glmfit_multilevel(subj_data_results.Yfit(dsnum99==5), subj_data_results.Y(dsnum99==5), [], 'plots'); % p = 0.00000**
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-2 8]); title('')

% Alc users
disp('=============================')
disp('Alc users')
disp('=============================')
stats_alcU = glmfit_multilevel(subj_data_results.Yfit(dsnum99==1), subj_data_results.Y(dsnum99==1), [], 'plots'); % p = 0.01407	* 
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-2 8]); title('')
% xlabel('Rating'); ylabel('Predicted rating')

% Coc users
disp('=============================')
disp('Coc users')
disp('=============================')
stats_cocU = glmfit_multilevel(subj_data_results.Yfit(dsnum99==2), subj_data_results.Y(dsnum99==2), [], 'plots'); % p = 0.00095*
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-2 8]); title('')

% Coc controls
disp('=============================')
disp('Coc non-users')
disp('=============================')
stats_cocN = glmfit_multilevel(subj_data_results.Yfit(dsnum99==3), subj_data_results.Y(dsnum99==3), [], 'plots'); % p = 0.00001**
set(gca, 'LineWidth', 3, 'FontSize', 20, 'XTick', 1:5)
xlim([0.5 5.5]); ylim([-2 8]); title('')

%% test effects of person-level covariates

l2cov = meancenter([scanner, user, sex, age, bmi, edu, motion]);

stats_all_cov = glmfit_multilevel(subj_data_results.Yfit, subj_data_results.Y, l2cov, 'beta_names', {'scanner' 'user' 'sex' 'age' 'bmi' 'edu' 'motion'}); % p < 0.0001 
% only effect for USER (users smaller effect than non-users), nothing else is significant



% plot NCS signal across studies
stats_alcall = glmfit_multilevel(subj_data_results.Yfit(dsnum99==1), subj_data_results.Y(dsnum99==1), [], 'plots'); % p = 0.00001**
stats_cocall = glmfit_multilevel(subj_data_results.Yfit(dsnum99>1 & dsnum99<4), subj_data_results.Y(dsnum99>1 & dsnum99<4), [], 'plots'); % p = 0.00001**
stats_cigall = glmfit_multilevel(subj_data_results.Yfit(dsnum99>3), subj_data_results.Y(dsnum99>3), [], 'plots'); % p = 0.00001**

%% test effects of WM/CSF signal

load ~/Dropbox/work/BOULDER_PROJECTS/15_Craving/Results/N99_controlanalyses/WM_ventricles_signal/WM_CSF_signal.mat

for s = 1:99
    sub_data.Y_WMCSF{s} = [subj_data_results.Y{s}, sub_data.WM_CSF_signal{s}];
end

wmcsf_stats = glmfit_multilevel(subj_data_results.Yfit, sub_data.Y_WMCSF, meancenter(user), 'names', {'Int' 'Rating' 'WM_CSF'}); % effect of rating on NCS is super sign (0.00000), but not WMCSF signal (0.156310

%% mean NCS and ratings for Food and Drug cues in each cohort (Suppl. Fig. 6)

figure;
subplot(2,2,1)
barplot_columns({pexp_D(dsnum99==4), pexp_D(dsnum99==1), pexp_D(dsnum99==2), pexp_D(dsnum99==5), pexp_D(dsnum99==3)}, ...
    'x', [1 2 3 5 6], 'nofig', 'nostars', 'noind', 'nobars', 'color', {[0.004, 0.616, 0.851]; [0.090, 0.812, 0.612]; [0.651, 0.757, 0.290]; [0.016, 0.376, 0.525]; [0.427, 0.498, 0.200]}); title('Drug cues')
ylim([-2 8]); ylabel('NCS response'); xlabel(''); set(gca, 'FontSize', 12, 'XTickLabel', {'1a Cig Smokers' '2 Alc Users' '3a Coc Users' '1b Non-smokers' '3b Non-users'}); xtickangle(45)

subplot(2,2,2)
barplot_columns({pexp_F(dsnum99==4), pexp_F(dsnum99==1), pexp_F(dsnum99==2), pexp_F(dsnum99==5), pexp_F(dsnum99==3)}, ...
    'x', [1 2 3 5 6], 'nofig', 'nostars', 'noind', 'nobars', 'MarkerSize', 3, 'color', {[0.004, 0.616, 0.851]; [0.090, 0.812, 0.612]; [0.651, 0.757, 0.290]; [0.016, 0.376, 0.525]; [0.427, 0.498, 0.200]}); title('Food cues')
ylim([-2 8]); ylabel('NCS response'); xlabel(''); set(gca, 'FontSize', 12, 'XTickLabel', {'1a Cig Smokers' '2 Alc Users' '3a Coc Users' '1b Non-smokers' '3b Non-users'}); xtickangle(45)

load('~/Dropbox/work/BOULDER_PROJECTS/15_Craving/Results/N99_controlanalyses/mean_ratings_behav_99.mat')

subplot(2,2,3) % behavior ratings
barplot_columns({mean(meanratings_cig(:,1:2),2), mean(meanratings_alc(:,1:2),2), mean(meanratings_coc(:,1:2),2), mean(meanratings_cig_hc(:,1:2),2), mean(meanratings_coc_hc(:,1:2),2)}, ...
    'x', [1 2 3 5 6], 'nofig', 'nostars', 'noind', 'nobars', 'color', {[0.004, 0.616, 0.851]; [0.090, 0.812, 0.612]; [0.651, 0.757, 0.290]; [0.016, 0.376, 0.525]; [0.427, 0.498, 0.200]}); title('Drug cues')
ylim([-2 8]); ylabel('Craving ratings'); xlabel(''); set(gca, 'FontSize', 12, 'XTickLabel', {'1a Cig Smokers' '2 Alc Users' '3a Coc Users' '1b Non-smokers' '3b Non-users'}); xtickangle(45)

subplot(2,2,4)
barplot_columns({mean(meanratings_cig(:,3:4),2), mean(meanratings_alc(:,3:4),2), mean(meanratings_coc(:,3:4),2), mean(meanratings_cig_hc(:,3:4),2), mean(meanratings_coc_hc(:,3:4),2)}, ...
    'x', [1 2 3 5 6], 'nofig', 'nostars', 'noind', 'nobars', 'MarkerSize', 3, 'color', {[0.004, 0.616, 0.851]; [0.090, 0.812, 0.612]; [0.651, 0.757, 0.290]; [0.016, 0.376, 0.525]; [0.427, 0.498, 0.200]}); title('Food cues')
ylim([-2 8]); ylabel('Craving ratings'); xlabel(''); set(gca, 'FontSize', 12, 'XTickLabel', {'1a Cig Smokers' '2 Alc Users' '3a Coc Users' '1b Non-smokers' '3b Non-users'}); xtickangle(45)

