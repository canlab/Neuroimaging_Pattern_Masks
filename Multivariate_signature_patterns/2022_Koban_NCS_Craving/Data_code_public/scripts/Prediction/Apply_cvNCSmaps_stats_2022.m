%% Analysis with pattern expression and LASSO weight map
% Leonie Koban, Oct 2018
% adapted Apr 2020
% updated Feb 2022 

clear all
close all
clc

basedir = '~/Dropbox/work/BOULDER_PROJECTS/15_Craving/A_Data_Code_FINAL_cleaned/';
cd (basedir)
%% need to use cv-fold  pattern expression values for unbiased prediction
% this is cv-pattern expression with intercept added

load([basedir, 'PatternExpressionData/pexp.mat'])


%% separate users from non-users (drug condition only)

pexp2_D_user_nuser = [pexp.D{1}; pexp.D{2}; pexp.D{4}; pexp.D{3}; pexp.D{5}];
use = [ones(17+21+21,1); -ones(18+22,1)];

figure('Color', [1 1 1], 'Name', 'Drug users vs non-users DRUG and FOOD Cues');
subplot(2,2,1);
USE_LASSO_ROC = roc_plot(pexp2_D_user_nuser, use==1, 'color', 'b');  set(gca, 'FontSize', 18);

subplot(2,2,2); 
barplot_colored([{[pexp.D{1}; pexp.D{2}; pexp.D{4}]}; {[pexp.D{3}; pexp.D{5}]}]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Drug users', 'Non-users'}, 'FontSize', 18); title('Pattern expression for DRUG cues')

% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	1.83	Sens:	 86% CI(77%-95%)	Spec:	 57% CI(41%-72%)	
% PPV:	 75% CI(65%-85%)	Nonparametric AUC:	0.76	Parametric d_a:	0.90	 
% Accuracy:	 75% +- 4.4% (SE), P = 0.002311
% in separating drug users from non-users for DRUG condition (significant)

%% separate users from non-users (FOOD condition only)

pexp2_F_user_nuser = [pexp.F{1}; pexp.F{2}; pexp.F{4}; pexp.F{3}; pexp.F{5}];
use = [ones(17+21+21,1); -ones(18+22,1)];

subplot(2,2,3);
USE_LASSO_ROC_F = roc_plot(pexp2_F_user_nuser, use==1, 'color', 'b');  set(gca, 'FontSize', 18);

subplot(2,2,4); 
barplot_colored([{[pexp.F{1}; pexp.F{2}; pexp.F{4}]}; {[pexp.F{3}; pexp.F{5}]}], 'MarkerSize', 4);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Drug users', 'Non-users'}, 'FontSize', 18); title('Pattern expression for FOOD cues')

% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-0.25	Sens:	100% CI(100%-100%)	Spec:	  0% CI(0%-0%)	PPV:	 60% CI(51%-70%)	
% Nonparametric AUC:	0.40	Parametric d_a:	-0.23	  Accuracy:	 60% +- 4.9% (SE), P = 1.000000

%% separate users from non-users (DRUG minus FOOD condition)

for ds = 1:5
    pexp.DminF{ds} = pexp.D{ds} - pexp.F{ds};
end

USE_LASSO_ROC_DmF = roc_plot([pexp.DminF{1}; pexp.DminF{2}; pexp.DminF{4}; pexp.DminF{3}; pexp.DminF{5}], use==1, 'color', 'b');  set(gca, 'FontSize', 18);
 
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-1.19	Sens:	 97% CI(91%-100%)	Spec:	 60% CI(44%-74%)	PPV:	 78% CI(67%-87%)	
% Nonparametric AUC:	0.87	Parametric d_a:	1.59	  Accuracy:	 82% +- 3.9% (SE), P = 0.000004

%% control for Age, Sex, Education

user = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'B3:B101');
scanner = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'F3:F101');
sex = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'G3:G101');
age = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'H3:H101');
bmi = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'I3:I101');
edu = xlsread([basedir, 'Metadata_Craving_20220913.xlsx'], 'Subjects_ROC', 'J3:J101');

[rage, page] = corr(cell2mat(pexp.D'), age); % r = -0.1024, p = 0.3131
[redu, pedu] = corr(cell2mat(pexp.D'), edu); % r =  -0.1054   0.2989
[rage2, page2] = corr(cell2mat(pexp.F'), age); % -0.2145, p = 0.0330
[redu2, pedu2] = corr(cell2mat(pexp.F'), edu); % r = 0. 0.1187    0.2419
[rage3, page3] = corr(cell2mat(pexp.DminF'), age); % r = 0.1041, p = 0.3051
[redu3, pedu3] = corr(cell2mat(pexp.DminF'), edu); % r = -0.2557    0.0106
[rbmi, pbmi] = corr(cell2mat(pexp.F'), bmi, 'Rows', 'pairwise'); % r = -0.2003    0.0629 (t)
[rbmi2, pbmi2] = corr(cell2mat(pexp.DminF'), bmi, 'Rows', 'pairwise'); % r = -0.0653    0.5478

% partial correlation to control for education
[rho, ppart] = partialcorr([user, cell2mat(pexp.DminF')], [age, sex, edu, scanner]); % association remains highly significant
[rho3, ppart3] = partialcorr([user, cell2mat(pexp.DminF')], [edu]);
[rho2, ppart2] = partialcorr([edu, cell2mat(pexp.DminF')], [age, sex, user, scanner]); % controlling for use, no sign. association between education and NCS

% regress out education from NCS pattern expression (Rev 2)
ncsdiff = cell2mat(pexp.DminF');
ncsD = cell2mat(pexp.D');
[b,dev,statsdiff] = glmfit([user./2, edu, age, sex], ncsdiff); %  use is significant, nothing else is sign. predictor of NCS-difference
[b2,dev2,statsD] = glmfit([user./2, edu, age, sex], ncsD); %  use is significant, nothing else is sign. predictor of NCS-difference

% mediation (education does not mediate group difference effect on NCS)
[tm, resmed] = mediation(meancenter(user), cell2mat(pexp.DminF'), edu);

%% three subplot version July 2021 (Figure 4)

pexp.us.F = [pexp.F{1}; pexp.F{2}; pexp.F{4}];
pexp.nu.F = [pexp.F{3}; pexp.F{5}];
pexp.us.D = [pexp.D{1}; pexp.D{2}; pexp.D{4}];
pexp.nu.D = [pexp.D{3}; pexp.D{5}];

colus = [.1 .1 .1];
colnu = [.5 .5 .5];

figure 
subplot(1,3,1)

plot(0.9, pexp.us.F, '.', 'Color', colus.*4); hold on
plot(1.9, pexp.us.D, '.', 'Color', colus.*4); hold on
plot(1.1, pexp.nu.F, '.', 'Color', colnu.*1.5); hold on
plot(2.1, pexp.nu.D, '.', 'Color', colnu.*1.5); hold on

line([0.95, 1.95], [nanmean(pexp.us.F), nanmean(pexp.us.D)], 'LineWidth', 2.5, 'Color', colus); hold on
line([1.05, 2.05], [nanmean(pexp.nu.F), nanmean(pexp.nu.D)], 'LineWidth', 2.5, 'Color', colnu, 'LineStyle', '--'); hold on
errorbar(0.95, nanmean(pexp.us.F), ste(pexp.us.F), 'LineWidth', 2.5, 'Color', colus); hold on
errorbar(1.95, nanmean(pexp.us.D), ste(pexp.us.D), 'LineWidth', 2.5, 'Color', colus); hold on
errorbar(1.05, nanmean(pexp.nu.F), ste(pexp.nu.F), 'LineWidth', 2.5, 'Color', colnu); hold on
errorbar(2.05, nanmean(pexp.nu.D), ste(pexp.nu.D), 'LineWidth', 2.5, 'Color', colnu); hold on
xlim([0.5 2.5]); ylim([1.5 4.5]); xlabel('Cue type', 'FontSize', 18); ylabel('NCS response', 'FontSize', 18); 
set(gca, 'XTick', 1:2, 'XTickLabel', {'Food', 'Drug'}, 'FontSize', 18, 'LineWidth', 2);  box('off')

subplot(1,3,2)
barplot_columns({pexp.us.D-pexp.us.F,  pexp.nu.D-pexp.nu.F}, 'nofig', 'nostars', 'nobars', 'MarkerSize', 3); 
set(gca, 'XTick', [1 2], 'XTickLabel', {'Users' 'Non-users'}, 'FontSize', 18, 'LineWidth', 2); box('off'); ylabel('NCS response Drug - Food'); xlabel('')
 
subplot(1,3,3)
USE_LASSO_ROC = roc_plot(pexp2_D_user_nuser, use==1, 'color', [.8 .2 .2]); hold on
USE_LASSO_ROC_F = roc_plot(pexp2_F_user_nuser, use==1, 'color', [0 .7 .85]); hold on
USE_LASSO_ROC_DmF = roc_plot([pexp.DminF{1}; pexp.DminF{2}; pexp.DminF{4}; pexp.DminF{3}; pexp.DminF{5}], use==1, 'color', [.5 0 .76]);  set(gca, 'FontSize', 18, 'LineWidth', 2); box('off')
 
%% check that prediction of USER is not driven by cocaine or cigarette or any other participants (new Suppl. Fig. 5)

figure;
% cigarette users only **********
use2 = [ones(21,1); -ones(22,1)];
subplot(1,2,1)
USE_LASSO_ROC_F2 = roc_plot([pexp.F{4}; pexp.F{5}], use2==1, 'color', [0 .7 .85]);  set(gca, 'FontSize', 18); hold on
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	4.07	Sens:	 33% CI(14%-55%)	Spec:	 86% CI(70%-100%)	PPV:	 70% CI(39%-100%)	Nonparametric AUC:	0.48	Parametric d_a:	0.06	  Accuracy:	 60% +- 7.5% (SE), P = 0.285442

USE_LASSO_ROC_D2 = roc_plot([pexp.D{4}; pexp.D{5}], use2==1, 'color', [.8 .2 .2]);  set(gca, 'FontSize', 18); hold on
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	1.62	Sens:	 90% CI(75%-100%)	Spec:	 50% CI(29%-71%)	PPV:	 63% CI(46%-80%)	Nonparametric AUC:	0.76	Parametric d_a:	0.91	  Accuracy:	 70% +- 7.0% (SE), P = 0.020718

USE_LASSO_ROC_DmF2 = roc_plot([pexp.DminF{4}; pexp.DminF{5}], use2==1, 'color', [.5 0 .76]);  set(gca, 'FontSize', 18);
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-1.12	Sens:	 95% CI(84%-100%)	Spec:	 68% CI(47%-86%)	PPV:	 74% CI(57%-89%)	Nonparametric AUC:	0.85	Parametric d_a:	1.33	  Accuracy:	 81% +- 5.9% (SE), P = 0.000079
set(gca, 'FontSize', 18, 'LineWidth', 2); box('off')

% cocaine users only **********
use3 = [ones(21,1); -ones(18,1)];
subplot(1,2,2)

USE_LASSO_ROC_F3 = roc_plot([pexp.F{2}; pexp.F{3}], use3==1, 'color', [0 .7 .85]);  set(gca, 'FontSize', 18); hold on
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	1.72	Sens:	 90% CI(77%-100%)	Spec:	 17% CI(0%-36%)	PPV:	 56% CI(40%-73%)	Nonparametric AUC:	0.31	Parametric d_a:	-0.61	  Accuracy:	 56% +- 7.9% (SE), P = 0.875879

USE_LASSO_ROC_D3 = roc_plot([pexp.D{2}; pexp.D{3}], use3==1, 'color', [.8 .2 .2]);  set(gca, 'FontSize', 18); hold on
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	2.26	Sens:	 86% CI(71%-100%)	Spec:	 72% CI(50%-92%)	PPV:	 78% CI(60%-95%)	Nonparametric AUC:	0.77	Parametric d_a:	0.77	  Accuracy:	 79% +- 6.5% (SE), P = 0.001621

USE_LASSO_ROC_DmF3 = roc_plot([pexp.DminF{2}; pexp.DminF{3}], use3==1, 'color', [.5 0 .76]);  set(gca, 'FontSize', 18);
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-0.59	Sens:	 90% CI(77%-100%)	Spec:	 78% CI(57%-95%)	PPV:	 83% CI(65%-96%)	Nonparametric AUC:	0.89	Parametric d_a:	1.59	  Accuracy:	 85% +- 5.8% (SE), P = 0.000105
set(gca, 'FontSize', 18, 'LineWidth', 2); box('off')


%%  ANOVA cv-NCS responses by group and condition

pexp2_vector = [pexp.DN{1}; pexp.DL{1}; pexp.FN{1}; pexp.FL{1}; ...
               pexp.DN{2}; pexp.DL{2}; pexp.FN{2}; pexp.FL{2}; ...
               pexp.DN{3}; pexp.DL{3}; pexp.FN{3}; pexp.FL{3}; ...
               pexp.DN{4}; pexp.DL{4}; pexp.FN{4}; pexp.FL{4}; ...
               pexp.DN{5}; pexp.DL{5}; pexp.FN{5}; pexp.FL{5}];
 
pexp2_vector_U = [pexp.DN{1}; pexp.DL{1}; pexp.FN{1}; pexp.FL{1}; ...
               pexp.DN{2}; pexp.DL{2}; pexp.FN{2}; pexp.FL{2}; ...
               pexp.DN{4}; pexp.DL{4}; pexp.FN{4}; pexp.FL{4}];
           
pexp2_vector_HC = [pexp.DN{3}; pexp.DL{3}; pexp.FN{3}; pexp.FL{3}; ...
               pexp.DN{5}; pexp.DL{5}; pexp.FN{5}; pexp.FL{5}];           

% [h, p, ci, stats] = ttest2(mean([pexp.DN{1} pexp.DL{1}; pexp.DN{2} pexp.DL{2}; pexp.DN{4} pexp.DL{4}], 2), mean([pexp.DN{3} pexp.DL{3}; pexp.DN{5} pexp.DL{5}], 2));
% [h, p, ci, stats] = ttest2(mean([pexp.FN{1} pexp.FL{1}; pexp.FN{2} pexp.FL{2}; pexp.FN{4} pexp.FL{4}], 2), mean([pexp.FN{3} pexp.FL{3}; pexp.FN{5} pexp.FL{5}], 2));
[h, p, ci, stats] = ttest2(mean([pexp.DN{1} pexp.FN{1}; pexp.DN{2} pexp.FN{2}; pexp.DN{4} pexp.FN{4}], 2), mean([pexp.DN{3} pexp.FN{3}; pexp.DN{5} pexp.FN{5}], 2));
[h, p, ci, stats] = ttest2(mean([pexp.DL{1} pexp.FL{1}; pexp.DL{2} pexp.FL{2}; pexp.DL{4} pexp.FL{4}], 2), mean([pexp.DL{3} pexp.FL{3}; pexp.DL{5} pexp.FL{5}], 2));
% [h, p, ci, stats] = ttest2(mean([pexp.DN{1} pexp.FN{1}; pexp.DN{2} pexp.FN{2}; pexp.DN{4} pexp.FN{4}], 2) - mean([pexp.DL{1} pexp.FL{1}; pexp.DL{2} pexp.FL{2}; pexp.DL{4} pexp.FL{4}], 2), mean([pexp.DN{3} pexp.FN{3}; pexp.DN{5} pexp.FN{5}], 2)- mean([pexp.DL{3} pexp.FL{3}; pexp.DL{5} pexp.FL{5}], 2));
% [h, p, ci, stats] = ttest(mean([pexp.DN{1} pexp.FN{1}; pexp.DN{2} pexp.FN{2}; pexp.DN{4} pexp.FN{4}], 2), mean([pexp.DL{1} pexp.FL{1}; pexp.DL{2} pexp.FL{2}; pexp.DL{4} pexp.FL{4}], 2));
% [h, p, ci, stats] = ttest(mean([pexp.DN{3} pexp.FN{3}; pexp.DN{5} pexp.FN{5}], 2), mean([pexp.DL{3} pexp.FL{3}; pexp.DL{5} pexp.FL{5}], 2));

xlsfile = 'NCS_weightmaps/CV_wmaps/Table_PatternExp_byGroupCondition_FINAL.xlsx';
[~, datname] = xlsread(xlsfile, 'All_byROCcondition', 'A2:A397');
datset = meancenter(xlsread(xlsfile, 'All_byROCcondition', 'F2:F397'));
subjects = meancenter(xlsread(xlsfile, 'All_byROCcondition', 'G2:G397'));
patients = xlsread(xlsfile, 'All_byROCcondition', 'H2:H397') ./2;
drug = xlsread(xlsfile, 'All_byROCcondition', 'I2:I397') ./2;
regul = xlsread(xlsfile, 'All_byROCcondition', 'J2:J397') ./2;
pexp2_val = pexp2_vector;

datpexp2 = table(datname, datset, subjects, patients, drug, regul, pexp2_val);

lme = fitlme(datpexp2, 'pexp2_val ~ 1 + patients + drug + regul + drug*regul + patients*drug + patients*regul + patients*regul*drug + (1 + drug + regul + drug*regul |subjects)', 'FitMethod', 'REML');
anova(lme)
cilme = coefCI(lme);

% Post-hoc tests: 
% compare drug users vs. non-users on Drug responses and Food Responses
[h, p, ci, stats] = ttest2([pexp.D{1}; pexp.D{2}; pexp.D{4}], [pexp.D{3}; pexp.D{5}]); % p = 2.8936e-05*
[h, p, ci, stats] = ttest2([pexp.F{1}; pexp.F{2}; pexp.F{4}], [pexp.F{3}; pexp.F{5}]); % p = 0.2597


%% add behavior, load and draw 2-dim errorbars, behavior on x-axis, pexp on y-axis (Fig. 3c)

load('BehavioralData/mean_ratings_behav_N99.mat')
dnames = {'meanratings_alc'; 'meanratings_coc'; 'meanratings_coc_hc'; 'meanratings_cig'; 'meanratings_cig_hc'};

for ds = 1:5
    dn{ds} = datpexp2.pexp2_val(datpexp2.datset==ds & datpexp2.drug==1 & datpexp2.regul==1);
    dl{ds} = datpexp2.pexp2_val(datpexp2.datset==ds & datpexp2.drug==1 & datpexp2.regul==-1);
    fn{ds} = datpexp2.pexp2_val(datpexp2.datset==ds & datpexp2.drug==-1 & datpexp2.regul==1);
    fl{ds} = datpexp2.pexp2_val(datpexp2.datset==ds & datpexp2.drug==-1 & datpexp2.regul==-1);    
end

%% cross-plots

cols = [1 .373 .475; .667 .047 .118; .463 .486 .965; .102 .098 .749];

figure
posi = [3 4 5 1 2]; 

for ds = 1:5
    
    beh_all_meanrate{ds} = eval(dnames{ds});
    
    subplot(1,5, posi(ds))

    % add subject-wise lines
    for s = 1:numel(pexp.DN{ds})
        scatter(beh_all_meanrate{ds}(s,:), [pexp.DN{ds}(s), pexp.DL{ds}(s), pexp.FN{ds}(s), pexp.FL{ds}(s)], 3, cols); refline; hold on
    end
    
    % DN condition
%     scatter(beh_all_meanrate{ds}(:,1), dn{ds}, 4, [.9 0 0], '*'); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,1)), mean(dn{ds}), ste(beh_all_meanrate{ds}(:,1)), 'horizontal', 'Color', cols(1,:), 'LineWidth', 3); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,1)), mean(dn{ds}), ste(dn{ds}), 'vertical', 'Color',  cols(1,:), 'LineWidth', 3); hold on
    
    % DL condition
%     scatter(beh_all_meanrate{ds}(:,2), dl{ds}, 4, [.7 0 .1], '*'); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,2)), mean(dl{ds}), ste(beh_all_meanrate{ds}(:,2)), 'horizontal', 'Color', cols(2,:), 'LineWidth', 3); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,2)), mean(dl{ds}), ste(dl{ds}), 'vertical', 'Color',   cols(2,:), 'LineWidth', 3); hold on

    % FN condition
%     scatter(beh_all_meanrate{ds}(:,3), fn{ds}, 4, [.5 0 .3]); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,3)), mean(fn{ds}), ste(beh_all_meanrate{ds}(:,3)), 'horizontal', 'Color', cols(3,:), 'LineWidth', 3); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,3)), mean(fn{ds}), ste(fn{ds}), 'vertical', 'Color',  cols(3,:), 'LineWidth', 3); hold on

    % FL condition
%     scatter(beh_all_meanrate{ds}(:,4), fl{ds}, 4, [.3 0 .5]); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,4)), mean(fl{ds}), ste(beh_all_meanrate{ds}(:,4)), 'horizontal', 'Color', cols(4,:), 'LineWidth', 3); hold on
    errorbar(mean(beh_all_meanrate{ds}(:,4)), mean(fl{ds}), ste(fl{ds}), 'vertical', 'Color',  cols(4,:), 'LineWidth', 3); hold on

%     title(stname{ds})
    set(gca, 'FontSize', 14, 'LineWidth', 2, 'box', 'off', 'XTick', 1:5); xlabel('Craving rating')
    if ds==4
        ylabel('Pattern expression'); 
    end
    xlim([.5 5.5]); ylim([-1 7])
    
end

%% ANOVA behavior (ratings)

beh_vector = [reshape(beh_all_meanrate{1}, [numel(beh_all_meanrate{1}),1]); ...
              reshape(beh_all_meanrate{2}, [numel(beh_all_meanrate{2}),1]); ...
              reshape(beh_all_meanrate{3}, [numel(beh_all_meanrate{3}),1]); ...
              reshape(beh_all_meanrate{4}, [numel(beh_all_meanrate{4}),1]); ...
              reshape(beh_all_meanrate{5}, [numel(beh_all_meanrate{5}),1])];
 
beh_vector_U = [reshape(beh_all_meanrate{1}, [numel(beh_all_meanrate{1}),1]); ...
               reshape(beh_all_meanrate{2}, [numel(beh_all_meanrate{2}),1]); ...
               reshape(beh_all_meanrate{4}, [numel(beh_all_meanrate{4}),1])];
           
beh_vector_HC = [reshape(beh_all_meanrate{3}, [numel(beh_all_meanrate{3}),1]); ...
                 reshape(beh_all_meanrate{5}, [numel(beh_all_meanrate{5}),1])];           

datpexp22 = table(datname, datset, subjects, patients, drug, regul, pexp2_val, beh_vector);

lme2 = fitlme(datpexp22, 'beh_vector ~ 1 + patients + drug + regul + drug*regul + patients*drug + patients*regul + patients*regul*drug + (1 + drug + regul + drug*regul |subjects)', 'FitMethod', 'REML');
anova(lme2)
cilme2 = coefCI(lme2);

% planned comparisons for behavior effects
[h, p, ci, stats] = ttest(mean([pexp.DN{1} pexp.FN{1}; pexp.DN{2} pexp.FN{2}; pexp.DN{4} pexp.FN{4}], 2), mean([pexp.DL{1} pexp.FL{1}; pexp.DL{2} pexp.FL{2}; pexp.DL{4} pexp.FL{4}], 2));

% in controls DN vs DL:
[h, p, ci, stats] = ttest([beh_all_meanrate{3}(:,1); beh_all_meanrate{5}(:,1)], [beh_all_meanrate{3}(:,2); beh_all_meanrate{5}(:,2)]);
% in controls FN vs FL:
[h, p, ci, stats] = ttest([beh_all_meanrate{3}(:,3); beh_all_meanrate{5}(:,3)], [beh_all_meanrate{3}(:,4); beh_all_meanrate{5}(:,4)]);
% regulation effect difference (D-F) in controls: p = 1.4117e-04*
[h, p, ci, stats] = ttest([beh_all_meanrate{3}(:,1)-beh_all_meanrate{3}(:,2); beh_all_meanrate{5}(:,1)-beh_all_meanrate{5}(:,2)], [beh_all_meanrate{3}(:,3)-beh_all_meanrate{3}(:,4); beh_all_meanrate{5}(:,3)-beh_all_meanrate{5}(:,4)]);

% in users DN vs DL:
[h, p, ci, stats] = ttest([beh_all_meanrate{1}(:,1); beh_all_meanrate{2}(:,1); beh_all_meanrate{4}(:,1)], [beh_all_meanrate{1}(:,2); beh_all_meanrate{2}(:,2); beh_all_meanrate{4}(:,2)]);
% in users FN vs FL:
[h, p, ci, stats] = ttest([beh_all_meanrate{1}(:,3); beh_all_meanrate{2}(:,3); beh_all_meanrate{4}(:,3)], [beh_all_meanrate{1}(:,4); beh_all_meanrate{2}(:,4); beh_all_meanrate{4}(:,4)]);
% regulation effect difference (D-F) in users: p =  0.0365* (in the other direction than non-users)
[h, p, ci, stats] = ttest([beh_all_meanrate{1}(:,1)-beh_all_meanrate{1}(:,2); beh_all_meanrate{2}(:,1)-beh_all_meanrate{2}(:,2); beh_all_meanrate{4}(:,1)-beh_all_meanrate{4}(:,2)], [beh_all_meanrate{1}(:,3)-beh_all_meanrate{1}(:,4); beh_all_meanrate{2}(:,3)-beh_all_meanrate{2}(:,4); beh_all_meanrate{4}(:,3)-beh_all_meanrate{4}(:,4)]);

% Post-hoc tests: 
% compare drug users vs. non-users on Drug responses and Food Responses
[h, p, ci, stats] = ttest2(mean([beh_all_meanrate{1}(:,1:2); beh_all_meanrate{2}(:,1:2); beh_all_meanrate{4}(:,1:2)],2), mean([beh_all_meanrate{3}(:,1:2); beh_all_meanrate{5}(:,1:2)],2)); % p = 5.9436e-26
[h, p, ci, stats] = ttest2(mean([beh_all_meanrate{1}(:,3:4); beh_all_meanrate{2}(:,3:4); beh_all_meanrate{4}(:,3:4)],2), mean([beh_all_meanrate{3}(:,3:4); beh_all_meanrate{5}(:,3:4)],2)); % p = 0.5038


%% only users

[~, datname] = xlsread(xlsfile, 'Users', 'A2:A245');
datset = xlsread(xlsfile, 'Users', 'F2:F245');
subjects = xlsread(xlsfile, 'Users', 'G2:G245');
drug = xlsread(xlsfile, 'Users', 'I2:I245');
regul = xlsread(xlsfile, 'Users', 'J2:J245');
pexp2_val_U = pexp2_vector_U;

datpexp2_U = table(datname, datset, subjects, drug, regul, pexp2_val_U);
glme_U = fitglme(datpexp2_U, 'pexp2_val_U ~ 1 + drug + regul + drug*regul + (1|subjects) + (1|datset)');


%% only NON-users

[~, datname] = xlsread(xlsfile, 'NonUsers', 'A2:A161');
datset = xlsread(xlsfile, 'NonUsers', 'F2:F161');
subjects = xlsread(xlsfile, 'NonUsers', 'G2:G161');
drug = xlsread(xlsfile, 'NonUsers', 'I2:I161');
regul = xlsread(xlsfile, 'NonUsers', 'J2:J161');
pexp2_val_C = pexp2_vector_HC;

datpexp2_C = table(datname, datset, subjects, drug, regul, pexp2_val_C);
glme_C = fitglme(datpexp2_C, 'pexp2_val_C ~ 1 + drug + regul + drug*regul + (1|subjects) + (1|datset)');
