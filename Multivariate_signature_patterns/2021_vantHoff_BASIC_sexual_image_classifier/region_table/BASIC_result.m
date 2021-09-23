%% Load data object
%DATA_OBJ = importdata(which('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/sophie/results/FINAL FOR PAPER PUBLISH/svm sex-pnn/forced choice/data_objects.mat'))
pos_stark = fmri_data(DATA_OBJ{1, 2});
neg_stark = fmri_data(DATA_OBJ{1, 3});
sex_stark = fmri_data(DATA_OBJ{1, 1});
neu_stark = fmri_data(DATA_OBJ{1, 4});

%% Load svm results

%svm_stats_results = importdata(which('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/sophie/results/FINAL FOR PAPER PUBLISH/svm sex-pnn/forced choice/svm_stats_results_contrasts_masked.mat'));
obj_sex_neu = svm_stats_results{1, 1}.weight_obj; 
obj_sex_pos = svm_stats_results{1, 2}.weight_obj;
obj_sex_neg = svm_stats_results{1, 3}.weight_obj;

%% Treshold images
% these weight maps are 0.05 unc by default
obj_sex_neu_fdr_05 = threshold(obj_sex_neu, 0.05, 'fdr');

obj_sex_pos_fdr_05 = threshold(obj_sex_pos, 0.05, 'fdr');

obj_sex_neg_fdr_05 = threshold(obj_sex_neg, 0.05, 'fdr');


%% Create image overlap three weight maps

% set colors
mycolors = colorcube_colors
color_sex_neu_stark = mycolors{3}
color_sex_pos_stark = mycolors{4};
color_sex_neg_stark = mycolors{54};

% unc 05
my_display_obj = canlab_results_fmridisplay(obj_sex_neu, 'compact2', 'color', color_sex_neu_stark,  'trans', 'transvalue', 0.6, 'nooutline');
my_display_obj = addblobs(my_display_obj, region(obj_sex_pos), 'color', color_sex_pos_stark, 'trans', 'transvalue', 0.6);
my_display_obj = addblobs(my_display_obj, region(obj_sex_neg), 'color', color_sex_neg_stark, 'trans', 'transvalue', 0.6, 'nooutline');
% legend({'Sex-Neu', 'Sex-Pos', 'Sex-Neg'}, 'Location','NorthEastOutside')
drawnow, snapnow

%% fdr 05
create_figure('Weight map overlap 3 contrasts')
% color_sex_neu = [  0.8889    0.9444    0.4000];
% color_sex_pos = [   0.3492    0.6746    0.4000];
% color_sex_neg = [   0    0.5000    0.6000];
my_display_obj = canlab_results_fmridisplay(obj_sex_neu_fdr_05, 'compact2', 'color', color_sex_neu_stark, 'nooutline', 'trans','transvalue', 0.6);
my_display_obj = addblobs(my_display_obj, region(obj_sex_pos_fdr_05), 'color', color_sex_pos_stark, 'nooutline', 'trans','transvalue', 0.6);
my_display_obj = addblobs(my_display_obj, region(obj_sex_neg_fdr_05), 'color', color_sex_neg_stark, 'nooutline', 'trans', 'transvalue', 0.6);
% legend({'Sex-Neu', 'Sex-Pos', 'Sex-Neg'}, 'Location','NorthEastOutside')
drawnow, snapnow
%% Create weight map of overlap unc 05

neu_pos = apply_mask(obj_sex_neu, obj_sex_pos);
BASIC = apply_mask(neu_pos, obj_sex_neg);

orthviews(BASIC)
drawnow, snapnow
create_figure('Overlap 3 SVM weight maps montage')
axis('off')
montage(BASIC)
drawnow, snapnow
%% table of results unc
r_BASIC = region(BASIC);    
keuken_atlas = load_atlas('keuken');
 
schaefer_atlas = load_atlas('schaefer400' );
shen_atlas = load_atlas('shen');
glasser_atlas = load_atlas('glasser')

%% Table of regions
r_BASIC = subdivide_by_atlas(r_BASIC);
 [rpos, rneg] = table(r_BASIC);       % add labels
    r_BASIC = [rpos rneg];               % re-concatenate labeled regions
drawnow, snapnow
%% Table of regions with Glasser atlas
      [r_BASIC, region_table, table_legend_text] = autolabel_regions_using_atlas(r_BASIC, glasser_atlas);
      
    [rpos, rneg] = table(r_BASIC);       % add labels
    r_BASIC = [rpos rneg];               % re-concatenate labeled regions
drawnow, snapnow

%% Write BASIC to nifti
fname = 'BASIC.nii';
write(BASIC, 'fname', fname, 'overwrite');
%% Visualization BASIC model

create_figure('BASIC on surface')
surface(BASIC);
drawnow, snapnow

create_figure('lateral surfaces');
surface_handles = surface(BASIC, 'foursurfaces');

render_on_surface(BASIC, surface_handles);
drawnow, snapnow

% full montage
my_display_obj = canlab_results_fmridisplay([], 'montagetype', 'full');
my_display_obj = addblobs(my_display_obj, region(BASIC));
title('BASIC')
drawnow, snapnow

create_figure('Montage BASIC unc 05')
axis off
my_display_obj = canlab_results_fmridisplay([], 'montagetype', 'compact');
my_display_obj = addblobs(my_display_obj, region(BASIC));
drawnow, snapnow

% 3d coronal slabs
create_figure('cutaways')
axis off
surface_handles = surface(BASIC, 'coronal_slabs');
drawnow, snapnow

create_figure('cutaways')
axis off
surface_handles = surface(BASIC, 'coronal_slabs_4');
drawnow, snapnow

%% create movie
% create_figure('lateral surfaces');
% surface_handles = surface(BASIC);
% 
% % Then, choose movie option
% % mov = movie_tools('rotate',270,30) % rotate
%   mov = movie_tools('batch360.1'); % rotate 360
% %  mov = movie_tools('zoom',.8,mov,1); % zooms out
% 
% %% Save video
%   
%    vid = VideoWriter('BASIC_unc05_maps', 'MPEG-4');
%   open(vid)
%   for i = 1:length(mov), writeVideo(vid, mov(i)); end
%   close(vid)
%% test SVM model Sex vs. Neutral on pos/neg/sex/neu affective images
% [dist_from_hyperplane, Y, svm_dist_pos_neg, svm_dist_pos_neg_matrix] = plugin_svm_contrasts_get_results_per_subject(DAT, svm_stats_results, DATA_OBJ)
    load('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/sophie/MATLAB/MASTER_FOLDER_SEXUALBIOMARKER/results/Y.mat');


pos_stark_pattern_expression=apply_mask(pos_stark,replace_empty(BASIC),  'pattern_expression', 'cosine_similarity');
neg_stark_pattern_expression=apply_mask(neg_stark,replace_empty(BASIC),  'pattern_expression', 'cosine_similarity');
sex_stark_pattern_expression=apply_mask(sex_stark,replace_empty(BASIC),  'pattern_expression', 'cosine_similarity');
neu_stark_pattern_expression=apply_mask(neu_stark,replace_empty(BASIC),  'pattern_expression', 'cosine_similarity');

%% All together in barplot
mycolors = seaborn_colors;
create_figure('Barplot Model Similarity to all conditions ')
barplot_columns([sex_stark_pattern_expression,pos_stark_pattern_expression, neg_stark_pattern_expression, neu_stark_pattern_expression],'dolines','names',{'Sexual','Positive', 'Negative','Neutral'},'colors',mycolors,'nofigure')
ylabel('Cosine Similarity')
drawnow, snapnow

%% ROC plot forced choice
create_figure('ROC plot')
ROC_sex_neu_forced_choice = roc_plot([sex_stark_pattern_expression; neu_stark_pattern_expression],[ones(97,1);zeros(97,1)],'twochoice','color',mycolors{1}, 'threshold_type', 'Optimal balanced error rate');
ROC_sex_neu_single_interval = roc_plot([sex_stark_pattern_expression; neu_stark_pattern_expression],[ones(97,1);zeros(97,1)],'color',mycolors{1}, 'threshold_type', 'Optimal balanced error rate');

ROC_sex_pos_forced_choice = roc_plot([sex_stark_pattern_expression; pos_stark_pattern_expression],[ones(97,1);zeros(97,1)],'twochoice','color',mycolors{2});
ROC_sex_pos_single_interval= roc_plot([sex_stark_pattern_expression; pos_stark_pattern_expression],[ones(97,1);zeros(97,1)],'color',mycolors{2},'threshold_type', 'Optimal balanced error rate');

ROC_sex_neg_forced_choice = roc_plot([sex_stark_pattern_expression; neg_stark_pattern_expression],[ones(97,1);zeros(97,1)],'twochoice','color',mycolors{3});
ROC_sex_neg_single_interval = roc_plot([sex_stark_pattern_expression; neg_stark_pattern_expression],[ones(97,1);zeros(97,1)],'color',mycolors{3}, 'threshold_type', 'Optimal balanced error rate');

%% Effect size

% Effect size, cross-validated, paired samples
dfun_paired = @(x, Y) mean(x(Y > 0) - x(Y < 0)) ./ std(x(Y > 0) - x(Y < 0))

d_paired_sex_neu = dfun_paired(svm_stats_results{1, 1}.dist_from_hyperplane_xval, Y{1,1});
    fprintf('Effect size sex-neu, cross-val: Forced choice: d = %3.2f\n\n', d_paired_sex_neu);

    d_paired_sex_pos = dfun_paired(svm_stats_results{1, 2}.dist_from_hyperplane_xval, Y{1,2});
    fprintf('Effect size sex-pos, cross-val: Forced choice: d = %3.2f\n\n', d_paired_sex_pos);
    
    d_paired_sex_neg = dfun_paired(svm_stats_results{1, 3}.dist_from_hyperplane_xval, Y{1,3});
    fprintf('Effect size sex-neg, cross-val: Forced choice: d = %3.2f\n\n', d_paired_sex_neg);

    

%% Compare to Phil's independent dataset

% Load images
    
    sex_kragel=fmri_data('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/neuroimaging data/IAPS_firstlevel_contrasts_sexual_positive_negative/IAPS_n18_sexual.nii');
    pos_kragel=fmri_data('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/neuroimaging data/IAPS_firstlevel_contrasts_sexual_positive_negative/IAPS_n18_positive.nii');
    neg_kragel=fmri_data('/Users/sophievanthof/Dropbox (Dartmouth College)/Cooperation Valencespace/neuroimaging data/IAPS_firstlevel_contrasts_sexual_positive_negative/IAPS_n18_negative.nii');


%% test SVM model Sex vs. Neutral on sex/pos/neg IAPS data
% BASIC.p=[];
% BASIC.ste=[];
% BASIC.sig=[];
sex_kragel_pattern_expression=apply_mask(sex_kragel,replace_empty(BASIC),'pattern_expression', 'cosine_similarity'); 
pos_kragel_pattern_expression=apply_mask(pos_kragel,replace_empty(BASIC),'pattern_expression', 'cosine_similarity');
neg_kragel_pattern_expression=apply_mask(neg_kragel,replace_empty(BASIC),'pattern_expression', 'cosine_similarity');

% cols=[0 114 198;217 84 26; 237 177 32; 114 14 138]./255;
% for c=1:length(cols); COLS{c}=cols(c,:); end; 
% COLS=COLS([4 1 2 3]);

mycolors = seaborn_colors;
cols=[0.5921    0.6418    0.1935; 0.3127    0.6929    0.1924; 0.2031    0.6881    0.5178];

for c=1:length(cols); COLS{c}=cols(c,:); end; 
% kragel_colors=COLS([1 2 3]);

create_figure('Barplot Model2 Similarity to Independent Dataset ')
barplot_columns([sex_kragel_pattern_expression,pos_kragel_pattern_expression,neg_kragel_pattern_expression],'dolines','names',{'Sexual','Positive','Negative'},'colors',COLS,'nofigure')
ylabel('Cosine similarity')
ymax( 0.5)

%% ROC plot
create_figure('ROC plot independent dataset')

ROC_sex_pos_kragel_forced_choice = roc_plot([sex_kragel_pattern_expression; pos_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'twochoice','color',mycolors{6});
ROC_sex_pos_kragel_single_interval = roc_plot([sex_kragel_pattern_expression; pos_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'color',mycolors{2}, 'threshold_type', 'Optimal balanced error rate');

ROC_sex_neg_kragel_forced_choice = roc_plot([sex_kragel_pattern_expression;neg_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'twochoice','color',mycolors{3});
ROC_sex_neg_kragel_single_interval = roc_plot([sex_kragel_pattern_expression;neg_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'color',mycolors{4}, 'threshold_type', 'Optimal balanced error rate');

plot([0 1],[0 1],'k','linestyle','-.')
set(gca,'Fontsize',11)
hl=findobj(gca,'type','line');
legend(hl([2 1]),{'Sex vs. Pos' 'Sex vs. Neg'})

%% ROC plot for sex-pos both dependent and independent dataset
create_figure('ROC plot sex-pos')
seaborn = seaborn_colors
scn = scn_standard_colors
colorcube = colorcube_colors
color_sex_pos_stark_FC = colorcube{4}
color_sex_pos_stark_SI = seaborn {20}

color_sex_pos_kragel_FC = scn{18}
color_sex_pos_kragel_SI = seaborn{12}

ROC_sex_pos_forced_choice = roc_plot([sex_stark_pattern_expression; pos_stark_pattern_expression],[ones(97,1);zeros(97,1)],'twochoice','color',color_sex_pos_stark_FC);
ROC_sex_pos_single_interval= roc_plot([sex_stark_pattern_expression; pos_stark_pattern_expression],[ones(97,1);zeros(97,1)],'color',color_sex_pos_stark_SI,'threshold_type', 'Optimal balanced error rate');
ROC_sex_pos_kragel_forced_choice = roc_plot([sex_kragel_pattern_expression; pos_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'twochoice','color',color_sex_pos_kragel_FC);
ROC_sex_pos_kragel_single_interval = roc_plot([sex_kragel_pattern_expression; pos_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'color',color_sex_pos_kragel_SI, 'threshold_type', 'Optimal balanced error rate');


plot([0 1],[0 1],'k','linestyle','-.')
set(gca,'Fontsize',22)

%% ROC plot for sex-neg both dependent and independent dataset
create_figure('ROC plot sex-neg')


seaborn = seaborn_colors
scn = scn_standard_colors
colorcube = colorcube_colors

color_sex_neg_stark_FC = colorcube{54}
color_sex_neg_stark_SI = colorcube{57}

color_sex_neg_kragel_FC = colorcube{23}
color_sex_neg_kragel_SI = seaborn{16}


ROC_sex_neg_forced_choice = roc_plot([sex_stark_pattern_expression; neg_stark_pattern_expression],[ones(97,1);zeros(97,1)],'twochoice','color',color_sex_neg_stark_FC);

ROC_sex_neg_single_interval = roc_plot([sex_stark_pattern_expression; neg_stark_pattern_expression],[ones(97,1);zeros(97,1)],'color',color_sex_neg_stark_SI, 'threshold_type', 'Optimal balanced error rate');
ROC_sex_neg_kragel_forced_choice = roc_plot([sex_kragel_pattern_expression;neg_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'twochoice','color',color_sex_neg_kragel_FC);
ROC_sex_neg_kragel_single_interval = roc_plot([sex_kragel_pattern_expression;neg_kragel_pattern_expression],[ones(18,1);zeros(18,1)],'color',color_sex_neg_kragel_SI, 'threshold_type', 'Optimal balanced error rate');

plot([0 1],[0 1],'k','linestyle','-.')
set(gca,'Fontsize',22)
