input_var = ST_cleaned.cPDM;
input_var_labels = ST_cleaned.cpdm_labels;

mainanalysislabel = 'cPDM';
partitionlabel = sprintf(mainanalysislabel, ' local patterns');

figsavedir = fullfile(pwd, 'figures');

pain2019_remove_pain_unrelated_nuisance % run this to get remove_nuisance function needed below

%% Correlation matrices, removing nuisance

X = scale(input_var, 1);
X = remove_nuisance(X);

p1 = ST_cleaned.gray_white_csf;
p2 = X;
p3 = [ST_cleaned.rel_temp ST_cleaned.pain_rating];
X = [p1 p2 p3];
Xlabels = [ST_cleaned.graywhitecsf_labels input_var_labels {'Temp' 'Pain'}];
Xpartitions = [ones(size(p1, 2), 1); 2 * ones(size(p2, 2), 1); 3 * ones(size(p3, 2), 1)]; 
partitionlabels = {'Nuisance' partitionlabel 'Temp Pain'};

figtitle = sprintf('%s_local_pattern_ST_connectivity_matrix_Zwithin_GWCSFcompRem', mainanalysislabel);
create_figure(figtitle, 1, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'p_thr', .001, 'nofigure', ...
    'var_names', Xlabels, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

% OUT = plot_correlation_matrix(X, 'dospearman', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
%     'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);


title('Full pairwise Spearman Zwithin GWCSFcomp removed');
drawnow

subplot(1, 2, 2);

OUT = plot_correlation_matrix(X, 'dospearman', true, 'dopartial', true, 'doimage', true, 'docircles', false, 'p_thr', .001, 'dofigure', false, ...
    'var_names', Xlabels, 'dotext', false, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

title('Partial r Spearman Zwithin GWCSFcomp removed');
drawnow

plugin_save_figure;

