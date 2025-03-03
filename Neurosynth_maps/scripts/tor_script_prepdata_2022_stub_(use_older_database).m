features = importdata('current_data/features.txt', '\t');
features.labels = features.textdata(1, 2:end)';

features.studies = features.textdata(2:end, 1);

% note: in latest version study ID is 1st numeric column of features.data

%% This stuff did not work properly:

% FID = fopen('current_data/database.txt')
% %d = textscan(FID, '%s', 'delimiter', '\t');
% d = textscan(FID, [repmat('%s', 1, 14)],  1, 'Delimiter', '\t');
% 
% %d = d{1};
% FID = fclose(FID);
% 
% names = d';
% for i = 1:length(names), names{i} = names{i}{1}; end
% 
% %%
% FID = fopen('current_data/database.txt');
% d = textscan(FID, '%s', 'delimiter', '\t');
% d = d{1};
% FID = fclose(FID);
% 
% % names = d';
% % for i = 1:length(names), names{i} = names{i}{1}; end
% 
% n = 5000;
% d2 = reshape(d(1:23*n), [23 n])'

%% This worked

%dbname = 'current_data/Neurosynth_4_8_22.txt';
dbname = 'current_data/database.txt';

% read_database;  % didn't work - spaces?

dat = readtable('current_data/database.txt', 'delimiter', '\t');
DB = table2struct(dat, 'ToScalar', true);

% rename
DB.CoordSys = DB.space;
DB = rmfield(DB, 'space');

%%
% other steps in read_database
% -----------------------------------------------------------------------------
% * convert from talairach to MNI, if we can
% -----------------------------------------------------------------------------
x = DB.x;
y = DB.y;
z = DB.z;

fprintf(1,'Converting T88 coordinates to MNI (M. Brett transform): ')
numt = 0;
for i = 1:length(DB.CoordSys)
    if strcmp(DB.CoordSys{i},'TAL')  % adapted this to code in Neurosynth
        XYZt = tal2mni([x(i) y(i) z(i)]);
        x(i) = XYZt(1); y(i) = XYZt(2); z(i) = XYZt(3);
        numt = numt + 1;
    end
end
fprintf(1,'%3.0f Transformed\n',numt)

XYZall = [x y z];

DB.xyz = XYZall;

DB.x_notransform = x;
DB.y_notransform = y;
DB.z_notransform = z;

DB.x = x;
DB.y = y;
DB.z = z;

%%
u_sort = unique(DB.table_id);
% Contrast = DB.table_id

u = unique(DB.table_id);
u = unique(DB.table_id, 'stable');

DB.Contrast = DB.table_id;

DB.N = ones(size(DB.x));
DB.study = DB.doi;
%%

DB = Meta_Setup(DB, 4);

disp('studies')
length(unique(DB.id))

disp('contrasts')
length(unique(DB.Contrast))

disp('points')
length(DB.x)

% ans =
% 
%         9721
% 
% contrasts
% 
% ans =
% 
%        18316
% 
% points
% 
% ans =
% 
%       347911

% Creates SETUP.mat

% Creates MC_Setup in MC_Info with data matrix
Meta_Activation_FWE('setup', DB);

%% Sort term matrix to be same as MC_Setup matrix

dois = DB.ID(DB.pointind);

% ia is for DB, ib is for features
% resort features by ib
% DB is already in order of ascending contrasts
[c, ia, ib] = intersect(dois, features.studies, 'stable');

features.data = features.data(ib, :);
features.studies = features.studies(ib, :);

DB.terms = features.data;

MC_Setup.terms = features.data;
MC_Setup.termlabels = features.labels;

% 5107 is missing from terms
wh = find(diff(ia) > 1) + 1;

MC_Setup.unweighted_study_data(:, wh) = [];;
MC_Setup.n(wh) = [];
MC_Setup.wts(wh) = [];
MC_Setup.cl(wh) = [];

save Wolfi_MC_Setup MC_Setup

%% Get clusters (BG Rois) and activation by cluster

cl = region('cluster_mask_7_pfc.nii', 'unique_mask_values');

[studybyroi,studybyset] = Meta_cluster_tools('getdata',region2struct(cl),MC_Setup.unweighted_study_data,MC_Setup.volInfo);

MC_Setup.bg_rois = cl;

MC_Setup.bg_data = studybyroi;

%% Average term profiles by cluster

yesno = MC_Setup.terms > .001;  % activation above thresh or not

% 3% of studies use a term > .001 freq on average
sum(MC_Setup.terms(:) > .001) ./ numel(MC_Setup.terms)

nterms = size(yesno, 2);

for i = 1:nterms
    
    %wh = MC_Setup.bg_data(:, i) == 1;
    t = yesno(:, i);
    
    % likelihood ratio.  
    % test = cluster activation (test), "disease" = term 
    % P(cluster | term) / P(cluster | no term)
    % P(term | cl) / P(term | no cl) -- flip it this way? no?
    
    h = MC_Setup.bg_data(t == 1, :);  % hits
    f = MC_Setup.bg_data(t == 0, :);  % false alarms
    
    hr = sum(h) ./ size(h, 1);
    far = sum(f) ./ size(f, 1);
    
    % lr matrix, terms x clusters
    lr(i, :) = hr ./ far;

end

% Interesting terms have high or low mean LR and high variance of LR across
% clusters
% For final analysis, we may want a priori groups of terms, but this is
% good to explore...
figure; plot_correlation_samefig(var(lr'), mean(lr'))
xlabel('variance of LR'); ylabel('Average LR');
interestingterms = var(lr') > .3;
interestingterms(415) = 0; % remove "reward" because everybody loads high on it
% tor: we can maybe deal with this in the plot stage?

colors = scn_standard_colors(7);

create_figure('LR');
for i = 1:7
    plot(lr(interestingterms, i), 1:sum(interestingterms), '-', 'Color', colors{i}, 'LineWidth', 2);
end
plot_vertical_line(1)

corrcoef(lr(interestingterms, :))

set(gca, 'YTicklabel', MC_Setup.termlabels(interestingterms), 'YTick', 1:sum(interestingterms))

%% Viz clustersColumnNames

create_figure('brains');
for i = 1:7
    phan(i) = imageCluster('cluster', region2struct(cl(i)), 'color', colors{i});
end

ph2 = addbrain('hires');
set(ph2, 'FaceAlpha', .1); axis image; axis tight; set(gca, 'YLim', [-40 40])

view(-156, 28);
lightRestoreSingle
material dull; lighting gouraud


%% RELOAD

load('Wolfi_MC_Setup.mat')


%% Regressions

% For each cluster, predict cluster activation *(yes/no) with combo of term maps
% * we may need to select terms of psychological interest
% among a large set...

MC_Setup.bg_data;  % 5808 studies x 7 BG clusters, 1/0 for activation in/out of cluster

MC_Setup.terms; % 5808 studies x 525 terms, frequency per 1000 words??
MC_Setup.termlabels; % the terms

cverr={};
stats={};
optout={};

% used this for one at a time
% parfor i = 1:7 % , index cluster...
% 
%     y = 2 * (double(MC_Setup.bg_data(:, i) == 1) - .5);
%     %y = double(MC_Setup.bg_data(:, i) == 1);
%     
%     X = MC_Setup.terms(:,interestingterms);
%     
%     dat = fmri_data;
%     dat.Y = y;
%     dat.dat = X';
% 
%     %%
% 
%     %%
%     % predict_test_suite(dat, 'nfolds', 5);
% 
%     %Elastic net with first 10 components:
%     %[cverr, stats, optout] = predict(dat_masked, 'algorithm_name', 'cv_lassopcrmatlab', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', .5); stats.pred_outcome_r
% 
%     % First try : standard SVM prediction, 4-fold xval -- 
%     [cverr{i}, stats{i}, optout{i}] = predict(dat, 'algorithm_name', 'cv_svm', 'nfolds', 4, 'error_type', 'mcr', 'C', 2, 'balanced_ridge', 8, 'rbf', 1);
% 
% end

% multiclass
y = 2 * (double(MC_Setup.bg_data) - .5);
    
X = MC_Setup.terms(:,interestingterms);
    
dat = fmri_data;
dat.Y = y;
dat.dat = X';

[cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'nfolds', 4, 'error_type', 'mcr', 'C', 2, 'balanced_ridge', 8, 'rbf', 1, 'MultiClass');

    
%%    
    
%save('svm_results.mat', 'cverr','stats','optout', '-v7.3'); % standard
%with 42 terms
%save('svm_results_1.mat', 'cverr','stats','optout', '-v7.3'); % 23 terms
%save('svm_results_2.mat', 'cverr','stats','optout', '-v7.3'); % 23 terms
%multiclass
load('svm_results_1.mat')

X = zeros(length(stats{1}.weight_obj.dat),7);
for i = 1:7
    X(:,i) = stats{i}.weight_obj.dat;
end

% Wolfi, I don't really understand what this code I commented out is doing.
% It produces a different solution/plot every time it is run (seems to be a
% stochastic element). 
% mns = repmat(mean(X), size(X,1), 1);
% sds = repmat(sqrt(var(X)), size(X,1), 1);
% sums = repmat(sum(X), size(X,1), 1);
% 
% %X = (X - mns)./sds;
% X = X./sums;
% 
% % some attempt at sorting the terms
% K = 8;
% kobj = kmeans(['k=', num2str(K)]);
% d=data('training data', X);
% [r,a] =train(kobj, d);

termlabels = MC_Setup.termlabels(interestingterms);
% X1 = [];
% termlabels1 = [];
% for k = 1:K
%     X1 = [X1; X(r.X == k,:)];
%     termlabels1 = [termlabels1; termlabels(r.X == k)];  
% end

X1 = X;
Xorig = X1;
X1 = scale(X1, 1);
X1 = scale(X1', 1)';

tor_polar_plot({X1}, colors, {termlabels});

% Barplot of accuracy
n = length(stats);
for i = 1:n
myacc(1, i) = 1 - (sum(stats{i}.err) ./ length(stats{i}.err));
end

create_figure('accuracy');
for i = 1:n
hh = bar(i, myacc(i));
set(hh, 'FaceColor', colors{i});
end
hh = plot_horizontal_line(.5); 
set(hh, 'LineStyle', '--');


% % Problem: always classifies as "1" -- need some work on this...
% [m,dprime,corr,far,missclass] = confusion_matrix(stats.Y,stats.yfit);
% 
% 
% svmobj = svm({'C=1', 'optimizer="andre"', 'balanced_ridge=8'});
% clear data
% sdata = data(dat.dat', dat.Y);
% %[svmobj, res] = train(sdata.....
    

%
% PRED = Meta_SVM_from_mkda(DB, fieldname, names_cell, 'balanced_ridge', balridgeval, 'holdout_meth', 'study');
 
% nfolds = 10;
% rbf_sigma = 15; % in mm
% balanced_ridge = 8; % arbitrary; 0 is no penalty for unbalanced classes
% holdout_meth = 'kfold'; % or 'study' for leave one study out
% test_method = 'vote';  % or 'distance'; method for multi-class classifications
% 
% % set up empty svm object for training
% % need spider package!
% % svmobj = svm({kernel('rbf', rbf_sigma), 'C=1', 'optimizer="andre"'});
% % balanced ridge=1 only shifts things a little bit...
% svmobj = svm({kernel('rbf', rbf_sigma), 'C=1', 'optimizer="andre"', ['balanced_ridge=' sprintf('%d', balanced_ridge)]});

%% MCA

% Y is 5808 studies x 7 clusters, activation
% X is 5808 studies x 24 (k) terms, term usage strength, mentions / 1000 words

Y = ((dat.Y + 1) ./ 2);
X = dat.dat';
%X = MC_Setup.terms;

% remove marginals for overall regional activation (Y) and overall term mentions (X)
Y = Y ./ repmat(sum(Y), size(Y, 1), 1);
X = X ./ repmat(sum(X), size(X, 1), 1);

% without centering rows, m becomes a matrix containing the average term mentions per activated
% study

% Y = Y';
% Y = Y ./ repmat(sum(Y), size(Y, 1), 1);
% Y = Y';
% Y(isnan(Y)) = 0;
% 
% X = X';
% X = X ./ repmat(sum(X), size(X, 1), 1);
% X = X';
% X(isnan(X)) = 0;

m = (Y' * X)';  % terms x regions, association matrix

% remove region marginals
m = m ./ repmat(sum(m), size(m, 1), 1);

% m = m';
% m = m';

% re-sort terms
c = clusterdata(m, 'linkage', 'average', 'maxclust', 7);
[c, ia] = sort(c);

m = m(ia, :);
myterms = termlabels(ia);

%create_figure('Counts_polarplot');
tor_polar_plot({m}, colors, {myterms});

figure; imagesc(m)
set(gca, 'YTicklabel', myterms, 'YTick', 1:sum(interestingterms))
colorbar

% NNMF http://cogsys.imm.dtu.dk/toolbox/nmf/index.html
% (this seems to be highly non-deterministic)
K = 7;
maxiter = 1000;
speak = true;
alg = 'mm'; % choices: mm, prob, cjlin, als, alsobs

[W,m]=nmf(X,K,alg,maxiter,speak);
m= m';

c = clusterdata(m, 'linkage', 'average', 'maxclust', 7);
[c, ia] = sort(c);

m = m(ia, :);
myterms = termlabels(ia);

%create_figure('Counts_polarplot');
tor_polar_plot({m}, colors, {myterms});

figure; imagesc(m)
set(gca, 'YTicklabel', myterms, 'YTick', 1:sum(interestingterms))
colorbar

