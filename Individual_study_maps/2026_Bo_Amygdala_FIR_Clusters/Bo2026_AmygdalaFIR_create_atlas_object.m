% Bo2026_AmygdalaFIR_create_atlas_object.m
%
% How the data-driven amygdala clusters in
%
%   Bo, K., Lindquist, M. A., Gianaros, P. J., & Wager, T. D. (2026).
%   Spatiotemporal dissociation of human amygdala response to negative
%   affect. bioRxiv. doi:10.64898/2026.08.18.745530
%
% were made, and the script that writes every file in this folder.
%
% PIPELINE (paper Methods, "Data-Driven Clustering of Amygdala Voxels")
%   1. Voxelwise FIR estimates (15 TRs, TR = 2 s) for Regulate Negative,
%      Look Negative and Look Neutral, N = 358 (AHAB n = 182, PIP n = 176),
%      in the 958 voxels of the union of the Amunts (CANlab 2018) and
%      Nacewicz amygdala masks.
%   2. Features: group-mean time course per voxel, concatenated across
%      conditions -> 958 voxels x 45 values.
%   3. Number of clusters: silhouette pseudo-Z permutation test
%      (clusterdata_permtest, k = 2..7, 1,000 permutations) -> k = 4.
%   4. k-means, k = 4, squared Euclidean distance.
%   5. Clusters mapped back to the brain and compared with the Nacewicz and
%      Amunts subregions: adjusted Rand index and Dice, each against
%      10,000 label permutations, Dice FDR-corrected (Benjamini-Hochberg,
%      q < .05) across each cluster-by-region matrix.
%
%   k-means starts from random centroids, so a new run can differ from
%   the paper's by a few voxels on cluster borders. The assignment used in
%   the paper is stored in Bo2026_AmygdalaFIR_cluster_idx.csv and is used
%   by default (USE_PAPER_ASSIGNMENT = true). With it, this script
%   reproduces every cluster size, ARI and Dice value in the paper.
%
% OUTPUTS (this folder)
%   Bo2026_AmygdalaFIR_4clusters_atlas_object.mat      CANlab atlas object
%   Bo2026_AmygdalaFIR_4clusters_atlas_regions.nii.gz  label image (1-4)
%   Bo2026_AmygdalaFIR_4clusters_atlas_regions.mat     region objects
%   Bo2026_AmygdalaFIR_C<k>_<name>.nii.gz              binary mask per cluster
%   Bo2026_AmygdalaFIR_cluster_timecourses.csv         Figure 5B data (mean, SEM)
%
% Steps 1-5 need the voxelwise FIR source data (to be released with the
% paper). Without it, set REBUILD_FROM_SOURCE = false to rebuild only the
% atlas object from the label image in this folder.
%
% Requirements: CanlabCore, SPM12, Statistics and Machine Learning
% Toolbox. kmeans_Matlab.m (a renamed copy of MATLAB's kmeans.m; CanlabCore
% shadows kmeans with an @kmeans class) is needed only for a new k-means run.

%% ------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------

REBUILD_FROM_SOURCE  = true;    % false = only build the atlas object from the label image
USE_PAPER_ASSIGNMENT = true;    % false = run k-means again (may differ by a few voxels)
RUN_K_SELECTION      = false;   % true = rerun the k = 2..7 permutation test (slow)
nPerm                = 10000;   % permutations for ARI / Dice significance (0 = skip)

data_root      = 'C:\Users\diban\Dartmouth College Dropbox\Ke Bo\2023_Ke_Bo_Amygdala_Timeseries';
voxFIR_file    = fullfile(data_root, 'All_FIR_Time_Sub_SingleVoxel_Shackman_Combined.mat');
combAtlas_file = fullfile(data_root, 'Ke_TimeSeries_Canlab2024_Amygdalatemplate', ...
                          'Amygdala_Timeseries_ShackmanTemplate', 'CombinedAtlas.mat');
nacewicz_file  = fullfile(data_root, 'ShackmanAmygdalObject.mat');

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end

atlas_name = 'Bo2026_AmygdalaFIR_4clusters';
label_file = fullfile(this_dir, [atlas_name '_atlas_regions.nii.gz']);
idx_file   = fullfile(this_dir, 'Bo2026_AmygdalaFIR_cluster_idx.csv');
tc_file    = fullfile(this_dir, 'Bo2026_AmygdalaFIR_cluster_timecourses.csv');

best_k = 4;
labels = {'AmyFIR_C1_early_transient', 'AmyFIR_C2_weak_lateral', ...
          'AmyFIR_C3_early_sustained', 'AmyFIR_C4_delayed_centromedial'};
labels_2 = {'BLBM_CoAHA_LB', 'LA_LB', 'CoAHA_SF', 'CeMe_CM'};   % best anatomical matches (Nacewicz / Amunts)
labels_3 = repmat({'Amygdala'}, 1, best_k);
label_descriptions = { ...
    'Cluster 1 (255 voxels): early, transient response (Look Negative peak TR 5, no late peak); overlaps Nacewicz BLBM and CoAHA, Amunts LB'; ...
    'Cluster 2 (474 voxels): weakest evoked response, modest negative vs neutral difference; overlaps Nacewicz LA, Amunts LB'; ...
    'Cluster 3 (140 voxels): robust early response (peak TR 5) that stays elevated; overlaps Nacewicz CoAHA (and CeMe), Amunts SF'; ...
    'Cluster 4 (89 voxels): most delayed response (peak TR 8, into the post-stimulus rating period); overlaps Nacewicz CeMe, Amunts CM'};
space_description = 'MNI152 2 mm (SPM12 normalization; 91 x 109 x 91 grid)';
references = ['Bo K, Lindquist MA, Gianaros PJ, Wager TD (2026). Spatiotemporal ' ...
              'dissociation of human amygdala response to negative affect. bioRxiv. ' ...
              'doi:10.64898/2026.08.18.745530'];
cluster_colors = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.89 0.10 0.11; 1.00 0.50 0.00];  % Figure 5

% Values reported in the paper (Results and Figure 5), used as checks
paper_sizes = [255 474 140 89];
paper_ari   = [0.276 0.161];                       % Nacewicz, Amunts
paper_dice  = { ...                                % {cluster, atlas (1 = Nacewicz, 2 = Amunts), region, Dice}
    1, 1, 'BLBM', 0.520;  1, 1, 'CoAHA', 0.505;  1, 2, 'LB', 0.465; ...
    2, 1, 'LA',   0.687;  2, 2, 'LB',    0.662; ...
    3, 1, 'CoAHA',0.517;  3, 1, 'CeMe',  0.289;  3, 2, 'SF', 0.389; ...
    4, 1, 'CeMe', 0.139;  4, 2, 'CM',    0.577};

if REBUILD_FROM_SOURCE

    %% --------------------------------------------------------------------
    % 1. LOAD SOURCE DATA
    % ---------------------------------------------------------------------

    load(voxFIR_file, 'H_All_Reg', 'H_All_Neg', 'H_All_Neu');          % 1 x 358 cells, voxels x TRs
    load(combAtlas_file, 'Combined_atlas', 'Combined_atlas_index');    % amygdala mask + voxel index
    load(nacewicz_file, 'atlas_separate_Value_All');                   % Nacewicz subregions

    nVox  = size(H_All_Neg{1}, 1);
    nTRs  = size(H_All_Neg{1}, 2);
    valid = find(cellfun(@(x) size(x, 1), H_All_Neg) == nVox);
    nSub  = numel(valid);
    fprintf('%d voxels, %d TRs, %d participants\n', nVox, nTRs, nSub);

    %% --------------------------------------------------------------------
    % 2. FEATURES: group-mean FIR, [Reg TR1..15 | Neg TR1..15 | Neu TR1..15]
    % ---------------------------------------------------------------------

    conds = {H_All_Reg, H_All_Neg, H_All_Neu};
    averageMatrix = zeros(nVox, nTRs, 3);
    for c = 1:3
        averageMatrix(:, :, c) = mean(cat(3, conds{c}{valid}), 3);
    end
    features = reshape(averageMatrix, nVox, nTRs * 3);
    features(isnan(features)) = mean(features(~isnan(features)));

    %% --------------------------------------------------------------------
    % 3. NUMBER OF CLUSTERS (paper result: k = 4)
    % ---------------------------------------------------------------------

    if RUN_K_SELECTION
        kstats = clusterdata_permtest(features, 'k', 2:7, 'nperm', 1000);
        fprintf('clusterdata_permtest best k = %d\n', kstats.best_k);
    end

    %% --------------------------------------------------------------------
    % 4. K-MEANS CLUSTERING
    % ---------------------------------------------------------------------

    % Anatomical reference maps, resampled to the cluster space as in the paper
    cluster_map = fmri_data(Combined_atlas);
    cluster_map.dat = zeros(size(cluster_map.dat, 1), 1);
    cluster_map.metadata_table = table();   % otherwise write() adds .metadata.csv files

    nac = resample_space(atlas_separate_Value_All, cluster_map);
    amu = resample_space(select_atlas_subset(load_atlas('canlab2018_2mm'), ...
                         {'Amygdala_CM_', 'Amygdala_SF_', 'Amygdala_LB_'}), cluster_map);
    anat_maps   = {nac.dat, amu.dat};
    anat_vals   = {[1 3 4 5], [1 2 3]};
    anat_labels = {{'BLBM', 'CeMe', 'CoAHA', 'LA'}, {'CM', 'SF', 'LB'}};
    atlas_names = {'Nacewicz', 'Amunts'};

    if USE_PAPER_ASSIGNMENT
        T = readtable(idx_file);
        assert(isequal(T.fir_row(:), (1:nVox)'), 'Row order in %s does not match the FIR data.', idx_file);
        idx = T.cluster(:);
        fprintf('Using the cluster assignment from the paper (%s)\n', idx_file);
    else
        idx_raw = kmeans_Matlab(features, best_k, 'Distance', 'sqeuclidean');
        % k-means numbers clusters arbitrarily: renumber to the paper's order
        % by anatomy (C4 = best Amunts CM, C2 = best Nacewicz LA,
        % C3 = best Amunts SF, C1 = the rest)
        tmp = cluster_map.dat;
        tmp(Combined_atlas_index) = idx_raw;
        bN = tmp > 0 & anat_maps{1} > 0;
        bA = tmp > 0 & anat_maps{2} > 0;
        Dn = local_dice(tmp(bN), anat_maps{1}(bN), anat_vals{1}, best_k);
        Da = local_dice(tmp(bA), anat_maps{2}(bA), anat_vals{2}, best_k);
        remaining = 1:best_k;
        new_of_old = zeros(1, best_k);
        [~, i] = max(Da(remaining, 1)); new_of_old(remaining(i)) = 4; remaining(i) = [];
        [~, i] = max(Dn(remaining, 4)); new_of_old(remaining(i)) = 2; remaining(i) = [];
        [~, i] = max(Da(remaining, 2)); new_of_old(remaining(i)) = 3; remaining(i) = [];
        new_of_old(remaining) = 1;
        idx = new_of_old(idx_raw)';
    end

    %% --------------------------------------------------------------------
    % 5. MAP TO BRAIN AND COMPARE WITH ANATOMY
    % ---------------------------------------------------------------------

    cluster_map.dat(Combined_atlas_index) = idx;

    sizes = arrayfun(@(k) sum(idx == k), 1:best_k);
    fprintf('\nCluster sizes: %s (paper: %s)\n', mat2str(sizes), mat2str(paper_sizes));

    rng(42);
    for a = 1:2
        both = (cluster_map.dat > 0) & (anat_maps{a} > 0);
        cl = cluster_map.dat(both);
        an = anat_maps{a}(both);

        ari_obs = local_ari(cl, an);
        D_obs   = local_dice(cl, an, anat_vals{a}, best_k);

        fprintf('\n%s: ARI = %.3f (paper %.3f)', atlas_names{a}, ari_obs, paper_ari(a));
        if nPerm > 0
            ari_null = zeros(nPerm, 1);
            D_ge     = zeros(size(D_obs));
            for p = 1:nPerm
                s = cl(randperm(numel(cl)));
                ari_null(p) = local_ari(s, an);
                D_ge = D_ge + (local_dice(s, an, anat_vals{a}, best_k) >= D_obs);
            end
            p_ari  = (sum(ari_null >= ari_obs) + 1) / (nPerm + 1);
            p_dice = (D_ge + 1) / (nPerm + 1);
            pID    = FDR(max(p_dice(:), eps), 0.05);
            if isempty(pID), pID = 0; end
            sig = p_dice <= pID & pID > 0;           % <= : pID is itself one of the p-values
            fprintf(', permutation p = %.4g', p_ari);
        else
            sig = false(size(D_obs));
        end
        fprintf('\n%10s', '');
        fprintf('%9s', anat_labels{a}{:});
        fprintf('\n');
        for k = 1:best_k
            fprintf('%10s', sprintf('C%d', k));
            for j = 1:numel(anat_vals{a})
                star = ' '; if sig(k, j), star = '*'; end
                fprintf('%8.3f%s', D_obs(k, j), star);
            end
            fprintf('\n');
        end

        % Check the Dice values reported in the paper
        for r = 1:size(paper_dice, 1)
            [k, aa, reg, val] = paper_dice{r, :};
            if aa ~= a, continue; end
            obs = D_obs(k, strcmp(anat_labels{a}, reg));
            if abs(obs - val) > 0.0015
                warning('C%d vs %s %s: Dice %.3f, paper %.3f', k, atlas_names{a}, reg, obs, val);
            end
        end
    end
    if ~isequal(sizes, paper_sizes)
        warning('Cluster sizes differ from the paper.');
    end

    %% --------------------------------------------------------------------
    % 6. WRITE IMAGES AND TABLES
    % ---------------------------------------------------------------------

    local_write_gz(cluster_map, label_file);
    for k = 1:best_k
        m = cluster_map;
        m.dat = double(cluster_map.dat == k);
        local_write_gz(m, fullfile(this_dir, sprintf('Bo2026_AmygdalaFIR_C%d_%s.nii.gz', ...
                       k, regexprep(labels{k}, '^AmyFIR_C\d_', ''))));
    end

    if ~USE_PAPER_ASSIGNMENT
        writetable(table((1:nVox)', idx, 'VariableNames', {'fir_row', 'cluster'}), ...
                   fullfile(this_dir, 'Bo2026_AmygdalaFIR_cluster_idx_new_run.csv'));
    end

    % Figure 5B data: per-participant cluster means, then group mean and SEM
    cond_names = {'Regulate_Negative', 'Look_Negative', 'Look_Neutral'};
    C = {}; L = {}; N = {}; B = []; TR = []; M = []; SE = [];
    for k = 1:best_k
        for c = 1:3
            tc = zeros(nSub, nTRs);
            for s = 1:nSub
                tc(s, :) = mean(conds{c}{valid(s)}(idx == k, :), 1);
            end
            C  = [C;  num2cell(repmat(k, nTRs, 1))];          %#ok<AGROW>
            L  = [L;  repmat(labels(k), nTRs, 1)];            %#ok<AGROW>
            N  = [N;  repmat(cond_names(c), nTRs, 1)];        %#ok<AGROW>
            B  = [B;  (1:nTRs)'];                             %#ok<AGROW>
            TR = [TR; (0:nTRs-1)' - 1];                       %#ok<AGROW> % TR 0 = cue onset
            M  = [M;  mean(tc, 1)'];                          %#ok<AGROW>
            SE = [SE; (std(tc, 0, 1) / sqrt(nSub))'];         %#ok<AGROW>
        end
    end
    writetable(table(cell2mat(C), L, N, B, TR, M, SE, repmat(nSub, numel(B), 1), ...
        'VariableNames', {'cluster', 'cluster_label', 'condition', 'fir_bin', ...
                          'tr_rel_cue', 'mean_beta', 'sem', 'n_subjects'}), tc_file);
end

%% ------------------------------------------------------------------------
% 7. ATLAS OBJECT (from the label image)
% -------------------------------------------------------------------------

% atlas() reads an unzipped copy if one exists, so remove any stale copy first
label_nii = regexprep(label_file, '\.gz$', '');
if exist(label_nii, 'file'), delete(label_nii); end

atlas_obj = atlas(label_file, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'labels_2', labels_2, ...
    'labels_3', labels_3, ...
    'label_descriptions', label_descriptions(:), ...
    'space_description', space_description, ...
    'references', references, ...
    'noverbose');

if exist(label_nii, 'file'), delete(label_nii); end   % tidy the copy atlas() unzipped

atlas_sizes = arrayfun(@(k) sum(atlas_obj.dat == k), 1:best_k);
assert(isequal(atlas_sizes, paper_sizes), 'Atlas cluster sizes %s differ from the paper.', mat2str(atlas_sizes));
save(fullfile(this_dir, [atlas_name '_atlas_object.mat']), 'atlas_obj');

r = atlas2region(atlas_obj);
R = struct('r', r, 'region_names', {labels}, 'cluster_colors', cluster_colors);
for k = 1:best_k
    R.(labels{k}) = r(k);
end
save(fullfile(this_dir, [atlas_name '_atlas_regions.mat']), '-struct', 'R');

fprintf('\nDone: atlas object, region file, label image and masks are in %s\n', this_dir);

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function D = local_dice(cl, an, vals, K)
% Dice between each cluster (1..K) and each anatomical label in vals.
D = zeros(K, numel(vals));
for k = 1:K
    a = (cl == k);
    for j = 1:numel(vals)
        b = (an == vals(j));
        if any(a) || any(b)
            D(k, j) = 2 * sum(a & b) / (sum(a) + sum(b));
        end
    end
end
end

function ari = local_ari(x, y)
% Adjusted Rand index (Hubert & Arabie 1985).
[~, ~, i1] = unique(x(:));
[~, ~, i2] = unique(y(:));
T = accumarray([i1 i2], 1);
n2 = @(v) v .* (v - 1) / 2;
s_ij = sum(n2(T(:)));
s_a  = sum(n2(sum(T, 2)));
s_b  = sum(n2(sum(T, 1)));
expected = s_a * s_b / n2(numel(i1));
mx = (s_a + s_b) / 2;
if mx == expected
    ari = 0;
else
    ari = (s_ij - expected) / (mx - expected);
end
end

function local_write_gz(obj, gz_name)
% Write an fmri_data object to <name>.nii.gz
nii_name = regexprep(gz_name, '\.gz$', '');
write(obj, 'fname', nii_name, 'overwrite');
gzip(nii_name);
delete(nii_name);
end
