function visualize_contents()
% visualize_contents  Render the Bo 2026 amygdala FIR temporal clusters
% into png_images/. See contents_description.md.
%
% Writes:
%   Bo2026_AmygdalaFIR_4clusters_montage.png      standard CANlab atlas montage
%   Bo2026_AmygdalaFIR_4clusters_isosurface.png   3-D isosurface
%   Bo2026_AmygdalaFIR_4clusters_fig5colors.png   axial/coronal/sagittal slices, paper colors
%   Bo2026_AmygdalaFIR_cluster_timecourses.png    group-mean FIR time course per cluster

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

atlas_file = fullfile(this_dir, 'Bo2026_AmygdalaFIR_4clusters_atlas_object.mat');

% Figure 5 colors: C1 blue, C2 green, C3 red, C4 orange
cluster_colors = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.89 0.10 0.11; 1.00 0.50 0.00];

% 1. Standard atlas montage + isosurface (shared CANlab helper)
canlab_render_atlas({'Bo2026_AmygdalaFIR_4clusters', atlas_file}, 'savedir', savedir);

% 2. Slices in the paper's colors (as in Figure 5A)
S = load(atlas_file, 'atlas_obj');
r = atlas2region(S.atlas_obj);

fig = figure('Name', 'Bo2026_AmygdalaFIR_4clusters_fig5colors', 'Color', 'w');
o = fmridisplay;
o = montage(o, 'axial',    'wh_slice', [-20 -16 -12], 'spacing', 6);
o = montage(o, 'coronal',  'wh_slice', [-6 -2 2],     'spacing', 6);
o = montage(o, 'saggital', 'wh_slice', [18 22 26],    'spacing', 6);
for k = 1:numel(r)
    o = addblobs(o, r(k), 'color', cluster_colors(k, :));
end
drawnow;
saveas(fig, fullfile(savedir, 'Bo2026_AmygdalaFIR_4clusters_fig5colors.png'));

% 3. Cluster time courses (Figure 5B) from the CSV
T = readtable(fullfile(this_dir, 'Bo2026_AmygdalaFIR_cluster_timecourses.csv'));
conds      = {'Look_Neutral', 'Look_Negative', 'Regulate_Negative'};
cond_names = {'Look neutral', 'Look negative', 'Regulate negative'};
cond_cols  = [0.30 0.65 0.30; 0.20 0.35 0.75; 0.75 0.20 0.20];

fig = figure('Name', 'Bo2026_AmygdalaFIR_cluster_timecourses', 'Color', 'w', ...
             'Position', [100 100 1500 380]);
for k = 1:4
    subplot(1, 4, k); hold on;
    h = gobjects(1, 3);
    for c = 1:3
        w = T.cluster == k & strcmp(T.condition, conds{c});
        h(c) = errorbar(T.tr_rel_cue(w), T.mean_beta(w), T.sem(w), ...
                        'Color', cond_cols(c, :), 'LineWidth', 1.8, 'CapSize', 3);
    end
    yline(0, '--k', 'Alpha', 0.3);
    ax = gca;
    ax.XColor = cluster_colors(k, :); ax.YColor = cluster_colors(k, :); ax.LineWidth = 2.5;
    title(sprintf('Cluster %d', k), 'Color', cluster_colors(k, :), 'FontSize', 13);
    xlabel('Time (TR; 0 = cue onset)'); ylabel('BOLD response (a.u.)');
    xlim([-1 13]); set(gca, 'XTick', -1:2:13, 'FontSize', 10);
    grid on;
    if k == 1, legend(h, cond_names, 'Location', 'northwest', 'Box', 'off'); end
end
drawnow;
saveas(fig, fullfile(savedir, 'Bo2026_AmygdalaFIR_cluster_timecourses.png'));

fprintf('\nDone. PNGs written to: %s\n', savedir);

end
