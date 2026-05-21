function visualize_contents()
% visualize_contents  Render Van Oudenhove 2020 Yeo-network coefficient bar chart.
% This classifier uses network-level betas rather than voxelwise weights,
% so we plot the betas instead of running surface/montage/isosurface.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

S = load(which('Visceral_vs_Somatic_betas_Yeo_Networks.mat'));
fields = fieldnames(S);

% Try to find the betas in a common-sense way.
betas = [];
labels = {};
for i = 1:numel(fields)
    val = S.(fields{i});
    if isnumeric(val) && isvector(val)
        betas = val(:);
        labels = arrayfun(@(k) sprintf('Yeo%d', k), 1:numel(betas), 'UniformOutput', false);
        break
    end
end
if isempty(betas)
    warning('visualize_contents:novector', 'Could not auto-find beta vector in .mat; dumping all fields instead.');
    fig = figure('Position',[100 100 720 360]); axis off;
    lines = {};
    for i = 1:numel(fields)
        v = S.(fields{i});
        lines{end+1} = sprintf('%s: %s [%s]', fields{i}, class(v), strjoin(string(size(v)),'x')); %#ok<AGROW>
    end
    text(0.05, 0.95, lines, 'FontSize', 11, 'VerticalAlignment', 'top', 'FontName', 'Courier');
    saveas(fig, fullfile(savedir, 'VanOudenhove2020_inventory.png'));
    return
end

fig = figure('Name','VanOudenhove2020_Yeo_betas','Position',[100 100 640 400]);
bar(betas);
set(gca, 'XTick', 1:numel(betas), 'XTickLabel', labels);
ylabel('SVC coefficient (visceral vs somatic)');
title('Van Oudenhove et al. 2020 — Yeo-network coefficients');
grid on;
saveas(fig, fullfile(savedir, 'VanOudenhove2020_Yeo_betas.png'));
fprintf('Done. PNG written to: %s\n', savedir);

end
