function visualize_contents()
% visualize_contents  Summarise Rosenberg 2016 saCPM connectivity model.
% This is a connectome-based predictive model, not a voxelwise pattern,
% so we render a summary figure of edge counts rather than the standard
% surface/montage/isosurface trio.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

cpm_path = which('saCPM.mat');
if isempty(cpm_path), cpm_path = fullfile(this_dir, 'saCPM.mat'); end
if ~exist(cpm_path, 'file')
    warning('visualize_contents:missing', 'saCPM.mat not found.');
    return
end

S = load(cpm_path);
fprintf('Loaded saCPM.mat — fields:\n');
disp(fieldnames(S));

fig = figure('Name', 'saCPM_summary', 'Position', [100 100 720 420]);

% Try to plot the edge sets if standard fields are present.
fields = fieldnames(S);
text_lines = {};
edge_counts = struct();
for i = 1:numel(fields)
    val = S.(fields{i});
    if isnumeric(val)
        if islogical(val) || all(val(:) == 0 | val(:) == 1)
            edge_counts.(fields{i}) = nnz(val);
        end
        text_lines{end+1} = sprintf('%s: %s [%s], %d nz', ...
            fields{i}, class(val), strjoin(string(size(val)), 'x'), nnz(val)); %#ok<AGROW>
    else
        text_lines{end+1} = sprintf('%s: %s', fields{i}, class(val)); %#ok<AGROW>
    end
end

axes('Position', [0.05 0.05 0.9 0.9]); axis off;
text(0, 1, ['saCPM.mat contents (Rosenberg et al. 2016)'], ...
    'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
text(0, 0.92, text_lines, 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'FontName', 'Courier');

saveas(fig, fullfile(savedir, 'saCPM_summary.png'));
fprintf('Done. Summary written to: %s\n', fullfile(savedir, 'saCPM_summary.png'));

end
