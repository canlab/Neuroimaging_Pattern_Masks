function visualize_contents()
% visualize_contents  Inventory panel for Matthewson/Woo 2019 SCR-pain folder.
% No NIfTI patterns are stored locally; render only a text summary.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

fig = figure('Name', 'Matthewson2019_inventory', 'Position', [100 100 720 360]);
axes('Position', [0.05 0.05 0.9 0.9]); axis off;

text(0, 1, 'Matthewson, Woo et al. 2019 — folder inventory', ...
    'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

d = dir(this_dir);
lines = {};
for i = 1:numel(d)
    if d(i).isdir, continue, end
    if d(i).name(1) == '.', continue, end
    lines{end+1} = sprintf('  %s  (%d bytes)', d(i).name, d(i).bytes); %#ok<AGROW>
end

text(0, 0.92, lines, 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'FontName', 'Courier');

saveas(fig, fullfile(savedir, 'Matthewson2019_inventory.png'));
fprintf('Done. Inventory PNG written to: %s\n', savedir);

end
