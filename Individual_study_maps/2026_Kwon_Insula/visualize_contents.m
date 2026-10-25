function visualize_contents()
% visualize_contents  Render Kwon 2026 insular zone masks to png_images/.
% See contents_description.md.
%
% The masks are small binary insular zones, so instead of the generic
% canlab_render_patterns triplet this draws slice montages through the
% insula, with each zone in the color used in the paper (Fig. 3).
% Requirements: CanlabCore and SPM12 on path.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

savedir = fullfile(this_dir, 'png_images');
if ~exist(savedir, 'dir'), mkdir(savedir); end

% label, file, color (paper colors)
zones = {
    'ins_gen',          'ins_gen.nii',          [0.55 0.25 0.75]   % domain-general, purple
    'ins_pain_sel',     'ins_pain_sel.nii',     [0.90 0.15 0.15]   % pain-selective, red
    'ins_posemo_sel',   'ins_posemo_sel.nii',   [1.00 0.80 0.00]   % appetitive-selective, yellow
    'ins_negemo_sel',   'ins_negemo_sel.nii',   [0.15 0.40 0.90]   % aversive-selective, blue
    'ins_cog_sel',      'ins_cog_sel.nii',      [0.15 0.70 0.25]   % cognitive control-selective, green
    };

% All five zones together
o = insula_montage();
for k = 1:size(zones, 1)
    o = addblobs(o, region(fmri_data(zones{k, 2}, 'noverbose')), 'color', zones{k, 3});
end
save_fig(fullfile(savedir, 'Kwon2026_insula_zones_montage.png'));

% One montage per zone
for k = 1:size(zones, 1)
    o = insula_montage();
    o = addblobs(o, region(fmri_data(zones{k, 2}, 'noverbose')), 'color', zones{k, 3}); %#ok<NASGU>
    save_fig(fullfile(savedir, [zones{k, 1} '_montage.png']));
end

fprintf('Done. PNGs written to: %s\n', savedir);

end


function o = insula_montage()
% Top row: sagittal slices through left and right insula.
% Bottom row: axial slices through the insula.
figure('Color', 'w', 'Position', [50 50 1600 650]);
xs = [-40 -36 -32 32 36 40];
zs = [-8 0 8 16];
ax_sag = gobjects(1, numel(xs));
ax_ax = gobjects(1, numel(zs));
for i = 1:numel(xs)
    ax_sag(i) = axes('Position', [(i - 1) / numel(xs), 0.5, 1 / numel(xs), 0.48]);
end
for i = 1:numel(zs)
    ax_ax(i) = axes('Position', [0.1 + (i - 1) * 0.2, 0.02, 0.2, 0.46]);
end
o = fmridisplay('overlay', which('keuken_2014_enhanced_for_underlay.img'));
o = montage(o, 'saggital', 'wh_slice', [xs' zeros(numel(xs), 2)], 'existing_axes', ax_sag);
o = montage(o, 'axial', 'wh_slice', [zeros(numel(zs), 2) zs'], 'existing_axes', ax_ax);
end


function save_fig(fname)
set(gcf, 'Color', 'w');
saveas(gcf, fname);
close(gcf);
end
