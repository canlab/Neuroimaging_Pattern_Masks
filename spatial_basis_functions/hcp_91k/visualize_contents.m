function visualize_contents()
% visualize_contents  Render selected HCP 91k spectral basis modes to png_images/.
% See contents_description.md.
%
% The dscalar file holds 200 eigenmodes on HCP 32k_fs_LR + 2mm subcortex.
% This script projects a small subset to a volumetric NIfTI on the fly
% (using wb_command if available) and renders surface/montage/isosurface
% triplets via canlab_render_patterns. Edit `mode_idx` below to pick
% which modes to render.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

cifti_file = fullfile(this_dir, 'spectral_bases_200.dscalar.nii');
if ~exist(cifti_file, 'file')
    warning('Spectral bases CIFTI not found: %s', cifti_file);
    return
end

% Edit this list to render additional modes (1..200).
mode_idx = [1 2 3 4 5 10 20 50];

vol_dir = fullfile(this_dir, 'derived_volumes');
if ~exist(vol_dir, 'dir'), mkdir(vol_dir); end

imgs = cell(numel(mode_idx), 2);
for i = 1:numel(mode_idx)
    k = mode_idx(i);
    out_nii = fullfile(vol_dir, sprintf('hcp91k_mode_%03d.nii.gz', k));
    if ~exist(out_nii, 'file')
        cmd = sprintf(['wb_command -cifti-separate "%s" COLUMN ' ...
                       '-volume-all "%s" -label '], cifti_file, out_nii);
        try
            system(cmd);
        catch
            warning('Could not run wb_command for mode %d; skipping.', k);
            continue
        end
    end
    imgs{i, 1} = sprintf('hcp91k_mode_%03d', k);
    imgs{i, 2} = out_nii;
end

% Drop empty rows in case wb_command was unavailable.
keep = ~cellfun(@isempty, imgs(:, 1));
imgs = imgs(keep, :);

if isempty(imgs)
    fprintf('No volumetric modes available to render. See README for Workbench instructions.\n');
    return
end

canlab_render_patterns(imgs);

end
