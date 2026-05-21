function visualize_contents()
% visualize_contents  Render Quah/Saggar 2025 RDoC factor maps into png_images/.
% See contents_description.md. Mirrors quah_factor_map_montages.m but uses
% the shared canlab_render_patterns helper to also produce isosurface PNGs.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

ff = dir(fullfile(this_dir, 'PA*.nii.gz'));
imgs = cell(numel(ff), 2);
for i = 1:numel(ff)
    label = regexprep(ff(i).name, '\.nii\.gz$', '');
    imgs{i, 1} = ['Quah2025_' label];
    imgs{i, 2} = ff(i).name;
end

canlab_render_patterns(imgs);

end
