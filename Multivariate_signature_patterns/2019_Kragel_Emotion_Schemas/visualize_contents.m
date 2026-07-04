function visualize_contents()
% visualize_contents  Render Kragel 2019 emotion-schema patterns into png_images/.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Find all PLS_betas_<emotion>.nii.gz files in this folder.
nii = dir(fullfile(this_dir, 'PLS_betas_*.nii.gz'));

imgs = cell(numel(nii), 2);
for i = 1:numel(nii)
    label = regexprep(nii(i).name, '^PLS_betas_', '');
    label = regexprep(label, '\.nii\.gz$', '');
    label = ['Kragel2019_' regexprep(label, '\s+', '_')]; %#ok<AGROW>
    imgs{i, 1} = label;
    imgs{i, 2} = nii(i).name;
end

canlab_render_patterns(imgs);

end
