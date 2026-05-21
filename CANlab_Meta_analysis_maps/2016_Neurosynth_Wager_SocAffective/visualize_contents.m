function visualize_contents()
% visualize_contents  Render CANlab-curated Neurosynth social-affective maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

ff = dir(fullfile(this_dir, '*_pFgA_z_FDR_0.01.nii.gz'));
imgs = cell(numel(ff), 2);
for i = 1:numel(ff)
    term = regexprep(ff(i).name, '_pFgA_z_FDR_0\.01\.nii\.gz$', '');
    imgs{i, 1} = ['NeurosynthCANlab_' term];
    imgs{i, 2} = ff(i).name;
end

canlab_render_patterns(imgs);

end
