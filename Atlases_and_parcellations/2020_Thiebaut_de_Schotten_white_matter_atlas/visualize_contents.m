function visualize_contents()
% visualize_contents  Optional renderer for the Thiebaut de Schotten
% functional-component white-matter atlas. The folder does NOT
% redistribute NIfTIs; download components from the NeuroVault
% collections listed in README.rst, drop them next to this script,
% then point the imgs cell at them. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Example: imgs = { 'component_01', 'TdS_component_01.nii.gz' };
imgs = {};

if isempty(imgs)
    fprintf(['No components configured. Download maps from NeuroVault ' ...
        '(see README.rst) and add entries to the imgs cell.\n']);
    return
end

canlab_render_patterns(imgs, 'kinds', {'montage', 'isosurface'});

end
