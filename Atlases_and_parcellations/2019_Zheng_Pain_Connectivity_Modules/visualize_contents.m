function visualize_contents()
% visualize_contents  Render Zheng 2019 pain-connectivity backbone +
% module label volumes into png_images/. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Render the Brainnetome backbone used by the paper as the parcellation
% reference. The module label volumes live in module/<condition>/ and
% can be rendered ad-hoc by passing their paths to canlab_render_patterns.
atlases = {
    'brainnetome_backbone',  'brainnetome'
    };

canlab_render_atlas(atlases);

end
