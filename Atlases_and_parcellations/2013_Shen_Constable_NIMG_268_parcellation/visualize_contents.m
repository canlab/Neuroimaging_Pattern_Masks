function visualize_contents()
% visualize_contents  Render Shen 268-node whole-brain functional
% parcellation into png_images/. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'Shen_MNI152NLin2009cAsym',  'shen_fmriprep20'
    'Shen_MNI152NLin6Asym',      'shen_fsl6'
    'Shen_Colin27v1998',         'shen'
    };

canlab_render_atlas(atlases);

end
