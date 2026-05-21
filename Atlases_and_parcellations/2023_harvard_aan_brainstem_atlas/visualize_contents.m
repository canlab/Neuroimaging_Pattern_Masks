function visualize_contents()
% visualize_contents  Render Harvard AAN brainstem atlas (v2.0) into
% png_images/. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'harvard_aan_v2_MNI152NLin2009cAsym',  'harvard_aan'
    'harvard_aan_v2_MNI152NLin6Asym',      'harvard_aan_fsl6'
    };

canlab_render_atlas(atlases);

end
