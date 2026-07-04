function visualize_contents()
% visualize_contents  Render CIT168 amygdala atlas v1.0.3 into
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
    'CIT168_Amygdala_MNI152NLin2009cAsym',  'cit168_amygdala_fmriprep20'
    'CIT168_Amygdala_MNI152NLin6Asym',      'cit168_amygdala_fsl6'
    };

canlab_render_atlas(atlases);

end
