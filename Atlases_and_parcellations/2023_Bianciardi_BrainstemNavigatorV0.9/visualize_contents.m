function visualize_contents()
% visualize_contents  Render Bianciardi Brainstem Navigator v0.9 atlas
% (4 space x resolution combinations) into png_images/. See
% contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'bianciardi_MNI152NLin2009cAsym',      'bianciardi'
    'bianciardi_MNI152NLin2009cAsym_2mm',  'bianciardi_2mm'
    'bianciardi_MNI152NLin6Asym',          'bianciardi_fsl6'
    'bianciardi_MNI152NLin6Asym_2mm',      'bianciardi_fsl6_2mm'
    };

canlab_render_atlas(atlases);

end
