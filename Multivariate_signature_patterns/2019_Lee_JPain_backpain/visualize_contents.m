function visualize_contents()
% visualize_contents  Render Lee 2019 chronic back-pain markers into png_images/.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Make subfolder .nii files findable.
addpath(genpath(this_dir));

imgs = {
    'Lee2019_PCASL_pairedSVM',  'nilearn_pairedSVM_W_PCASL53.nii.gz'
    'Lee2019_S1_pairedSVM',     'nilearn_pairedSVM_W_S1conn53.nii'
    'Lee2019_S1back_ROIs',      'S1back_Lee_2018coords.nii'
    };

canlab_render_patterns(imgs);

end
