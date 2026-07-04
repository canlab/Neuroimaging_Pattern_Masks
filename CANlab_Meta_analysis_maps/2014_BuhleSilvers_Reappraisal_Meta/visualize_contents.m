function visualize_contents()
% visualize_contents  Render Buhle & Silvers 2014 reappraisal meta maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'BuhleSilvers2014_EmoReg_Meta',         'Buhle_Silvers_2014_Emotion_Regulation_Meta.nii.gz'
    'BuhleSilvers2014_EmoReg_Meta_thresh',  'Buhle_Silvers_2014_Emotion_Regulation_Meta_thresh.hdr'
    'BuhleSilvers2014_EmoReg_z_rescaled',   'Buhle_Silvers_2014_Emotion_Regulation_Meta_z_rescaled_.hdr'
    'BuhleSilvers2014_Activation_FWE_height','Activation_FWE_height.hdr'
    'BuhleSilvers2014_Activation_FWE_extent','Activation_FWE_extent.hdr'
};

canlab_render_patterns(imgs);

end
