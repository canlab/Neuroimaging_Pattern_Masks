function visualize_contents()
% visualize_contents  Render Miao 2026 Social/ToM Bayes-factor maps to png_images/.
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
    't_social_with_tom_audio',        't_social_with_tom_audio.nii'
    't_social_with_tom_text',         't_social_with_tom_text.nii'
    't_social_with_tom_conjunction',  't_social_with_tom_conjunction.nii'
    't_tom_with_social_audio',        't_tom_with_social_audio.nii'
    't_tom_with_social_text',         't_tom_with_social_text.nii'
    't_tom_with_social_conjunction',  't_tom_with_social_conjunction.nii'
    'BF_class_audio',                 'BF_class_audio.nii'
    'BF_class_text',                  'BF_class_text.nii'
    'BF_class_conjunction',           'BF_class_conjunction.nii'
    };

canlab_render_patterns(imgs);

end
