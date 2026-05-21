function visualize_contents()
% visualize_contents  Render Bo 2024 emotion-regulation Bayes-factor maps.
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
    'Bo2024_CommonAppraisal',       'Common Appraisal.nii'
    'Bo2024_ModifiableEmotion',     'Modifiable Emotion.nii'
    'Bo2024_NonmodifiableEmotion',  'Non-modifiable Emotion.nii'
    'Bo2024_ReappraisalOnly',       'Reappraisal Only.nii'
    };

canlab_render_patterns(imgs);

end
