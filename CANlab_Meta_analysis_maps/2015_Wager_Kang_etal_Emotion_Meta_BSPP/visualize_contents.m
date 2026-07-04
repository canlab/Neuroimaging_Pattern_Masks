function visualize_contents()
% visualize_contents  Render Wager/Kang 2015 BSPP 5-emotion meta maps.
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
    'WagerKang2015_Anger',   'Wager_Kang_PlosCB_emometa_2015_anger.nii.gz'
    'WagerKang2015_Disgust', 'Wager_Kang_PlosCB_emometa_2015_disgust.nii.gz'
    'WagerKang2015_Fear',    'Wager_Kang_PlosCB_emometa_2015_fear.nii.gz'
    'WagerKang2015_Happy',   'Wager_Kang_PlosCB_emometa_2015_happy.nii.gz'
    'WagerKang2015_Sad',     'Wager_Kang_PlosCB_emometa_2015_sad.nii.gz'
};

canlab_render_patterns(imgs);

end
