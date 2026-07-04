function visualize_contents()
% visualize_contents  Render Silvestrini/Rainville 2020 dACC + Stroop patterns.
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
    'Silvestrini2020_dACC_pain',  'dACC_pain_pattern_wani_121416.nii'
    'Silvestrini2020_dACC_task',  'dACC_task_pattern_wani_121416.nii'
    'Silvestrini2020_Stroop_WB',  'stroop_pattern_wani_121416.nii'
    };

canlab_render_patterns(imgs);

end
