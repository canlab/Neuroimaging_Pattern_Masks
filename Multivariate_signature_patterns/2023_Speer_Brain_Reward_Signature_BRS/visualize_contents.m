function visualize_contents()
% visualize_contents  Render Speer 2023 BRS into png_images/.
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
    'Speer2023_BRS', 'Reward_Signature_bootstrapped_0.5.nii.gz'
    };

canlab_render_patterns(imgs);

end
