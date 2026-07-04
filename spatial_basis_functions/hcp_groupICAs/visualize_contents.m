function visualize_contents()
% visualize_contents  Render HCP PTN1200 group-ICA maps into png_images/.
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
    'hcp_groupICA_d15', 'hcp_d15_ICs.nii.gz'
    'hcp_groupICA_d25', 'hcp_d25_ICs.nii.gz'
    'hcp_groupICA_d50', 'hcp_d50_ICs.nii.gz'
    };

canlab_render_patterns(imgs);

end
