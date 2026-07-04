function visualize_contents()
% visualize_contents  Render Woo 2014 rejection signature into png_images/.
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
    'Rejection_dpsp_weights',  'dpsp_rejection_vs_others_weights_final.nii'
    'Heat_dpsp_weights',       'dpsp_hot_vs_others_weights_final.nii.gz'
    'dACC_hw_searchlight',     'dACC_hw_pattern_sl6mm.nii.gz'
    'dACC_rf_searchlight',     'dACC_rf_pattern_sl6mm.nii.gz'
    };

canlab_render_patterns(imgs);

end
