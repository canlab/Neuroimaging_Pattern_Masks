function visualize_contents()
% visualize_contents  Render Murillo 2026 PiFoneM patterns into png_images/.
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
    'Murillo2026_PiFoneM_unthresh',       'PiFoneM_unthresholded.nii'
    'Murillo2026_SVM_RFE_importance',     'svm_rfe_importance_pcr_5000by2500_k10.nii'
    'Murillo2026_LASSOPCR_FDR05_fear',    'cvlassopcr_boots_thresholdedFDR05_k10_fear.nii'
    };

canlab_render_patterns(imgs);

end
