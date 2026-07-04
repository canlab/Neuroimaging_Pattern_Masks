function visualize_contents()
% visualize_contents  Render Geuter 2020 cPDM + 10-PDM stack into png_images/.
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
    'Geuter2020_cPDM_combined',   'Geuter_2020_cPDM_combined_pain_map.nii'
    'Geuter2020_PDM10_unthresh',  'All_PDM10_unthresholded.nii.gz'
    'Geuter2020_PDM10_FDR',       'All_PDM10_FDR_thresholded.nii.gz'
    };

canlab_render_patterns(imgs);

end
