function visualize_contents()
% visualize_contents  Render PINES (Chang 2015) into png_images/.
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
    'PINES_Rating_Weights_LOSO',           'Rating_Weights_LOSO_2.nii'
    'PINES_Rating_Boot5000_001unc',        'Rating_LASSO_PCR_Boot5000_2_001_unc.nii.gz'
    'PINES_Rating_Boot5000_fdr05_k10',     'Rating_LASSO_PCR_boot5000_fdr05_k10_2.nii.gz'
    };

canlab_render_patterns(imgs);

end
