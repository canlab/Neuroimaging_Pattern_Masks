function visualize_contents()
% visualize_contents  Render Zunhammer 2021 N=603 pain-placebo meta maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

pain_re = fullfile(this_dir, 'pain_effect_study_as_random_effect');
pla_re  = fullfile(this_dir, 'placebo_effect_study_as_random_effect');
rcorr_re= fullfile(this_dir, 'placebo_rating_correlation_study_as_random_effect');

imgs = {
    'Zunhammer2021_Pain_g_unthresh',          fullfile(pain_re, 'full_pain_g_unthresh.nii.gz')
    'Zunhammer2021_Pain_g_pperm_FWE05',       fullfile(pain_re, 'full_pain_g_pperm_FWE05.nii.gz')
    'Zunhammer2021_Placebo_g_unthresh',       fullfile(pla_re,  'full_pla_g_unthresh.nii.gz')
    'Zunhammer2021_Placebo_g_pperm_tfceFWE05',fullfile(pla_re,  'full_pla_g_pperm_tfce_FWE05.nii.gz')
    'Zunhammer2021_PlaceboRatingCorr_unthresh',fullfile(rcorr_re, 'full_pla_rrating_unthresh.nii.gz')
    'Zunhammer2021_PlaceboRatingCorr_FWE05',  fullfile(rcorr_re, 'full_pla_rrating_pperm_FWE05.nii.gz')
};

canlab_render_patterns(imgs);

end
