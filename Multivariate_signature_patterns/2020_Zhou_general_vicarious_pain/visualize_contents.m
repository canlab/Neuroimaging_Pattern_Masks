function visualize_contents()
% visualize_contents  Render Zhou 2020 vicarious-pain patterns into png_images/.
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
    'Zhou2020_General_unthresh',  'General_vicarious_pain_pattern_unthresholded.nii'
    'Zhou2020_General_FDR05',     'General_vicarious_pain_pattern_FDR05_Boot10000.nii'
    'Zhou2020_NS_unthresh',       'NS_vicarious_pain_pattern_unthresholded.nii'
    'Zhou2020_NS_FDR05',          'NS_vicarious_pain_pattern_FDR05.nii'
    'Zhou2020_FE_unthresh',       'FE_vicarious_pain_pattern_unthresholded.nii'
    'Zhou2020_FE_FDR05',          'FE_vicarious_pain_pattern_FDR05.nii'
    };

canlab_render_patterns(imgs);

end
