function visualize_contents()
% visualize_contents  Render Açıl 2026 mentalizing signatures into png_images/.
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
    'Acil2026_MS_unthresh',        'Mentalization_Boot_Unthr_11-Jun-2024.nii'
    'Acil2026_MS_Self_unthresh',   'Self_Boot_Unthr_11-Jun-2024.nii'
    'Acil2026_MS_Other_unthresh',  'Other_Boot_Unthr_11-Jun-2024.nii'
    'Acil2026_MS_SvO_unthresh',    'SvO_Boot_Unthr_11-Jun-2024.nii'
    };

canlab_render_patterns(imgs);

end
