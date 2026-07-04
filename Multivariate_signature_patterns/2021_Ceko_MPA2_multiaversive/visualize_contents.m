function visualize_contents()
% visualize_contents  Render Ceko 2022 MPA2 patterns into png_images/.
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
    'MPA2_General',     'General_bplsF_unthr.nii'
    'MPA2_Mechanical',  'Mechanical_bplsF_unthr.nii'
    'MPA2_Thermal',     'Thermal_bplsF_unthr.nii'
    'MPA2_Sound',       'Sound_bplsF_unthr.nii'
    'MPA2_Visual',      'Visual_bplsF_unthr.nii'
    };

canlab_render_patterns(imgs);

end
