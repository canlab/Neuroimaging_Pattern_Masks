function visualize_contents()
% visualize_contents  Render Pourmajidian 2026 mitochondrial
% energetic-capacity maps (six continuous patterns) into png_images/.
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
    'mito_CI',     'CI.nii'
    'mito_CII',    'CII.nii'
    'mito_CIV',    'CIV.nii'
    'mito_MRC',    'MRC.nii'
    'mito_TRC',    'TRC.nii'
    'mito_MitoD',  'MitoD.nii'
    };

canlab_render_patterns(imgs);

end
