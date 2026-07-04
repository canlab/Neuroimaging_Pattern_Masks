function visualize_contents()
% visualize_contents  Render Mosharov/Picard 2025 mitochondrial maps to png_images/.
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
    'Mito_CI',    'CI.nii.gz'
    'Mito_CII',   'CII.nii.gz'
    'Mito_CIV',   'CIV.nii.gz'
    'Mito_MitoD', 'MitoD.nii.gz'
    'Mito_MRC',   'MRC.nii.gz'
    'Mito_TRC',   'TRC.nii.gz'
    };

canlab_render_patterns(imgs);

end
