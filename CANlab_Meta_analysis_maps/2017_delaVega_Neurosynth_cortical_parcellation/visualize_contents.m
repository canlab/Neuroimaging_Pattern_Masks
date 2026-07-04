function visualize_contents()
% visualize_contents  Render de la Vega 2017 cortical parcellation.
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
    'delaVega2017_WholeBrain_70',  'wb_70.nii.gz'
    'delaVega2017_LFC_70',         'lfc_70.nii.gz'
    'delaVega2017_LFC_5',          'lfc_5.nii.gz'
    'delaVega2017_Atlas_regions',  'delaVega2017_neurosynth_atlas_regions.hdr'
};

% Integer-coded parcellations -- use atlas-style isosurface threshold.
canlab_render_patterns(imgs, 'thresh', 'fdr05');

end
