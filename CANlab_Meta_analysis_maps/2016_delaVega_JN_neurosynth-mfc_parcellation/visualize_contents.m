function visualize_contents()
% visualize_contents  Render de la Vega 2016 Neurosynth MFC parcellation.
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
    'delaVega2016_MFC_kmeans12', 'delaVega_JN_2016_neurosynth-mfc_kmeans_12.nii.gz'
};

% Parcellation -- use an atlas-style isosurface threshold.
canlab_render_patterns(imgs, 'thresh', 'fdr05');

end
