function visualize_contents()
% visualize_contents  Render Pauli 2016 striatal parcellation + cortical maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Make NIfTIs in subfolders findable by canlab_render_patterns.
addpath(fullfile(this_dir, 'striatum'));
addpath(fullfile(this_dir, 'cortical'));

imgs = {
    'Pauli2016_striatum_5cluster',  'Pauli_bg_cluster_mask_5.nii'
    'Pauli2016_striatum_17cluster', 'Pauli_bg_cluster_mask_17.nii.gz'
    'Pauli2016_cortex_Ca',          'Pauli_bg_nb_param_rank_fst_Ca.nii.gz'
    'Pauli2016_cortex_Cp',          'Pauli_bg_nb_param_rank_fst_Cp.nii.gz'
    'Pauli2016_cortex_Pa',          'Pauli_bg_nb_param_rank_fst_Pa.nii.gz'
    'Pauli2016_cortex_Pp',          'Pauli_bg_nb_param_rank_fst_Pp.nii.gz'
    'Pauli2016_cortex_VS',          'Pauli_bg_nb_param_rank_fst_VS.nii.gz'
};

canlab_render_patterns(imgs);

end
