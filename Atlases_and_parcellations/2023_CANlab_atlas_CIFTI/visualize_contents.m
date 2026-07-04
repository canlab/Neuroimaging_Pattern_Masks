function visualize_contents()
% visualize_contents  Renderer for the CIFTI version of CANLab2023.
% This folder holds CIFTI artefacts (when present); the volumetric
% versions used by canlab_render_atlas live in
% ../2023_CANLab_atlas. We render those here for visual reference;
% if the CIFTI .dlabel.nii files are present locally they can be
% inspected in HCP Workbench. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'canlab2023_coarse_fmriprep20',  'canlab2023'
    'canlab2023_fine_fmriprep20',    'canlab2023_fine'
    };

canlab_render_atlas(atlases);

end
