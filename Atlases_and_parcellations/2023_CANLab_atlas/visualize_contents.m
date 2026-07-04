function visualize_contents()
% visualize_contents  Render CANLab2023 combined atlas across all
% 8 build combos into png_images/. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'canlab2023_coarse_fmriprep20_2mm',  'canlab2023'
    'canlab2023_coarse_fmriprep20_1mm',  'canlab2023_1mm'
    'canlab2023_coarse_fsl6_2mm',        'canlab2023_fsl6'
    'canlab2023_coarse_fsl6_1mm',        'canlab2023_fsl6_1mm'
    'canlab2023_fine_fmriprep20_2mm',    'canlab2023_fine'
    'canlab2023_fine_fmriprep20_1mm',    'canlab2023_fine_1mm'
    'canlab2023_fine_fsl6_2mm',          'canlab2023_fine_fsl6'
    'canlab2023_fine_fsl6_1mm',          'canlab2023_fine_fsl6_1mm'
    };

canlab_render_atlas(atlases);

end
