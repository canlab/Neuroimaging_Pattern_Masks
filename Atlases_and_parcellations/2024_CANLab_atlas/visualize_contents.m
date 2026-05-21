function visualize_contents()
% visualize_contents  Render CANLab2024 (and openCANLab2024) combined
% atlas across the canonical build combos into png_images/, plus the
% pain_pathways2024 subset. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'canlab2024_coarse_fmriprep20_2mm',      'canlab2024'
    'canlab2024_coarse_fsl6_2mm',            'canlab2024_fsl6'
    'canlab2024_fine_fmriprep20_2mm',        'canlab2024_fine'
    'canlab2024_fine_fsl6_2mm',              'canlab2024_fine_fsl6'
    'opencanlab2024_coarse_fmriprep20_2mm',  'opencanlab2024'
    'opencanlab2024_fine_fmriprep20_2mm',    'opencanlab2024_fine'
    'painpathways2024',                      'painpathways2024'
    'painpathways2024_finegrained',          'painpathways2024_finegrained'
    };

canlab_render_atlas(atlases);

end
