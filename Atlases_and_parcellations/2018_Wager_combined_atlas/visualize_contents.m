function visualize_contents()
% visualize_contents  Render CANlab 2018 combined whole-brain atlas
% (plus its main sub-atlases) into png_images/. See
% contents_description.md.
%
% Note: pre-existing PNGs in png_images/ that match the labels below
% will be overwritten by canlab_render_atlas. Author-curated reference
% figures should be moved out of png_images/ first if they must be
% preserved.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'CANlab_2018_combined',        'canlab2018'
    'CANlab_2018_combined_2mm',    'canlab2018_2mm'
    'Thalamus_combined',           'thalamus'
    'BasalGanglia_combined',       'Basal_ganglia_combined_atlas_object.mat'
    'Brainstem_combined',          'brainstem_combined_atlas_object.mat'
    };

canlab_render_atlas(atlases);

end
