function visualize_contents()
% visualize_contents  Render Schaefer-Yeo 400-parcel / 17-network
% cortical parcellation into png_images/. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_atlas.m'), 'file')
    addpath(helper_dir);
end

atlases = {
    'Schaefer2018_400parcels',     'Schaefer2018Cortex_atlas_object.mat'
    'Schaefer2018_17networks',     'yeo17networks'
    };

canlab_render_atlas(atlases);

end
