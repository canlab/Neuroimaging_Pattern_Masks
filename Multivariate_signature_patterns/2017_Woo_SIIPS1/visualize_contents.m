function visualize_contents()
% visualize_contents  Render SIIPS1 figures into png_images/.
%
% Woo et al. 2017, Nature Communications. See contents_description.md for
% loading instructions and citations. Delegates rendering to
% docs/canlab_render_patterns.m (shared helper).

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

% Make the shared helper resolvable.
helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'SIIPS1_weighted_mean',       'nonnoc_v11_4_137subjmap_weighted_mean.nii'
    'SIIPS1_subcluster_pattern',  'nonnoc_v11_4_subcluster_maps_fdr05_pattern_wttest.nii'
    };

canlab_render_patterns(imgs);

end
