function visualize_contents()
% visualize_contents  Render Wager 2004 attention-switching meta map.
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
    'Wager2004_AttentionSwitching_clusters', 'switch_clusters_mix_final.hdr'
};

canlab_render_patterns(imgs);

end
