function visualize_contents()
% visualize_contents  Render Etkin & Wager 2007 anxiety meta maps.
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
    'Etkin2007_AnxietyDisorders_chi2_p005', 'chi2p_map_UNC_p005.hdr'
    'Etkin2007_FearConditioning_p005',      'Fear_conditioning_p005.hdr'
};

canlab_render_patterns(imgs);

end
