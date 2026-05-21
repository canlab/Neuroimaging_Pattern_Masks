function visualize_contents()
% visualize_contents  Render Wager 2003 emotion meta-analysis density maps.
% See contents_description.md. Writes surface / montage / isosurface PNGs
% into png_images/ using the shared canlab_render_patterns helper.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'Wager2003_Approach_density',          'dens_ap.hdr'
    'Wager2003_Avoidance_density',         'dens_av.hdr'
    'Wager2003_Approach_gt_Avoidance',     'ap-av_enl.hdr'
    'Wager2003_Avoidance_gt_Approach',     'av-ap_enl.hdr'
};

canlab_render_patterns(imgs);

end
