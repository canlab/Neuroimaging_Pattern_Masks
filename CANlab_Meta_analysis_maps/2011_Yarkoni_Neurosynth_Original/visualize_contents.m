function visualize_contents()
% visualize_contents  Render Yarkoni 2011 Neurosynth original pain maps.
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
    'Yarkoni2011_Pain_z',             'pain_2s_z_val.hdr'
    'Yarkoni2011_Pain_z_FDR05',       'pain_2s_z_val_FDR_05.hdr'
    'Yarkoni2011_PainMinusEmo_z',     'pain-emotion_2s_z_val.hdr'
    'Yarkoni2011_PainMinusEmo_FDR05', 'pain-emotion_2s_z_val_FDR_05.hdr'
};

canlab_render_patterns(imgs);

end
