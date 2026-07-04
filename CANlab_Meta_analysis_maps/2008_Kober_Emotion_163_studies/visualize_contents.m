function visualize_contents()
% visualize_contents  Render Kober 2008 emotion meta-analysis FWE maps.
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
    'Kober2008_EmoActivation_FWE_all',         'EMOActivation_FWE_all.hdr'
    'Kober2008_ValenceNeg_gt_Pos_FWE_all',     'Valence_Neg-Pos_Neg_FWE_all.hdr'
    'Kober2008_ValencePos_gt_Neg_FWE_all',     'Valence_Neg-Pos_Pos_FWE_all.hdr'
    'Kober2008_ExpPer_NegEmotion_FWE_all',     'Mode_Exp-Per_Neg_FWE_all.hdr'
    'Kober2008_ExpPer_PosEmotion_FWE_all',     'Mode_Exp-Per_Pos_FWE_all.hdr'
};

canlab_render_patterns(imgs);

end
