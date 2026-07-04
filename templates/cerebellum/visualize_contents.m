function visualize_contents()
% visualize_contents  Render cerebellar template / mask to png_images/.
% See contents_description.md.
%
% The folder is mostly SPM warp / segmentation files rather than
% patterns. We render the isolated cerebellum and SUIT reference for
% visual sanity-checking of the inter-template alignment.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

xfm_dir = fullfile(this_dir, 'transforms', ...
    'MNI152NLin6Asym_T1_1mm_cerebellum_to_SUIT');

imgs = {
    'SUIT_template',                          fullfile(xfm_dir, 'SUIT.nii')
    'MNI152NLin6Asym_cerebellum_isolated',    fullfile(xfm_dir, 'c_MNI152NLin6Asym_T1_1mm.nii')
    'MNI152NLin6Asym_cerebellum_pmask',       fullfile(xfm_dir, 'c_MNI152NLin6Asym_T1_1mm_pcereb.nii')
    };

canlab_render_patterns(imgs, 'kinds', {'montage', 'isosurface'});

end
