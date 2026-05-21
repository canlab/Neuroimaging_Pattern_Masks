function visualize_contents()
% visualize_contents  Render Margulies 2016 principal gradient into png_images/.
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
    'Margulies_MNI152NLin2009cAsym_grad1', 'MNI152NLin2009cAsym_margulies_grad1.nii.gz'
    'Margulies_MNI152NLin6Asym_grad1',     'MNI152NLin6Asym_margulies_grad1.nii.gz'
    };

canlab_render_patterns(imgs);

end
