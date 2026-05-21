function visualize_contents()
% visualize_contents  Render Krishnan 2016 VPS into png_images/.
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
    'VPS_unthresholded',    'bmrk4_VPS_unthresholded.nii'
    'VPS_fdr05',            'bmrk4_VPS_fdr05.img.gz'
    'VPS_nooccipital',      'Krishnan_2016_VPS_bmrk4_Without_Occipital_Lobe.nii'
    };

canlab_render_patterns(imgs);

end
