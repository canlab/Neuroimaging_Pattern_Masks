function visualize_contents()
% visualize_contents  Render Reddan 2018 CS+ vs CS- SVM pattern into png_images/.
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
    'Reddan2018_CSplus_unthresh', 'IE_ImEx_Acq_Threat_SVM_nothresh.nii'
    'Reddan2018_CSplus_01unc',    'IE_ImEx_Acq_Threat_SVM_01thresh.nii.gz'
    'Reddan2018_CSplus_05FDR',    'IE_ImEx_Acq_Threat_SVM_05FDR.nii.gz'
    };

canlab_render_patterns(imgs);

end
