function visualize_contents()
% visualize_contents  Render Ashar 2017 empathic care/distress into png_images/.
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
    'Ashar2017_care_unthresh',     'Ashar_2017_empathic_care_marker.nii'
    'Ashar2017_care_fdr05',        'Ashar_2017_empathic_care_markerFDR_05.nii.gz'
    'Ashar2017_distress_unthresh', 'Ashar_2017_empathic_distress_marker.nii'
    'Ashar2017_distress_fdr05',    'Ashar_2017_empathic_distress_markerFDR_05.nii.gz'
    };

canlab_render_patterns(imgs);

end
