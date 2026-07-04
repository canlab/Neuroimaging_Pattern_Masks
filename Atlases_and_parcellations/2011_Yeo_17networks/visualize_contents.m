function visualize_contents()
% visualize_contents  Render Yeo 17-network cortical parcellation into
% png_images/. See contents_description.md.
%
% The file shipped here is a hard-labelled integer NIfTI rather than a
% CANlab atlas object, so we render via canlab_render_patterns().

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'Yeo2011_17Networks',  'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm.nii.gz'
    };

canlab_render_patterns(imgs);

end
