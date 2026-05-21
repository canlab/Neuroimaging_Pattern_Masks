function visualize_contents()
% visualize_contents  Render Eisenbarth 2016 GSR/HR patterns into png_images/.
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
    'Eisenbarth2016_GSR_unthresh',  'ANS_Eisenbarth_JN_2016_GSR_pattern.img'
    'Eisenbarth2016_GSR_p005thr',   'ANS_Eisenbarth_JN_2016_GSR_pattern_p005thr.img.gz'
    'Eisenbarth2016_HR_unthresh',   'ANS_Eisenbarth_JN_2016_HR_pattern.img'
    'Eisenbarth2016_HR_p005thr',    'ANS_Eisenbarth_JN_2016_HR_pattern_p005thr.img.gz'
    };

canlab_render_patterns(imgs);

end
