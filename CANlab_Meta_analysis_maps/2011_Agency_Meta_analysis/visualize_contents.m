function visualize_contents()
% visualize_contents  Render Miele 2011 agency MKDA meta maps.
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
    'Miele2011_Agency_Activation_proportion', 'Activation_proportion.hdr'
    'Miele2011_Agency_Activation_FWE_height', 'Activation_FWE_height.hdr'
    'Miele2011_Agency_Activation_FWE_extent', 'Activation_FWE_extent.hdr'
    'Miele2011_Agency_Activation_FWE_all',    'Activation_FWE_all.hdr'
};

canlab_render_patterns(imgs);

end
