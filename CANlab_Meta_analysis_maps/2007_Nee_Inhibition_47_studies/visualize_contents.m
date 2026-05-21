function visualize_contents()
% visualize_contents  Render Nee 2007 inhibition meta maps.
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
    'Nee2007_Inhibition_combined', 'inhib_idens15_enl.hdr'
    'Nee2007_Stroop',              'stroop15.hdr'
    'Nee2007_Flanker',             'flanker15.hdr'
    'Nee2007_GoNoGo',              'gng15.hdr'
    'Nee2007_SRC',                 'src15.hdr'
};

canlab_render_patterns(imgs);

end
