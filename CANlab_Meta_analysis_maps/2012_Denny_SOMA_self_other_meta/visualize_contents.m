function visualize_contents()
% visualize_contents  Render Denny 2012 SOMA self/other meta maps.
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
    'Denny2012_Self',                  'self.hdr'
    'Denny2012_AnyOther',              'anyother.hdr'
    'Denny2012_Other_minus_Self',      'SelfOthe_Other-Self.hdr'
    'Denny2012_Activation_proportion', 'Activation_proportion.hdr'
};

canlab_render_patterns(imgs);

end
