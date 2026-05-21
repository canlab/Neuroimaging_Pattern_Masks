function visualize_contents()
% visualize_contents  Render Wager & Smith 2003 working-memory meta maps.
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
    'WagerSmith2003_NonExec_density',           'dens_noexec.hdr'
    'WagerSmith2003_Exec_minus_NonExec',        'dens_exec-noexec.hdr'
    'WagerSmith2003_Exec_minus_NonExec_thresh', 'dens_exec-noexec_thr.hdr'
};

canlab_render_patterns(imgs);

end
