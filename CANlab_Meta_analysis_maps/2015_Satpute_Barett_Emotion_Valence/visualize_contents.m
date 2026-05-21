function visualize_contents()
% visualize_contents  Stub renderer for Satpute/Barrett 2015 valence DB.
% No NIfTI result maps are bundled here -- the folder ships only the
% MKDA SETUP.mat and a coordinate database. Re-run the MKDA pipeline
% with SETUP.mat to produce voxelwise maps. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = cell(0, 2);   % no NIfTIs bundled in this folder

if isempty(imgs)
    warning('visualize_contents:nothing_to_render', ...
        ['No NIfTI maps are bundled in this folder. ', ...
         'Re-run the MKDA pipeline against SETUP.mat to generate maps.']);
    return
end

canlab_render_patterns(imgs);

end
