function visualize_contents()
% visualize_contents  Render one montage / isosurface per Pandora
% white-matter bundle method into png_images/. The Pandora atlases are
% 4-D probabilistic NIfTIs rather than labelled atlas objects, so we
% defer to canlab_render_patterns. See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

imgs = {
    'AFQ',         fullfile('AFQ',         'AFQ.nii.gz')
    'AFQclipped',  fullfile('AFQclipped',  'AFQclipped.nii.gz')
    'Recobundles', fullfile('Recobundles', 'Recobundles.nii.gz')
    'TractSeg',    fullfile('TractSeg',    'TractSeg.nii.gz')
    'Tracula',     fullfile('Tracula',     'Tracula.nii.gz')
    'Xtract',      fullfile('Xtract',      'Xtract.nii.gz')
    };

canlab_render_patterns(imgs, 'kinds', {'montage', 'isosurface'});

end
