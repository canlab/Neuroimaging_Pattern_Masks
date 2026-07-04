function visualize_contents()
% visualize_contents  Render Kragel & LaBar 2015 emotion BPLS into png_images/.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

emos = {'amused','angry','content','fearful','neutral','sad','surprised'};
imgs = cell(numel(emos), 2);
for i = 1:numel(emos)
    imgs{i, 1} = ['Kragel2015_BPLS_' emos{i}];
    imgs{i, 2} = sprintf('mean_3comp_%s_group_emotion_PLS_beta_BSz_10000it.img', emos{i});
end

canlab_render_patterns(imgs);

end
