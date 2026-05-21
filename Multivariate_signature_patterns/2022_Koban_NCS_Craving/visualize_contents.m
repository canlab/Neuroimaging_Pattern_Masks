function visualize_contents()
% visualize_contents  Render Koban 2023 NCS craving pattern into png_images/.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Make NCS_weightmaps subfolder findable for the drug/food maps.
addpath(genpath(this_dir));

imgs = {
    'Koban2023_NCS_general', 'craving_wmapN99_boot10K_02-May-2022.img'
    };

canlab_render_patterns(imgs);

end
