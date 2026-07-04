function visualize_contents()
% visualize_contents  Render Ashar 2017 placebo-review meta maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

% Put the MKDA subfolders on path so canlab_render_patterns can resolve.
sub_inc = fullfile(this_dir, 'Placebo_increases');
sub_dec = fullfile(this_dir, 'Placebo_decreases');
sub_val = fullfile(this_dir, 'Experience_imagery_recall_Neg-Pos');
addpath(sub_inc); addpath(sub_dec); addpath(sub_val);

imgs = {
    'Ashar2017_CoreAppraisalNetwork',  'core_appraisal_network.nii'
    'Ashar2017_PlaceboIncreases_prop', fullfile(sub_inc, 'Activation_proportion.img')
    'Ashar2017_PlaceboDecreases_prop', fullfile(sub_dec, 'Activation_proportion.img')
    'Ashar2017_NegVsPos_Valence',      fullfile(sub_val, 'Valence_Neg-Pos.img')
};

canlab_render_patterns(imgs);

end
