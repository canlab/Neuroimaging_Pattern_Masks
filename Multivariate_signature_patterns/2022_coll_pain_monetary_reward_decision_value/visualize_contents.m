function visualize_contents()
% visualize_contents  Render Coll 2022 pain/money/shock patterns into png_images/.
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
    'Coll2022_painvalue_unthresh',         'painvalue_weights.nii.gz'
    'Coll2022_painvalue_bootz',            'painvalue_bootz.nii.gz'
    'Coll2022_painvalue_fdr05',            'painvalue_weights_fdr05.nii.gz'
    'Coll2022_moneyvalue_unthresh',        'moneyvalue_weights.nii.gz'
    'Coll2022_moneyvalue_bootz',           'moneyvalue_bootz.nii.gz'
    'Coll2022_moneyvalue_fdr05',           'moneyvalue_weights_fdr05.nii.gz'
    'Coll2022_shockintensity_unthresh',    'shockintensity_weights.nii.gz'
    'Coll2022_shockintensity_bootz',       'shockintensity_bootz.nii.gz'
    'Coll2022_shockintensity_fdr05',       'shockintensity_weights_fdr05.nii.gz'
    };

canlab_render_patterns(imgs);

end
