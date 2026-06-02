function visualize_contents()
% visualize_contents  Render Lopez-Sola 2017 fibromyalgia patterns into png_images/.
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
    'LopezSola2017_FM_pain',           'FM_pain_wholebrain.nii.gz'
    'LopezSola2017_FM_multisensory',   'FM_Multisensory_wholebrain.nii.gz'
    'LopezSola2017_NPSp',              'NPSp_Lopez-Sola_2017_PAIN.img.gz'
    'LopezSola2017_NPSn',              'NPSn_Lopez-Sola_2017_PAIN.img.gz'
    'LopezSola2017_rNPS_fdr_pospeaks', 'rNPS_fdr_pospeaks_smoothed.img.gz'
    };

canlab_render_patterns(imgs);

end
