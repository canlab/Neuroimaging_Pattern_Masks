function visualize_contents()
% visualize_contents  Render Sha 2018 common-networks-of-psychopathology maps.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

sbfc = fullfile(this_dir, 'MKDA_Maps', 'SB_FC');
vbm  = fullfile(this_dir, 'MKDA_Maps', 'VBM');

imgs = {
    'Sha2018_DMN_Pool_FWEheight', fullfile(sbfc, 'DMN', 'DMN_Pool', 'DMN_Pool_FWE_height.nii.gz')
    'Sha2018_DMN_Inc_FWEheight',  fullfile(sbfc, 'DMN', 'DMN_Inc',  'DMN_Inc_FWE_height.nii.gz')
    'Sha2018_DMN_Dec_FWEheight',  fullfile(sbfc, 'DMN', 'DMN_Dec',  'DMN_Dec_FWE_height.nii.gz')
    'Sha2018_FPN_Pool_FWEheight', fullfile(sbfc, 'FPN', 'FPN_Pool', 'FPN_Pool_FWE_height.nii.gz')
    'Sha2018_SN_Pool_FWEheight',  fullfile(sbfc, 'SN',  'SN_Pool',  'SN_Pool_FWE_height.nii.gz')
    'Sha2018_VBM_FWEheight',      fullfile(vbm,  'VBM_FWE_height.nii.gz')
};

canlab_render_patterns(imgs);

end
