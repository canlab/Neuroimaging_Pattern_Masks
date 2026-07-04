function visualize_contents()
% visualize_contents  Render Kragel 2018 MFC bPLS patterns into png_images/.
% See contents_description.md.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
old_pwd = pwd; cd(this_dir);
cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>

helper_dir = fullfile(this_dir, '..', '..', 'docs');
if exist(fullfile(helper_dir, 'canlab_render_patterns.m'), 'file')
    addpath(helper_dir);
end

regions = {'Wholebrain','MFC','dMFC','aMCC','pMCC','pACC','sgACC','vmPFC'};
conds   = {'Pain','Negative_Emotion','Cognitive_Control'};

imgs = {};
for r = 1:numel(regions)
    for c = 1:numel(conds)
        label = sprintf('Kragel2018_bPLS_%s_%s', regions{r}, conds{c});
        fname = sprintf('bPLS_%s_%s.nii.gz', regions{r}, conds{c});
        imgs(end+1, :) = {label, fname}; %#ok<AGROW>
    end
end

canlab_render_patterns(imgs);

end
