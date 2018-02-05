% create_brainstem_atlas_group
% 
% Create a brain atlas with anatomically defined region groups, from
% various atlases/papers.  Uses canlab_load_ROI
%
% Notes: functional atlases probably do not [yet] have very good
% subdivisions, and there is a clear demarcation of functions, inputs, and
% outputs by anatomical subnuclei.

% Define: 1 mm space by default, based on HCP image
% This initial image covers the whole space

bstem_atlas = atlas(which('canlab_brainstem.img'));
bstem_atlas = threshold(bstem_atlas, .2);

bstemimg = fmri_data(which('brainstem_mask_tight_2018.img'));
bstem_atlas = apply_mask(bstem_atlas, bstemimg);

% this has some other regions in it (e.g., mammillary bodies), so would be better to use the
% 'noreplace' option when merging it with other atlases.

orthviews(bstem_atlas);

%see also: bstemimg.fullpath = fullfile(pwd, 'brainstem_mask_tight_2018.img');

%% add other regions
% ...replacing voxels where new one overlaps

regionnames = {'pag' 'PBP' 'sn' 'VTA' 'rn' 'pbn' 'lc' 'rvm' 'nts' 'drn' 'mrn' 'sc' 'ic'};
% NEW ONES TOO

for i = 1:length(regionnames)
    regionname = regionnames{i};
    
    [~, roi_atlas] = canlab_load_ROI(regionname);
    orthviews(roi_atlas);

    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

