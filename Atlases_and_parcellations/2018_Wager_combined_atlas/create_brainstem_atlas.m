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

%% add other regions not in canlab_load_ROI
% Load morel, and select more region groups

morelfile = which('Morel_thalamus_atlas_object.mat');

morel = load(morelfile); morel = morel.atlas_obj;
%%
group_codes = {'Pu' 'VA' 'VM' 'AM' 'AV' 'LD'};  % each is a group to load and add

% AD is tiny and next to DM

for i = 1:length(group_codes)
    
roi_atlas = select_atlas_subset(morel, group_codes(i), 'flatten'); orthviews(roi_atlas)

% a = input('press...');

thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);

end

% view stuff left out so far
% orthviews(select_atlas_subset(thalamus_atlas, 1))

roi_atlas = select_atlas_subset(morel,{'VL'}, 'flatten'); orthviews(roi_atlas)
thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);

roi_atlas = select_atlas_subset(morel,{'LP'}, 'flatten'); orthviews(roi_atlas)
thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);

thalamus_atlas.atlas_name = 'Morel_groups_combined';
thalamus_atlas.labels{1} = 'Other_thalamus';

%% Save

cd('/Users/tor/Documents/Code_Repositories/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2018_Wager_combined_atlas');

savefile = 'Thalamus_atlas_combined_Morel.mat';
save(savefile, 'thalamus_atlas');


%% Load thalamic regions from Brainnetome
% better maybe for functional divisions of anterior nuc.
% but maybe not...

% bn = load(bnfile); bn = bn.atlas_obj;
% bnthal = select_atlas_subset(bn, {'Tha'});
