% create_thalamus_atlas_group
% 
% Create a thalamus atlas with anatomically defined region groups, coarser
% than Morel but still containing groups that match different functions.
%
% Notes: functional atlases probably do not [yet] have very good
% subdivisions, and there is a clear demarcation of functions, inputs, and
% outputs by anatomical subnuclei, which can be identified histologically.
% So the anatomical atlases (Morel) are preferred

% 'cm'      Centromedian thalamus   Morel thalamus atlas, Krauth 2010
% 'md'      Mediodorsal thalamus    Morel thalamus atlas, Krauth 2010
% 'lgn'     Lateral geniculate nuc  Morel thalamus atlas, Krauth 2010
% 'mgn'     Medial geniculate nuc   Morel thalamus atlas, Krauth 2010
% 'VPthal'  Ventral posterior thal  Morel thalamus atlas, Krauth 2010
% 'intralaminar_thal' Intralaminar  Morel thalamus atlas, Krauth 2010

% Load selected thalamic regions from Morel

% Define: 1 mm space by default, based on HCP image
% This initial image covers the whole space

[~, thalamus_atlas] = canlab_load_ROI('thalamus');

%% add other regions
% ...replacing voxels where new one overlaps

regionnames = {'cm' 'md' 'lgn' 'mgn' 'VPthal' 'intralaminar_thal'} ;

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
