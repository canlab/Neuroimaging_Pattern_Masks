% Add pain_pathways_canlab2018indices to painpathways atlas object .mat file
% ----------------------------------------------------------------------
% This is a subset of the canlab_2018 combined atlas
% Designed so that it can be pulled out in a single step from the combined
% atlas, using original parcel definitions in the combined atlas.
% % Some labels will be different from painpathways, as in painpathways:
% (a) some regions are combined (e.g., thal), esp. in the main 24-region object
% (b) regions are separated into L/R subclusters in painpathways (thal, amy)
% (c) some use different definitions from other atlases besides canlab2018 combined, like Hy 

% 67 defined regions known to be involved in nociceptive transmission in animals 
% - does not include regions important for pain/pain regulation that do not
% receive nociceptive affererents. These regions may be of great importance
% for pain, but they are not included in this set of nociceptive brain
% pathways.
%
% ----------------------------------------------------------------------

load pain_pathways_atlas_obj

dosave = true;

atlas_obj = load_atlas('canlab2018_2mm');
[thal, wh_labels] = select_atlas_subset(atlas_obj, {'Thal'});

% select three zones that comprise pain-related pathways in the thalamus
% ----------------------------------------------------------------------

% sensory thalamus is comprised of VPL and VPM zones.
[vplm, wh1] = select_atlas_subset(atlas_obj, [{'VPL' 'VPM'}]);

% intralaminar thalamus is comprised of intrlaminar and midline groups
[ilthal, wh2] = select_atlas_subset(atlas_obj, [{'Intralam' 'Midline'}]);

% mediodorsal nucleus is a large 'association' nucleus in the thalamus,
% connecting to multiple limbic areas (e.g., amygdala)

[md, wh3] = select_atlas_subset(atlas_obj, {'MD'});

thal = [vplm ilthal md];

% Build NUMERIC index of wh_labels to avoid ambiguity if text labels select
% some wrong regions
wh_labels = wh_labels & (wh1 | wh2 | wh3);


%% select hypothalamus, masking out other regions to clean up (VTA/midbrain)

% hy = select_atlas_subset(atlas_obj, {'Hy'});
[hy, wh1] = select_atlas_subset(atlas_obj, {'Hy'});

wh_labels = wh_labels | wh1;

% test = select_atlas_subset(atlas_obj, find(wh_labels)); test.labels, orthviews(test)

%% Add other regions

[other_regions, wh1] = select_atlas_subset(atlas_obj, {'pbn' 'PAG' 'rvm' 'Amy' '_OP1' 'Ctx_RI' '_OP2' 'Ctx_Ig' '_OP4' 'Ctx_PFcm' 'Ctx_PoI2' 'Ctx_FOP2' 'Ctx_FOP3' 'Ctx_MI_' 'Ctx_45' 'Ctx_FOP4' 'Ctx_AVI'  'Ctx_FOP5'});

wh_labels = wh_labels | wh1;

%% Add cingulate/MPFC regions

obj = select_regions_near_crosshairs(atlas_obj, 'coords', [0 15 42], 'thresh', 20);

% remove 8 BM - this is medial area 8, not of interest here
wh = strfind(obj.labels, 'Ctx_8BM');
wh = cellfun(@isempty, wh);
mylabels = obj.labels(wh);

[mpfc_finegrained, wh1] = select_atlas_subset(atlas_obj, mylabels);

wh_labels = wh_labels | wh1;


%% add S1 regions

% foot areas
obj = select_regions_near_crosshairs(atlas_obj, 'coords', [0 -34 61], 'thresh', 15);
orthviews(obj)
wh = strfind(obj.labels, 'Ctx_PCV');
wh = cellfun(@isempty, wh);
mylabels = obj.labels(wh);  % exclude the one matched above
[s1_foot, wh1] = select_atlas_subset(atlas_obj, [mylabels 'Ctx_5mv_L']); % fine-grained

wh_labels = wh_labels | wh1;


[s1_handplus, wh1] = select_atlas_subset(atlas_obj, {'Ctx_1_' 'Ctx_2_' 'Ctx_3a_' 'Ctx_3b_'});  % fine-grained

wh_labels = wh_labels | wh1;

%% add back S1 region from Lee/Napadow 2018

obj = select_regions_near_crosshairs(atlas_obj, 'coords', [18 -38 72], 'thresh', 15);

[obj, wh1] = select_atlas_subset(atlas_obj, {'Ctx_5L_L' 'Ctx_5L_R'});  % fine-grained

wh_labels = wh_labels | wh1;

% test = select_atlas_subset(atlas_obj, find(wh_labels)); test.labels, orthviews(test)

%% compile into overall atlases

pain_pathways_canlab2018_indx = wh_labels;

[pain_pathways_canlab2018indices] = select_atlas_subset(atlas_obj, find(pain_pathways_canlab2018_indx));


%% save

if dosave
    
    save pain_pathways_atlas_obj -append pain_pathways_canlab2018indices pain_pathways_canlab2018_indx
    
end

%% Plot these and save figures

mkdir figures

montage(pain_pathways_canlab2018indices, 'regioncenters');

if dosave
    saveas(gcf, fullfile('figures', 'pain_pathways_canlab2018indices_regioncenters.png'));
end

montage(pain_pathways_canlab2018indices, 'regioncenters', 'heatmap');

if dosave
    saveas(gcf, fullfile('figures', 'pain_pathways_canlab2018indices2.png'));
end

o2 = montage(pain_pathways_canlab2018indices, 'outline', 'trans');

if dosave
    saveas(gcf, fullfile('figures', 'pain_pathways_canlab2018indices_outlines.png'));
end


% --------------------------------------------------------------
% --------------------------------------------------------------
%  Main atlas construction is done at this point
% --------------------------------------------------------------
% --------------------------------------------------------------

