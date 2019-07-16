atlas_obj = load_atlas('canlab2018_2mm');
thal = select_atlas_subset(atlas_obj, {'Thal'});

% select three zones that comprise pain-related pathways in the thalamus
% ----------------------------------------------------------------------

% sensory thalamus is comprised of VPL and VPM zones.
vplm = select_atlas_subset(thal, [{'VPL' 'VPM'}], 'flatten');

% intralaminar thalamus is comprised of intrlaminar and midline groups
ilthal = select_atlas_subset(thal, [{'Intralam' 'Midline'}], 'flatten');

% mediodorsal nucleus is a large 'association' nucleus in the thalamus,
% connecting to multiple limbic areas (e.g., amygdala)

md = select_atlas_subset(thal, {'MD'});

thal = [vplm ilthal md];

% split lateral thalamus into L and R
thal = split_atlas_into_contiguous_regions(thal);

thal_finegrained = select_atlas_subset(thal, [{'VPL' 'VPM' 'Intralam' 'Midline' 'MD'}]);

%% select hypothalamus, masking out other regions to clean up (VTA/midbrain)

% hy = select_atlas_subset(atlas_obj, {'Hy'});
% don't like this quite as much, it didn't reslice so nicely. use Morel:

thal2 = load_atlas('thalamus');
hy = select_atlas_subset(thal2, {'Hythal'});

%% Add brainstem regions

pbn = select_atlas_subset(atlas_obj, {'pbn'});
pag = select_atlas_subset(atlas_obj, {'PAG'});
rvm = select_atlas_subset(atlas_obj, {'rvm'});

%% Add amygdala regions, dividing into L and R hemisphere for each

amy_finegrained = select_atlas_subset(atlas_obj, {'Amy'});
amy_finegrained = split_atlas_into_contiguous_regions(amamy_finegrainedy);

amy = select_atlas_subset(atlas_obj, {'Amy'}, 'flatten');
amy = split_atlas_into_contiguous_regions(amy);

%% Add dpINS / S2 and other insula regions

% dpins zones
% op1 = select_atlas_subset(atlas_obj, {'_OP1'});
% ri = select_atlas_subset(atlas_obj, {'Ctx_RI'});
% op2 = select_atlas_subset(atlas_obj, {'_OP2'});     % medial posterior operculum
% ig = select_atlas_subset(atlas_obj, {'Ctx_Ig'});

dpins_finegrained = select_atlas_subset(atlas_obj, {'_OP1' 'Ctx_RI' '_OP2' 'Ctx_Ig'}); % fine-grained regions

dpins = select_atlas_subset(atlas_obj, {'_OP1' 'Ctx_RI' '_OP2' 'Ctx_Ig'}, 'flatten'); % larger regions more suitable for identifying patterns within
dpins = split_atlas_into_contiguous_regions(dpins);

% S2 zones
% op4 = select_atlas_subset(atlas_obj, {'_OP4'});
% pfcm = select_atlas_subset(atlas_obj, {'Ctx_PFcm'});

s2_finegrained = select_atlas_subset(atlas_obj, {'_OP4' 'Ctx_PFcm'});

s2 = select_atlas_subset(atlas_obj, {'_OP4' 'Ctx_PFcm'}, 'flatten');
s2 = split_atlas_into_contiguous_regions(s2);

% mid-insula zones
% poi2 = select_atlas_subset(atlas_obj, {'Ctx_PoI2'}); % posterior/ventral. extends into mid-insula, and large. keep separate from post insula.
% fop2 = select_atlas_subset(atlas_obj, {'Ctx_FOP2'});
% fop3 = select_atlas_subset(atlas_obj, {'Ctx_FOP3'});
% mi = select_atlas_subset(atlas_obj, {'Ctx_MI'});     % anterior mid-insula

midins_finegrained = select_atlas_subset(atlas_obj, {'Ctx_PoI2' 'Ctx_FOP2' 'Ctx_FOP3' 'Ctx_MI'});

midins = select_atlas_subset(atlas_obj, {'Ctx_PoI2' 'Ctx_FOP2' 'Ctx_FOP3' 'Ctx_MI'}, 'flatten');
midins = split_atlas_into_contiguous_regions(midins);

% auditory/temporal regions
% lbelt = select_atlas_subset(atlas_obj, {'Ctx_LBelt'}); % nope
% five2 = select_atlas_subset(atlas_obj, {'Ctx_52'});
% poi1 = select_atlas_subset(atlas_obj, {'Ctx_PoI1'});
% ccc = select_atlas_subset(atlas_obj, {'Ctx_MBelt'}); orthviews(ccc)
% ccc = select_atlas_subset(atlas_obj, {'Ctx_LBelt'}); orthviews(ccc)

% still missing some coverage after just OP
% use cluster_graphic_select to find names

% Anterior insula regions
% obj_within_15mm = select_regions_near_crosshairs(atlas_obj);

ains_finegrained = select_atlas_subset(atlas_obj, {'Ctx_45' 'Ctx_FOP4' 'Ctx_AVI'  'Ctx_FOP5'});

ains = select_atlas_subset(atlas_obj, {'Ctx_45' 'Ctx_FOP4' 'Ctx_AVI'  'Ctx_FOP5'}, 'flatten');
ains = split_atlas_into_contiguous_regions(ains);

% 'Ctx_AAIC' very ventral anterior insula
% 'Ctx_45' VLPFC
% 'Ctx_47' Multiple zones in posterior OFC

%% Add cingulate/MPFC regions

obj = select_regions_near_crosshairs(atlas_obj, 'coords', [0 15 42], 'thresh', 20);

% remove 8 BM - this is medial area 8, not of interest here
wh = strfind(obj.labels, 'Ctx_8BM');
wh = cellfun(@isempty, wh);
mylabels = obj.labels(wh);

mpfc_finegrained = select_atlas_subset(atlas_obj, mylabels);

mpfc = select_atlas_subset(atlas_obj, mylabels, 'flatten');

%% add S1 regions

% foot areas
obj = select_regions_near_crosshairs(atlas_obj, 'coords', [0 -34 61], 'thresh', 15);
orthviews(obj)
wh = strfind(obj.labels, 'Ctx_PCV');
wh = cellfun(@isempty, wh);
mylabels = obj.labels(wh);  % exclude the one matched above
s1_foot = select_atlas_subset(atlas_obj, [mylabels 'Ctx_5mv_L']); % fine-grained

s1_foot_L = select_atlas_subset(s1_foot, {'_L'}, 'flatten');
s1_foot_R = select_atlas_subset(s1_foot, {'_R'}, 'flatten');
s1_foot_L.labels = {'s1_foot_L'};
s1_foot_R.labels = {'s1_foot_R'};

s1_handplus = select_atlas_subset(atlas_obj, {'Ctx_1_' 'Ctx_2_' 'Ctx_3a_' 'Ctx_3b_'});  % fine-grained

s1_handplus_L = select_atlas_subset(s1_handplus, {'_L'}, 'flatten');
s1_handplus_R = select_atlas_subset(s1_handplus, {'_R'}, 'flatten');
s1_handplus_L.labels = {'s1_handplus_L'};
s1_handplus_R.labels = {'s1_handplus_R'};


%% compile into overall atlases

pain_pathways_finegrained = [thal_finegrained hy pbn pag rvm amy_finegrained dpins_finegrained s2_finegrained midins_finegrained ains_finegrained mpfc_finegrained s1_foot s1_handplus];

pain_pathways = [thal hy pbn pag rvm amy dpins_finegrained s2 midins ains mpfc s1_foot_L s1_foot_R s1_handplus_L s1_handplus_R];

%% save

save pain_pathways_atlas_obj pain_pathways_finegrained pain_pathways

%% Get NPS patterns within each region, in a series of clusters

%% Mask PDM1 with these regions
% Extract local PDM1 patterns within each region, save in .vals field

[pdm, pdmnames] = load_image_set('pain_pdm');
pdm1 = get_wh_image(pdm, 1);

create_figure('montage'); axis off
o2 = montage(pdm1);

pain_regions = atlas2region(pain_pathways);

pain_regions_pdm1 = extract_data(pain_regions, pdm1);
% pain_regions_pdm1(1).all_data -> weights are stored in in all_data
% save in .val field, which extract_data will use to extract
for i = 1:length(pain_regions_pdm1)
    pain_regions_pdm1(i).val = pain_regions_pdm1(i).all_data';
    pain_regions_pdm1(i).Z = pain_regions_pdm1(i).all_data;
end
k = length(pain_regions_pdm1);
is_empty = false;
for i = 1:k, is_empty(i) = all(abs(pain_regions_pdm1(i).val) < 1000*eps | isnan( pain_regions_pdm1(i).val)); end
pain_regions_pdm1(is_empty) = [];

o2 = montage(pain_regions_pdm1, 'colormap');  % PDM1 weights within anatomically defined nociceptive pathways
o2 = legend(o2);
o2 = title_montage(o2, 5, 'PDM1 weights within anatomically defined nociceptive pathways');

%% Extract NPS patterns within each region, save in .vals field

[nps, npsnames] = load_image_set('npsplus');
siips = get_wh_image(nps, 4);
nps = get_wh_image(nps, 1);

pain_regions_nps = extract_data(pain_regions, nps);
for i = 1:length(pain_regions_nps)
    pain_regions_nps(i).val = pain_regions_nps(i).all_data';
    pain_regions_nps(i).Z = pain_regions_nps(i).all_data;
end
is_empty = false;
k = length(pain_regions_nps);
for i = 1:k, is_empty(i) = all(abs(pain_regions_nps(i).val) < 1000*eps | isnan( pain_regions_nps(i).val)); end
pain_regions_nps(is_empty) = [];

pain_regions_siips = extract_data(pain_regions, siips);
for i = 1:length(pain_regions_siips)
    pain_regions_siips(i).val = pain_regions_siips(i).all_data';
    pain_regions_siips(i).Z = pain_regions_siips(i).all_data;
end
is_empty = false;
k = length(pain_regions_siips);
for i = 1:k, is_empty(i) = all(abs(pain_regions_siips(i).val) < 1000*eps | isnan( pain_regions_siips(i).val)); end
pain_regions_siips(is_empty) = [];

o2 = montage(pain_regions_nps, 'colormap');  % PDM1 weights within anatomically defined nociceptive pathways
o2 = legend(o2);
o2 = title_montage(o2, 5, 'NPS weights within anatomically defined nociceptive pathways');

o2 = montage(pain_regions_siips, 'colormap');  % PDM1 weights within anatomically defined nociceptive pathways
o2 = legend(o2);
o2 = title_montage(o2, 5, 'SIIPS weights within anatomically defined nociceptive pathways');


%%
save pain_pathways_region_obj_with_local_patterns pain_regions_pdm1 pain_regions_nps pain_regions_siips

%% extract data from MV-mediation

% - method to apply local signatures stored within regions
% - region(x).vals contains the pattern weights
test_dat = load_image_set('emotionreg');
[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_pdm1, test_dat);

k = length(pain_regions_pdm1);
rr = [];
for i = 1:k
    rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
    % rr(i) = rvals(1, 2);
end





