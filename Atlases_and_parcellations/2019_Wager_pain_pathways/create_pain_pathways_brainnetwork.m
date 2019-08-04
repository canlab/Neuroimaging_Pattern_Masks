
% ----------------------------------------------------------------------
% 68 defined regions known to be involved in nociceptive transmission in animals ('finegrained')
% - does not include regions important for pain/pain regulation that do not
% receive nociceptive affererents. These regions may be of great importance
% for pain, but they are not included in this set of nociceptive brain
% pathways.
% 
% grouped into 24 larger zones for application of local patterns and
% pathway identification
%
% - 5 in thalamus and hypothalamus (sensory thalamus x 2, intralaminar,
% mediodorsal, and hypothalamus)
% - 4 in brainstem (PBN x 2, PAG, RVM)
% - Amygdala (2)
% - Cortex (13): Bilateral dpINS, S2, mIns, aIns, S1lowerlimb S1 upperlimb, and aMCC/MPFC.
% ----------------------------------------------------------------------

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

% re-do names for ease of visualization
thal.labels = {'Thal_VPLM_R' 'Thal_VPLM_L' 'Thal_IL' 'Thal_MD'};


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
amy_finegrained = split_atlas_into_contiguous_regions(amy_finegrained);

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

% clean up - remove small 3rd region in R angular gyrus
dpins = select_atlas_subset(dpins, 1:2);

% rename for ease of viz
dpins.labels = {'dpIns_L' 'dpIns_R'};

% S2 zones
% op4 = select_atlas_subset(atlas_obj, {'_OP4'});
% pfcm = select_atlas_subset(atlas_obj, {'Ctx_PFcm'});

s2_finegrained = select_atlas_subset(atlas_obj, {'_OP4' 'Ctx_PFcm'});

s2 = select_atlas_subset(atlas_obj, {'_OP4' 'Ctx_PFcm'}, 'flatten');
s2 = split_atlas_into_contiguous_regions(s2);

s2.labels = {'S2_L' 'S2_R'};

% mid-insula zones
% poi2 = select_atlas_subset(atlas_obj, {'Ctx_PoI2'}); % posterior/ventral. extends into mid-insula, and large. keep separate from post insula.
% fop2 = select_atlas_subset(atlas_obj, {'Ctx_FOP2'});
% fop3 = select_atlas_subset(atlas_obj, {'Ctx_FOP3'});
% mi = select_atlas_subset(atlas_obj, {'Ctx_MI'});     % anterior mid-insula

midins_finegrained = select_atlas_subset(atlas_obj, {'Ctx_PoI2' 'Ctx_FOP2' 'Ctx_FOP3' 'Ctx_MI_'});

midins = select_atlas_subset(atlas_obj, {'Ctx_PoI2' 'Ctx_FOP2' 'Ctx_FOP3' 'Ctx_MI_'}, 'flatten');
midins = split_atlas_into_contiguous_regions(midins);
midins.labels = {'mIns_L' 'mIns_R'};

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
ains.labels = {'aIns_L' 'aIns_R'};

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

mpfc.labels = {'aMCC_MPFC'};

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

pain_pathways = [thal hy pbn pag rvm amy dpins s2 midins ains mpfc s1_foot_L s1_foot_R s1_handplus_L s1_handplus_R];

%% save

save pain_pathways_atlas_obj pain_pathways_finegrained pain_pathways

%% Plot these and save figures

mkdir figures

montage(pain_pathways, 'regioncenters');

saveas(gcf, fullfile('figures', 'pain_pathways_regioncenters.png'));

montage(pain_pathways_finegrained, 'regioncenters');

saveas(gcf, fullfile('figures', 'pain_pathways_finegrained.png'));

montage(pain_pathways_finegrained, 'regioncenters', 'heatmap');

saveas(gcf, fullfile('figures', 'pain_pathways_finegrained2.png'));

o2 = montage(pain_pathways_finegrained, 'outline', 'trans');

saveas(gcf, fullfile('figures', 'pain_pathways_finegrained_outlines.png'));


%% Get pain signatures patterns within each region, in a series of clusters
% -----------------------------------------------------------------------------

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

saveas(gcf, fullfile('figures', 'pain_pathways_PDM1.png'));

%% Extract NPS patterns within each region, save in .vals field

% SEE ALSO:
% [parcel_means, parcel_pattern_expression, parcel_valence] = extract_data(pain_pathways, test_data_obj, 'pattern_expression', nps, 'cosine_similarity');

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

saveas(gcf, fullfile('figures', 'pain_pathways_NPS.png'));

o2 = montage(pain_regions_siips, 'colormap');  % PDM1 weights within anatomically defined nociceptive pathways
o2 = legend(o2);
o2 = title_montage(o2, 5, 'SIIPS weights within anatomically defined nociceptive pathways');


saveas(gcf, fullfile('figures', 'pain_pathways_SIIPS.png'));
%% Save

save pain_pathways_region_obj_with_local_patterns pain_regions_pdm1 pain_regions_nps pain_regions_siips


%% Add cPDM (done later)
% ---------------------------------------------------------------

% load pain_pathways_region_obj_with_local_patterns  % has local weights extracted
% load pain_pathways_atlas_obj

cPDM = load('CombinedPDM.mat');
cPDM = cPDM.CPDM2;

pain_regions = atlas2region(pain_pathways);

pain_regions_cpdm = extract_data(pain_regions, cPDM);

% pain_regions_pdm1(1).all_data -> weights are stored in in all_data
% save in .val field, which extract_data will use to extract
for i = 1:length(pain_regions_cpdm)
    pain_regions_cpdm(i).val = pain_regions_cpdm(i).all_data';
    pain_regions_cpdm(i).Z = pain_regions_cpdm(i).all_data;
end
k = length(pain_regions_cpdm);
is_empty = false;
for i = 1:k, is_empty(i) = all(abs(pain_regions_cpdm(i).val) < 1000*eps | isnan( pain_regions_cpdm(i).val)); end
pain_regions_cpdm(is_empty) = [];

montage(cPDM, 'colormap')
saveas(gcf, fullfile('figures', 'cPDM_full_brain_weights.png'));

montage(pain_regions_cpdm, 'colormap')
saveas(gcf, fullfile('figures', 'cPDM_painpathways_weights.png'));

o2 = montage(pain_regions_cpdm, 'colormap');
% add outlines:
for i = 1:length(pain_regions_cpdm) 
    o2 = addblobs(o2, pain_regions_cpdm(i), 'outline', 'color', [.2 .2 .2]); 
    %set(o2.activation_maps{i + 1}.blobhandles, 'LineWidth', 1);  % one map aleady registered, so add 1 here
end
saveas(gcf, fullfile('figures', 'cPDM_painpathways_weights_and_regions.png'));

o2 = removeblobs(o2);
% add outlines only:
for i = 1:length(pain_regions_cpdm) 
    o2 = addblobs(o2, pain_regions_cpdm(i), 'outline', 'color', [.2 .2 .2]); 
    %set(o2.activation_maps{i + 1}.blobhandles, 'LineWidth', 1);  % one map aleady registered, so add 1 here
end
saveas(gcf, fullfile('figures', 'Painpathways_region_outlines.png'));

% regioncenters
montage(pain_regions_cpdm, 'colormap', 'regioncenters');
saveas(gcf, fullfile('figures', 'cPDM_painpathways_regioncenters.png'));

save pain_pathways_region_obj_with_local_patterns -append pain_regions_cpdm


%% 

% --------------------------------------------------------------
% --------------------------------------------------------------
%  Data extraction and saving is done at this point
% --------------------------------------------------------------
% --------------------------------------------------------------


%% extract data from test pain dataset (BMRK3) and plot

% Load test dataset
test_dat = load_image_set('pain');  % bmrk3 data

% Load regions with local patterns stored in them
% Contains pain_regions_nps pain_regions_siips pain_regions_pdm1
load pain_pathways_region_obj_with_local_patterns % in Neuroimaging_Pattern_Masks

% NPS
% ---------------------------------------------------------------
[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_nps, test_dat);

k = length(pain_regions_nps);
rr = [];
for i = 1:k
    rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
end

names = strrep({pain_regions_nps.shorttitle}, '_', ' ');
create_figure('Barplot_correlation_with_region_avgs: NPS');
bar(rr);
ylabel('Pearson''s r');
title('Correlation between local pattern response and region average');
set(gca, 'XTick', 1:length(rr), 'XLim', [0 length(rr)+1], 'XTickLabel', names, 'XTickLabelRotation', 45);
saveas(gcf, fullfile('figures', 'r_local_pattern_region_average_nps.png'));

% SIIPS
% ---------------------------------------------------------------

[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_siips, test_dat);

k = length(pain_regions_siips);
rr = [];
for i = 1:k
    rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
    % rr(i) = rvals(1, 2);
end

names = strrep({pain_regions_siips.shorttitle}, '_', ' ');

create_figure('Barplot_correlation_with_region_avgs: SIIPS');
bar(rr);
ylabel('Pearson''s r');
title('Correlation between local pattern response and region average');
set(gca, 'XTick', 1:length(rr), 'XLim', [0 length(rr)+1], 'XTickLabel', names, 'XTickLabelRotation', 45);
saveas(gcf, fullfile('figures', 'r_local_pattern_region_average_siips.png'));

% PDM1
% ---------------------------------------------------------------

[regions_with_testdata, local_pattern_responses] = extract_data(pain_regions_pdm1, test_dat);

k = length(pain_regions_pdm1);
rr = [];
for i = 1:k
    rr(i) = corr(regions_with_testdata(i).dat, local_pattern_responses{i});
    % rr(i) = rvals(1, 2);
end

names = strrep({pain_regions_pdm1.shorttitle}, '_', ' ');

create_figure('Barplot_correlation_with_region_avgs: PDM1');
bar(rr);
ylabel('Pearson''s r');
title('Correlation between local pattern response and region average');
set(gca, 'XTick', 1:length(rr), 'XLim', [0 length(rr)+1], 'XTickLabel', names, 'XTickLabelRotation', 45);
saveas(gcf, fullfile('figures', 'r_local_pattern_region_average_pdm1.png'));





%% Extract data from single_trial_database

run('/Users/torwager/Google Drive/Wagerlab_Single_Trial_Pain_Analyses/Individual_Diffs/scripts/a_set_up_paths_always_run_first.m')
run('/Users/torwager/Google Drive/Wagerlab_Single_Trial_Pain_Analyses/Individual_Diffs/scripts/b_load_check_prep_for_analysis.m')

for i = 1:nstudies
    
    matfile = get_var(study_canlab_dataset{i}, 'ST_matfile');
    
end

% See script: prep1_extrract_single_trial_pain_data.m

% and:
% prep1_extract_single_trial_pain_data_cpdm_only

load data/pain_pathways_single_trial_data_2019 Subject_Table ST_*

% Cleaned data
load data/pain_pathways_single_trial_data_2019 ST_cleaned

%% Examine associations with key brainstem regions - fine-grained

%rmat = corr(ST_cleaned.fine_regions);
rmat = partialcorr(ST_cleaned.fine_regions, 'type', 'Spearman');

labels = pain_pathways_finegrained.labels;

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)


% pain_pathways_finegrained.labels

% max row values: most connected regions

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% pre-windsorizing (z-score only) standard Pearson's r
% VPLVPM_R connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_Ig_R Ctx_PoI2_R Ctx_MI_R
% VPLVPM_L connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Ctx_Ig_L Ctx_MI_L
% IntralamMidline_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Ctx_a24pr_L Ctx_a32pr_R
% Thal_MD_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_a24pr_L Ctx_a24pr_R Ctx_a32pr_R
% Hythal connected to  Hythal Bstem_PAG Amygdala_CM__L Amygdala_CM__R Amygdala_SF__R Amygdala_SF__L Amygdala_LB__L
% pbn_R connected to  VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R Ctx_PoI2_R
% pbn_L connected to  VPLVPM_R VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG Ctx_PoI2_R
% Bstem_PAG connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Hythal pbn_R pbn_L Bstem_PAG
% rvm_R connected to  pbn_R pbn_L rvm_R Ctx_PoI2_R Ctx_MI_L Ctx_MI_R Ctx_p32pr_R

% post-windsorizing  standard Pearson's r
% VPLVPM_R connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_Ig_R Ctx_PoI2_R Ctx_MI_R
% VPLVPM_L connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Amygdala_AStr__L Ctx_Ig_L
% IntralamMidline_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Bstem_PAG Ctx_a24pr_L Ctx_a32pr_R
% Thal_MD_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Ctx_a24pr_L Ctx_a24pr_R Ctx_a32pr_R
% Hythal connected to  Hythal Bstem_PAG Amygdala_CM__L Amygdala_CM__R Amygdala_SF__R Amygdala_SF__L Amygdala_LB__L
% pbn_R connected to  Hythal pbn_R pbn_L Bstem_PAG rvm_R Ctx_PoI2_R Ctx_MI_R
% pbn_L connected to  VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R Ctx_PoI2_R
% Bstem_PAG connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Hythal pbn_R pbn_L Bstem_PAG
% rvm_R connected to  pbn_R pbn_L rvm_R Ctx_PoI2_R Ctx_MI_L Ctx_MI_R Ctx_p32pr_R

% VPLVPM to ipsiateral Ctx_Ig_L Ctx_MI_L maybe Ctx_PoI2_R

% IntralamMidline_M to Ctx_a24pr_L Ctx_a32pr_R
% Thal_MD_M to Ctx_a24pr_L Ctx_a24pr_R Ctx_a32pr_R

% Hythal to Bstem_PAG 
% PBN to Hythal Bstem_PAG Ctx_PoI2_R (R from both L and R)
% rvm to Ctx_PoI2_R Ctx_MI_L Ctx_MI_R Ctx_p32pr_R

% Spearman's 
% VPLVPM_R connected to  VPLVPM_R VPLVPM_L IntralamMidline_M pbn_R rvm_R Amygdala_CM__R Ctx_RI_R Ctx_Ig_R Ctx_PoI2_R Ctx_5mv_R Ctx_2_R
% VPLVPM_L connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Amygdala_CM__L Amygdala_AStr__L Ctx_RI_L Ctx_Ig_L Ctx_p24pr_L Ctx_2_L
% IntralamMidline_M connected to  VPLVPM_R VPLVPM_L IntralamMidline_M Thal_MD_M Hythal Bstem_PAG Ctx_RI_R Ctx_Ig_L Ctx_Ig_R Ctx_AVI_R Ctx_1_R Ctx_3a_L
% Thal_MD_M connected to  IntralamMidline_M Thal_MD_M Amygdala_CM__L Amygdala_CM__R Amygdala_SF__L Ctx_OP1_R Ctx_a24pr_L Ctx_a32pr_R
% Hythal connected to  IntralamMidline_M Hythal pbn_R pbn_L Bstem_PAG Amygdala_SF__R Amygdala_SF__L Ctx_45_L Ctx_AVI_L Ctx_a24pr_L Ctx_a32pr_R
% pbn_R connected to  VPLVPM_R VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R Amygdala_AStr__R Ctx_PFcm_R Ctx_PoI2_R Ctx_MI_L Ctx_AVI_R Ctx_p32pr_L
% pbn_L connected to  VPLVPM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R Amygdala_LB__L Ctx_PoI2_L Ctx_FOP3_R Ctx_AVI_L Ctx_a24pr_L Ctx_1_L
% Bstem_PAG connected to  IntralamMidline_M Hythal pbn_R pbn_L Bstem_PAG Amygdala_CM__R Amygdala_AStr__L Ctx_PFcm_L Ctx_FOP2_L Ctx_33pr_R Ctx_a24pr_R
% rvm_R connected to  VPLVPM_R VPLVPM_L pbn_R pbn_L rvm_R Ctx_OP1_L Ctx_RI_L Ctx_PoI2_R Ctx_FOP5_L Ctx_p32pr_R Ctx_5mv_R

%% Examine associations with key brainstem regions: PDM1

% rmat = corr([ST_cleaned.pdm1 ST_cleaned.pain_rating]);
% labels = format_strings_for_legend({pain_regions_pdm1.shorttitle ['pain_rating']});

%rmat = corr(ST_cleaned.pdm1);
rmat = partialcorr(ST_cleaned.pdm1, 'type', 'Spearman');
labels = format_strings_for_legend({pain_regions_pdm1.shorttitle});

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('PDM1 inter-region single-trial correlations');
saveas(gcf, 'figures/PDM1_local_pattern_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% Thal VPLM R connected to  Thal VPLM R Thal VPLM L Thal IL dpIns R
% Thal VPLM L connected to  Thal VPLM R Thal VPLM L Thal IL dpIns L
% Thal IL connected to  Thal VPLM R Thal VPLM L Thal IL Thal MD Hythal Bstem PAG aMCC MPFC
% Thal MD connected to  Thal IL Thal MD Amy R
% Hythal connected to  Thal IL Hythal pbn R Bstem PAG mIns R
% pbn R connected to  Thal VPLM R Hythal pbn R pbn L Bstem PAG rvm R aIns R
% pbn L connected to  Thal VPLM L Hythal pbn R pbn L Bstem PAG rvm R aIns L
% Bstem PAG connected to  Thal IL Hythal pbn R pbn L Bstem PAG aMCC MPFC
% rvm R connected to  Thal VPLM R Thal VPLM L Hythal pbn R pbn L rvm R s1 foot R


%% Examine associations with key brainstem regions: SIIPS

% rmat = corr(ST_cleaned.SIIPS);
rmat = partialcorr(ST_cleaned.SIIPS, 'type', 'Spearman');
labels = format_strings_for_legend({pain_regions_siips.shorttitle});

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('SIIPS inter-region single-trial correlations');
saveas(gcf, 'figures/SIIPS_local_pattern_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% Thal VPLM R connected to  Thal VPLM R Thal VPLM L Thal MD dpIns R
% Thal VPLM L connected to  Thal VPLM R Thal VPLM L Thal MD Bstem PAG aIns L
% Thal IL connected to  Thal IL Thal MD Hythal s1 handplus R
% Thal MD connected to  Thal VPLM R Thal VPLM L Thal IL Thal MD Bstem PAG Amy R
% Hythal connected to  Thal IL Hythal Amy L
% pbn R connected to  Thal VPLM R pbn R pbn L rvm R mIns R
% pbn L connected to  Thal VPLM R Thal VPLM L pbn R pbn L rvm R s1 handplus R
% Bstem PAG connected to  Thal VPLM R Thal MD Bstem PAG aMCC MPFC
% rvm R connected to  pbn R pbn L rvm R Amy L

%% Examine associations with key brainstem regions: NPS

%rmat = corr(ST_cleaned.NPS);
rmat = partialcorr(ST_cleaned.NPS, 'type', 'Spearman');

labels = format_strings_for_legend({pain_regions_nps.shorttitle});

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('NPS inter-region single-trial correlations');
saveas(gcf, 'figures/NPS_local_pattern_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 
    
    fprintf('%s connected to ', labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% Partial Spearman corr:
% Thal VPLM R connected to  Thal VPLM R Thal VPLM L Thal IL Thal MD dpIns R
% Thal VPLM L connected to  Thal VPLM R Thal VPLM L Thal IL dpIns L
% Thal IL connected to  Thal VPLM R Thal VPLM L Thal IL Thal MD Hythal pbn R Bstem PAG Amy R aIns R
% Thal MD connected to  Thal VPLM R Thal IL Thal MD aMCC MPFC
% Hythal connected to  Hythal Amy L aIns L
% pbn R connected to  Thal VPLM R Thal VPLM L Thal IL pbn R Bstem PAG Amy R mIns R
% Bstem PAG connected to  Thal IL pbn R Bstem PAG aMCC MPFC
% Amy R connected to  Amy R mIns R
% Amy L connected to  Hythal Amy L aIns L

%% Examine associations with key brainstem regions: regions

rmat = corr(ST_cleaned.big_regions);
labels = format_strings_for_legend({pain_regions_siips.shorttitle});

create_figure('r');
imagesc(rmat, [-1 1])
axis tight
set(gca, 'YDir', 'reverse', 'YTick', 1:length(rmat),  'YTickLabel', labels, 'XTick', 1:length(rmat), 'XTickLabel', labels, 'XTickLabelRotation', 45);
colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)
title('Region averages inter-region single-trial correlations');
saveas(gcf, 'figures/Region_averages_ST_connectivity_matrix.png');

for i = 1:9 % do for brainstem regions 9 = rvm
    
    fprintf('%s connected to ', pain_pathways.labels{i});
    wh = rmat(i, :) > prctile(rmat(i, 10:end), 90); lab = pain_pathways.labels(wh);
    for j = 1:length(lab), fprintf(' %s', lab{j}); end
    fprintf('\n');
    
end

% Thal_VPLM_R connected to  Thal_VPLM_R Thal_VPLM_L Thal_IL Thal_MD *mIns_R
% Thal_VPLM_L connected to  Thal_VPLM_R Thal_VPLM_L Thal_IL Thal_MD Bstem_PAG *aMCC_MPFC
% Thal_IL connected to  Thal_VPLM_R Thal_VPLM_L Thal_IL Thal_MD Bstem_PAG *aMCC_MPFC
% Thal_MD connected to  Thal_VPLM_R Thal_VPLM_L Thal_IL Thal_MD *aMCC_MPFC
% Hythal connected to  Hythal Bstem_PAG *Amy_L
% pbn_R connected to  Thal_VPLM_R Thal_VPLM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R *mIns_R
% pbn_L connected to  Thal_VPLM_R Thal_VPLM_L Hythal pbn_R pbn_L Bstem_PAG rvm_R *mIns_R
% Bstem_PAG connected to  Thal_VPLM_R Thal_VPLM_L Thal_IL Thal_MD Hythal pbn_R pbn_L Bstem_PAG *aMCC_MPFC
% rvm_R connected to  pbn_R pbn_L rvm_R *mIns_R

%% Notes

% Notes on analysis so far:
% Fine-grained regions seem to work fairly consistently in terms of
% sensible thalamocortical pathways (spinothalamic) and lateral symmetry.
% Pathways for patterns (PDM1, siips, nps) are less bilaterally symmetric
% and make less sense. 
% Thresholding methods are arbitrary (so far), so we need better connectivity metrics I think. 

% ****NOTE: could do average within-person r and t-test on fisher's z
% instead. ***
% This does not allow for good covariate modeling and/or stable estimation
% of partial correlations. We have a lot of trials (27174) so we can do
% this well, but we probably want to make more assumptions and not allow
% every subject to have a random intercept and random slope.  Actually,
% data are z-scored so there is already a random subject intercept removed.





