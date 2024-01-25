% this script rereates the canlab2018 atlas except it constrains certain
% structures to more accurately abide by CIFTI volumetric boundaries. For
% example, CA1-3, DG and Subiculum are restricted to the hippocampal
% volume, and other structures are excluded from it, even though the fornix 
% approaches the thalamus, and the=% entorhinal cortex potentially bleeds 
% into this volume in MNI space.

% labels 1 - fine scale
% labels 2 - coarse scale
% labels 3 - larger structures (put, caud, 2nd level glasser areas, e.g.)
% labels 4 - large structural divisions (large cortical areas (glasser colors), thalamus, brainstem, etc.)
% labels 5 - source atlas

% uncomment these lines to run as a standalone script

clear all; close all;
SPACE = 'MNI152NLin2009cAsym';

LIB = '/dartfs-hpc/rc/home/m/f0042vm/software';
%LIB = '/home/bogdan/.matlab';
ROOT = [LIB, '/canlab/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2023_CANLab_atlas/'];

addpath([LIB, '/spm12']);
%addpath([LIB, '/spm/spm12']);

addpath(genpath([LIB, '/canlab/CanlabCore']))
addpath(genpath([LIB, '/canlab/Neuroimaging_Pattern_Masks']))
addpath(genpath([LIB, '/canlab/MasksPrivate']))


if isempty(gcp('nocreate'))
    parpool(32);
end

switch SPACE
    case 'MNI152NLin2009cAsym'
        ALIAS = 'fmriprep20';
        OVERLAY = which('fmriprep20_template.nii.gz');
        TEMPLATE = which('MNI152NLin2009cAsym_T1_1mm.nii.gz');


        labels = fmri_data(sprintf('%s/src/hcp_cifti_subctx_labels_%s.nii.gz',ROOT,SPACE));
        labels_txt = textscan(fopen([ROOT, 'src/hcp_cifti_subctx_labels.txt']),'%s');
    case 'MNI152NLin6Asym'
        ALIAS = 'fsl6';
        OVERLAY = which('fsl6_hcp_template.nii.gz');
        TEMPLATE = which('MNI152NLin6Asym_T1_1mm.nii');
        
        labels = fmri_data([ROOT, 'src/hcp_cifti_subctx_labels.nii']);
        labels_txt = textscan(fopen([ROOT, 'src/hcp_cifti_subctx_labels.txt']),'%s');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only the above code is needed to run child scripts (I think) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cifti_atlas = atlas(labels);
cifti_atlas.labels = labels_txt{1}';


% this is the qsiprep reference volume. The 1mm template is designed for
% use with qsiprep, so let's pull it
%{
qsiprep_template = which('MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz');
if isempty(qsiprep_template)
    websave([ROOT, 'MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz'], 'https://github.com/PennLINC/qsiprep/raw/master/qsiprep/data/mni_1mm_t1w_lps.nii.gz');
end
ref = fmri_data(which('MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz'));
%}
ref = fmri_data(TEMPLATE);

% create L/R masks
brain_mask = atlas(fmri_mask_image(fmri_data().resample_space(TEMPLATE)),...
    'labels',{'brain'},...
    'label_descriptions',{'brain'});
hemi_masks = lateralize(brain_mask);
hemi_L = fmri_mask_image(hemi_masks.select_atlas_subset({'brain_L'}));
hemi_R = fmri_mask_image(hemi_masks.select_atlas_subset({'brain_R'}));

%% create basal ganglia
% produces atlas object bg_dil 
create_bg2023_atlas

%% Cerebellum
if strcmp(SPACE, 'MNI152NLin6Asym')
    suitfile = which('SUIT_Cerebellum_MNI_atlas_object.mat');
else
    suitfile = which(sprintf('SUIT_Cerebellum_%s_atlas_object.mat',SPACE));
end

cerebellum = load(suitfile); cerebellum = cerebellum.atlas_obj.resample_space(ref);
cerebellum.labels_2 = cerebellum.labels;
cerebellum.labels_3 = cerebellum.labels;
labels_4 = {};
for i = 1:num_regions(cerebellum)
    if contains(cerebellum.labels{i},'_L')
        labels_4{end+1} = 'cerebellar_cortex_L';
    elseif contains(cerebellum.labels{i},'_R')
        labels_4{end+1} = 'cerebellar_cortex_R';
    elseif contains(cerebellum.labels{i},'Vermis')
        labels_4{end+1} = 'vermis';
    else
        error('Unexpected cerebellar lobule');
    end
end
cerebellum.labels_4 = labels_4;
cerebellum.labels_5 = repmat({'SUIT/Diedrichsen'},1,num_regions(cerebellum));

% remove white matter structures that aren't well represented in
% grayordinate space
cerebellum = cerebellum.select_atlas_subset(find(~ismember(cerebellum.labels,...
    {'Cblm_Dentate_L', 'Cblm_Dentate_R', 'Cblm_Interposed_L', 'Cblm_Interposed_R', ...
    'Cblm_Fastigial_L', 'Cblm_Fastigial_R', 'Cblm_Vermis_CrusI'})));

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'cerebellum'))));

tic
cerebellum_dil = dilate(cerebellum, cifti_mask);
clear cerebellum
toc

% renormalize the ~52vx p>1
total_p = sum(cerebellum_dil.probability_maps,2);
renorm = total_p > 1;
cerebellum_dil.probability_maps(renorm,:) = cerebellum_dil.probability_maps(renorm,:)./total_p(renorm);


%% Hippocampus and Hippocampus
julich = load_atlas(sprintf('julich_%s',ALIAS)).resample_space(ref);

% SF - superficial amygdala
% LB - laterobasal amygdala
% CM - centromedian amygdala
% CA - cornu ammonis
% DG - dentate gyrus
% HATA - hippocampal-amygdaloid transition area
group_names = {'Subic', {'CA', 'DG'}};
group_labels = {'Subic', 'Hippocampus'};

% note that {'VTM','MF','IF'} are listed as fiber bundles, so we want to
% exclude them from the dilution above
% VTM - Tractus ventromedialis (Stria terminalis) of Brockhaus
% MF - medial fiber bundles
% IF - intermediate fiber bundles
hipp = julich.select_atlas_subset(cat(2,group_names{:}));
hipp.label_descriptions = {};
hipp.labels_2 = {};
hipp.labels_3 = {};
hipp.labels_4 = {};
% add label description and group regions
for i = 1:num_regions(hipp)
    if contains(hipp.labels{i},'L_')
        hipp.label_descriptions{end+1} = 'Left';
    elseif contains(hipp.labels{i},'R_')
        hipp.label_descriptions{end+1} = 'Right';
    end

    if contains(hipp.labels{i},'CA')
        hipp.label_descriptions{end} = [hipp.label_descriptions{end}, ...
            regexprep(hipp.labels{i},'[LR]_CA([0-9])',' hippocampal cornu ammonis \1')];
        hipp.labels_2{end+1} = regexprep(hipp.labels{i},'([LR])_.*','$1_CA');
        hipp.labels_3{end+1} = regexprep(hipp.labels{i},'([LR])_.*','$1_Hipp');
    elseif contains(hipp.labels{i},'DG')
        hipp.label_descriptions{end} = [hipp.label_descriptions{end}, ' hippocampal dentate gyrus'];
        hipp.labels_2{end+1} = hipp.labels{i};
        hipp.labels_3{end+1} = regexprep(hipp.labels{i},'([LR])_.*','$11_Hipp');
    elseif contains(hipp.labels{i},'Subiculum')
        hipp.label_descriptions{end} = [hipp.label_descriptions{end}, ' hippocampal subiculum'];
        hipp.labels_2{end+1} = hipp.labels{i};
        hipp.labels_3{end+1} = regexprep(hipp.labels{i},'([LR])_.*','$1_Subiculum');
    end
    hipp.labels_4{end+1} = regexprep(hipp.labels{i},'([LR])_.*','$1_Hippocampal_Formation');
end
hipp.labels_5 = repmat({'Julich/Amunts'},1,num_regions(hipp));

hipp.references = char([{'Amunts, K., Kedo, O., Kindler, M., Pieperhoff, P., Mohlberg, H., Shah, N. J., Habel, U., Schneider, F., & Zilles, K. (2005). Cytoarchitectonic mapping of the human amygdala, hippocampal region and entorhinal cortex: intersubject variability and probability maps. Anatomy and Embryology, 210(5–6), 343–352. https://doi.org/10.1007/s00429-005-0025-5 DOI: 10.1007/s00429-005-0025-5'}; ...
    {'Kedo, O., Zilles, K., Palomero-Gallagher, N., Schleicher, A., Mohlberg, H., Bludau, S., & Amunts, K. (2017). Receptor-driven, multimodal mapping of the human amygdala. Brain Structure and Function. https://doi.org/10.1007/s00429-017-1577-x DOI: 10.1007/s00429-017-1577-x'}; ...
    {'Amunts, K., Mohlberg, H., Bludau, S., & Zilles, K. (2020). Julich-Brain: A 3D probabilistic atlas of the human brain s cytoarchitecture. Science, 369(6506), 988–992. https://doi.org/10.1126/science.abb4588 DOI: 10.1126/science.abb4588'}; ...
    {'Amunts et al (2023) [Dataset v3.0.3] DOI:10.25493/56EM-75H'}]);

% renormalize the 1 or 2vx p>1
total_p = sum(hipp.probability_maps,2);
renorm = total_p > 1;
hipp.probability_maps(renorm,:) = hipp.probability_maps(renorm,:)./total_p(renorm);



amyg = load_atlas(sprintf('cit168_amygdala_%s',ALIAS));
amyg = lateralize(amyg);
amyg.labels_4 = amyg.labels_3;
amyg.labels_3 = amyg.labels_2;
amyg.labels_2 = amyg.labels;
amyg.labels_5 = repmat({'CIT168 amydala v1.0.3'},1,num_regions(amyg));

hipp = hipp.merge_atlases(amyg);

% renormalize ~186vx voxels shared between hipp and amygdala parcellations
total_p = sum(hipp.probability_maps,2);
renorm = total_p > 1;
hipp.probability_maps(renorm,:) = hipp.probability_maps(renorm,:)./total_p(renorm);

% get cifti hippocampal mask and for each voxel find its nearest
% hippocampal structure from hipp.
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,{'hipp','amygdala'}))));
cifti_mask = cifti_mask.replace_empty();
hipp_amyg_dil = dilate(hipp, cifti_mask);

hipp_amyg_dil = hipp_amyg_dil.apply_mask(cifti_mask);
clear julitch hipp amyg

% add area prefix
for i = 1:length(hipp_amyg_dil.labels)
    hipp_amyg_dil.labels{i} = [ 'MTL_' hipp_amyg_dil.labels{i}]; 
end

%% thalamus and brainstem
% returns atlas object in variable called thalamus_atlas
create_thalamus2023_atlas
% make LGN/MGN label_4 into metathalamus instead of thalamus

% returns atlas object in variable called bstem_atlas 
create_brainstem2023_atlas_unrestricted
%}
% Add labels to make more consistent with other atlases

for i = 1:length(bstem_atlas.labels)
    bstem_atlas.labels{i} = [ 'Bstem_' bstem_atlas.labels{i}]; 
end

thal_bstem = thalamus_atlas.merge_atlases(bstem_atlas,'noreplace');

% add prefixes for the more esoteric structures
thal_bstem.labels{contains(thal_bstem.labels,'Mamm_Nuc_L')} = 'Hythal_Mamm_Nuc_L';
thal_bstem.labels{contains(thal_bstem.labels,'Mamm_Nuc_R')} = 'Hythal_Mamm_Nuc_R';

thal_bstem.labels{contains(thal_bstem.labels,'Haben_L')} = 'Epithal_Haben_L';
thal_bstem.labels{contains(thal_bstem.labels,'Haben_R')} = 'Epithal_Haben_R';

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,{'thalamus','brain_stem','diencephalon'}))));

thal_bstem = thal_bstem.apply_mask(cifti_mask);

%% combinate atlases

atlas_obj = hipp_amyg_dil.merge_atlases(cerebellum_dil).merge_atlases(bg_dil).merge_atlases(thal_bstem,'noreplace');

% renorm
total_p = sum(atlas_obj.probability_maps,2);
renorm = total_p > 1;
atlas_obj.probability_maps(renorm,:) = atlas_obj.probability_maps(renorm,:)./total_p(renorm);


atlas_obj.probability_maps = sparse(atlas_obj.probability_maps);

clear thal_bstem hipp_amyg_dil bg_dil cerebellum_dil

%% create a full brain NIFTI atlas
% combine subctx and canlab20187 cortical parcels, but resort the cortical
% parcels to match the surface based Glasser parcels we intend to combine
% with the subctx volumes in our CIFTI volume.
% this is somewhat redundant with canlab2018 but there are a couple of
% differences:
% - subcortical structures are bilateral.
% - canlab2018 structures missing from canlab2023: Cau_L, Cau_R, BST_SLEA,
%   Cblm_Interposed_R, Cblm_Fastigial_L, Cblm_Fastigial_R

glasser = load_atlas(sprintf('glasser_%s',ALIAS)).apply_mask(fmri_mask_image(atlas_obj), 'invert');
% get rid of regions dismembered by cifti_atlas masking
regions_to_fix = find(contains(glasser.labels, {'PHA','EC','PreS'}));
resid_regions = find(~ismember(1:num_regions(glasser), regions_to_fix));

fixed_regions = glasser.select_atlas_subset(regions_to_fix).threshold(0.05,'k',5,'remove_parcel_fragments');
glasser = glasser.select_atlas_subset(resid_regions).merge_atlases(fixed_regions);
[~,reorder] = sort([resid_regions, regions_to_fix]);
glasser = glasser.reorder_atlas_regions(reorder);

glasser.labels_5 = repmat({'Glasser2016+Petre2023_vol_proj'},1,num_regions(glasser));

% remove hippocampal ROI because it's redundant with volumes
canlab_mask = fmri_mask_image(glasser);
ctx_regions_less_hipp = find(~ismember(1:num_regions(glasser),find(contains(glasser.labels,'_H'))));
glasser = glasser.select_atlas_subset(ctx_regions_less_hipp);
glasser.probability_maps = single(glasser.probability_maps);
%delete(gcp('nocreate')); parpool(2); % this is memory intensive, but a short process, so let's limit number of threads to avoid out of memory errors
%glasser = dilate(glasser, canlab_mask);
%delete(gcp('nocreate'))

glasser.labels_5 = repmat({'Glasser2016+Petre2023VolProj'},1,num_regions(glasser));
glasser.labels_4 = glasser.labels_3;
glasser.labels_3 = glasser.labels_2;
glasser.labels_2 = glasser.labels;

%% add macro structural glasser labels

glasser = glasser.resample_space(ref);
glasser_L = glasser.select_atlas_subset(1:num_regions(glasser)/2);
glasser_R = glasser.select_atlas_subset(num_regions(glasser)/2+1:num_regions(glasser));
clear glasser
% make sure left and right ctx labels are in the same order as in the
% Glasser atlas

% note these txt files are generated by
% script_2023_wagerlab_combined_atlas_prep.sh
fid = fopen('lctx_labels.txt');
lbls_L = textscan(fid,'%s');
lbls_L = lbls_L{1}(1:6:end);
ctx_regions_less_hipp = ~ismember(1:length(lbls_L),find(contains(lbls_L,'_H')));
lbls_L = lbls_L(ctx_regions_less_hipp);
fclose(fid);

fid = fopen('rctx_labels.txt');
lbls_R = textscan(fid,'%s');
lbls_R = lbls_R{1}(1:6:end);
ctx_regions_less_hipp = ~ismember(1:length(lbls_R),find(contains(lbls_R,'_H')));
lbls_R = lbls_R(ctx_regions_less_hipp);
fclose(fid);

labels = glasser_L.labels;
is_equal = [];
for i = 1:length(labels)
    is_equal(end+1) = strcmp(['L_', labels{i}(5:end-2), '_ROI'], strrep(lbls_L{i},'-','_'));
end
assert(all(is_equal))

labels = glasser_R.labels;
is_equal = [];
for i = 1:length(labels)
    is_equal(end+1) = strcmp(['R_', labels{i}(5:end-2), '_ROI'], strrep(lbls_R{i},'-','_'));
end
assert(all(is_equal))

canlab = glasser_L.merge_atlases(glasser_R).merge_atlases(atlas_obj);
clear glasser_L glasser_R

canlab.references = unique(canlab.references,'rows');

% reformat laterality labels
for fnames = {'labels','labels_2','labels_3','labels_4'}
    this_lbl = fnames{1};
    canlab.(this_lbl) = cellfun(@(x1)(regexprep(x1,'_([LR])_(.*)','_$2_$1')), canlab.(this_lbl), 'UniformOutput', false);
    canlab.(this_lbl) = cellfun(@(x1)(regexprep(x1,'_rh$','_R')), canlab.(this_lbl), 'UniformOutput', false);
    canlab.(this_lbl) = cellfun(@(x1)(regexprep(x1,'_lh$','_L')), canlab.(this_lbl), 'UniformOutput', false);
    canlab.(this_lbl) = cellfun(@(x1)(regexprep(x1,'^([LR])_(.*)','$2_$1')), canlab.(this_lbl), 'UniformOutput', false);
end

atlas_name = sprintf('CANLab2023_%s', SPACE);
canlab.atlas_name = atlas_name;

canlab.probability_maps = sparse(canlab.probability_maps);
save(sprintf('%s_scaffold.mat',atlas_name), 'canlab');  

% we need the memory for this
%clear all;
%load(sprintf('%s_scaffold.mat',atlas_name))

nii = fmri_data(canlab);
nii.dat = single(full(canlab.probability_maps));
nii.fullpath = sprintf('%s_scaffold.nii', canlab.atlas_name);
nii.write('overwrite');
gzip(nii.fullpath);

labels = table(canlab.labels', canlab.label_descriptions, canlab.labels_2', canlab.labels_3', canlab.labels_4', canlab.labels_5', ...
    'VariableNames', {'labels', 'label_descriptions', 'labels_2', 'labels_3', 'labels_4', 'labels_5'});
writetable(labels,[canlab.atlas_name, '_labels.csv']);
writetable(table(canlab.references,'VariableNames',{'references'}),[canlab.atlas_name, '_ref.txt']);
