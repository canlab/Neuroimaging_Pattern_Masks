% this script rereates the canlab2018 atlas except it constrains certain
% structures to more accurately abide by CIFTI volumetric boundaries. For
% example, CA1-3, DG and Subiculum are restricted to the hippocampal
% volume, and other structures are excluded from it, even though the fornix 
% approaches the thalamus, and the=% entorhinal cortex potentially bleeds 
% into this volume in MNI space.

% ToDo:
% Deal with overlapping Glasser vs. volumetric hippocampal segmentations.
%  The hippocampus and presubiculum surface ROIs overlap the hippocampal 
%  volume. Glasser calls these 
%    {'Ctx_PreS_L'  }
%    {'Ctx_H_L'     }
%  Question: are grayordinates redundant, or are our volumetric probability
%  maps confusing the situation?
% If redundant: Merge ROI indices if possible. The same area should not
%  have two different labels
% If it's a probability map issue set MTL volume probabilities outside of
%   grayordinate space to zero. Cede it all to Glasser.

clear all; close all;

%LIB = '/dartfs-hpc/rc/home/m/f0042vm/software';
LIB = '/home/bogdan/.matlab';
ROOT = [LIB, '/canlab/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2023_CANLab_atlas/'];

%addpath([LIB, '/spm12']);
addpath([LIB, '/spm/spm12']);

addpath(genpath([LIB, '/canlab/CanlabCore']))
addpath(genpath([LIB, '/canlab/Neuroimaging_Pattern_Masks']))
addpath(genpath([LIB, '/canlab/MasksPrivate']))

if isempty(gcp('nocreate'))
    parpool(12);
end

SPACE = 'MNI152NLin2009cAsym';
SCALE = 'fine';

switch SPACE
    case 'MNI152NLin2009cAsym'
        ALIAS = 'fmriprep20';
        OVERLAY = which('fmriprep20_template.nii.gz');
        TEMPLATE = which('MNI152NLin2009cAsym_T1_1mm.nii.gz');


        labels = fmri_data(sprintf('%s/hcp_cifti_subctx_labels_%s.nii.gz',ROOT,SPACE));
        labels_txt = textscan(fopen([ROOT, 'hcp_cifti_subctx_labels.txt']),'%s');
    case 'MNI152NLin6Asym'
        ALIAS = 'fsl6';
        OVERLAY = which('fsl6_hcp_template.nii.gz');
        TEMPLATE = which('MNI152NLin6Asym_T1_1mm.nii.gz');
        
        labels = fmri_data([ROOT, 'hcp_cifti_subctx_labels.nii']);
        labels_txt = textscan(fopen([ROOT, 'hcp_cifti_subctx_labels.txt']),'%s');
end


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
cerebellum.labels_5 = repmat({'SUIT/Diedrichsen'},1,num_regions(cerebellum));

% remove white matter structures that aren't well represented in
% grayordinate space
cerebellum = cerebellum.select_atlas_subset(find(~ismember(cerebellum.labels,...
    {'Cblm_Dentate_L', 'Cblm_Dentate_R', 'Cblm_Interposed_L', 'Cblm_Interposed_R', 'Cblm_Fastigial_L', 'Cblm_Fastigial_R'})));

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'cerebellum'))));

tic
cerebellum_dil = dilate(cerebellum, cifti_mask);
toc


%% Hippocampus and Hippocampus
julich = load_atlas(sprintf('julich_%s',ALIAS)).resample_space(ref);

% SF - superficial amygdala
% LB - laterobasal amygdala
% CM - centromedian amygdala
% CA - cornu ammonis
% DG - dentate gyrus
% HATA - hippocampal-amygdaloid transition area
hipp_amyg = julich.select_atlas_subset({'Subic', 'HATA', 'CA', 'DG', 'Transsubic','L_SF','R_SF','LB','CM'}).replace_empty();
hipp_amyg.references = char([{'Amunts, K., Kedo, O., Kindler, M., Pieperhoff, P., Mohlberg, H., Shah, N. J., Habel, U., Schneider, F., & Zilles, K. (2005). Cytoarchitectonic mapping of the human amygdala, hippocampal region and entorhinal cortex: intersubject variability and probability maps. Anatomy and Embryology, 210(5–6), 343–352. https://doi.org/10.1007/s00429-005-0025-5 DOI: 10.1007/s00429-005-0025-5'}; ...
    {'Kedo, O., Zilles, K., Palomero-Gallagher, N., Schleicher, A., Mohlberg, H., Bludau, S., & Amunts, K. (2017). Receptor-driven, multimodal mapping of the human amygdala. Brain Structure and Function. https://doi.org/10.1007/s00429-017-1577-x DOI: 10.1007/s00429-017-1577-x'}; ...
    {'Amunts, K., Mohlberg, H., Bludau, S., & Zilles, K. (2020). Julich-Brain: A 3D probabilistic atlas of the human brain s cytoarchitecture. Science, 369(6506), 988–992. https://doi.org/10.1126/science.abb4588 DOI: 10.1126/science.abb4588'}; ...
    {'Amunts et al (2023) [Dataset v3.0.3] DOI:10.25493/56EM-75H'}]);
hipp_amyg.labels_5 = repmat({'Julich/Amunts'},1,num_regions(hipp_amyg));

% get cifti hippocampal mask and for each voxel find its nearest
% hippocampal structure from hipp.
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,{'hipp','amygdala'}))));
cifti_mask = cifti_mask.replace_empty();
hipp_amyg_dil = dilate(hipp_amyg, cifti_mask);

% note that {'VTM','MF','IF'} are listed as fiber bundles, so we want to
% exclude them from the dilution above
% VTM - Tractus ventromedialis (Stria terminalis) of Brockhaus
% MF - medial fiber bundles
% IF - intermediate fiber bundles
amygdalar_fiber_structures = julich.select_atlas_subset({'VTM','MF','IF'});
amygdalar_fiber_structures.labels_5 = repmat({'Julich/Kedo'},1,num_regions(amygdalar_fiber_structures));
hipp_amyg_dil = hipp_amyg_dil.merge_atlases(amygdalar_fiber_structures, 'always_replace').apply_mask(cifti_mask);


%% create other areas

thalamusfile = which(sprintf('thalamus2023_%s_combined_atlas_object.mat',SPACE));
brainstemfile = which(sprintf('brainstem2023_%s_%s_combined_atlas_object.mat', SCALE, SPACE));

% These are used indirectly and not loaded:
shenfile = which(sprintf('Shen_%s_atlas_object.mat',ALIAS));
if strcmp(ALIAS, 'MNI152NLin6Asym')
    morelfile = which('Morel_thalamus_atlas_object.mat');
else
    morelfile = which(sprintf('Morel_thalamus_%s_atlas_object.mat',ALIAS));
end

%% thalamus and brainstem
% this block's dilation takes abount 40 minutes to run
% Map in Morel atlas-based combined atlas regions for thalamus

% we lateralize some areas and exclude others, like Hb, because we already
% have habenula from the CIT atlas
thal_atlas = load(thalamusfile).thalamus_atlas;
lat_thal = lateralize(thal_atlas.select_atlas_subset(find(~ismember(thal_atlas.labels,{'Midline','Hb','Hythal'}))));
bilat_thal = thal_atlas.select_atlas_subset(find(ismember(thal_atlas.labels,'Midline')));
thal_atlas = lat_thal.merge_atlases(bilat_thal);

% Add labels to make more consistent with other atlases
% create probability maps too and default to 50% so that any conflicts with
% non-thalamic regions that are better than not to be correct are ceded to 
% said regions
pmap = zeros(size(thal_atlas.dat,1), num_regions(thal_atlas));
for i = 1:num_regions(thal_atlas)
    thal_atlas.labels{i} = [ 'Thal_' thal_atlas.labels{i}]; 
    pmap(thal_atlas.dat == i,i) = 0.5;
end
thal_atlas.probability_maps = pmap;
thal_atlas.labels_2 = {};
thal_atlas.labels_3 = {};
thal_atlas.labels_4 = repmat({'Restricted (contact the secretariat of the Computer Vision Lab, ETH Zürich to confirm permission for use/distribution)'}, 1, num_regions(thal_atlas));
thal_atlas.labels_5 = repmat({'Morel'}, 1, num_regions(thal_atlas));


bstem = load(brainstemfile);
bstem = bstem.atlas_obj.resample_space(ref);

% Add labels to make more consistent with other atlases
for i = 1:length(bstem.labels)
    bstem.labels{i} = [ 'Bstem_' bstem.labels{i}]; 
end

cit = load_atlas(citfile);
cit = lateralize(cit.select_atlas_subset({'Haben','Hythal'}));
cit.labels_5 = repmat({'CIT168'}, 1, num_regions(cit));

thal_bstem = bstem.merge_atlases(cit).merge_atlases(thal_atlas);

% these have had Thal_ or bstem_ prefixed onto them. Let's get rid of that
thal_bstem.labels{contains(thal_bstem.labels,'Haben_L')} = 'Haben_L';
thal_bstem.labels{contains(thal_bstem.labels,'Haben_R')} = 'Haben_R';
thal_bstem.labels{contains(thal_bstem.labels,'Mamm_Nuc_L')} = 'Mamm_Nuc_L';
thal_bstem.labels{contains(thal_bstem.labels,'Mamm_Nuc_R')} = 'Mamm_Nuc_R';
thal_bstem.labels{contains(thal_bstem.labels,'Hythal_L')} = 'Hythal_L';
thal_bstem.labels{contains(thal_bstem.labels,'Hythal_R')} = 'Hythal_R';


%% combinate atlases

atlas_obj = thal_bstem.merge_atlases(cerebellum_dil).merge_atlases(bg_dil).merge_atlases(hipp_amyg_dil);
atlas_obj.probability_maps = sparse(atlas_obj.probability_maps);

%% create label file
% we only use the subcortical volumes for CIFTI file creation. We handle
% NIFTI at the full brain level below. Meanwhile we only need 2mm data for
% CIFTI, so we make a downsampled copy here and save that.

cmap = round(255*colormap('lines'));
cmap = [cmap;cmap;cmap;cmap];

atlas_obj_ds = atlas_obj.resample_space(cifti_atlas);

atlas_obj_ds.fullpath = sprintf('%s/subctx_atlas_%s_%s.nii', ROOT, SCALE, SPACE);
atlas_obj_ds.write('overwrite')
gzip(atlas_obj_ds.fullpath);
delete(atlas_obj_ds.fullpath);

n_roi = length(atlas_obj_ds.labels);

fid = fopen(sprintf('%s/subctx_atlas_%s_%s.txt', ROOT, SCALE, SPACE),'w+');
[right_ind, left_ind] = deal(0);
bilateral_ind = size(cmap,1);
for i = 1:n_roi
    fprintf(fid, [atlas_obj_ds.labels{i}, '\n']);
    % keep colors on left and right sides in sync. Color bilateral areas
    % independently
    if strcmp(atlas_obj_ds.labels{i}(end-1:end),'_R')
        right_ind = right_ind + 1;
        fprintf(fid, [int2str(i), ' ' num2str(cmap(right_ind,1)), ' ', num2str(cmap(right_ind,2)), ' ', num2str(cmap(right_ind,3)), ' 255\n']);
    elseif strcmp(atlas_obj_ds.labels{i}(end-1:end),'_L')
        left_ind = left_ind + 1;
        fprintf(fid, [int2str(i), ' ' num2str(cmap(left_ind,1)), ' ', num2str(cmap(left_ind,2)), ' ', num2str(cmap(left_ind,3)), ' 255\n']);
    else
        fprintf(fid, [int2str(i), ' ' num2str(cmap(bilateral_ind,1)), ' ', num2str(cmap(bilateral_ind,2)), ' ', num2str(cmap(bilateral_ind,3)), ' 255\n']);
        bilateral_ind = bilateral_ind - 1;
    end
end
fclose(fid);

%% create a full brain NIFTI atlas
% combine subctx and canlab20187 cortical parcels, but resort the cortical
% parcels to match the surface based Glasser parcels we intend to combine
% with the subctx volumes in our CIFTI volume.
% this is somewhat redundant with canlab2018 but there are a couple of
% differences:
% - subcortical structures are bilateral.
% - canlab2018 structures missing from canlab2023: Cau_L, Cau_R, BST_SLEA,
%   Cblm_Interposed_R, Cblm_Fastigial_L, Cblm_Fastigial_R

glasser = load_atlas('glasser_fmriprep20').apply_mask(fmri_mask_image(atlas_obj), 'invert');
% get rid of regions dismembered by cifti_atlas masking
glasser = glasser.threshold(0.05,'remove_parcel_fragments');
glasser.labels_5 = repmat({'Glasser with registration fusion volume projection'},1,num_regions(glasser));

% remove hippocampal ROI because it's redundant with volumes
canlab_mask = fmri_mask_image(glasser);
ctx_regions_less_hipp = find(~ismember(1:num_regions(glasser),find(contains(glasser.labels,'_H'))));
glasser = glasser.select_atlas_subset(ctx_regions_less_hipp);
delete(gcp('nocreate')); parpool(2); % this is memory intensive, but a short process, so let's limit number of threads to avoid out of memory errors
glasser = dilate(glasser, canlab_mask);
delete(gcp('nocreate'))

glasser = glasser.resample_space(ref);
glasser_L = glasser.select_atlas_subset(1:num_regions(glasser)/2);
glasser_R = glasser.select_atlas_subset(num_regions(glasser)/2+1:num_regions(glasser));
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

canlab.references = unique(canlab.references,'rows');
canlab.atlas_name = sprintf('CANLab_2023_%s_%s', SCALE, SPACE);
canlab.fullpath = sprintf('%s_1mm.nii', canlab.atlas_name);
canlab.write();
gzip(sprintf('%s_1mm.nii', canlab.atlas_name))
delete(sprintf('%s_1mm.nii', canlab.atlas_name));

canlab_ds = canlab.resample_space(cifti_atlas);
canlab.probability_maps = sparse(canlab.probability_maps);
save(sprintf('%s_1mm.mat', canlab.atlas_name)','canlab');

%% make 2mm version
canlab_ds.atlas_name = sprintf('CANLab_2023_%s_%s_2mm', SCALE, SPACE);

canlab_ds.fullpath = sprintf('%s.nii', canlab_ds.atlas_name);
canlab_ds.write('overwrite');
gzip(sprintf('%s.nii', canlab_ds.atlas_name));
delete(sprintf('%s.nii', canlab_ds.atlas_name))

save(sprintf('%s.mat', canlab_ds.atlas_name)','canlab_ds');