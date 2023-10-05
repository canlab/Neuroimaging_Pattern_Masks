% this script rereates the canlab2018 atlas except it constrains certain
% structures to more accurately abide by CIFTI volumetric boundaries. For
% example, CA1-3, DG and Subiculum are restricted to the hippocampal
% volume, and other structures are excluded from it, even though the fornix 
% approaches the thalamus, and the% entorhinal cortex potentially bleeds 
% into this volume in MNI space.

% this script needs to be modified to respect lateralization

% There are several ROIs we'll want to drop for hyperalignment because
% they're too fragmented. We can use them for post-hoc BSC, but we don't
% want to restirct hyperalignment to their boundaries. These could be
% 'VeP' - split across CIFTI pallidum, putamen and NAc, and quite 
%   fragmented in some cases
% 'BST_SLEA' - part of the extended amgydala, dorsal to pallidum and 
%   putamen, but much of it is in white matter excluded from the
%   grayordinate space.

clear all; close all;

addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/spm12'));

addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/CanlabCore'))
addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/Neuroimaging_Pattern_Masks'))
addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/MasksPrivate'))

ROOT = '/dartfs-hpc/rc/lab/C/CANlab/labdata/projects/bogdan_hcp_hyperalignment/sandbox/2023_wagerlab_combined_atlas_cifti/';

labels = fmri_data([ROOT, 'hcp_cifti_subctx_labels.nii']);
labels_txt = textscan(fopen([ROOT, 'hcp_cifti_subctx_labels.txt']),'%s');

cifti_atlas = atlas(labels);
cifti_atlas.labels = labels_txt{1};

% these imports are copied from 
% Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2018_Wager_combined_atlas/script_2018_Wager_combined_atlas.m
bgfile = which('Basal_ganglia_combined_atlas_object.mat'); % see create_basal_ganglia_atlas.m
citfile = which('CIT168_MNI_subcortical_atlas_object.mat');
suitfile = which('SUIT_Cerebellum_MNI_atlas_object.mat');
thalamusfile = which('Thalamus_combined_atlas_object.mat');
brainstemfile = which('brainstem_combined_atlas_object.mat');
spmanatfile = which('SPMAnatomy22c_atlas_object.mat');

% These are used indirectly and not loaded:
shenfile = which('Shen_atlas_object.mat');
morelfile = which('Morel_thalamus_atlas_object.mat');
paulifile = which('Pauli2016_striatum_atlas_object.mat');
bnfile = which('Brainnetome_atlas_object.mat');

bg = load(bgfile); bg = bg.atlas_obj;
atlas_obj = load(citfile); atlas_obj = atlas_obj.atlas_obj;
wh_regions = {'BST_SLEA' 'Haben' 'Mamm_Nuc'};

bg.merge_atlases(atlas_obj.select_atlas_subset(find(contains(atlas_obj.labels, 'BST_SLEA'))));

%% Caudate
% some of V_Striatum is in the CIFTI accumben (13%), but most is in the
% CIFTI caudate mask (45.9%). This is a 70-30 split, and it's hard to
% assign it to one over the other like with entorhinal cortex below. We'll
% give it to both.
% We could also include 'Cau' here, but it's very sparse within the CIFTI
% mask regions and is more likely to cause problems down the line than not
% (e.g. selecting it as an ROI for seed based analysis or alignment results
% in a nonsensical fragmented region of caudate), so we drop it.
caud = bg.select_atlas_subset(find(contains(bg.labels,{'Caudate','V_Striatum', 'BST_SLEA'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'caudate'))));
cifti_mask = cifti_mask.replace_empty();

tic
caud_dil = dilate(caud, cifti_mask);
toc

%% Putamen
% It's not clear why Caudate_Ca overlaps with the putamen ROI, but it's
% a major hunk, so we include it here.
%put = bg.select_atlas_subset(find(contains(bg.labels,{'V_Striatum','Putamen','VeP', 'BST_SLEA','Caudate_Ca_L', 'Caudate_Ca_R'}))).remove_empty();
put = bg.select_atlas_subset(find(contains(bg.labels,{'V_Striatum','Putamen','VeP', 'BST_SLEA'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'putamen'))));
cifti_mask = cifti_mask.replace_empty();

tic
put_dil = dilate(put, cifti_mask);
toc


%% Pallidum
pal = bg.select_atlas_subset(find(contains(bg.labels,{'GP','VeP', 'BST_SLEA'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'pallidum'))));
cifti_mask = cifti_mask.replace_empty();

tic
pal_dil = dilate(pal, cifti_mask);
toc



%% Accumbens
% We're dropping VeP because it's split across multiple structures (Putamen
% and NAc) and neither has much of it.
accumbens = bg.select_atlas_subset(find(contains(bg.labels,{'V_Striatum','NAC','VeP', 'BST_SLEA'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'accumbens'))));
cifti_mask = cifti_mask.replace_empty();

tic
accumbens_dil = dilate(accumbens, cifti_mask);
toc

%% bg atlas
bg_dil = caud_dil.merge_atlases(put_dil).merge_atlases(pal_dil).merge_atlases(accumbens_dil);
% deal with redundant labels
uniq_labels = unique(bg_dil.labels);
remove = [];
for i = 1:length(uniq_labels)
    this_lbl = uniq_labels(i);
    this_ind = find(contains(bg_dil.labels, this_lbl));
    for j = 2:length(this_ind)
        bg_dil.dat(bg_dil.dat == this_ind(j)) = this_ind(1);
        remove(end+1) = this_ind(j);
    end
end
[~,~,dat] = unique(bg_dil.dat);
bg_dil.dat = dat - 1;
bg_dil.probability_maps = [];
bg_dil.labels(remove) = [];

%% Cerebellum
cerebellum = load(suitfile); cerebellum = cerebellum.atlas_obj;

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'cerebellum'))));

tic
cerebellum_dil = dilate(cerebellum, cifti_mask);
toc

%% thalamus and brainstem
% this block's dilation takes abount 40 minutes to run
% Map in Morel atlas-based combined atlas regions for thalamus
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,{'thalamus','dienc','brain_stem'}))));

% Add in smaller areas we don't want to dilute: mammilary nucleus and habenula
atlas_obj = load(citfile); atlas_obj = atlas_obj.atlas_obj;
wh_regions = {'Haben' 'Mamm_Nuc'};
atlas_obj = atlas_obj.select_atlas_subset(find(contains(atlas_obj.labels, wh_regions)));

atlas_obj = merge_atlases(atlas_obj, load(thalamusfile).thalamus_atlas, 'noreplace');

% Add labels to make more consistent with other atlases
for i = 3:length(atlas_obj.labels)
    atlas_obj.labels{i} = [ 'Thal_' atlas_obj.labels{i}]; 
end


bstem = load(brainstemfile);
bstem = bstem.atlas_obj;

% Add labels to make more consistent with other atlases
for i = 1:length(atlas_obj.labels)
    bstem.labels{i} = [ 'Bstem_' bstem.labels{i}]; 
end

% merge, using 'noreplace' to avoid overwriting 
atlas_obj = merge_atlases(atlas_obj, bstem, 'noreplace');

tic
thal_bstem_dil = dilate(atlas_obj, cifti_mask);
toc

% these have had Thal_ or bstem_ prefixed onto them. Let's get rid of that
thal_bstem_dil.labels{contains(thal_bstem_dil.labels,'Haben')} = 'Haben';
thal_bstem_dil.labels{contains(thal_bstem_dil.labels,'Mamm_Nuc')} = 'Mamm_Nuc';
thal_bstem_dil.labels{contains(thal_bstem_dil.labels,'Dorsal_raphe_DR')} = 'Dorsal_raphe_DR';
thal_bstem_dil.labels{contains(thal_bstem_dil.labels,'IC_L')} = 'IC_L';

%% Hippocampus
spmanat = load_atlas(spmanatfile);

% 12.8% of the entorhinal cortex volume is in the cifti hippocampal volume
% mask, but this should be a cortical structure and most is subsumed by the 
% glasser atlas, so we only include CA1-3, DG and Subiculum here
hipp = spmanat.select_atlas_subset(find(contains(spmanat.labels,{'Hipp', 'Subic'}))).replace_empty();

% make values lateralized
xyz = hipp.volInfo.mat*[hipp.volInfo.xyzlist'; ones(1,length(hipp.volInfo.xyzlist))];
hipp.dat(xyz(1,:)' > 0 & hipp.dat ~= 0) = hipp.dat(xyz(1,:)' > 0 & hipp.dat ~= 0) + length(hipp.labels);
hipp = hipp.remove_empty();

n_labels = length(hipp.labels);
for i = 1:length(hipp.labels)
    hipp.labels{i + n_labels} = [hipp.labels{i}, 'R'];
    hipp.labels{i} = [hipp.labels{i}, 'L'];
end

% get cifti hippocampal mask and for each voxel find its nearest
% hippocampal structure from hipp.
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'hipp'))));
cifti_mask = cifti_mask.replace_empty();

hipp_dil = dilate(hipp, cifti_mask);


%% Amygdala

% 12.8% of the entorhinal cortex volume is in the cifti hippocampal volume
% mask, but this should be a cortical structure and most is subsumed by the 
% glasser atlas, so we only include CA1-3, DG and Subiculum here
amyg = spmanat.select_atlas_subset(find(contains(spmanat.labels,{'Amy'}))).replace_empty();

% make values lateralized
xyz = amyg.volInfo.mat*[hipp.volInfo.xyzlist'; ones(1,length(amyg.volInfo.xyzlist))];
amyg.dat(xyz(1,:)' > 0 & amyg.dat ~= 0) = amyg.dat(xyz(1,:)' > 0 & amyg.dat ~= 0) + length(amyg.labels);
amyg = amyg.remove_empty();

n_labels = length(amyg.labels);
for i = 1:length(amyg.labels)
    amyg.labels{i + n_labels} = [amyg.labels{i}, 'R'];
    amyg.labels{i} = [amyg.labels{i}, 'L'];
end 

% get cifti hippocampal mask and for each voxel find its nearest
% hippocampal structure from hipp.
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'amygdala'))));
cifti_mask = cifti_mask.replace_empty();

amyg_dil = dilate(amyg, cifti_mask);

%% combinate atlases

atlas_obj = bg.select_atlas_subset(find(contains(bg.labels,'STN'))).merge_atlases(thal_bstem_dil,'noreplace').merge_atlases(cerebellum_dil).merge_atlases(bg_dil).merge_atlases(hipp_dil).merge_atlases(amyg_dil);

%% create label file

cmap = round(255*colormap('lines'));

atlas_obj.fullpath = [ROOT, 'subctx_atlas.nii'];
atlas_obj.write('overwrite')
n_roi = length(atlas_obj.labels);

fid = fopen([ROOT, 'subctx_atlas.txt'],'w+');
for i = 1:n_roi
    fprintf(fid, [atlas_obj.labels{i}, '\n']);
    fprintf(fid, [int2str(i), ' ' num2str(cmap(i,1)), ' ', num2str(cmap(i,2)), ' ', num2str(cmap(i,3)), ' 255\n']);
end
fclose(fid);

%% create a full brain NIFTI atals
% combine subctx and canlab20187 cortical parcels, but resort the cortical
% parcels to match the surface based Glasser parcels we intend to combine
% with the subctx volumes in our CIFTI volume.
% this is somewhat redundant with canlab2018 but there are a couple of
% differences:
% - subcortical structures are bilateral.
% - canlab2018 structures missing from canlab2023: Cau_L, Cau_R, BST_SLEA,
%   Cblm_Interposed_R, Cblm_Fastigial_L, Cblm_Fastigial_R

canlab = load_atlas('canlab2018_2mm');
canlab_ctx_L = canlab.select_atlas_subset(1:2:360);
canlab_ctx_R = canlab.select_atlas_subset(2:2:360);
% make sure left and right ctx labels are in the same order as in the
% Glasser atlas

% note these txt files are generated by
% script_2023_wagerlab_combined_atlas_prep.sh
fid = fopen('lctx_labels.txt');
lbls_L = textscan(fid,'%s');
lbls_L = lbls_L{1}(1:6:end);
fclose(fid);

fid = fopen('rctx_labels.txt');
lbls_R = textscan(fid,'%s');
lbls_R = lbls_R{1}(1:6:end);
fclose(fid);

labels = canlab_ctx_L.labels;
is_equal = [];
for i = 1:length(labels)
    is_equal(end+1) = strcmp(['L_', labels{i}(5:end-2), '_ROI'], strrep(lbls_L{i},'-','_'));
end
assert(all(is_equal))

labels = canlab_ctx_R.labels;
is_equal = [];
for i = 1:length(labels)
    is_equal(end+1) = strcmp(['R_', labels{i}(5:end-2), '_ROI'], strrep(lbls_R{i},'-','_'));
end
assert(all(is_equal))

reordered_canlab = canlab_ctx_L.merge_atlases(canlab_ctx_R).merge_atlases(atlas_obj);

% except for hippocampus and amygdala, which have been lateralized in this
% new atlas, all labels should be the same as in the old atlas, so let's
% just import the corresponding label descriptions.
reordered_canlab.label_descriptions = {};
reordered_canlab.labels_2 = {};
reordered_canlab.labels_3 = {};
reordered_canlab.labels_4 = {};
reordered_canlab.labels_5 = {};

for i = 1:length(reordered_canlab.labels)
    if ~contains(reordered_canlab.labels{i}, {'Hippocampus','Amygdala','Subiculum'})
        ind2018 = find(ismember(canlab.labels, reordered_canlab.labels{i}));
        for fname = {'label_descriptions','labels_2','labels_3'}
            reordered_canlab.(fname{1}){end+1} = canlab.(fname{1}){ind2018};
        end
    else
        % deal with the hippocampi and amygdalae
        ind2018 = find(ismember(canlab.labels, reordered_canlab.labels{i}(1:end-1)));
        for fname = {'label_descriptions','labels_2','labels_3'}
            reordered_canlab.(fname{1}){end+1} = canlab.(fname{1}){ind2018};
        end
    end
end

reordered_canlab.references = unique(reordered_canlab.references,'rows');
reordered_canlab.atlas_name = 'CANlab_2023_combined';
reordered_canlab.fullpath = 'canlab_2023_2mm.nii';
reordered_canlab.write();
system('gzip canlab_2023_2mm.nii')

save('canlab_2023_2mm.mat','reordered_canlab');