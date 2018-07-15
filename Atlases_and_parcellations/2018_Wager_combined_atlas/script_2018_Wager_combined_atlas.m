cd('/Users/torwager/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2018_Wager_combined_atlas')

savedir = '/Users/torwager/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2018_Wager_combined_atlas';

% Cortex from Glasser
% Striatum, STN, VTA, pallidum, BST from Pauli 2018  [alt striatum: Brainnetome] [alt: Keuken RN SN STN Pallidum]
% Cerebellum from Diedrichsen
% Striatal subregions from Pauli 2016
% Thalamus from Morel  [alt: Brainnetome]
% Basal forebrain, hipp/MTL, amy from SPM Anatomy 
% Brainstem from Shen
% Specific brainstem nuclei from Keren, Brooks, Tor SPM5 PAG [redo], 

% SPM Anatomy
% Septum/DBB, BNM, Astr?, 
% CA1, CA2, CA3, Subiculum, DG
% BL, SF, CM, HATA

% Glasser: 180 regions
% CIT168 : 'NAC' 'BST_SLEA' 'GPe' 'GPi' 'SNc' 'RN' 'SNr' 'PBP' 'VTA' 'VeP' 'Haben' 'Hythal' 'Mamm_Nuc' 'STN'


% Brainstem mask: 'canlab_brainstem.img'.  CIT168 brainstem masked with
% SPM8 canonical brainstem segmentation, then manually edited/cleaned.

% Atlases done:
% Keuken_create_atlas_object
% CIT168_create_atlas_object
% GlasserHCP_create_atlas_object
% SUIT_create_atlas_object
% Brainnetome_create_atlas_object
% SPManatomy_create_atlas_object

glasserfile = which('Glasser2016HCP_atlas_object.mat');
citfile = which('CIT168_MNI_subcortical_atlas_object.mat');
suitfile = which('SUIT_Cerebellum_MNI_atlas_object.mat');
shenfile = which('Shen_atlas_object.mat');
morelfile = which('Morel_thalamus_atlas_object.mat');
spmanatfile = which('SPMAnatomy22c_atlas_object.mat');
paulifile = which('Pauli2016_striatum_atlas_object.mat');
bnfile = which('Brainnetome_atlas_object.mat');

%% Start with Glasser cortex

glasser = load(glasserfile, 'atlas_obj'); glasser = glasser.atlas_obj;


%% Map in CIT168 regions (BG, some subcortex)

% ** may need 'split into contiguous' function **

cit = load(citfile); cit = cit.atlas_obj;
wh_regions = cit.labels(3:end);

%subregions = select_atlas_subset(cit, wh_regions); orthviews(subregions); n = num_regions(subregions)

atlas_obj = merge_atlases(glasser, cit, wh_regions);

atlas_obj.atlas_name = 'CANlab_2018_combined';

clear cit glasser

%% Map in SUIT MNI regions for cerebellum

suit = load(suitfile); suit = suit.atlas_obj;

%wh_regions = cit.labels(3:end);
%subregions = select_atlas_subset(cit, wh_regions); orthviews(subregions); n = num_regions(subregions)

atlas_obj = merge_atlases(atlas_obj, suit);

clear suit

%% Map in Morel combined atlas regions for thalamus

%savefile = which('Thalamus_atlas_combined_Morel.mat');

%load(savefile, 'thalamus_atlas');
% thalamus_atlas = load(morelfile);
% thalamus_atlas = thalamus_atlas.thalamus_atlas;

thalamus_atlas = load(which('Thalamus_combined_atlas_object.mat'));
thalamus_atlas = thalamus_atlas.thalamus_atlas;

% merge, using 'noreplace' to avoid overwriting habenula region already in.
atlas_obj = merge_atlases(atlas_obj, thalamus_atlas, 'noreplace');
%% Map in brainstem regions from canlab_load_ROI
% Combined brainstem atlas

% Load thalamic regions from Brainnetome
% better maybe for functional divisions of anterior nuc.

% bn = load(bnfile); 
% bn = bn.atlas_obj;
% bnthal = select_atlas_subset(bn, {'Tha'});

bstem = load(which('brainstem_combined_atlas_object.mat'));
bstem = bstem.atlas_obj;

% merge, using 'noreplace' to avoid overwriting 
atlas_obj = merge_atlases(atlas_obj, bstem, 'noreplace');

%% Map in Pauli regions for striatum

pauli = load(paulifile); pauli = pauli.atlas_obj;

% Threshold some to clean up and avoid bleed-over
pauli = threshold(pauli, .6);

%wh_regions = cit.labels(3:end);
%subregions = select_atlas_subset(cit, wh_regions); orthviews(subregions); n = num_regions(subregions)

atlas_obj = merge_atlases(atlas_obj, pauli, 'noreplace');

clear pauli


%% Save temporary

save(fullfile(savedir, 'CANlab_combined_atlas_object_2018.mat'), 'atlas_obj');



%% SPM Anatomy

% ****labels are wrong

spmanat = load(spmanatfile); spmanat = spmanat.atlas_obj;

amy_plus = select_atlas_subset(spmanat, {'Amy' 'Hipp' 'Subic'});

orthviews(amy_plus);

atlas_obj = merge_atlases(atlas_obj, amy_plus, 'noreplace');

%orthviews(atlas_obj)

%% SAVE combined atlas

save(fullfile(savedir, 'CANlab_combined_atlas_object_2018.mat'), 'atlas_obj');

%% Downsample for space efficiency and re-save

mask = fmri_data(which('brainmask.nii'));

atlas_obj = resample_space(atlas_obj, mask);

save(fullfile(savedir, 'CANlab_combined_atlas_object_2018_2mm.mat'), 'atlas_obj');

%% Add brainnetome without replacement
% Glasser is too sparse in some areas of cortex.  Fill in missing areas
% with brainnetome.
% NOT DONE - NOT SURE WE SHOULD DO THIS
% 
% bn = load(bnfile); bn = bn.atlas_obj;
% 
% %wh_regions = cit.labels(3:end);
% %subregions = select_atlas_subset(cit, wh_regions); orthviews(subregions); n = num_regions(subregions)
% 
% atlas_obj = merge_atlases(atlas_obj, bn, 'noreplace');


%% Map in Morel thalamus regions, without replacing existing
% NOT NEEDED - COMBINED THAL ATLAS NOW
% morel = load(morelfile); morel = morel.atlas_obj;
% 
% %wh_regions = cit.labels(3:end);
% %subregions = select_atlas_subset(cit, wh_regions); orthviews(subregions); n = num_regions(subregions)
% 
% atlas_obj = merge_atlases(atlas_obj, morel, 'noreplace');
% 
% % *** seems to break after adding thal in orthviews
% 
% % clear morel
% 
% orthviews(atlas_obj)


