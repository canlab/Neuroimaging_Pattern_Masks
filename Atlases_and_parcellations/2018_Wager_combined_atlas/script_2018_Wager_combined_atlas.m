
%% Save dir

savedir = what('2018_Wager_combined_atlas');
savedir = savedir.path;

cd(savedir)

%%

% NOTE: Some files/atlases are in MasksPrivate.  

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

% These are used directly:
glasserfile = which('Glasser2016HCP_atlas_object.mat');
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

% Scripts to recreate required files and sub-atlases:
% GlasserHCP_create_atlas_object
% CIT168_create_atlas_object
% pauli2016_create_atlas_object
% create_basal_ganglia_atlas
% SUIT_create_atlas_object

% Create_Morel_atlas_object.m
% create_thalamus_atlas.m

% create_brainstem_atlas

%% Start with Glasser cortex

glasser = load(glasserfile, 'atlas_obj'); glasser = glasser.atlas_obj;


%% Map in combined BG atlas

bg = load(bgfile); bg = bg.atlas_obj;

atlas_obj = merge_atlases(glasser, bg);

atlas_obj.atlas_name = 'CANlab_2018_combined';

clear bg glasser

%% Map in CIT168 regions (only those not in other atlases)

cit = load(citfile); cit = cit.atlas_obj;
wh_regions = {'BST_SLEA' 'Haben' 'Mamm_Nuc'}; % some regions already in BG atlas.  cit.labels(3:end);  'Hythal' 

atlas_obj = merge_atlases(atlas_obj, cit, wh_regions);

clear cit

%% Map in SUIT MNI regions for cerebellum

suit = load(suitfile); suit = suit.atlas_obj;

atlas_obj = merge_atlases(atlas_obj, suit);

clear suit

%% Map in Morel atlas-based combined atlas regions for thalamus

%savefile = which('Thalamus_atlas_combined_Morel.mat');

%load(savefile, 'thalamus_atlas');
% thalamus_atlas = load(morelfile);
% thalamus_atlas = thalamus_atlas.thalamus_atlas;

thalamus_atlas = load(thalamusfile);
thalamus_atlas = thalamus_atlas.thalamus_atlas;

% Add labels to make more consistent with other atlases
for i = 1:length(thalamus_atlas.labels)
    thalamus_atlas.labels{i} = [ 'Thal_' thalamus_atlas.labels{i}]; 
end

% Note: July 2018 - this fails for some reason (incorrect separation)
% needs debugging.  this work ok on other atlases... 
% thalamus_atlas = split_atlas_by_hemisphere(thalamus_atlas);

% merge, using 'noreplace' to avoid overwriting habenula region already in.
atlas_obj = merge_atlases(atlas_obj, thalamus_atlas, 'noreplace');

%% Map in brainstem regions from canlab_load_ROI
% Combined brainstem atlas

% Load thalamic regions from Brainnetome
% better maybe for functional divisions of anterior nuc.

% bn = load(bnfile); 
% bn = bn.atlas_obj;
% bnthal = select_atlas_subset(bn, {'Tha'});

bstem = load(brainstemfile);
bstem = bstem.atlas_obj;

% Add labels to make more consistent with other atlases
for i = 1:length(thalamus_atlas.labels)
    bstem.labels{i} = [ 'Bstem_' bstem.labels{i}]; 
end

% merge, using 'noreplace' to avoid overwriting 
atlas_obj = merge_atlases(atlas_obj, bstem, 'noreplace');


%% Save temporary
% 
% save(fullfile(savedir, 'CANlab_combined_atlas_object_2018.mat'), 'atlas_obj');
% 


%% SPM Anatomy for amygdala and hippocampus

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

%%  Add larger units to labels_2
% Also save and write .nii images

plugin_canlab_atlas_2018_relabel_larger_units

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

 %% save figure

 atlas_name = atlas_obj.atlas_name;
 
 r = atlas2region(atlas_obj);
 
if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1);
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end
 
%% Isosurface

if dosave
    
    figure; han = isosurface(atlas_obj);
    
    set(han,'FaceAlpha', .5)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
