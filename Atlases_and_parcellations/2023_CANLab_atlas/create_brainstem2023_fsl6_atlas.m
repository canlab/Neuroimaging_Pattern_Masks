% create_brainstem_atlas_group
% 
% Create a brain atlas with anatomically defined region groups, from
% various atlases/papers.  Uses canlab_load_ROI
%
% Notes: functional atlases probably do not [yet] have very good
% subdivisions, and there is a clear demarcation of functions, inputs, and
% outputs by anatomical subnuclei.

addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'));
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'));

dosave = true;
space = 'fsl6';
scale = 'coarse';

template = 'fmriprep20_template.nii.gz';
if strcmp(space,'fsl6')
    template = 'fsl6_hcp_template.nii.gz';
end

% for Diedrichsen mask, which is publically available. I don't know if
% there are restrictions on the rest of the repo, but the mask at least
% could be moved into a public repo and is only in MasksPrivate for
% organization reasons (to keep it with the rest of the Diedrichsen atlas
% which was already there).
addpath(genpath('/home/bogdan/.matlab/canlab/MasksPrivate'));

% Define: 1 mm space by default. Depends on diedrichsen
% brainstem-cerebellum mask in MasksPrivate
bstem_mask = atlas(which('canlab_brainstem_dil.nii.gz'));
bstem_mask = threshold(bstem_mask, .2);

bstemimg = fmri_mask_image(which('tpl-MNI152NLin6AsymC_desc-pcereb_to_MNI152NLin2009cAsym.nii.gz'));
bstem_mask = apply_mask(bstem_mask, bstemimg);

% this has some other regions in it (e.g., mammillary bodies), so would be better to use the
% 'noreplace' option when merging it with other atlases.

orthviews(bstem_mask, 'overlay', which(template));

%see also: bstemimg.fullpath = fullfile(pwd, 'brainstem_mask_tight_2018.img');

%% Shen regions: Fill in parcels not assigned to a named region

shen_file = sprintf('shen_filler_%s_atlas_object.mat', space);

if isempty(which(shen_file))
    shen = load_atlas(sprintf('shen_%s',space));
    
    bstem_mask = bstem_mask.resample_space(shen);
    
    % select brainstem regions and give them descriptive labels
    % r/c - rostral/caudal
    % d/v - dorsal/ventral
    % R/L - right/left
    shen = shen.select_atlas_subset(shen.apply_mask(bstem_mask).labels);
    shenR = lateralize(shen.select_atlas_subset([10,12:16])).select_atlas_subset('_R');
    shenL = lateralize(shen.select_atlas_subset([26,28:31])).select_atlas_subset('_L');
    shen = shenR.merge_atlases(shenL);
    
    % Ponscd becomes Ponscd_R in canlab2023
    labels = {'Shen_Midb_Rrd', 'Shen_Med_R', 'Shen_Pons_R', 'Shen_Pons_Rcv', ...
        'Shen_Midb_Rd', 'Shen_Pons_Rcd', 'Shen_Midb_Lrd', 'Shen_Midb_Ld', ...
        'Shen_Med_L', 'Shen_Pons_Lcd', 'Shen_Pons_Lcv'};
    
    shen.labels_3 = shen.labels_2;
    shen.labels_2 = shen.labels;
    shen.labels = labels;
    
    shen = dilate(shen, bstem_mask);
    
    % to find clusters by hand:
    % [~,wh] = find_closest_cluster(r, spm_orthviews('pos')) 
    
    % Reorder Shen regions
    
    [~, left] = select_atlas_subset(shen, {'_L'});
    [~, right] = select_atlas_subset(shen, {'_R'});
    
    [~, wh_midb] = select_atlas_subset(shen, {'Midb'});
    wh_midb = [find([wh_midb & left]) find([wh_midb & right])];
    
    [~, wh_pons] = select_atlas_subset(shen, {'Pons'});
    wh_pons = [find([wh_pons & left]) find([wh_pons & right])];
    
    [~, wh_med] = select_atlas_subset(shen, {'Med'});
    wh_med = [find([wh_med & left]) find([wh_med & right])];
    
    
    wh_order = [wh_midb wh_pons wh_med];
    
    shen = reorder_atlas_regions(shen, wh_order);
    
    
    %% Add brainstem nuclei, replacing old ones
    
    bstem_atlas = shen.remove_empty;
    
    % add probability map to privilege other regions over vague functionally
    % defined parcels
    uniq_roi = unique(bstem_atlas.remove_empty.dat);
    bstem_atlas = bstem_atlas.replace_empty;
    pmap = zeros(size(bstem_atlas.dat,1), length(uniq_roi));
    for i = 1:length(uniq_roi)
        pmap(bstem_atlas.dat == i, i) = 0.35;
    end
    bstem_atlas.probability_maps = sparse(double(pmap));
    
    bstem_atlas.labels_5 = repmat({'Shen 268 parcellation (dilated)'}, 1, num_regions(bstem_atlas));
    
    save(shen_file, 'bstem_atlas');
    shen_references = bstem_atlas.references;
else
    load(which(shen_file), 'bstem_atlas');
    shen_references = bstem_atlas.references;
end

switch space
    case 'fsl6'
        cifti_mask = fmri_mask_image('hcp_cifti_subctx_labels.nii');
    case 'fmriprep20'
        cifti_mask = fmri_mask_image('hcp_cifti_subctx_labels_MNI152NLin2009cAsym.nii');
    otherwise 
        error('No cifti atlas available in space %s', space);
end
bstem_atlas = bstem_atlas.apply_mask(cifti_mask);

bstem_atlas.labels_2 = repmat({'Brainstem'},1,num_regions(bstem_atlas));
diencephalic_ind = find(contains(bstem_atlas.labels, {'Midb_Lrd','Midb_Rrd'}));
bstem_atlas.labels_2(diencephalic_ind) = repmat({'Diencephalic'},1,length(diencephalic_ind));

%% add other regions

biancia = load_atlas(sprintf('bianciardi_%s_%s', scale, space));
biancia.labels_4 = repmat({'Restricted (see BrainstemNavigator v0.9 license)'}, 1, num_regions(biancia));
biancia.labels_5 = repmat({'Bianciardi brainstem navigator v.0.9'}, 1, num_regions(biancia));

% note that the locus coerulues has more rigorous segmentations based on 
% T1-TSE sequences, but they produce regions that overlap very well with 
% the bianciardi atlas' segmentation at this resolution.
% MnR replaces NCS
% VSM replaces dmnx_nts
% RMg replaces ncs_B6_B8
% NRP (nucleus raphe pontis) - dropped in CANLab2023
% nuc_ambiguus - dropped in CANLab2023
% medullary_raphe - replaced by RPa and ROb
% spinal_trigeminal - temporarily excluded. Can be added later with
%   transformation from SUIT brainstem space (transforms into MNI available
%   from Diedrichsen lab)
% CLi_RLi - new in CANLab2023
% iMRt, isRt, mRt, PCRtA, PnO_PnC, sMRt - reticular formations, new in CANLab2023
% SOC, ION - superior and inferior olives, new in CANLab2023
% MiTg_PBG
% PDT_CGPn
% PTg
% RN - red nucleus new in CANLab2023
% SubC - subcoeruleus new in CANLab2023
% Ve - vestibular nucle icomplex
% STh - Subthalamic nucelus new in CANLab2023
biancia_regions = {'PAG', 'SC' 'IC', 'DR', 'MnR', 'SN', 'LC', 'RMg', 'VSM', ...
    'CnF', 'RPa', 'ROb', 'CLi_RLi', 'Rt', 'ION', 'SOC', 'LDTg_CGPn', ...
    'MiTg_PBG', 'PnO_PnC', 'PTg', 'RN', 'SubC', 'Ve', 'STh', 'LPB', 'MPB'};
% there's a clear right side bias in the alignment
% MnR seems to be a single voxel. Look more carefully at this
% 
biancia = biancia.select_atlas_subset(biancia_regions);


cit = load_atlas(sprintf('cit168_%s', space));

cit_regions = {'PBP', 'VTA','Mamm'};

cit = lateralize(cit.select_atlas_subset(cit_regions));
cit.labels_5 = repmat({'CIT168'}, 1, num_regions(cit));

%biancia = biancia.merge_atlases(cit.select_atlas_subset(cit_regions),'noreplace');

% also include regions in other atlases that we want to remove here - so
% that we remove these voxels

% to-do: 'pbn' 'nts'

regionnames = {'rvm', 'spinal_trigeminal'};
lateralized = {'spinal_trigeminal'};

% NEW ONES TOO

% we put these into cit, since there's minimal overlap.

for i = 1:length(regionnames)
    regionname = regionnames{i};
    
    [~, roi_atlas] = canlab_load_ROI(regionname);

    % since these ROIs are poorly specified we cede ground to competition
    % anytime odds are better than not that the competition is right by
    % setting these ROIs probabilities to 50%.
    roi_atlas = roi_atlas.replace_empty;
    pmap = zeros(size(roi_atlas.dat,1),num_regions(roi_atlas));
    for j = 1:num_regions(roi_atlas)
        pmap(roi_atlas.dat == j,j) = 0.5;
    end
    roi_atlas.probability_maps = pmap;
    roi_atlas.labels_5 = repmat({'Manually drawn coordinate based ROI/Nash'}, 1, num_regions(roi_atlas));
    roi_atlas.labels_2 = repmat({'Brainstem'},1,num_regions(roi_atlas));

    orthviews(roi_atlas, 'overlay', which(template));

    cit = merge_atlases(cit, roi_atlas);
end


% the order here is important. The most sensitive atlases (i.e. those with
% the smallest regions) go first so that we resample to their spaces.
atlas_obj = biancia.merge_atlases(cit).merge_atlases(bstem_atlas);

%atlas_obj = atlas_obj.apply_mask(bstem_mask);
%atlas_obj = atlas_obj.threshold(0.2,'k',10);


%% Adjust labels
% make more consistent with other atlases
% relabel L and R

pat = 'Reg_1';
atlas_obj.labels = regexprep(atlas_obj.labels, pat, 'other');

pat = 'Shen_';
atlas_obj.labels = regexprep(atlas_obj.labels, pat, '');

atlas_obj = atlas_add_L_R_to_labels(atlas_obj);



%% Add references

references = {biancia.references, ...
    cit.references, ...
    shen_references,...
    'Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.'};

atlas_obj.references = unique(char(references),'rows');

% 'PAG'        Periaqueductal gray
% 'SC'         Superior colliculus
% 'IC'         Inferior colliculus
% 'DR'         Dorsal raphe nucleus
% 'MnR'        Median raphe nucleus (replaces NCS in CANLab2018)
% 'PBP'        Parabrachial pigmented nuc.      % Pauli 2018 CIT168 subcortical atlas
% 'SN'         Substantia Nigra (SN1 and SN2 are approximately the pars compacta and pars reticulata
% 'VTA'        Ventral tegmental area           % Pauli 2018 CIT168 subcortical atlas
% 'RN'         Red nucleus; 
% 'LPB'        Lateral Parabrachial
% 'MPB'        Medial Parabrachial
% 'LC'         Locus coeruleus
% 'SubC'       Subcoeruleus
% 'rvm_old'    Hand-drawn rostral ventral medulla (Tor) in anatomical rvm
% 'rvm'        Rostral ventral medulla from Brooks et al. 2016(??)
% 'nts'        Nuc. tractus solitarius (rough; hand-drawn, Tor)
% 'olive'      Inferior olive; MISSING
% 'RMg'        Nuc. raphe magnus (replaces previously redundant ncs_B6_B8 in canlab2018)
% 'VSM'        Visero-sensory-motor nucleu, including nucleus ambiguus and dorsal motor nucleus of the vagus (DMNX)
% 'Cnf'        Nuc. cuneiformis
% 'RPa'        Raphe pallidus (replaces medullary raphe in canlab2018)
% 'ROb'        Raphe obscurus (replaces medullary raphe in canlab2018)
% 'CLi_RLi'    Caudal-rostral linear raphe (new in CANLab2023)
% 'iMRt[l|m]'  Inferior medullary reticular formation (lateral|medial) (new in CANLab2023)
% 'isRt'       Isthmic reticular formation (new in CANLab2023)
% 'mRt[a|d|l]' Mesemcephalic reticular formation (anterior|dorsal|lateral) (new in CANLab2023)
% 'PCRtA'      Parvicellular reticular nucleus alpha part (new in CANLab2023)
% 'sMRt[l|m]'  Superior medullary reticular formation (lateral|medial) (new in CANLab2023)
% 'PnO_PnC'    Pontine reticular nucleus, oral and caudal parts (pontis oralis and caudalis) (new in CANLab2023)
% 'ION'        Inferior olivary nucleus (new in CANLab2023)
% 'SOC'        Superior olivary complex (new in CANLab2023)
% 'LDTg_CGPn'  Lateraodorsal tegemental nucleus - central gray of the rhomboencephalon (new in CANLab2023)
% 'MiTg_PBG'   Microcellular tegmental nucleus - prabigeminal nucleus
% 'PTg'        Pedunculotegmental nucleus (aka pedunculopontine nucleus)
% 'Ve'         Vestibular nuclei complex (new in CANLab2023)
% 'STh'        Subthalamic nucleus (new in CANLab2023)
% 'spinal_trigeminal' and 'RVM',
%              Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.

%% REMOVE general region - not needed

try
    atlas_obj = remove_atlas_region(atlas_obj, {'other'});
end

%% Save dir

savedir = pwd;

cd(savedir)

%% save object

atlas_name = sprintf('brainstem2023_%s_%s_combined', scale, space);

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save(savename, 'atlas_obj');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

r = atlas2region(atlas_obj);
labels = atlas_obj.labels;

for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    cellfun(@(x1)set(x1,'FaceAlpha', .5), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off

    savedir = fullfile(pwd, 'png_images');
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end

%% save figure

if dosave
    o2 = canlab_results_fmridisplay([], 'full2', 'overlay', which(template));
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

end
 