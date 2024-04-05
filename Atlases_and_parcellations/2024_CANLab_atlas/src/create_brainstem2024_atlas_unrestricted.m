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

% for Diedrichsen mask, which is publically available. I don't know if
% there are restrictions on the rest of the repo, but the mask at least
% could be moved into a public repo and is only in MasksPrivate for
% organization reasons (to keep it with the rest of the Diedrichsen atlas
% which was already there).
addpath(genpath('/home/bogdan/.matlab/canlab/MasksPrivate'));

% Define: 1 mm space by default. Depends on diedrichsen
% brainstem-cerebellum mask in MasksPrivate
bstem_mask = atlas(which('canlab_brainstem_dil.nii.gz'));
bstem_mask = fmri_mask_image(threshold(bstem_mask, .2).flip('mirror'));

bstemimg = fmri_mask_image(which('tpl-MNI152NLin6AsymC_desc-pcereb_to_MNI152NLin2009cAsym.nii.gz'));
bstem_mask = apply_mask(bstem_mask, bstemimg);

% this has some other regions in it (e.g., mammillary bodies), so would be better to use the
% 'noreplace' option when merging it with other atlases.

%see also: bstemimg.fullpath = fullfile(pwd, 'brainstem_mask_tight_2018.img');

%% Shen regions: Fill in parcels not assigned to a named region

shen_file = sprintf('shen_filler_%s_atlas_object.mat', ALIAS);

if isempty(which(shen_file))
    shen = load_atlas(sprintf('shen_%s',ALIAS));

    bstem_mask = bstem_mask.resample_space(shen);
    
    % select brainstem regions and give them descriptive labels
    % r/c - rostral/caudal
    % d/v - dorsal/ventral
    % R/L - right/left
    shen = shen.select_atlas_subset(shen.apply_mask(bstem_mask).labels);
    %{
    shenR = lateralize(shen.select_atlas_subset([10,12:16])).select_atlas_subset({'_R'});
    shenL = lateralize(shen.select_atlas_subset([26,28:31])).select_atlas_subset({'_L'});
    shen = shenR.merge_atlases(shenL);
    
    % Ponscd becomes Ponscd_R in canlab2023
    labels = {'Shen_Midb_Rrd', 'Shen_Med_R', 'Shen_Pons_R', 'Shen_Pons_Rcv', ...
        'Shen_Midb_Rd', 'Shen_Pons_Rcd', 'Shen_Midb_Lrd', 'Shen_Midb_Ld', ...
        'Shen_Med_L', 'Shen_Pons_Lcd', 'Shen_Pons_Lcv'};
    %}
    shenR = lateralize(shen.select_atlas_subset({'R_126_Network_4',...
        'R_129_Network_4', ...
        'R_130_Network_4',...
        'R_131_Network_4',...
        'R_132_Network_4',...
        'R_133_Network_4'})).select_atlas_subset({'_R'}).apply_mask(hemi_R);
    shenL = lateralize(shen.select_atlas_subset({'R_251_Network_4',...
        'R_262_Network_4',...
        'R_265_Network_4',...
        'R_266_Network_4'....
        'R_267_Network_4',...
        'R_268_Network_4'})).select_atlas_subset({'_L'}).apply_mask(hemi_L);
    labels = {'Shen_Midb_Rrd','Shen_Med_R','Shen_Pons_Rrv','Shen_Pons_Rcv',...
        'Shen_Midb_Rcd','Shen_Pons_Rcd','Shen_Pons_Lrd','Shen_Midb_Lrd',...
        'Shen_Midb_Lc','Shen_Med_L','Shen_Pons_Lcd','Shen_Pons_Lv'};
    label_descript = {'Midbrain right rostral dorsal',...
        'Medulla right',...
        'Pons right rostral ventral',...
        'Pons right caudal ventral',...
        'Midbrain right caudal dorsal'...
        'Pons right caudal dorsal',...
        'Pons left rostral dorsal',...
        'Midbrain left rostral dorsal',...
        'Midbrain left caudal',...
        'Medulla left',...
        'Pons left cudal dorsal',...
        'Pons left ventral'};
    labels_4 = {'Midbrain_R','Medulla_R','Pons_R','Pons_R','Midbrain_R','Pons_R','Pons_L',...
        'Midbrain_L','Midbrain_L','Medulla_L','Pons_L','Pons_L'};bstem_atlas = atlas_add_L_R_to_labels(bstem_atlas);
    shenL = dilate(shenL, bstem_mask.apply_mask(hemi_L));
    shenR = dilate(shenR, bstem_mask.apply_mask(hemi_R));
    shen = shenR.merge_atlases(shenL).apply_mask(bstem_mask);


    shen.labels = labels;
    shen.label_descriptions = label_descript;
    shen.labels_2 = labels;
    shen.labels_3 = labels;
    shen.labels_4 = labels_4;
    shen.labels_5 = repmat({'Shen2013'},1,num_regions(shen));
    
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
    
    bstem_atlas.labels_5 = repmat({'Shen268'}, 1, num_regions(bstem_atlas));
    
    save(shen_file, 'bstem_atlas');
    shen_references = bstem_atlas.references;
else
    load(which(shen_file), 'bstem_atlas');
    shen_references = bstem_atlas.references;
end

diencephalic_ind_L = find(contains(bstem_atlas.labels, {'Midb_Lrd'}));
bstem_atlas.labels_4(diencephalic_ind_L) = repmat({'Midbrain_L'},1,length(diencephalic_ind_L));
diencephalic_ind_R = find(contains(bstem_atlas.labels, {'Midb_Rrd'}));
bstem_atlas.labels_4(diencephalic_ind_R) = repmat({'Midbrain_R'},1,length(diencephalic_ind_R));

bstem_atlas = atlas_add_L_R_to_labels(bstem_atlas);

%% add CIT regions

cit = load_atlas(sprintf('cit168_%s', ALIAS));

% bianciardi atlas also has SNc, SNr, RN, STH, and even subdivides the RN
% and STH into subparcels, but it's not open source so we prefer these
% alternatives. The subdivisions of the RN and especially STh are also 
% non-contiguous, which makes using them awkward.
cit_regions = {'PBP', 'VTA', 'SNc', 'SNr', 'RN', 'STH'};

cit = lateralize(cit.select_atlas_subset(cit_regions));
cit.labels_2 = {};
cit.labels_3 = {};
cit.labels_4 = {};
% add label description and group regions
for i = 1:num_regions(cit)
    if contains(cit.labels{i},'SN')
        cit.labels_2{end+1} = regexprep(cit.labels{i},'.*_([LR])','SN_$1');
        cit.labels_3{end+1} = regexprep(cit.labels{i},'.*_([LR])','VTA_PBP_SN_$1');
    elseif contains(cit.labels{i}, {'VTA','PBP'})
        cit.labels_2{end+1} = regexprep(cit.labels{i},'.*_([LR])','VTA_PBP_$1');
        cit.labels_3{end+1} = regexprep(cit.labels{i},'.*_([LR])','VTA_PBP_SN_$1');
    else
        cit.labels_2{end+1} = cit.labels{i};
        cit.labels_3{end+1} = cit.labels{i};
    end
    cit.labels_4{end+1} = regexprep(cit.labels{i},'.*_([LR])','Midbrain_$1');
end
cit.labels_5 = repmat({'CIT168 subcortical v1.1.0'}, 1, num_regions(cit));

% the order here is important. The most sensitive atlases (i.e. those with
% the smallest regions) go first so that we resample to their spaces.
bstem_atlas = cit.merge_atlases(bstem_atlas);

%% Add Levinson-Bari regions

lb_atlas = load_atlas(sprintf('limbic_brainstem_atlas_%s', ALIAS));

% bianciardi atlas also has SNc, SNr, RN, STH, and even subdivides the RN
% and STH into subparcels, but it's not open source so we prefer these
% alternatives. The subdivisions of the RN and especially STh are also 
% non-contiguous, which makes using them awkward.
lb_regions = {'DRN','LC','NTS'};

lb_atlas = lb_atlas.select_atlas_subset(lb_regions);
lb_atlas.labels_2 = {};
lb_atlas.labels_3 = {};
lb_atlas.labels_4 = {};
% add label description and group regions
for i = 1:num_regions(lb_atlas)
    if contains(lb_atlas.labels{i},'DR')
        lb_atlas.labels{i} = 'DR_B7'; % rename for consistency with Bianciardi and CANlab2023
        lb_atlas.labels_2{end+1} = 'DR_B7';
        lb_atlas.labels_3{end+1} = 'Rostral Raphe (Serotonergic)';
        lb_atlas.labels_4{end+1} = 'Midbrain';
    elseif contains(lb_atlas.labels{i}, 'LC')
        lb_atlas.labels_2{end+1} = lb_atlas.labels{i};
        lb_atlas.labels_3{end+1} = regexprep(lb_atlas.labels{i},'.*_([LR])','LC+_$1');
        lb_atlas.labels_4{end+1} = regexprep(lb_atlas.labels{i},'.*_([LR])','Pons_$1');
    elseif contains(lb_atlas.labels{i}, 'NTS')
        lb_atlas.labels_2{end+1} = lb_atlas.labels{i};
        lb_atlas.labels_3{end+1} = regexprep(lb_atlas.labels{i},'.*_([LR])','Cranial_nuclei_$1');
        lb_atlas.labels_4{end+1} = regexprep(lb_atlas.labels{i},'.*_([LR])','Medulla_$1');
    else
        error('Unrecognized Levinson-Bari atlas region');
    end
end
lb_atlas.labels_5 = repmat({'Levinson-Bari Limbic Brainstem Atlas'}, 1, num_regions(lb_atlas));
lb_atlas = lb_atlas.threshold(0.01);

% the order here is important. The most sensitive atlases (i.e. those with
% the smallest regions) go first so that we resample to their spaces.
bstem_atlas = lb_atlas.merge_atlases(bstem_atlas);


%% Add Harvard AAN regions
% these are non-probablistic and will be replaced in the non-open atlas by
% Bianciardi regions

aan_atlas = load_atlas(sprintf('harvard_aan_%s', ALIAS));

% bianciardi atlas also has SNc, SNr, RN, STH, and even subdivides the RN
% and STH into subparcels, but it's not open source so we prefer these
% alternatives. The subdivisions of the RN and especially STh are also 
% non-contiguous, which makes using them awkward.
aan_regions = {'LDTg','PBC','PTg','PnO','mRt','MnR'};
aan_regions_to_dilate = {'LDTg','MPB_LPB','PnO'}; % dilate these a bit because they're very small and unlikely to generalize well otherwise

aan_atlas = aan_atlas.select_atlas_subset(aan_regions);
aan_atlas.labels_2 = {};
aan_atlas.labels_3 = {};
aan_atlas.labels_4 = {};
% add label description and group regions
for i = 1:num_regions(aan_atlas)
    if contains(aan_atlas.labels{i},'mRt')
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'^([LR])_mRt','isRt_$1'); % rename for consistency with Bianciardi labels
        aan_atlas.labels_2{end+1} = aan_atlas.labels{i};
        aan_atlas.labels_3{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Rostral reticular formation_$1');
        aan_atlas.labels_4{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Midbrain_$1');
    elseif contains(aan_atlas.labels{i},'MnR')
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'MnR','MnR_B6_B8'); % rename for consistency with Bianciardi labels
        aan_atlas.labels_2{end+1} = aan_atlas.labels{i};
        aan_atlas.labels_3{end+1} = 'Rostral Raphe (Serotonergic)';
        aan_atlas.labels_4{end+1} = 'Midbrain';
    elseif contains(aan_atlas.labels{i}, 'LDTg')
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'^([LR])_(.*)','$2_$1');
        aan_atlas.labels_2{end+1} = aan_atlas.labels{i};
        aan_atlas.labels_3{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Cholinergic nuclei_$1');
        aan_atlas.labels_4{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Pons_$1');
    elseif contains(aan_atlas.labels{i}, 'PBC')0.
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'^([LR    ])_(.*)','MPB_LPB_$1'); % rename for consistency with Bianciardi labels
        aan_atlas.labels_2{end+1} = aan_atlas.labels{i};
        aan_atlas.labels_3{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Parabrachial nuclei_$1');
        aan_atlas.labels_4{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Pons_$1');
    elseif contains(aan_atlas.labels{i}, 'PTg')
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'^([LR])_(.*)','$2_$1');
        aan_atlas.labels_2{end+1} = aan_atlas.labels{i};
        aan_atlas.labels_3{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Cholinergic nuclei_nuclei_$1');
        aan_atlas.labels_4{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Pons_$1');
    elseif contains(aan_atlas.labels{i}, 'PnO')
        aan_atlas.labels{i} = regexprep(aan_atlas.labels{i},'^([LR])_(.*)','$2_$1');
        [aan_atlas.labels_2{end+1}, aan_atlas.labels_3{end+1}] = deal(aan_atlas.labels{i});
        aan_atlas.labels_4{end+1} = regexprep(aan_atlas.labels{i},'.*_([LR])','Pons_$1');
    else
        error('Unrecognized Harvard AAN region');
    end
end
aan_atlas.labels_5 = repmat({'Harvard Ascending Arousal Network'}, 1, num_regions(aan_atlas));

% add dummy probabilities
max_p = 0.8;
pmap = zeros(size(aan_atlas.dat));
for i = 1:num_regions(aan_atlas)
    pmap(aan_atlas.dat == i,i) = max_p;
end
aan_atlas.probability_maps = pmap;
aan_atlas = aan_atlas.probability_maps_to_region_index;

% dilate the smallest ROIs to make them usable
aan_dilated = aan_atlas.select_atlas_subset(aan_regions_to_dilate);
for i = 1:num_regions(aan_dilated)
    this_region = aan_dilated.select_atlas_subset(i).probability_maps;
    this_region(this_region > 0) = max_p;
    aan_atlas_ind = find(strcmp(aan_atlas.labels, aan_dilated.labels{i}));
    these_p = iimg_smooth_3d(this_region, aan_dilated.volInfo, 3);

    % renormalize so that max = 0.8
    these_p = these_p*max_p/max(these_p);

    aan_atlas.probability_maps(:,aan_atlas_ind) = these_p;
end
aan_atlas = aan_atlas.probability_maps_to_region_index;

aan_atlas = aan_atlas.threshold(0.01);

% diagnostics code to compare with bianciardi after dilation
%{
emap = {{{'LDTg_L'},{'L_LDTg_CGPn'}},...
    {{'LDTg_R'},{'R_LDTg_CGPn'}},...
    {{'MPB_LPB_L'},{'L_MPB','L_LPB'}},...
    {{'MPB_LPB_R'},{'R_MPB','R_LPB'}},...
    {{'PTg_L'},{'L_PTg'}},...
    {{'PTg_R'},{'R_PTg'}},...
    {{'PnO_L'},{'L_PnO_PnC_B5'}},...
    {{'PnO_R'},{'R_PnO_PnC_B5'}},...
    {{'isRt_L'},{'L_isRt'}},...
    {{'isRt_R'},{'R_isRt'}},...
    {{'MnR_B6_B8'},{'MnR_B6_B8','PMnR_B6_B8'}}};
biancia = load_atlas(sprintf('bianciardi_%s',ALIAS));
[newAtlas1, newAtlas2] = deal({});
for i = 1:length(emap)
    newAtlas1{end+1} = aan_atlas.select_atlas_subset(emap{i}{1},'flatten');
    newAtlas2{end+1} = biancia.select_atlas_subset(emap{i}{2},'flatten');
end
aanAtlas = [newAtlas1{:}];
bianciaAtlas = [newAtlas2{:}];


thisAtlas = aanAtlas.threshold(0.2); 
thisAtlas.probability_maps = [];
thisLeadsAtlas = bianciaAtlas.threshold(0.2);
thisLeadsAtlas.probability_maps = [];
for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('nofigure','transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            leads_roi = thisLeadsAtlas.select_atlas_subset(i);
        
            if num_regions(leads_roi) == 1
                o3 = o2;
                o3.activation_maps = o2.activation_maps(i);
                o3.montage = o2.montage(i);
                leads_roi.montage(o3,'existing_figure','existing_axes', o2.montage{i}.axis_handles,'outline','color',[0,0,0]);
            end
        end
    end
    set(gcf,'Tag',orientation{1})
    drawnow()
end
%}

% the order here is important. The most sensitive atlases (i.e. those with
% the smallest regions) go first so that we resample to their spaces.
bstem_atlas = aan_atlas.merge_atlases(bstem_atlas);

%% Adjust labels
% make more consistent with other atlases
% relabel L and R

pat = 'Reg_1';
bstem_atlas.labels = regexprep(bstem_atlas.labels, pat, 'other');

pat = 'Shen_';
bstem_atlas.labels = regexprep(bstem_atlas.labels, pat, '');

%% Add PAG
% we'll start with the MNI152NLin6Asym (original) space data regardless of
% our target data and project to the target space as needed in matlab.

% note the gunzipping/gzipping here shouldn't be necessary but there's a
% strange bug happening for me right now where
% exist('kragelpag_MNI152NLin6Asym.nii') returns true even when no such
% file is present, which breaks automated gunzipping. We do it manually as
% a temp workaround.
gunzip('kragelpag_MNI152NLin6Asym.nii.gz');
kragelmasks = fmri_data('kragelpag_MNI152NLin6Asym.nii');
gzip('kragelpag_MNI152NLin6Asym.nii')
delete('kragelpag_MNI152NLin6Asym.nii')
kragelpmaps = kragelmasks.mean();

% the aqueduct masks will be used to edit the Shen filler regions
kragelmasks_aqueduct = fmri_data('kragelaqueduct_MNI152NLin6Asym.nii.gz');
kragelpmaps_aqueduct = kragelmasks_aqueduct.mean();

% we get this to label columns
kragelPAG = load_atlas('Kragel2019PAG_atlas_object.mat');
kragelPAG = kragelPAG.resample_space(kragelpmaps);

% expand columns to fill our new probablistic PAG area using nearest
% neighbor labeling
kragelPAG_dil = dilate(kragelPAG, fmri_mask_image(kragelpmaps));

% split probability map into subregions and asign labels from column map
% above
pmap = zeros(size(kragelmasks.dat,1),num_regions(kragelPAG_dil));
for i = 1:num_regions(kragelPAG_dil)
    ind = kragelPAG_dil.dat == i;
    pmap(ind,i) = kragelmasks.dat(ind);
end
kragelpmaps.dat = pmap;

if strcmp(SPACE,'MNI152NLin2009cAsym')
    kragelpmaps.fullpath = '/tmp/kragelpmaps.nii';
    kragelpmaps.write('overwrite');
    apply_spm_warp(kragelpmaps.fullpath, which(sprintf('%s_T1_1mm.nii.gz',SPACE)), ...
        which('00_fsl_to_fmriprep_subctx_AffineTransform.csv'), ...
        which('y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii'),...
        [],...
        '/tmp/kragelpmaps_MNI152NLin2009cAsym.nii', ...
        1)
    delete(kragelpmaps.fullpath);
    kragelpmaps = fmri_data('/tmp/kragelpmaps_MNI152NLin2009cAsym.nii');
    delete('/tmp/kragelpmaps_MNI152NLin2009cAsym.nii');


    kragelpmaps_aqueduct.fullpath = '/tmp/kragelpmaps_aqueduct.nii';
    kragelpmaps_aqueduct.write('overwrite');
    apply_spm_warp(kragelpmaps_aqueduct.fullpath, which(sprintf('%s_T1_1mm.nii.gz',SPACE)), ...
        which('00_fsl_to_fmriprep_subctx_AffineTransform.csv'), ...
        which('y_01_fsl_to_fmriprep_subctx_DisplacementFieldTransform.nii'),...
        [],...
        '/tmp/kragelpmaps_aqueduct_MNI152NLin2009cAsym.nii', ...
        1)
    delete(kragelpmaps_aqueduct.fullpath);
    kragelpmaps_aqueduct = fmri_data('/tmp/kragelpmaps_aqueduct_MNI152NLin2009cAsym.nii');
    delete('/tmp/kragelpmaps_MNI152NLin2009cAsym.nii');
end

atlas_obj = atlas(kragelpmaps, ...
    'labels',kragelPAG_dil.labels,...
    'label_descriptions',kragelPAG_dil.label_descriptions,...
    'labels_2',repmat({'PAG'},1,num_regions(kragelPAG_dil)),...
    'labels_3',repmat({'PAG'},1,num_regions(kragelPAG_dil)),...
    'labels_4',repmat({'Midbrain'},1,num_regions(kragelPAG_dil)),...
    'labels_5',repmat({'Kragel2019'},1,num_regions(kragelPAG_dil)),...
    'space_description',SPACE);


% edit shen to erode PAG and aqueductal regions
shen_ind = find(contains(bstem_atlas.labels_5, 'Shen'));
kragelpmaps_aqueduct = kragelpmaps_aqueduct.resample_space(bstem_atlas);
prob_pag = min(1,kragelpmaps_aqueduct.dat + sum(atlas_obj.probability_maps,2));
pag_mask = prob_pag > 0;
prob_pag(pag_mask) = prob_pag(pag_mask)+0.1;
for i=1:length(shen_ind)
    ind = shen_ind(i);
    bstem_atlas.probability_maps(pag_mask,ind) = min(bstem_atlas.probability_maps(pag_mask,ind), 1-prob_pag(pag_mask));
end

bstem_atlas = atlas_obj.merge_atlases(bstem_atlas);

bstem_atlas.labels = cellfun(@(x1)['BStem_' x1], bstem_atlas.labels, 'UniformOutput', false);
bstem_atlas.labels{contains(bstem_atlas.labels,'STH_L')} = 'STH_L';
bstem_atlas.labels{contains(bstem_atlas.labels,'STH_L')} = 'STH_R';

bstem_atlas.labels_2 = cellfun(@(x1)['BStem_' x1], bstem_atlas.labels_2, 'UniformOutput', false);
bstem_atlas.labels_2{contains(bstem_atlas.labels_2,'STH_L')} = 'STH_L';
bstem_atlas.labels_2{contains(bstem_atlas.labels_2,'STH_L')} = 'STH_R';

%% Add references
%{
references = {cit.references, ...
    shen_references,...
    'Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.'};
%}

references = {cit.references, ...
    kragelPAG.references, ...
    'Shen X, Tokoglu F, Papademetris X, Constable R. Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. Neuroimage 82, 403-415, 2013.'};

bstem_atlas.references = unique(char(references),'rows');

%% REMOVE general region - not needed

try
    bstem_atlas = remove_atlas_region(bstem_atlas, {'other'});
end

clear shen cit