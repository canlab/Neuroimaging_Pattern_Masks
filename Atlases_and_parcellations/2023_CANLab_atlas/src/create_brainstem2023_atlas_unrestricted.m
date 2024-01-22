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
bstem_mask = threshold(bstem_mask, .2);

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
    shenR = lateralize(shen.select_atlas_subset([10,12:16])).select_atlas_subset({'_R'});
    shenL = lateralize(shen.select_atlas_subset([26,28:31])).select_atlas_subset({'_L'});
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

bstem_atlas.labels_4 = repmat({'Brainstem'},1,num_regions(bstem_atlas));
diencephalic_ind = find(contains(bstem_atlas.labels, {'Midb_Lrd','Midb_Rrd'}));
bstem_atlas.labels_4(diencephalic_ind) = repmat({'Diencephalon'},1,length(diencephalic_ind));

%% add other regions

cit = load_atlas(sprintf('cit168_%s', ALIAS));

% bianciardi atlas also has SNc, SNr, RN, STH, and even subdivides the RN
% and STH into subparcels, but it's not open source so we prefer these
% alternatives. The subdivisions of the RN and especially STh are also 
% non-contiguous, which makes using them awkward.
cit_regions = {'PBP', 'VTA', 'Mamm', 'SNc', 'SNr', 'RN', 'STH'};

cit = lateralize(cit.select_atlas_subset(cit_regions));
cit.labels_2 = {};
cit.labels_3 = {};
cit.labels_4 = {};
% add label description and group regions
for i = 1:num_regions(cit)
    if contains(cit.labels{i},'SN')
        cit.labels_2{end+1} = cit.labels{i};
        cit.labels_3{end+1} = regexprep(cit.labels{i},'.*_([LR])','SN_$1');
    else
        cit.labels_2{end+1} = cit.labels{i};
        cit.labels_3{end+1} = cit.labels{i};
    end
    cit.labels_4{end+1} = regexprep(cit.labels{i},'.*_([LR])','Diencephalic_nuclei_$1');
end
cit.labels_5 = repmat({'CIT168 v1.1.0 subcortical'}, 1, num_regions(cit));

% the order here is important. The most sensitive atlases (i.e. those with
% the smallest regions) go first so that we resample to their spaces.
bstem_atlas = cit.merge_atlases(bstem_atlas);


%% Adjust labels
% make more consistent with other atlases
% relabel L and R

pat = 'Reg_1';
bstem_atlas.labels = regexprep(bstem_atlas.labels, pat, 'other');

pat = 'Shen_';
bstem_atlas.labels = regexprep(bstem_atlas.labels, pat, '');

bstem_atlas = atlas_add_L_R_to_labels(bstem_atlas);



%% Add references
%{
references = {cit.references, ...
    shen_references,...
    'Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.'};
%}

references = {cit.references, ...
    'Shen X, Tokoglu F, Papademetris X, Constable R. Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. Neuroimage 82, 403-415, 2013.'};

bstem_atlas.references = unique(char(references),'rows');

%% REMOVE general region - not needed

try
    bstem_atlas = remove_atlas_region(bstem_atlas, {'other'});
end

clear shen cit