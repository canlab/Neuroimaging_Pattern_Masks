% this script is meant to be invoked by create_CANLab202mid3_unrestricted.m
% it's only saved separately for readability purposes

bg = load_atlas(sprintf('tian_3t_%s',ALIAS));
bg = bg.resample_space(ref,'nearest');
bg.labels_5 = repmat({'Tian2020'},1,num_regions(bg));

% replace rh/lh with R/L
fnames = {'labels','labels_2','labels_3','labels_4'};
for i = 1:length(fnames)
    bg.(fnames{i}) = cellfun(@(x1)strrep(x1,'_rh','_R'), bg.(fnames{i}), 'UniformOutput', false);
    bg.(fnames{i}) = cellfun(@(x1)strrep(x1,'_lh','_L'), bg.(fnames{i}), 'UniformOutput', false);
end

%% incorporate BST_SLEA into BG and make it lateralized
% we include parts that overlap with acumbens, GP and caudate, but not
% putamen because it's a trivial segment (4vxls) and not worth the mess,
% and not GP because it's a small segment and fragments the ROI which is
% mostly in caudate.

citfile = which(sprintf('CIT168_%s_subcortical_v1.1.0_atlas_object.mat',SPACE));
atlas_obj = load(citfile); atlas_obj = atlas_obj.atlas_obj.resample_space(ref);

BST = atlas_obj.select_atlas_subset(find(contains(atlas_obj.labels,{'BST_SLEA'}))).replace_empty();
BST.labels_2 = repmat({'BST_SLEA'}, 1, num_regions(BST));
[BST.labels_3, BST.labels_4] = deal(repmat({'VStriatum'}, 1, num_regions(BST)));
BST = lateralize(BST);
BST.labels_5 = repmat({'CIT168 subcortical v1.1.0'},1,num_regions(BST));

%% Pallidum
% tian atlas doesn't distinguish interior/exterior GP, but CIT168 does, so 
% we use that instead.
pal = load_atlas(sprintf('cit168_%s',ALIAS)).select_atlas_subset({'GP','VeP'}).remove_empty();
pal = lateralize(pal);

pal.labels_2 = pal.labels;
pal.labels_3 = pal.labels_2;
pal.labels_4 = {'GP_L','GP_L','GP_L','GP_R','GP_R','GP_R'};
pal.labels_5 = repmat({'CIT168 subcortical v1.1.0'},1,num_regions(pal));

%% expand CIFTI mask for VeP and BST
% These structures extend somewhat outside the CIFTI mask. Let's expand it
% to include most of these

cifti_atlas = cifti_atlas.replace_empty();
BST_VeP = BST.merge_atlases(pal.select_atlas_subset('VeP')).resample_space(cifti_atlas);
BST_VeP_mask = fmri_mask_image(BST_VeP.threshold(0.05,'spin_off_parcel_fragments'));
vx_ind = find(BST_VeP_mask.dat);
target_mm = cifti_atlas.volInfo.mat*[cifti_atlas.volInfo.xyzlist(vx_ind,:)';ones(1,length(vx_ind))];
target_mm = target_mm(1:3,:);
region = zeros(length(vx_ind),1);
cifti_atlas = cifti_atlas.remove_empty();
for i = 1:length(vx_ind)
    region(i) = cifti_atlas.find_closest_region(target_mm(:,i)).region_number;
end
cifti_atlas = cifti_atlas.replace_empty();
cifti_atlas.dat(vx_ind) = region;


%% finish up the pallidal region
cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'pallidum'))));
cifti_mask = cifti_mask.replace_empty();

tic
pal_dil = dilate(pal, cifti_mask);
clear pal;
toc


%% Accumbens
accumbens = load_atlas(sprintf('cartmell_NAc_%s',ALIAS));
% the right core is split with a couple of voxels ventral to shell, which
% isn't right. Let's touch it up.
accumbens = accumbens.threshold(0,'spin_off_parcel_fragments');
L_nac = accumbens.select_atlas_subset({'_L'});
R_nac_core = accumbens.select_atlas_subset({'NAc_core_R_frag1'});
R_nac_core.labels = {'NAc_core_R'};
R_nac_shell = accumbens.select_atlas_subset({'NAc_shell_R','NAc_core_R_frag2'}, 'flatten');
R_nac_shell.labels = {'NAc_shell_R'};
R_nac_shell.label_descriptions = {'NucleusAccumbens, putative shell (right)'};
R_nac_shell.labels_2 = {'NAc_R'};
accumbens = [L_nac, R_nac_core, R_nac_shell];
[accumbens.labels_3, accumbens.labels_4] = deal(cellfun(@(x1)(strrep(x1,'NAc','VStriatum')),accumbens.labels_2,'UniformOutput',false));
accumbens.labels_5 = repmat({'Cartmell2019'},1,num_regions(accumbens));

accumbens = accumbens.remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'accumbens'))));
cifti_mask = cifti_mask.replace_empty();

% dilate to fill the accumbens CIFTI mask without masking areas in CIFTI
% caudate mask
tic
accumbens_dil = dilate(accumbens, cifti_mask).replace_empty;
accumbens_dil.probability_maps = max(cat(3,accumbens_dil.probability_maps, accumbens.replace_empty.probability_maps),[],3);
accumbens_dil = accumbens_dil.probability_maps_to_region_index;
%[accumbens_dil.labels_3, accumbens_dil.labels_4] = deal(cellfun(@(x1)(strrep(x1,'NAc','VStriatum')),accumbens_dil.labels_2,'UniformOutput',false));
%accumbens_dil.labels_5 = repmat({'Cartmell2019'},1,num_regions(accumbens_dil));
accumbens_dil = BST.apply_mask(cifti_mask).merge_atlases(accumbens_dil,'noreplace');
toc

%% Caudate
% some of the ventral striatum is in the CIFTI accumbens (13%), but most is 
% in the CIFTI caudate mask (45.9%). This is a 70-30 split, and it's hard to
% assign it to one over the other like with entorhinal cortex below. We'll
% give it to both.
% We could also include 'Cau' here, but it's very sparse within the CIFTI
% mask regions and is more likely to cause problems down the line than not
% (e.g. selecting it as an ROI for seed based analysis or alignment results
% in a nonsensical fragmented region of caudate), so we drop it.
caud = bg.select_atlas_subset(find(contains(bg.labels,{'CAU'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'caudate'))));
cifti_mask = cifti_mask.replace_empty();

tic
caud_dil = dilate(caud, cifti_mask);
caud_dil = accumbens.apply_mask(cifti_mask).merge_atlases(caud_dil,'noreplace');
clear accumbens
caud_dil = BST.apply_mask(cifti_mask).merge_atlases(caud_dil,'noreplace');
caud_dil.labels_3 = deal(cellfun(@(x1)(strrep(x1,'NAc_core','VStriatum')),caud_dil.labels_3,'UniformOutput',false));
caud_dil.labels_3 = deal(cellfun(@(x1)(strrep(x1,'NAc_shell','VStriatum')),caud_dil.labels_3,'UniformOutput',false));
caud_dil.labels_4 = deal(cellfun(@(x1)(strrep(x1,'NAc','VStriatum')),caud_dil.labels_4,'UniformOutput',false));
clear caud
toc

%% Putamen
put = bg.select_atlas_subset(find(contains(bg.labels,{'PUT'}))).remove_empty();
%put = bg.select_atlas_subset(find(contains(bg.labels,{'V_Striatum','Putamen','VeP', 'BST_SLEA'}))).remove_empty();

cifti_mask = fmri_mask_image(cifti_atlas.select_atlas_subset(find(contains(cifti_atlas.labels,'putamen'))));
cifti_mask = cifti_mask.replace_empty();

tic
put_dil = dilate(put, cifti_mask);
toc

% trim dilated left and right margins of dilated putamen to avoid dilated parts bleeding over into insula
for i = 1:5
    put_xyz = put_dil.volInfo.xyzlist;
    trim_x = [min(put_xyz(put_dil.dat > 0,1)), max(put_xyz(put_dil.dat > 0,1))];
    trim_ind = ismember(put_xyz(:,1), trim_x);
    put_dil.probability_maps(trim_ind,:) = 0;
    put_dil = put_dil.probability_maps_to_region_index();
end
% let side looks like it needs more trimming than the right
for i = 1:2
    put_xyz = put_dil.volInfo.xyzlist;
    trim_x = min(put_xyz(put_dil.dat > 0,1));
    trim_ind = ismember(put_xyz(:,1), trim_x);
    put_dil.probability_maps(trim_ind,:) = 0;
    put_dil = put_dil.probability_maps_to_region_index();
end

% add back in any original voxels we'd removed (we only want to trim from
% the newly dilated voxels)
put = put.replace_empty;
put_dil.probability_maps(put.dat > 0,:) = put.probability_maps(put.dat > 0,:);
clear put

%put_dil = bg.select_atlas_subset('NAc').apply_mask(cifti_mask).merge_atlases(put_dil,'noreplace');

%% bg atlas
bg_dil = caud_dil.merge_atlases(put_dil).merge_atlases(pal_dil).merge_atlases(accumbens_dil);
%clear caud_dil put_dil pal_dil accumbens_dil

% deal with redundant labels
n_labels = length(bg_dil.labels);
uniq_labels = unique(bg_dil.labels);
remove_lbl = [];
for i = 1:length(uniq_labels)
    this_lbl = uniq_labels(i);
    this_ind = find(contains(bg_dil.labels, this_lbl));
    for j = 2:length(this_ind)
        bg_dil.dat(bg_dil.dat == this_ind(j)) = this_ind(1);
        remove_lbl(end+1) = this_ind(j);
        bg_dil.probability_maps(:,this_ind(1)) = max(bg_dil.probability_maps(:,this_ind),[],2);
    end
end
[~,~,dat] = unique(bg_dil.dat);
bg_dil.dat = dat - 1;
fnames = {'labels','label_descriptions','labels_2','labels_3','labels_4','labels_5'};
for i = 1:length(fnames)
    if length(bg_dil.(fnames{i})) == n_labels
        bg_dil.(fnames{i})(remove_lbl) = [];
    end
end
bg_dil.probability_maps(:,remove_lbl) = [];


% add BG prephix

for i = 1:length(bg_dil.labels)
    bg_dil.labels{i} = [ 'BG_' bg_dil.labels{i}]; 
end

%% eval NAc overlap with original
%{
thisAtlas = bg_dil; 
thisAtlas.probability_maps = [];
thisAtlas = thisAtlas.select_atlas_subset('NAc');
thisAtlas.labels = cellfun(@(x1)(strrep(x1,'BG_','')), thisAtlas.labels, 'UniformOutput', false);
thisAtlas.probability_maps = [];
thisLeadsAtlas = load_atlas(sprintf('cartmell_nac_%s',ALIAS));
thisLeadsAtlas.probability_maps = [];
for orientation = {'saggital','coronal','axial'}
    %%
    o2 = thisAtlas.montage('nofigure','transvalue',0.5,'regioncenters',orientation{1});
    for i = 1:num_regions(thisAtlas)
        try
            leads_roi = thisLeadsAtlas.select_atlas_subset(thisAtlas.labels(i),'exact');
        
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