
thalamus_atlas = load_atlas(sprintf('iglesias_thal_%s', ALIAS));
thalamus_atlas.labels_5 = repmat({'Iglesias2018'},1,num_regions(thalamus_atlas));

% Add labels to make more consistent with other atlases
% create probability maps too and default to 50% so that any conflicts with
% non-thalamic regions that are better than not to be correct are ceded to 
% said regions
for i = 1:num_regions(thalamus_atlas)
    for fname = {'labels','labels_2','labels_3','labels_4'}
        thalamus_atlas.(fname{1}){i} = regexprep(thalamus_atlas.(fname{1}){i},'([LR])_(.*)','Thal_$2_$1');
    end
end

%% add dummy probabilities
thalamus_atlas.probability_maps = sparse(double(thalamus_atlas.probability_maps));
thalamus_atlas.references = unique(thalamus_atlas.references,'rows');

%% Load CIT atlas and add regions

cit168 = load_atlas(sprintf('CIT168_%s', ALIAS));

group_codes = {'Haben'};
group_descript = {'Habenula'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i});
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i));
    end
    roi_atlas.label_descriptions = group_descript(i);
    [roi_atlas.labels_2, roi_atlas.labels_3, roi_atlas.labels_4] = deal(roi_atlas.labels);
    % this right here introduces potential whacky behavior when
    % downsampling because the rest of the posterior thalamus comes from
    % iglesias, not CIT like the habenula. It should be ok conceptually,
    % because higher order labels are just for identifying sources, not
    % actually meant to be used as an atlas.
    roi_atlas.labels_4 = repmat({'Thal_Posterior'},1,num_regions(roi_atlas));
    roi_atlas = lateralize(roi_atlas);
    roi_atlas.labels_5 = repmat({'CIT168 v1.1.0 subcortical'},1,num_regions(roi_atlas));
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

%% Load Hypothalamic atlas and add regions

hypothal_atlas = load_atlas(sprintf('iglesias_hypothal_%s', ALIAS));

group_codes = {'hypothalamus'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(hypothal_atlas, group_codes{i}, 'labels_3');
    else
        roi_atlas = select_atlas_subset(hypothal_atlas, group_codes(i), 'labels_3');
    end
    roi_atlas.labels_4 = roi_atlas.labels_3;
    roi_atlas.labels_5 = repmat({'Billot/Iglesias2020'},1,num_regions(roi_atlas));
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
end

% consolidate optic nuclei for coarser parcellations
thalamus_atlas.labels_2(ismember(thalamus_atlas.labels_2,{'L_hypothalamus_anterior_inferior','L_hypothalamus_tubular_inferior'})) = {'L_hypothalamus_anterior_and_tubular_inferior'};
thalamus_atlas.labels_2(ismember(thalamus_atlas.labels_2,{'R_hypothalamus_anterior_inferior','R_hypothalamus_tubular_inferior'})) = {'R_hypothalamus_anterior_and_tubular_inferior'};

thalamus_atlas.labels_2(ismember(thalamus_atlas.labels_2,{'L_hypothalamus_anterior_superior','L_hypothalamus_tubular_superior'})) = {'L_hypothalamus_anterior_and_tubular_superior'};
thalamus_atlas.labels_2(ismember(thalamus_atlas.labels_2,{'R_hypothalamus_anterior_superior','R_hypothalamus_tubular_superior'})) = {'R_hypothalamus_anterior_and_tubular_superior'};

%% dilate the cifti atlas to include the entire hypothalamus
% otherwise this trunctates the chiasmatic nuclei

cifti_atlas = cifti_atlas.replace_empty();
hypothal_atlas = hypothal_atlas.resample_space(cifti_atlas);
hypothal_atlas = fmri_mask_image(hypothal_atlas.threshold(0.05,'spin_off_parcel_fragments'));
vx_ind = find(hypothal_atlas.dat);
target_mm = cifti_atlas.volInfo.mat*[cifti_atlas.volInfo.xyzlist(vx_ind,:)';ones(1,length(vx_ind))];
target_mm = target_mm(1:3,:);
region = zeros(length(vx_ind),1);
cifti_atlas = cifti_atlas.remove_empty();
for i = 1:length(vx_ind)
    region(i) = cifti_atlas.find_closest_region(target_mm(:,i)).region_number;
end
cifti_atlas = cifti_atlas.replace_empty();
cifti_atlas.dat(vx_ind) = region;

%% renormalize probabilities
total_p = sum(thalamus_atlas.probability_maps,2);
renorm = total_p > 1;
thalamus_atlas.probability_maps(renorm,:) = thalamus_atlas.probability_maps(renorm,:)./total_p(renorm);