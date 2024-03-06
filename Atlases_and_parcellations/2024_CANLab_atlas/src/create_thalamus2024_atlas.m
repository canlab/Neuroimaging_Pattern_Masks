
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

group_codes = {'Hythal','Haben'}; % we'll add some more with the brainstem
group_descript = {'Hypothalamus', 'Habenula'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i});
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i));
    end
    roi_atlas.label_descriptions = group_descript(i);
    [roi_atlas.labels_2, roi_atlas.labels_3, roi_atlas.labels_4] = deal(roi_atlas.labels);
    roi_atlas.labels_4 = repmat({'Midbrain'},1,num_regions(roi_atlas));
    roi_atlas = lateralize(roi_atlas);
    roi_atlas.labels_5 = repmat({'CIT168 v1.1.0 subcortical'},1,num_regions(roi_atlas));
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

%% renormalize probabilities
total_p = sum(thalamus_atlas.probability_maps,2);
renorm = total_p > 1;
thalamus_atlas.probability_maps(renorm,:) = thalamus_atlas.probability_maps(renorm,:)./total_p(renorm);

clear morel