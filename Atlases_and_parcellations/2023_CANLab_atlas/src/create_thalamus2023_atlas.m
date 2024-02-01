
morel = load_atlas(sprintf('morel_%s', ALIAS));

%% Thalamus 
group_codes = {{'PuM', 'PuI', 'PuL', 'PuA'}, ...
    {'LGNmc', 'LGNpc'}, ...
    {'MGN','SG','Li','Po'}, ...
    {'CL', 'CeM', 'MV' 'CM' 'Pf'}, ...
    {'AD','AM','AV'}, ...
    'LP', ...
    'LD', ...
    {'MDmc', 'MDpc'}, ...
    {'VPLa', 'VPLp','VPM','VPI'}...
    {'VAmc','VApc'},...
    {'VLa','VLpd','VLpv','VM'} };  % each is a group to load and add

group_labels2 = {'Pulvinar' 'Lateral Geniculate Nucleus', 'Medial Geniculate Nucleus', 'Intralaminar', ...
    'Anterior', 'Lateral Posterior', 'Lateral dorsal', 'Mediodorsal', 'Ventral posterior', 'Ventral anterior', 'Ventral lateral-medial'};

group_labels3 = {'Pulvinar' 'Geniculate Nuclei', 'Geniculate Nuclei', 'Intralaminar', ...
    'Anterior', 'Lateral', 'Lateral', 'Mediodorsal', 'Ventral', 'Anterior', 'Ventral'};

label_descriptions = {'Anterior pulvinar', 'Inferior pulvinar',...
            'Lateral pulvinar', 'Medial pulvinar', 'Lateral geniculate nucleus (magnocellular)', ...
            'Lateral geniculate nucleus (parvocellular)', 'Limitans nucleus', 'Medial geniculate nucleus','Posterior nucleus', ...
            'Suprageniculate nucleus','Central lateral nucleus','Centre median nucleus', ...
            'Central medial nucleus', 'Medioventral nucleus', ...
            'Parafascicular nucleus', 'Anterior dorsal nucleus', 'Anterior medial nucleus',...
            'Anterior ventral nucleus', 'Lateral posterior nucleus', 'Lateral dorsal nucleus',...
            'Medial dorsal nucleus (magnocellular)','Medial dorsal nucleus (parvocellular)',...
            'Ventral posterior inferior nucleus', ...
            'Ventral posterior lateral nucleus (anterior part)', ...
            'Ventral posterior lateral nucleus (posterior part)', ...
            'Ventral posterior medial nucleus', 'Ventral anterior nucleus (magnocellular)',...
            'Ventral anterior nucleus (parvocellular)', ...
            'Ventral lateral anterior nucleus',...
            'Ventral lateral posterior nucleus (drosal part)',...
            'Ventral lateral posterior nucleus (ventral part)',...
            'Ventral medial nucleus'};

group_labels2 = cellfun(@(x1)regexprep(x1,'[ -]','_'),group_labels2, 'UniformOutput', false);
group_labels3 = cellfun(@(x1)regexprep(x1,'[ -]','_'),group_labels3, 'UniformOutput', false);

%group_codes = [group_codes {'mtt'}];
%group_labels2 = [group_labels2 {'Mammillothalamic tract'}];
%label_descriptions = [label_descriptions {'Mammillothalamic tract'}];

new_desc = {};
for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(morel, group_codes{i});

        ind = find(contains(cat(2,group_codes{:}), group_codes{i}));
    else
        roi_atlas = select_atlas_subset(morel, group_codes(i));

        ind = find(contains(cat(2,group_codes{:}), group_codes(i)));
    end

    new_desc0 = label_descriptions(ind(1));
    for j = ind(2:end)
        new_desc0{end+1} = label_descriptions{j};
    end
    
    new_desc_L = cellfun(@(x1)([x1 ' (left)']),new_desc0,'UniformOutput',false);
    new_desc_R = cellfun(@(x1)([x1 ' (right)']),new_desc0,'UniformOutput',false);
    new_desc = [new_desc new_desc_L new_desc_R];
    
    roi_atlas.labels_2 = repelem({[group_labels2{i}, '_L'], [group_labels2{i}, '_R']}, 1, num_regions(roi_atlas)/2);
    %error('redo label_3 so that it''s a coarser parcellation than label_2');
    roi_atlas.labels_3 = repelem({[group_labels3{i}, '_L'], [group_labels3{i}, '_R']}, 1, num_regions(roi_atlas)/2);
    roi_atlas.labels_4 = repelem({'Thalamus_L', 'Thalamus_R'}, 1, num_regions(roi_atlas)/2);
    roi_atlas.labels_5 = repmat({'Morel2010'},1,num_regions(roi_atlas));

    if i == 1
        thalamus_atlas = roi_atlas;
    else
        thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    end
    
end

thalamus_atlas.label_descriptions = new_desc(:);

% Add labels to make more consistent with other atlases
% create probability maps too and default to 50% so that any conflicts with
% non-thalamic regions that are better than not to be correct are ceded to 
% said regions
for i = 1:num_regions(thalamus_atlas)
    for fname = {'labels','labels_2','labels_3'}
        thalamus_atlas.(fname{1}){i} = [ 'Thal_' thalamus_atlas.(fname{1}){i}]; 
    end
end

%% add dummy probabilities
uniq_roi = unique(thalamus_atlas.labels);
pmap = zeros(size(thalamus_atlas.dat,1), num_regions(thalamus_atlas));
for i = 1:num_regions(thalamus_atlas)
    pmap(thalamus_atlas.dat == i,i) = 0.5;
end

thalamus_atlas.probability_maps = sparse(pmap);
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