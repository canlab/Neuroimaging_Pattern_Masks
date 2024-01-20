
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

group_labels2 = {'Pulvinar' 'Lateral Geniculate Nucleus', 'Auditory Geniculate Nuclei', 'Intralaminar', ...
    'Anterior', 'Lateral Posterior', 'Lateral dorsal', 'Mediodorsal', 'Ventral posterior', 'Ventral anterior', 'Ventral lateral-medial'};

label_descriptions = {'Anterior pulvinar', 'Inferior pulvinar',...
            'Lateral pulvinar', 'Medial pulvinar', 'Lateral geniculate nucleus (magnocellular)', ...
            'Lateral geniculate nucleus (parvocellular)', 'Limitans nucleus', 'Medial geniculate nucleus','Posterior nucleus', ...
            'Suprageniculate nucleus','Central medial nucleus','Central lateral nucleus', ...
            'Centre median nucleus', 'Medioventral nucleus', ...
            'Parafascicular nucleus', 'Anterior dorsal nucleus', 'Anterior medial nucleus',...
            'Anterior ventral nucleus', 'Lateral posterior nucleus', 'Lateral dorsl nucleus',...
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
switch SCALE
    case 'coarse'
        % See Create_Morel_thalamus_atlas_object, bottom
        group_names = {'Pulv' 'LGN', 'MGN_SG_Li_Po', 'Intralam', 'An', 'LP', 'LD', 'MD', 'VP','VA', 'VL_VM'};
        
        % VL: Motor/Emo, dentate/cerebellum
        % VA, CM: Connections with BG, pallidum
        % MD: ACC, insula
        %https://www.dartmouth.edu/~rswenson/NeuroSci/chapter_10.html
        % https://www.ncbi.nlm.nih.gov/pubmed/18273888
        
        % AD is tiny and next to DM
        
        new_desc = {};
        for i = 1:length(group_codes)
            if iscell(group_codes{i})
                roi_atlas = select_atlas_subset(morel, group_codes{i}, 'flatten');
        
                ind = find(contains(cat(2,group_codes{:}), group_codes{i}));
                new_desc{end+1} = label_descriptions{ind(1)};
                for j = ind(2:end)
                    new_desc{end} = [new_desc{end} ', ' label_descriptions{j}];
                end

            else
                roi_atlas = select_atlas_subset(morel, group_codes(i), 'flatten');
                
                ind = find(contains(cat(2,group_codes{:}), group_codes(i)));
                new_desc{end+1} = label_descriptions{ind(1)};
                for j = ind(2:end)
                    new_desc{end} = [new_desc{end} ', ' label_descriptions{j}];
                end
            end
            
            if i == 1
                thalamus_atlas = roi_atlas;
            else
                thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
            end
            
        end
                
        thalamus_atlas.labels = group_names;
        thalamus_atlas.labels_2 = group_labels2;
        thalamus_atlas.label_descriptions = new_desc(:);

        thalamus_atlas = lateralize(thalamus_atlas);

    case 'fine'
        group_codes = [group_codes {'mtt'}];
        group_labels2 = [group_labels2 {'Mammillothalamic tract'}];
        label_descriptions = [label_descriptions {'Mammillothalamic tract'}];
        
        new_desc = {};
        for i = 1:length(group_codes)
            
            if iscell(group_codes{i})
                roi_atlas = select_atlas_subset(morel, group_codes{i});

                ind = find(contains(cat(2,group_codes{:}), group_codes{i}));
                new_desc0 = label_descriptions(ind(1));
                for j = ind(2:end)
                    new_desc0{end+1} = label_descriptions{j};
                end
                
                new_desc_L = cellfun(@(x1)([x1 ' (left)']),new_desc0,'UniformOutput',false);
                new_desc_R = cellfun(@(x1)([x1 ' (right)']),new_desc0,'UniformOutput',false);
                new_desc = [new_desc new_desc_L new_desc_R];
            else
                roi_atlas = select_atlas_subset(morel, group_codes(i));

                ind = find(contains(cat(2,group_codes{:}), group_codes(i)));
                new_desc0 = label_descriptions(ind(1));
                for j = ind(2:end)
                    new_desc0{end+1} = label_descriptions{j};
                end
                new_desc_L = cellfun(@(x1)([x1 ' (left)']),new_desc0,'UniformOutput',false);
                new_desc_R = cellfun(@(x1)([x1 ' (right)']),new_desc0,'UniformOutput',false);
                new_desc = [new_desc new_desc_L new_desc_R];
            end
            
            roi_atlas.labels_2 = repmat(group_labels2(i),1,length(group_codes{i}));

            if i == 1
                thalamus_atlas = roi_atlas;
            else
                thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
            end
            
        end

        thalamus_atlas.label_descriptions = new_desc(:);
    otherwise 
        error('Did not understand scale %s', SCALE);
end

% Add labels to make more consistent with other atlases
% create probability maps too and default to 50% so that any conflicts with
% non-thalamic regions that are better than not to be correct are ceded to 
% said regions
for i = 1:num_regions(thalamus_atlas)
    thalamus_atlas.labels{i} = [ 'Thal_' thalamus_atlas.labels{i}]; 
end

%% Load CIT atlas and add regions

cit168 = load_atlas(sprintf('CIT168_%s', ALIAS));

group_codes = {'Hythal','Haben'};  % 'STH' will be in BG
group_descript = {'Hypothalamus', 'Habenula'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i});
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i));
    end
    roi_atlas.label_descriptions = group_descript(i);
    roi_atlas.labels_2 = {''};
    roi_atlas = lateralize(roi_atlas);
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

%% add dummy probabilities
uniq_roi = unique(thalamus_atlas.labels);
pmap = zeros(size(thalamus_atlas.dat,1), num_regions(thalamus_atlas));
for i = 1:num_regions(thalamus_atlas)
    pmap(thalamus_atlas.dat == i,i) = 0.5;
end

thalamus_atlas.probability_maps = sparse(pmap);
thalamus_atlas.labels_4 = repmat({'Restricted (contact the secretariat of the Computer Vision Lab, ETH Zürich to confirm permission for use/distribution)'}, 1, num_regions(thalamus_atlas));
thalamus_atlas.labels_5 = repmat({'Morel'}, 1, num_regions(thalamus_atlas));

thalamus_atlas.references = unique(thalamus_atlas.references,'rows');
