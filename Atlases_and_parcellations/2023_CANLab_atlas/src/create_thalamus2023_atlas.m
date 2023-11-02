
morel = load_atlas(sprintf('morel_%s', ALIAS));

%% Thalamus 
group_codes = {{'PuM', 'PuI', 'PuL', 'PuA'}, ...
    'LGN', ...
    {'MGN','SG','Li','Po'}, ...
    {'CL', 'CeM', 'MV' 'CM' 'sPf','Pf'}, ...
    {'AD','AM','AV'}, ...
    'LP', ...
    'LD', ...
    {'MDmc', 'MDpc'}, ...
    {'VPLa', 'VPLp','VPM','VPI'}...
    {'VAmc','VApc'},...
    {'VLa','VLpd','VLpv','VM'} };  % each is a group to load and add

group_labels2 = {'Pulvinar' 'Lateral Geniculate Nucleus', 'Auditory Geniculate Nuclei', 'Intralaminar', ...
    'Anterior', 'Lateral Posterior', 'Lateral dorsal', 'Mediodorsal', 'Ventral posterior', 'Ventral anterior', 'Ventral lateral-medial'};

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
        
        for i = 1:length(group_codes)
            
            if iscell(group_codes{i})
                roi_atlas = select_atlas_subset(morel, group_codes{i}, 'flatten');
            else
                roi_atlas = select_atlas_subset(morel, group_codes(i), 'flatten');
            end
            
            if i == 1
                thalamus_atlas = roi_atlas;
            else
                thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
            end
            
        end
                
        thalamus_atlas.labels = group_names;
        thalamus_atlas.labels_2 = group_labels2;

        thalamus_atlas = lateralize(thalamus_atlas);

    case 'fine'
        
        for i = 1:length(group_codes)
            
            if iscell(group_codes{i})
                roi_atlas = select_atlas_subset(morel, group_codes{i});
            else
                roi_atlas = select_atlas_subset(morel, group_codes(i));
            end
            
            roi_atlas.labels_2 = repmat(group_labels2(i),1,length(group_codes{i}));

            if i == 1
                thalamus_atlas = roi_atlas;
            else
                thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
            end
            
        end

        thalamus_atlas.label_descriptions = repmat({'Medial pulvinar', 'Inferior pulvinar',...
            'Lateral pulvinar', 'Anterior pulvinar', 'Lateral geniculate nucleus', ...
            'Medial geniculate nucleus','Suprageniculate nucleus', 'Limitans nucleus', ...
            'Posterior nucleus', 'Central lateral nucleus', 'Central medial nucleus', ...
            'Medioventral nucleus', 'Centre median nucleus', 'Subparafascicular nucleus', ...
            'Parafascicular nucleus', 'Anterior dorsal nucleus', 'Anterior medial nucleus',...
            'Anterior ventral nucleus', 'Lateral posterior nucleus', 'Lateral dorsl nucleus',...
            'Medial dorsal nucleus (magnocellular)','Medial dorsal nucleus (parvocellular)',...
            'Ventral posterior lateral nucleus (anterior part)', ...
            'Ventral posterior lateral nucleus (posterior part)', ...
            'Ventral posterior medial nucleus', 'Ventral posterior inferior nucleus', ...
            'Ventral anterior nucleus (magnocellular)',...
            'Ventral anterior nucleus (parvocellular)', ...
            'Ventral lateral anterior nucleus',...
            'Ventral lateral posterior nucleus (drosal part)',...
            'Ventral lateral posterior nucleus (ventral part)',...
            'Ventral medial nucleus'},1,2)';
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

group_codes = {'Hythal','Haben'};  % 'STN' will be in BG
group_descript = {'Hypothalamus', 'Habenula'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i});
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i));
    end
    roi_atlas.label_descriptions = group_descript(i);
    roi_atlas.label_descriptions = {''};
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