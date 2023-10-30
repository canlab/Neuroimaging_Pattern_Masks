
morel = load_atlas(sprintf('morel_%s', ALIAS));

%% Thalamus 
switch SCALE
    case 'coarse'
        % See Create_Morel_thalamus_atlas_object, bottom
        
        group_codes = {'Pu' 'LGN', 'MGN', 'VPL', 'VPM', {'CL' 'CeM' 'CM' 'Pf'}, {'Pv' 'SPf'}, 'LD' 'VL', 'LP', 'VA' 'VM' 'MD' 'AM' 'AV' 'Hb'};  % each is a group to load and add
        group_names = {'Pulv' 'LGN', 'MGN', 'VPL', 'VPM', 'Intralam', 'Midline' 'LD', 'VL', 'LP' 'VA' 'VM', 'MD' 'AM', 'AV' 'Hb'};
        
        group_labels2 = {'Pulv' 'Genic', 'Genic', 'VPgroup', 'VPgroup', 'LD', 'Intralam', 'Midline' 'Lat_group' 'Lat_group' 'Ant_group' 'Ant_group' 'MD_group' 'Ant_group' 'Ant_group' 'Habenula'};
        
        % VL: Motor/Emo, dentate/cerebellum
        % VA, CM: Connections with BG, pallidum
        % MD: ACC, insula
        %https://www.dartmouth.edu/~rswenson/NeuroSci/chapter_10.html
        % https://www.ncbi.nlm.nih.gov/pubmed/18273888
        
        % AD is tiny and next to DM
        
        for i = 1:length(group_codes)
            
            if iscell(group_codes{i})
                roi_atlas = select_atlas_subset(morel, group_codes{i}, 'flatten'); orthviews(roi_atlas)
            else
                roi_atlas = select_atlas_subset(morel, group_codes(i), 'flatten'); orthviews(roi_atlas)
            end
            
            if i == 1
                thalamus_atlas = roi_atlas;
            else
                thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
            end
            
        end
        
        thalamus_atlas.labels_2 = group_labels2;
        
        thalamus_atlas.atlas_name = 'thalamus2023_fmriprep20_combined';

    switch 'fine'
        exclude = {'mtt', ... % mammilothalamic tract
            'RN', ... % red nucleus, we get it from bianciardi
            'SubThalamic nucleus', ... % we get it from bianciardi
            'AD', ... % anterior dorsal nucleus, 9 and 13 voxels, too small
            'MV', ... %
        
        % AD - anterior dorsal nucleus
        % Li - limitans nucleus
        % Po - Posterior nucleus
        % MV - medioventral nucleus
        % PO - posterior nucleus
        % SG - suprageniculate nucleus
        % VPi - Ventral posterior inferior nucleus
        thalamus_atlas = morel.select_atlas_subset(find(~contains()));
%% Re-label with group names

thalamus_atlas.labels = group_names;

thalamus_atlas.labels_2 = group_labels2;
thalamus_atlas.labels_2{end+1} = 'Hythal';

%% Load CIT atlas and add regions

cit168 = load_atlas(sprintf('CIT168_%s', ALIAS));

group_codes = {'Hythal'};  % 'STN' will be in BG

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i}, 'flatten'); orthviews(roi_atlas)
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i), 'flatten'); orthviews(roi_atlas)
    end
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

%% 
thalamus_atlas = lateralize(thalamus_atlas);