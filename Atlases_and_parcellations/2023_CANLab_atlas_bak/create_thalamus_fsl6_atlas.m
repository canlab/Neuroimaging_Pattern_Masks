%% Group atlas

% Combine atlases to create a thalamus, hypothalamus, epithalamus atlas set
% Remove the "other thalamus"
% Add CIT168 STN and Hypothal regions

morel = load_atlas('morel_fsl6');

%% Thalamus 
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

thalamus_atlas.atlas_name = 'thalamus2023_fsl6_combined';

%% Re-label with group names

thalamus_atlas.labels = group_names;

thalamus_atlas.labels_2 = group_labels2;
thalamus_atlas.labels_2{end+1} = 'Hythal';

%% Load CIT atlas and add regions

cit168 = load_atlas('CIT168_fsl6');

group_codes = {'Hythal'};  % 'STN' will be in BG

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i}, 'flatten'); orthviews(roi_atlas)
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i), 'flatten'); orthviews(roi_atlas)
    end
    
    thalamus_atlas = merge_atlases(thalamus_atlas, roi_atlas);
    
end

%% Resample to standard space - 1 x 1 x 1 mm

% ref obj in 1 mm space
ref_obj = fmri_data(which('MNI152NLin6Asym_T1_1mm.nii.gz'));

thalamus_atlas = resample_space(thalamus_atlas, ref_obj);

%% Enforce some var types and compress

thalamus_atlas = check_properties(thalamus_atlas);
thalamus_atlas = remove_empty(thalamus_atlas);


%% Save

savefile = 'thalamus2023_fsl6_combined_atlas_object.mat';
save(savefile, 'thalamus_atlas');
