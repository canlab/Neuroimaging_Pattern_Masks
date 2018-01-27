% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI


% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the full probabilistic atlas file:
parcellation_file = 'CIT168toMNI152_prob_atlas_bilat_1mm.nii';  

labels = {'Put' 'Cau' 'NAC' 'BST_SLEA' 'GPe' 'GPi' 'SNc' 'RN' 'SNr' 'PBP' 'VTA' 'VeP' 'Haben' 'Hythal' 'Mamm_Nuc' 'STN'};

atlas_obj = atlas(which(parcellation_file), ...
    'atlas_name', 'CIT168', ...
    'labels', labels, ...
    'space_description', 'MNI152 space', ...
    'references', 'Pauli 2018 Bioarxiv: CIT168 from Human Connectome Project data', 'noverbose');

% Threshold at probability 0.2 or greater and k = 3 voxels or greater
atlas_obj = threshold(atlas_obj, .2, 'k', 3);

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique');
 
% Convert to regions:
 r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
 montage(r);
 
%% save object

if dosave
    
    save CIT168_MNI_subcortical_atlas_object atlas_obj
    
end

%% write - this writes only the label image

if dosave
    
    atlas_obj.fullpath = fullfile(pwd, 'CIT168_MNI_subcortical_atlas_regions.img');
    write(atlas_obj);
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([labels{i} ' = r;']);
    
    region_names{i} = r(i).shorttitle;
    
end

save('CIT168_MNI_atlas_regions.mat', 'r', 'region_names', labels{:});

%% Older code - handles contiguous region separation
% to do this need to update atlas2regions

% region_names = repmat(names(1), 1, length(r));
% 
% for i = 2:size(atlas_obj.atlas_obj, 2)
%     
%     dat_wh = get_wh_image(atlas_obj, i); 
%     
%     % threshold at 20% probability, get rid of stray voxels
%     dat_wh = threshold(dat_wh, [.2 Inf], 'raw-between', 'k', 3);
% 
%     r_wh = region(dat_wh);
%     
%     r = [r r_wh];
%     
%     region_names = [region_names repmat(names(i), 1, length(r_wh))];
%     
%     % save regions with separate names, for loading later
%     eval([names{i} ' = r_wh;']);
%     
% end

% save('CIT168_MNI_atlas_regions.mat', 'r', 'region_names', names{:});
