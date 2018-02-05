atlas_obj = load_atlas('pauli_bg');

atlas_obj = check_properties(atlas_obj);

%% add regions from CIT168

cit168 = load_atlas('CIT168');

group_codes = {'NAC' 'GPe' 'GPi' 'STN'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i}, 'flatten'); orthviews(roi_atlas)
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i), 'flatten'); orthviews(roi_atlas)
    end
    
    atlas_obj = merge_atlases(atlas_obj, roi_atlas, 'always_replace');
    
end

%% add caudate tail from CIT168 without replacing existing regions

roi_atlas = select_atlas_subset(cit168, {'Cau'}, 'flatten'); orthviews(roi_atlas)

% cut off: keep only behind caudate tail
% ******

atlas_obj = merge_atlases(atlas_obj, roi_atlas, 'noreplace');

%% Enforce some var types and compress

atlas_obj = check_properties(atlas_obj);
atlas_obj = remove_empty(atlas_obj);


%% Save

cd('/Users/tor/Documents/Code_Repositories/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2018_Wager_combined_atlas');

savefile = 'Basal_ganglia_combined_atlas_object.mat';
save(savefile, 'atlas_obj');

