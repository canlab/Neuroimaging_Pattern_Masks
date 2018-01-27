function [output_obj, indx, output_names] = get_brainnetome_regions_by_name(string_to_find)
% [output_obj, indx] = get_brainnetome_regions_by_name(string_to_find_in_names)
% indx: index values
% output_obj: fmri_data object with mask
%
% Amygdala from Brainnetome atlas, Fan 2016
% parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');
% parcel_obj = fmri_data(parcellation_file);
% output_obj = select_voxels_by_value(parcel_obj, [211:214]); % Amygdala
% output_obj = select_voxels_by_value(parcel_obj, [223:224]); % NAC
% output_obj = select_voxels_by_value(parcel_obj, [219   220   227   228]); % Caudate
%
% [output_obj, indx, output_names] = get_brainnetome_regions_by_name('NAC');  % bad stuff unusable!
% [output_obj, indx, output_names] = get_brainnetome_regions_by_name('Ca');   % pretty good
% [output_obj, indx, output_names] = get_brainnetome_regions_by_name('Pu');   % ok but some problems (overlap with Glo)
% [output_obj, indx, output_names] = get_brainnetome_regions_by_name('Tha');  % ok but too big in posterior
% [output_obj, indx, output_names] = get_brainnetome_regions_by_name('INS');  % pretty good

parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');

names_file = fullfile(fileparts(parcellation_file), 'cluster_names.mat');
load(names_file, 'names');

% Find which names match
wh = ~cellfun(@isempty, strfind(names, string_to_find));

% Index values to look for
indx = find(wh)';

output_names = names(indx);

%%

parcel_obj = fmri_data(parcellation_file);

%%
output_obj = select_voxels_by_value(parcel_obj, indx);


