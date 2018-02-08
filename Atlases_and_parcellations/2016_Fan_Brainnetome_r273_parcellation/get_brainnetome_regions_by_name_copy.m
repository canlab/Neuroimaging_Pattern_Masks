function [output_obj, indx] = brainnetome_get_roi(string_to_find)
% [output_obj, indx] = brainnetome_get_roi(string_to_find_in_names)
% indx: index values
% output_obj: fmri_data object with mask
%
% Amygdala from Brainnetome atlas, Fan 2016
% parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');
% parcel_obj = fmri_data(parcellation_file);
% output_obj = select_voxels_by_value(parcel_obj, [211:214]);

%string_to_find = 'Amyg';

parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');

names_file = fullfile(fileparts(parcellation_file), 'cluster_names.mat');
load(names_file, 'names');

% Find which names match
wh = ~cellfun(@isempty, strfind(names, string_to_find));

% Index values to look for
indx = find(wh)';

%%

parcel_obj = fmri_data(parcellation_file);

%%
output_obj = select_voxels_by_value(parcel_obj, indx);


