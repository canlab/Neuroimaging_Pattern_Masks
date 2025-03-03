dat = readtable(dbname);

% Format for MKDA setup

% D = table2dataset(dat);
% D = table2struct(dat);

% Convert table to structure
% -----------------------------------

N = dat.Properties.VariableNames;

S = struct;

for i = 1:length(N)
    
    S.(N{i}) = dat.(N{i});
    
end

% Create MKDA required fields
% -----------------------------------

% Create .study cell vector of strings that MKDA Meta_Setup
u = unique(S.id);
S.study = {};
for i = 1:length(u), S.study(S.id == u(i)) = {sprintf('Study_%d', u(i))}; end
S.study = S.study';

% Coordinate space: use field names MKDA Meta_Setup will recognize
S.CoordSys = S.space;
S.CoordSys = strrep(S.CoordSys, 'TAL', 'T88');

% Contrast: Assign sequential integers for each unique table (contrast)
u = unique(S.table_id);
S.Contrast = zeros(size(S.x));
for i = 1:length(u), S.Contrast(S.table_id == u(i)) = i; end

% Sample size not available, so create a placeholder that will weight all studies equally
S.Subjects = ones(size(S.x));

% Remove empty and duplicate fields

S = rmfield(S, 'id');
S = rmfield(S, 'space');
S = rmfield(S, 'table_id');
S = rmfield(S, 'doi');

% Subjects          : Sample size of the study to which the coordinate belongs
%   FixedRandom       : Study used fixed or random effects. 
%                       Values should be Fixed or Random.
%                       Fixed effects coordinates will be automatically
%                       downweighted
%   SubjectiveWeights : A coordinate or contrast weighting vector based on FixedRandom 
%                       and whatever else you want to weight by; e.g., study reporting threshold
%                       The default is to use FixedRandom only if available
%   x, y, z           : X, Y, and Z coordinates
%   study             : name of study
%   Contrast          : unique indices (e.g., 1:k) for each independent
%                       contrast. This is a required variable!
%                       All rows belonging to the same contrast should
%                       (almost) always have the same values for every
%                       variable other than x, y, and z.
%   CoordSys          : Values should be MNI or T88, for MNI space or Talairach space
%                       Talairach coordinates will be converted to MNI using
%                       Matthew Brett's transform

%% Run Meta_Setup and prepare matrix of contrasts x voxels

DB = Meta_Setup(S, 4);  % Prepare to run with 4 mm radius

% Print a brief summary
fprintf('This database contains %s studies, %s unique contrasts (tables), and %s activation coordinates\n', ...
    length(unique(DB.study)), length(unique(DB.Contrast)), length(DB.x));

%% Run Meta_Activation_FWE
% Convolve contrast maps and create data matrix
% This can be used for co-activation, etc.

MC_Setup = Meta_Activation_FWE('setup', DB);

% Creates MC_Info.mat with MC_Setup variable

%% Add mean activation data to an fmri_data image

obj = fmri_data(MC_Setup.volInfo.fname);
obj.volInfo = MC_Setup.volInfo;
obj.dat = activation_proportions;
obj.volInfo.cluster = ones(obj.volInfo.n_inmask, 1);

obj.dat = MC_Setup.unweighted_study_data;

% Notes: can't use single format (sparse), can't take mean (too large!)

%% Extract average data from each contrast for 500+ parcels

% parcel_means = apply_parcellation(obj, atlas_obj); % doesn't work - memory/array size problems

% resample atlas first
atlas_obj = resample_space(atlas_obj, obj);

%% do it manually for matched spaces
k = max(atlas_obj.dat);
n = size(obj.dat, 2);

parcel_means = NaN * ones(n, k);

for i = 1:k
    
    wh = atlas_obj.dat == i;
    
    parcel_means(:, i) = mean(obj.dat(wh, :))';
    
end

save canlab_combined_atlas_2018_parcel_means_neurosynth parcel_means atlas_obj

% these are means over different numbers of voxels for each parcel.
% normalize.

%% Create brainpathway object
% ---------------------------------------------------------------------

atlas_obj = load_atlas('canlab2018_2mm');
b = brainpathway(atlas_obj); % Construct a brainpathway object from an atlas object
b.region_dat = parcel_means;
  
% Reorder regions by larger network units to plot blocks
% ---------------------------------------------------------------------


network_names = unique(b.region_atlas.labels_2, 'stable');
[indic,network_names] = string2indicator(b.region_atlas.labels_2', network_names);

% new order for regions
new_order = [];
for i = 1:size(indic, 2)
    new_order = [new_order; find(indic(:, i))]; 
end

b.region_atlas = reorder_atlas_regions(b.region_atlas, new_order);

network_names = unique(b.region_atlas.labels_2, 'stable');
[indic,network_names] = string2indicator(b.region_atlas.labels_2', network_names);
condf = indic2condf(indic);

b.region_dat = parcel_means(:, new_order);
b.node_clusters = condf;
b.node_cluster_labels = network_names;

plot_connectivity(b, 'partitions', b.node_clusters, 'partitionlabels', b.node_cluster_labels);

b = degree_calc(b);

% save canlab_combined_atlas_2018_resorted_brainpathway_obj b



