%% Cluster Def Mode A parcels based on connectivity with all other regions

% obj_subset = cluster_region_subset_by_connectivity(b, 'Cortex_Default_ModeA', 3);

[obj_subset, region_objects, region_indices] = cluster_region_subset_by_connectivity(b, 'Cortex_Default_ModeA', 3);

n_clusters = max(obj_subset.node_clusters);

cluster_connect_maps = cell(1, n_clusters);

%% Show regions in each cluster in a unique color

clear connection_vals wh_connected cluster_connect_b
colors = scn_standard_colors(n_clusters);

create_figure('fmridisplay'); axis off
o2 = canlab_results_fmridisplay;

for i = 1:n_clusters
    
    % regions in this cluster (cols)
    wh = obj_subset.node_clusters == i;
    
    my_regions = select_atlas_subset(obj_subset.region_atlas, find(wh));
    my_regions = atlas2region(my_regions);
    
    o2 = addblobs(o2, my_regions, 'color', colors{i});
end

drawnow, snapnow

%% Plot connectivity across clusters

obj_subset = reorder_regions_by_node_cluster(obj_subset);

plot_connectivity(obj_subset, 'partitions', obj_subset.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');


%% Identify regions in whole brain most connected with regions in each cluster

o2 = removeblobs(o2);
disp(' ');

for i = 1:n_clusters
    
    regions_in_cluster = obj_subset.node_labels(obj_subset.node_clusters == i);
    
    fprintf('Regions in cluster %d:\n', i);
    disp(regions_in_cluster);
    disp(' ');
    
    % Get connectivity between seed regions in cluster and ALL regions
    % Use region names to index these, as index numbers are different
    % between the subset object and full-brain object
    [fmri_dat_connectivity_maps, rr] = seed_connectivity(b, regions_in_cluster);
    
    % connectivity of all regions in 'connected' subset (rows) with regions in this cluster (cols)
    connection_vals{i} = (obj_subset.connectivity.regions.r(:, region_indices{i}));
    
    % Threshold each in-cluster region connectivity map
    % Get 95th percentile of unique region values
    clear thr
    rr(rr > .999) = NaN; % exclude seed region
    for j = 1:size(rr, 2)
        thr(j) = prctile(rr(:, j), 95);
        
        mydat = fmri_dat_connectivity_maps.dat(:, j);
        mydat(mydat < thr(j)) = 0;
        
        fmri_dat_connectivity_maps.dat(:, j) = mydat;
    end
    
    % Boil down to a single region
    % threshold based on max - any region above threshold
    for j = 1:size(rr, 2)
        in_set(:, j) = fmri_dat_connectivity_maps.dat(:, j) > thr(j);
    end
    in_set = double(max(in_set, [], 2));
    
    fmri_dat_connectivity_maps.dat = max(fmri_dat_connectivity_maps.dat, [], 2);
    fmri_dat_connectivity_maps.dat = fmri_dat_connectivity_maps.dat .* in_set;
    
    o2 = removeblobs(o2);
    
    my_regions = region(fmri_dat_connectivity_maps);
    o2 = addblobs(o2, my_regions, 'color', colors{i}, 'trans');
    %
    drawnow, snapnow;
    
end


    
%%
