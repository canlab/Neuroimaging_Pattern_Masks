% Default Mode A regions
% ---------------------------------------------------------------------
% A: vmPFC, PCC, TPJ

b2 = select_atlas_subset(b, {'Default_ModeA'}, 'labels_2');

plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');

montage(b2.region_atlas);

%% Default Mode B regions
% ---------------------------------------------------------------------
% A: dmPFC, aINS, vlPFC

b2 = select_atlas_subset(b, {'Default_ModeB'}, 'labels_2');

plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');

montage(b2.region_atlas);

%% Default Mode C regions
% ---------------------------------------------------------------------
% C: PCC/rs, MTL/hipp, R vTC/Occ/pTPJ

b2 = select_atlas_subset(b, {'Default_ModeC'}, 'labels_2');

plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');

montage(b2.region_atlas);
%% Limbic regions
% ---------------------------------------------------------------------
% Limbic: mOFC, sgACC, TP/ventral TC, pst mod-lat OFC, cortex around NAC/BF  

b2 = select_atlas_subset(b, {'Limbic'}, 'labels_2');

plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');

montage(b2.region_atlas);

%% Fronto-parietal B

b2 = select_atlas_subset(b, {'Cortex_Fronto_ParietalB'}, 'labels_2'); orthviews(b2.region_atlas); spm_orthviews('Reposition', [0 35 -12]);

% NOTE: searching for OFC regions missing from other networks - found in: Cortex_Fronto_ParietalB


%% Full network - hubs
% Graph of most connected regions

[provincial_regions, connector_regions] = montage_tables_provincial_connector_hubs(b);


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
    
    % mean connectivity across regions in this cluster
    %     mymean = mean(connection_vals{i}, 2);
    %     thr = prctile(mymean, 90);
    %
    %     wh_connected{i} = find(mymean > thr);
    %
    %     % brainpathway object with regions connected to this cluster
    %     cluster_connect_b{i} = select_atlas_subset(obj_subset, wh_connected{i});
    %
    %     my_regions = atlas2region(cluster_connect_b{i}.region_atlas);
    %
    o2 = removeblobs(o2);
    
    my_regions = region(fmri_dat_connectivity_maps);
    o2 = addblobs(o2, my_regions, 'color', colors{i}, 'trans');
    %
    drawnow, snapnow;
    
end


    
%%

obj_subset = reorder_regions_by_node_cluster(obj_subset);

plot_connectivity(obj_subset, 'partitions', obj_subset.node_clusters, 'partitionlabels', obj_subset.node_cluster_labels, 'nonumbers');



%% Cortical networks
% ------------------------------------------------------------------
% Connector Hubs: Most frequently connected Def Mode A regions with others
%
% connnector hubs in red: includes vmPFC
% provincial hubs in blue: includes pCC and parts of lateral OFC

b2 = select_atlas_subset(b, {'Def' 'Limbic' 'Cortex_Fronto_ParietalB'}, 'labels_2');

[provincial_regions, connector_regions] = montage_tables_provincial_connector_hubs(b2);

%% Default-mode <-> Brainstem+
% ------------------------------------------------------------------
% Connector Hubs: Most frequently connected Def Mode A regions with others
%
% connnector hubs in red: includes vmPFC
% provincial hubs in blue: includes pCC and parts of lateral OFC

% 'Cortex_Fronto_ParietalB' has OFC, so add that
b2 = select_atlas_subset(b, {'Cortex_Default_ModeA' 'Brainstem' 'Basal_forebrain' 'Diencephalon' 'Amygdala'}, 'labels_2');

[provincial_regions, connector_regions] = montage_tables_provincial_connector_hubs(b2);

% Hubs change depending on which networks are selected, but they should -
% whether a region is broadly connected depends on what it is connected to.


%%
% Identify regions most strongly connected with other-network regions


b_defA = select_atlas_subset(b2, find(wh_vec));
orthviews(b_defA.region_atlas);


b_defA = select_atlas_subset(b2, find(wh_vec));
orthviews(b_defA.region_atlas);




% 
% function [provincial_regions, connector_regions] = montage_tables_provincial_connector_hubs(b2)
% 
% b2 = degree_calc(b2);
% plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', unique(b2.region_atlas.labels_2));
% 
% o2 = montage(b2.region_atlas, 'color', [.5 1 .5], 'trans');
% drawnow
% 
% hubs = find(b2.graph_properties.regions.between_hubs);
% b_connector = select_atlas_subset(b2, hubs);
% 
% connector_regions = region(b_connector.region_atlas);
% o2 = addblobs(o2, connector_regions, 'color', [1 .4 .3]);
% 
% hubs = find(b2.graph_properties.regions.within_hubs);
% b_provincial = select_atlas_subset(b2, hubs);
% 
% provincial_regions = region(b_provincial.region_atlas);
% o2 = addblobs(o2, provincial_regions, 'color', [.3 .4 1]);
% 
% disp('Connector Hubs (red):')
% table(atlas2region(b_connector.region_atlas), 'nolegend');
% 
% disp('Provincial Hubs (blue):')
% table(atlas2region(b_provincial.region_atlas), 'nolegend');
% 
% end
