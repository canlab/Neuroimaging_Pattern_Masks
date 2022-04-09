function [provincial_regions, connector_regions] = montage_tables_provincial_connector_hubs(b2)

b2 = degree_calc(b2);
plot_connectivity(b2, 'partitions', b2.node_clusters, 'partitionlabels', b2.node_cluster_labels);

o2 = montage(b2.region_atlas, 'color', [.5 1 .5], 'trans');
drawnow, snapnow

hubs = find(b2.graph_properties.regions.between_hubs);
b_connector = select_atlas_subset(b2, hubs);

o2 = removeblobs(o2);

connector_regions = region(b_connector.region_atlas);
o2 = addblobs(o2, connector_regions, 'color', [1 .4 .3]);

hubs = find(b2.graph_properties.regions.within_hubs);
b_provincial = select_atlas_subset(b2, hubs);

provincial_regions = region(b_provincial.region_atlas);
o2 = addblobs(o2, provincial_regions, 'color', [.3 .4 1]);
drawnow, snapnow

disp('Connector Hubs (red):')
table(atlas2region(b_connector.region_atlas), 'nolegend');

disp('Provincial Hubs (blue):')
table(atlas2region(b_provincial.region_atlas), 'nolegend');

end
