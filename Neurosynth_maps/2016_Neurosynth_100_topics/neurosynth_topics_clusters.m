% Multidimensional scaling 
% ------------------------------------------------
mdsinfo_struct = [];
mdsinfo_struct.names = topic_obj_forwardinference.metadata_table.Topic;
mdsinfo_struct.names = format_strings_for_legend(mdsinfo_struct.names);

mdsinfo_struct.r = corr(topic_obj_forwardinference.dat);

% Correlation distance
% figure; imagesc(mdsinfo_struct.r)
mdsinfo_struct.D = (1 - mdsinfo_struct.r) ./ 2;
% figure; imagesc(mdsinfo_struct.D)

[mdsinfo_struct.MDScoords,mdsinfo_struct.obs,mdsinfo_struct.implied_dissim] = shepardplot(mdsinfo_struct.D,[]);
% 10

create_figure('mds plot'); 
point_handles = plot(mdsinfo_struct.MDScoords(:, 1), mdsinfo_struct.MDScoords(:, 2), 'o', 'MarkerFaceColor', [.2 .2 .2]);

% for lines - 1/0 matrix of which lines to draw
mdsinfo_struct.sig = sign(mdsinfo_struct.r) .* abs(mdsinfo_struct.r) > 0.4;

line_handles =  nmdsfig_tools('drawlines',mdsinfo_struct.MDScoords, mdsinfo_struct.sig);
set(line_handles.hhp, 'LineWidth', .5, 'Color', [.7 .7 .7])

%% Clustering
% ------------------------------------------------

% Convert to condensed vector form needed by linkage()
distvec = squareform(mdsinfo_struct.D, 'tovector');

% Perform hierarchical clustering using Ward linkage
linkage_tree = linkage(distvec, 'average');

k = 12;  % specify the number of desired clusters
mdsinfo_struct.cluster_labels = cluster(linkage_tree, 'maxclust', k);

% Optional: Plot dendrogram
create_figure('dendrogram');
[~, wh_order] = dendrogram(gca, linkage_tree, 0, 'Labels', mdsinfo_struct.names, 'ClusterIndices', mdsinfo_struct.cluster_labels, 'ShowCut', true);
title('Hierarchical Clustering Dendrogram');
ylabel('Distance');
set(gca, 'FontSize', 18)


% nmdsfig_plot showing MDS space and clusters
% ------------------------------------------------
% for compatibility with nmdsfig_plot
mdsinfo_struct.GroupSpace = mdsinfo_struct.MDScoords;
mdsinfo_struct.ClusterSolution.classes =  mdsinfo_struct.cluster_labels;
mdsinfo_struct.STATS.sigmat = mdsinfo_struct.sig;
mdsinfo_struct.STATS.sigmat2 = [];
mdsinfo_struct.colors = seaborn_colors(k);

out = nmdsfig_tools('nmdsfig_plot',mdsinfo_struct, 0, 0, 'fill');

%%
T = addvars(topic_obj_forwardinference.metadata_table, mdsinfo_struct.cluster_labels, 'Before', 2, 'NewVariableNames', {'Cluster_number'})

topic_obj_forwardinference.metadata_table = addvars(topic_obj_forwardinference.metadata_table, mdsinfo_struct.cluster_labels, 'Before', 2, 'NewVariableNames', {'Cluster_number'});

topic_group = {'Audition' 'Speech' 'Action, Perception, Body' 'Language, memory, control' 'Events' 'Reasoning' 'Degeneration' ...
    'Social cog, emotion, motivation' 'Pain' 'Alcohol' 'Somatomotor' 'Acupuncture'};

Cluster_label = {};
for i = 1:length(topic_group)

    wh = mdsinfo_struct.cluster_labels == i;

    Cluster_label(wh) = topic_group(i);

end

Cluster_label = Cluster_label'

topic_obj_forwardinference.metadata_table = addvars(topic_obj_forwardinference.metadata_table, Cluster_label, 'Before', 3, 'NewVariableNames', {'Cluster_label'});


%%



