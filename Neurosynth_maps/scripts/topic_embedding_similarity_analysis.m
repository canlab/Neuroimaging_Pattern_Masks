% Analyze document embeddings exported by generate_topic_embeddings.py.

%%  RELOAD HERE

cd /Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/Neurosynth_maps/

load neurosynth_data_obj.mat


%% Save matlab file with embeddings

repo_root = fileparts(fileparts(mfilename('fullpath')));
csv_path = fullfile(repo_root, 'topic_embeddings_text-embedding-3-large.csv');
mat_path = fullfile(repo_root, 'topic_embeddings_text-embedding-3-large.mat');

embedding_table = readtable(csv_path, 'VariableNamingRule', 'preserve');
% save(mat_path, 'embedding_table'); % run once, not on reload

embedding_matrix = table2array(embedding_table);
vector_norms = vecnorm(embedding_matrix, 2, 1);
vector_norms(vector_norms == 0) = 1;
normalized_embeddings = embedding_matrix ./ vector_norms;
cosine_similarity = normalized_embeddings' * normalized_embeddings;

figure('Name', 'Topic Embedding Matrices', 'Color', 'w');

subplot(1, 2, 1);
imagesc(embedding_matrix);
axis tight;
set(gca, 'XTick', 1:width(embedding_table), ...
    'XTickLabel', embedding_table.Properties.VariableNames, ...
    'XTickLabelRotation', 90);
title('Embedding Matrix');
xlabel('Document');
ylabel('Embedding Dimension');
colorbar;

subplot(1, 2, 2);
imagesc(cosine_similarity, [-1 1]);
axis square;
set(gca, 'XTick', 1:width(embedding_table), ...
    'XTickLabel', embedding_table.Properties.VariableNames, ...
    'XTickLabelRotation', 90, ...
    'YTick', 1:width(embedding_table), ...
    'YTickLabel', embedding_table.Properties.VariableNames);
title('Cosine Similarity');
xlabel('Document');
ylabel('Document');
colorbar;

%% Plot sorted based on hierarchical clustering

[r, indx] = canlab_sort_distance_matrix(cosine_similarity, 'correlation_matrix', true);

names = format_strings_for_legend(embedding_table.Properties.VariableNames(indx));

create_figure('cos sim sorted')
imagesc(r, [0 1]);
axis square;
set(gca, 'XTick', 1:width(embedding_table), ...
    'XTickLabel', names, ...
    'XTickLabelRotation', 90, ...
    'YTick', 1:width(embedding_table), ...
    'YTickLabel', names);
title('Cosine Similarity');
xlabel('Document');
ylabel('Document');
colorbar;

set(gca, 'YDir', 'Reverse', 'FontSize', 10); axis tight;

cm = colormap_tor([0 0 1], [1 0 0], [0 .5 .5], [.5 0 .5], [1 .5 0]);
colormap(cm)

%% Brain pattern similarity - plot in same order

obj = topic_obj_forwardinference;
obj = get_wh_image(obj, indx); % reorder maps

create_figure('r_brain'); 
plot_correlation_matrix(obj.dat, 'nofigure', 'docircles', false, 'image', 'names', names);
colormap(cm)
set(gca, 'CLim', [0 1])
set(gca, 'YDir', 'Reverse', 'FontSize', 10); axis tight;
title('Forward inference brain similarity (r)')

%%
obj = topic_obj_reverseinference;
obj = get_wh_image(obj, indx); % reorder maps

create_figure('r_brain'); 
plot_correlation_matrix(obj.dat, 'nofigure', 'docircles', false, 'image', 'names', names);
colormap(cm)
set(gca, 'CLim', [-1 1])
set(gca, 'YDir', 'Reverse', 'FontSize', 10); axis tight;
title('Reverse inference brain similarity (r)')

%% Similarity across spaces
% note: used corr for brain, but did not try cosine sim for brain

obj = topic_obj_forwardinference;
obj = get_wh_image(obj, indx); % reorder maps
obj = remove_empty(obj);
r_fwd = corr(obj.dat);
sim_fwd = squareform((1 - r_fwd) ./ 2)';
sim_fwd = double(sim_fwd);

obj = topic_obj_reverseinference;
obj = get_wh_image(obj, indx); % reorder maps
obj = remove_empty(obj);
r_rev = corr(obj.dat);
sim_rev = squareform((1 - r_rev) ./ 2)';
sim_rev = double(sim_rev);

c_sem = (1 - cosine_similarity) ./ 2;
c_sem = c_sem .* (1 - eye(size(c_sem, 1))); %  enforce zero diagonal (rounding error issue)
sim_semantic = squareform(c_sem)';

figure; plotmatrix([sim_semantic sim_fwd sim_rev])
corr([sim_semantic sim_fwd sim_rev])

% Semantic similarity does not correlate with brain similarity

%% Plot forward inf in order sorted based on hierarchical clustering

obj = topic_obj_forwardinference;
obj = remove_empty(obj);

[r, indx] = canlab_sort_distance_matrix(corr(obj.dat), 'correlation_matrix', true);

names = format_strings_for_legend(obj.metadata_table.("Topic name_Summary")(indx));

create_figure('fwd sorted')
imagesc(r, [0 1]);
axis square;
set(gca, 'XTick', 0:length(names)+1, ...
    'XTickLabel', names, ...
    'XTickLabelRotation', 90, ...
    'YTick', 0:length(names)+1, ...
    'YTickLabel', names);
title('Forward inf');
xlabel('Topic');
ylabel('Topic');
colorbar;

set(gca, 'YDir', 'Reverse', 'FontSize', 10); axis tight;

cm = colormap_tor([0 0 1], [1 0 0], [0 .5 .5], [.5 0 .5], [1 .5 0]);
colormap(cm)

%% Plot reverse inf in order sorted based on hierarchical clustering

obj = topic_obj_reverseinference;
obj = remove_empty(obj);

[r, indx] = canlab_sort_distance_matrix(corr(obj.dat), 'correlation_matrix', true);

names = format_strings_for_legend(obj.metadata_table.("Topic name_Summary")(indx));

create_figure('rev sorted')
imagesc(r, [-1 1]);
axis square;
set(gca, 'XTick', 0:length(names)+1, ...
    'XTickLabel', names, ...
    'XTickLabelRotation', 90, ...
    'YTick', 0:length(names)+1, ...
    'YTickLabel', names);
title('Reverse inf');
xlabel('Topic');
ylabel('Topic');
colorbar;

set(gca, 'YDir', 'Reverse', 'FontSize', 10); axis tight;

cm = colormap_tor([0 0 1], [1 0 0], [0 .5 .5], [.5 0 .5], [1 .5 0]);
colormap(cm)

%% Find maximal cliques

obj = topic_obj_forwardinference;
obj = remove_empty(obj);
names = format_strings_for_legend(obj.metadata_table.("Topic name_Summary"));

r = corr(obj.dat);

thr = 0.65;
A = r > thr; % adjacency matrix

% remove diagonal
A = A .* (1 - eye(width(A)));

% find cliques
cliques = maximalCliques(A);

% Display
for i = 1:length(cliques)
    fprintf('Clique %d: ', i);
    disp(cliques{i});
    disp(names(cliques{i}))
end

[unique_cliques, out_obj, membership, stats] = assign_unique_cliques_from_maximal_cliques(obj, r, cliques, thr, names);
stats
out_obj.metadata_table

obj_forwardinference_grouped = out_obj;

%%

obj = topic_obj_reverseinference;
obj = remove_empty(obj);
names = format_strings_for_legend(obj.metadata_table.("Topic name_Summary"));

r = corr(obj.dat);

thr = 0.4;
A = r > thr; % adjacency matrix

% remove diagonal
A = A .* (1 - eye(width(A)));

% find cliques
cliques = maximalCliques(A);

% Display
for i = 1:length(cliques)
    fprintf('Clique %d: ', i);
    disp(cliques{i});
    disp(names(cliques{i}))
end

[unique_cliques, out_obj, membership, stats] = assign_unique_cliques_from_maximal_cliques(obj, r, cliques, thr, names);
stats
out_obj.metadata_table

obj_reverseinference_grouped = out_obj;

%%

% save neurosynth_topics_grouped_by_clique obj_forwardinference_grouped obj_reverseinference_grouped


