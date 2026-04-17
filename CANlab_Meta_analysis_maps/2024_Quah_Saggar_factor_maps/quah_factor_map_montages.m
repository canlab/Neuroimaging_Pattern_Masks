% Get a list of all factor maps (PA*.nii.gz)
ff = filenames('PA*gz')

factor_obj = {};

% For each, load and plot a montage of the map

% Load the image
for i = 1:size(ff, 1)

factor_obj{i} = fmri_data(ff{i});

% Strip out the name
name = regexp(ff{i}, 'PA\d{1,2}', 'match', 'once');

create_figure([name 'surf']); axis off; surface(factor_obj{i}, 'foursurfaces_hcp', 'nolegend', 'nofigure');
subplot(2, 2, 1);
title(name, 'FontSize', 24)
saveas(gcf, fullfile('figures', [name '_surface.png']))

% Create montage
create_figure(name); axis off; o2 = montage(factor_obj{i}, 'nolegend');
title_montage(o2, 5, name)

% Save
% mkdir('figures')
saveas(gcf, fullfile('figures', [name '_montage.png']))

drawnow

end % i


%% Combine them

factor_obj_combined = factor_obj{1};
for i = 2:8
    factor_obj_combined = cat(factor_obj_combined, factor_obj{i}); 
end

%% Plot correlations

names = cellstr(factor_obj_combined.image_names)';
names = cellfun(@(x) strrep(x, '.nii', ''), names, 'UniformOutput', false);

plot_correlation_matrix(factor_obj_combined.dat, 'names', names);

title('Spatial correlations among factor maps')


%% Report neurosynth associations

top_terms = {};

for i = 1:size(ff, 1)

    disp(names{i})

    % Topics
    [~, t] = neurosynth_feature_labels(factor_obj{i}, 'noverbose', 'topics_ri');

    NS_RI_topic_label{i, 1} = t{1}.Term_or_Topic_highest{1};

    [~, terms] = neurosynth_feature_labels(factor_obj{i}, 'noverbose');
    top_terms{i} = terms{1}.Term_or_Topic_highest;

end

imageName = names';

% read original names
orig_names = readtable('factor_names.txt', 'NumHeaderLines', 0);
orig_topic_label = orig_names.FactorName(1:8);

factor_obj_combined.metadata_table = table(imageName, orig_topic_label, NS_RI_topic_label);
factor_obj_combined.metadata_table

%% save

save quah_factor_obj_combined factor_obj_combined
