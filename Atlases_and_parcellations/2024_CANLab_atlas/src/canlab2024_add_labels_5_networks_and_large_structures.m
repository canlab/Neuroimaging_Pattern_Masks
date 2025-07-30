%canlab2024_add_labels_5_networks_and_large_structures
%
% Load the canlab2024 atlas and build a new label field in the canlab2024
% atlas with large structures, including Schafer/Yeo resting-state networks for cortical regions.
% add this as labels_5, and save the atlas_obj back in the associated file
% for future use.

atlas_obj = load_atlas('canlab2024');
yeo = load_atlas('yeo17networks');

[matchingLabels, diceMatrix, maxdice, has_match, atlas_subset_match] = match_atlas_labels(atlas_obj, yeo, 'dice_threshold', 0.01);

% Make a contingency table of the match between labels and cortical regions
contingencyTable = @(vec1, vec2)  [sum(vec1 & vec2), sum(vec1 & ~vec2); sum(~vec1 & vec2), sum(~vec1 & ~vec2)];
logicalVector = contains(atlas_obj.labels, 'Ctx')';
contingencyTable(logicalVector, has_match)

% mask out matches for non-cortical regions
matchingLabels(~logicalVector) = {'nomatch'};

% fill in missing labels from other large-scale structures
for labels = {'Thal' 'Cblm' 'BG' 'hypothalamus' 'AMY' 'BStem' 'MTL_Hipp'}
    logicalVector = contains(atlas_obj.labels, labels{1})';
    matchingLabels(logicalVector) = {labels{1}};
end

% check for any non-matched labels
atlas_obj.labels(contains(matchingLabels, 'nomatch'))'

% assign matchingLabels to labels_5 in atlas_obj 
atlas_obj.labels_5 = matchingLabels';

% save('CANLab2024_MNI152NLin2009cAsym_coarse_2mm_atlas_object', '-append', 'atlas_obj')
