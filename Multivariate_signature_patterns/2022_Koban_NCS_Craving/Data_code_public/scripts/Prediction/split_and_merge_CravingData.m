% save as two fmri_data objects due to size limits on Github, merge again
% Leonie Koban, Sept 2022


%% separate into two objects

Craving_5levels_part1 = get_wh_image(Craving_5levels, 1:250);
Craving_5levels_part2 = get_wh_image(Craving_5levels, 251:469);

%% merge again

Craving_5levels = cat(Craving_5levels_part1, Craving_5levels_part2);