

load('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/neurosynth/neurosynth_data_obj.mat')
load(which("neurosynth_data_obj.mat")) % github/canlab/Neuroimaging_Pattern_Masks/Neurosynth_maps/neurosynth_data_obj.mat
neurosynth_data % !
m = mean(neurosynth_data);

%%

create_figure('overall p activation montage'); axis off;
o2 = montage(m);
% Mean across 27072 contrast maps with diverse functional tasks, from Neurosynth.org, downloaded April 2022 (latest version)
toc
figure;
subplot(1, 2, 1); surface(m, 'left_insula_slab', 'colormap', 'summer');
subplot(1, 2, 2); surface(m, 'right_insula_slab');

%%
m = fmri_data('Activation_proportion.img');
d = descriptives(m);
m = threshold(m, [d.prctile_vals(8) Inf], 'raw-between'); % 95th percentile

create_figure('overall p activation montage'); axis off; h = addbrain('insula surfaces');
% surface(m, 'surface_handles', h, 'colormap', 'hsv');

%%
cm = colormap_tor([0 0 1], [1 1 0], [0 1 1], [1 0 1], [1 0 0]);
create_figure('overall p activation montage'); axis off; h = addbrain('insula surfaces');
surface(m, 'surface_handles', h, 'colormap', cm);

%%
m = fmri_data('Activation_proportion.img');
d = descriptives(m);
m = threshold(m, [d.prctile_vals(9) Inf], 'raw-between'); % 99th percentile

%% by parcel


%% neurosynth topic count
% each topic map is FDR corrected q < 0.05

% Load 100-topic Neurosynth maps (forward + reverse inference)
% load('/Users/f003vz1/Documents/GitHub/Neuroimaging_Pattern_Masks/Atlases_and_parcellations/2016_Neurosynth_100_topics/neurosynth_topics_v4.mat')
% canlab/Neuroimaging_Pattern_Masks/Neurosynth_maps/2016_Neurosynth_100_topics_Oldfiles/neurosynth_topics_v4.mat
load(which("neurosynth_topics_v4.mat")) 

%% forward inference
% The forward‑inference topic maps are used to examine where each topic 
% tends to activate (diversity of activation across topics).

% Start from the mean of all forward-inference topic maps
% (this produces a fmri_data object with average forward-inference strength)
count_obj = mean(topic_obj_forwardinference);
% Now *replace* the data (.dat) with a "topic count" measure:
%   - For each voxel, count how many topics have a positive value (dat > 0)
%   - Then divide by the total number of topics to get a fraction
%   - So count_obj.dat(i) = fraction of topics that activate voxel i
count_obj.dat = sum(topic_obj_forwardinference.dat > 0, 2) ./ size(topic_obj_forwardinference.dat, 2);
figure; montage(count_obj)

% Get descriptive stats on the fraction-of-topics map
d = descriptives(count_obj);
% Threshold at the 99th percentile:
% keep voxels that are activated (forward inference) by *many* topics
count_obj_thr = threshold(count_obj, [d.prctile_vals(9) Inf], 'raw-between'); % 99th percentile
% View thresholded map as slices and surfaces
figure; montage(count_obj_thr)
figure; surface(count_obj_thr, 'foursurfaces');

%% reverse inference
% Areas activated by everything should tend not to show these
% associations as strongly (i.e., reverse inference should be weaker there)

count_obj = mean(topic_obj_reverseinference);
count_obj.dat = sum(topic_obj_reverseinference.dat > 0, 2) ./ size(topic_obj_reverseinference.dat, 2);
figure; montage(count_obj)

d = descriptives(count_obj);
count_obj_thr = threshold(count_obj, [d.prctile_vals(8) Inf], 'raw-between'); % 95th percentile
figure; montage(count_obj_thr)
figure; surface(count_obj_thr, 'foursurfaces');