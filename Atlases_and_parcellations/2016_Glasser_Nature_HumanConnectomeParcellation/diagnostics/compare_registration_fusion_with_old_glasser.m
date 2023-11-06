old = load_atlas('glasser');
new = load_atlas('glasser_fmriprep20');


new = new.resample_space(old,'nearest');
new_reduced = new.threshold(0.2);
new_reduced.probability_maps = []; % otherwise this regenerates regions when we extract them

new_reduced = new_reduced.replace_empty();
old = old.replace_empty();


assert(all(ismember(old.labels, new_reduced.labels)))
labels = old.labels;
dicecoef = zeros(size(labels));
parfor i = 1:length(labels)
    old_roi = old.select_atlas_subset(labels(i));
    new_roi = new_reduced.select_atlas_subset(labels(i));
    dicecoef(i) = dice(old_roi.dat, new_roi.dat);
end

figure(1)
hist(dicecoef);
title({'Old (Horn 2016) Glasser vs.','New Registration Fusion Glasser'})
xlabel('Dice Coefficient ')
ylabel('Number of Parcels')
xlim([0.2,1])
set(gca,'FontSize',13)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),388,330])
saveas(gcf,'dice_coef_old_vs_new_glasser.png')

old = old.remove_empty;
volmap = old;
volmap.dat = dicecoef(old.dat)';

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_old_vs_new_glasser.png')

% old atlas is much thinner than new one, so let's mask to have a fair
% comparison
new_masked = new.apply_mask(fmri_mask_image(old));
new_masked = new_masked.threshold(0.2);
new_masked.probability_maps = []; % otherwise this regenerates regions when we extract them

new_masked = new_masked.replace_empty();
old = old.replace_empty();


assert(all(ismember(old.labels, new_masked.labels)))
labels = old.labels;
dicecoef = zeros(size(labels));
parfor i = 1:length(labels)
    old_roi = old.select_atlas_subset(labels(i));
    new_roi = new_masked.select_atlas_subset(labels(i));
    dicecoef(i) = dice(old_roi.dat, new_roi.dat);
end

figure(11);
hist(dicecoef);
title({'Old (Horn 2016) Glasser vs.','New registration fusion Glasser','after masking latter by former'})
xlabel('Dice Coefficient');
set(gca,'FontSize',13)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),388,330])
ylabel('Number of parcels')
saveas(gcf,'dice_hist_old_vs_new_glasser_masked_by_old.png')

old = old.remove_empty;
volmap = old;
volmap.dat = dicecoef(old.dat)';

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_old_vs_new_glasser_masked_by_old.png');


%% let's also compare the two different glasser spaces in the new atlases to get a sense of perspective for the above
fsl6 = load_atlas('glasser_fsl6');
fmriprep20 = load_atlas('glasser_fmriprep20');

fsl6 = fsl6.resample_space(fmriprep20,'nearest');
fsl6 = fsl6.threshold(0.2);
fmriprep20 = fmriprep20.threshold(0.2);
fsl6.probability_maps = [];
fmriprep20.probability_maps = [];

assert(all(ismember(fsl6.labels, fmriprep20.labels)))
labels = fsl6.labels;
dicecoef = zeros(size(labels));
parfor i = 1:length(labels)
    fsl_roi = fsl6.select_atlas_subset(labels(i));
    fmriprep_roi = fmriprep20.select_atlas_subset(labels(i));
    dicecoef(i) = dice(logical(fsl_roi.dat), logical(fmriprep_roi.dat));
end


figure(21);
hist(dicecoef);
title({'MNI152NLin2009cAsym (fmriprep) vs.','MNI152NLin6Asym (fsl)'})
xlabel('Dice Coefficient');
ylabel('Number of parcels');
xlim([0.2,1])
set(gca,'FontSize',13)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),388,330])
saveas(gcf,'dice_hist_fmriprep_vs_fsl_glassers.png')

fmriprep20 = fmriprep20.remove_empty;
volmap = fmriprep20;
volmap.dat = dicecoef(fmriprep20.dat)';

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_fmriprep_vs_fsl_glassers.png');

%% let's also compare between studies within template

addpath('/dartfs-hpc/rc/home/m/f0042vm/software/spm12')
addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/CanlabCore/'))

SPACE = 'MNI152NLin2009cAsym';

lh_lbls = readtable('../lctx_labels.txt');
rh_lbls = readtable('../rctx_labels.txt');
lbls = [lh_lbls{:,1}; rh_lbls{:,1}];

glasser = cell(2,1);
studies = {'paingen','bmrk5'};
for i = 1:length(glasser)
    study = studies{i}
    glasser{i} = atlas(sprintf('../%s_%s_atlas.nii.gz', study, SPACE));
    glasser{i}.labels = lbls';
end

% save individual atlases for inspection
glasser{2} = glasser{2}.resample_space(glasser{1}, 'nearest');


assert(all(ismember(glasser{1}.labels, glasser{2}.labels)))
labels = glasser{1}.labels;
dicecoef = zeros(size(labels));
parfor i = 1:length(labels)
    bmrk5_roi = glasser{1}.select_atlas_subset(labels(i));
    paingen_roi = glasser{2}.select_atlas_subset(labels(i));
    dicecoef(i) = dice(logical(bmrk5_roi.dat), logical(paingen_roi.dat));
end


figure(21);
hist(dicecoef);
title({'Paingen vs.','Bmrk5 parcels'});
xlabel('Dice Coefficient');
ylabel('Number of parcels')
xlim([0.2,1])
set(gca,'FontSize',13)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),388,330])
saveas(gcf,'dice_hist_paingen_vs_bmrk5_glassers_MNI152NLin2009cAsym.png')

fmriprep20 = fmriprep20.remove_empty;
volmap = fmriprep20;
volmap.dat = dicecoef(fmriprep20.dat)';

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_paingen_vs_bmrk5_glassers_MNI152NLin2009cAsym.png');