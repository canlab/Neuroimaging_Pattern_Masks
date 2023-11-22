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

% this is not cross validation, it's not an estimate of how similar
% the glasser atlas will be to atlases derived from future studies, it's a
% comparison of how homogeneous the studies we used to generate the atlas
% are with respect to one another.

addpath('/dartfs-hpc/rc/home/m/f0042vm/software/spm12')
addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/canlab/CanlabCore/'))

SPACE = 'MNI152NLin2009cAsym';

lh_lbls = readtable('../lctx_labels.txt');
rh_lbls = readtable('../rctx_labels.txt');
lbls = [lh_lbls{:,1}; rh_lbls{:,1}];

glasser = cell(3,1);
studies = {'paingen','bmrk5','spacetop'};
for i = 1:length(glasser)
    study = studies{i}
    glasser{i} = atlas(sprintf('../%s_%s_atlas.nii.gz', study, SPACE));
    glasser{i}.labels = lbls';
end

glasser{4} = glasser{2};
glasser{5} = glasser{3}.resample_space(glasser{2}, 'nearest');
assert(all(ismember(glasser{4}.labels, glasser{5}.labels)))

glasser{2} = glasser{2}.resample_space(glasser{1}, 'nearest');
assert(all(ismember(glasser{1}.labels, glasser{2}.labels)))

glasser{3} = glasser{3}.resample_space(glasser{1}, 'nearest');
assert(all(ismember(glasser{1}.labels, glasser{3}.labels)))

labels = glasser{1}.labels;
dicecoef = zeros(size(labels,2),3);
parfor i = 1:length(labels)
    spacetop_roi = glasser{3}.select_atlas_subset(labels(i));
    bmrk5_roi = glasser{2}.select_atlas_subset(labels(i));
    paingen_roi = glasser{1}.select_atlas_subset(labels(i));
    pg_bk = dice(logical(bmrk5_roi.dat), logical(paingen_roi.dat));
    pg_st = dice(logical(spacetop_roi.dat), logical(paingen_roi.dat));

    bmrk5_roi =  glasser{4}.select_atlas_subset(labels(i));
    spacetop_roi = glasser{5}.select_atlas_subset(labels(i));
    bk_st = dice(logical(spacetop_roi.dat), logical(bmrk5_roi.dat));

    dicecoef(i,:) = [pg_bk, pg_st, bk_st];
end

figure(21);
subplot(1,3,1);
hist(dicecoef(:,1),[0.2:0.05:1]);
title({'PainGen vs.','BMRK5 parcels'});
xlabel('Dice Coefficient');
ylabel('Number of parcels')
xlim([0.2,1])

subplot(1,3,2);
hist(dicecoef(:,2),[0.2:0.05:1]);
title({'Paingen vs.','Spactop parcels'});
xlabel('Dice Coefficient');
ylabel('Number of parcels')
xlim([0.2,1])

subplot(1,3,3);
hist(dicecoef(:,3),[0.2:0.05:1]);
title({'BMRK5 vs.','SpaceTop parcels'});
xlabel('Dice Coefficient');
ylabel('Number of parcels')
xlim([0.2,1])

for i = 1:3
    subplot(1,3,i);
    set(gca,'FontSize',13)
end

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),1200,330])
saveas(gcf,'dice_hist_interstudy_glassers_MNI152NLin2009cAsym.png')


fmriprep20 = fmriprep20.remove_empty;
volmap = fmriprep20;
volmap.dat = mean(dicecoef(fmriprep20.dat,:),2);

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_mean_bt_study_glassers_MNI152NLin2009cAsym.png');

%% cross validated out of study dice coefficients
% this tells us how similar the glasser atlas is likely to be with a random
% new study
fmriprep20 = load_atlas('glasser_fmriprep20').threshold(0.2);

SPACE = 'MNI152NLin2009cAsym';

lh_lbls = readtable('../lctx_labels.txt');
rh_lbls = readtable('../rctx_labels.txt');
lbls = [lh_lbls{:,1}; rh_lbls{:,1}];

[glasser, glasser_wh] = deal(cell(3,1));
studies = {'paingen','bmrk5','spacetop'};
for i = 1:length(glasser_wh)
    study = studies{i};
    glasser_wh{i} = atlas(sprintf('../no_%s_%s_atlas.nii.gz', study, SPACE));
    glasser_wh{i}.labels = lbls';
    glasser{i} = atlas(sprintf('../%s_%s_atlas.nii.gz', study, SPACE));
    glasser{i}.labels = lbls';
end

labels = glasser{1}.labels;
dicecoef = zeros(size(labels,2),3);
for k = 1:3
    parfor i = 1:length(labels)
        roi_test = glasser{k}.select_atlas_subset(labels(i));
        roi_train = glasser_wh{k}.select_atlas_subset(labels(i));
        dicecoef(i,k) = dice(logical(roi_test.dat), logical(roi_train.dat));
    end
end


figure(31);
hist(mean(dicecoef,2),[0.2:0.05:1]);
title({'Study similarity to template','Out-of-study 3-fold CV Est.'});
xlabel('Dice Coefficient');
ylabel('Number of parcels')
xlim([0.2,1])

set(gca,'FontSize',13)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),400,330])
saveas(gcf,'dice_hist_out_of_study_CV_glassers_MNI152NLin2009cAsym.png')


fmriprep20 = fmriprep20.remove_empty;
volmap = fmriprep20;
volmap.dat = mean(dicecoef(fmriprep20.dat,:),2);

volmap = fmri_data(volmap);

figure;
volmap.montage('full2')

saveas(gcf,'dice_map_mean_out_of_study_CV_glassers_MNI152NLin2009cAsym.png');