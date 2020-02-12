cd('/Users/tor/Documents/Code_Repositories/MasksPrivate/Masks_private/Wager_Kang_etal_Emotion_Meta_2015_BSPP')
dat = fmri_data(filenames('Wager*nii'))


cd('/Users/tor/Documents/Code_Repositories/MasksPrivate/Masks_private/Kragel_emotionClassificationBPLS_2015')
dat2 = fmri_data(filenames('mean*img'))
dat2 = resample_space(dat2, dat);


wh_exclude = [mean(dat.dat, 2) == 0 mean(dat2.dat, 2) == 0];
sum(wh_exclude)
datavals = [dat.dat dat2.dat];
datavals(any(wh_exclude, 2),  :) = []; % exclude values missing in either set


%% names


names = cellstr(dat.image_names);
names = strrep(names, 'Wager_Kang_PlosCB_emometa_2015_', '');

names2 = cellstr(dat2.image_names);
names2 = strrep(names2, 'mean_3comp_', '');
names2 = strrep(names2, 'group_emotion_PLS_beta_BSz_10000it.img', '');

names = strrep(names, '.nii', '');

names2 = strrep(names2, '_', '');

%% plot

r = corr(datavals)

figure; imagesc(r); colorbar

set(gca, 'XTick', [1:12], 'YTick', 1:12, 'XTickLabel', [names' names2'], 'YTickLabel', [names' names2'])

%%


r = r(1:5, 7:end)

figure; imagesc(r); colorbar


set(gca, 'XTick', [1:7], 'YTick', 1:5, 'XTickLabel', [names2'], 'YTickLabel', [names'], 'FontSize', 24)
xlabel('Kragel');
ylabel('Wager meta');

scn_export_papersetup(400); saveas(gcf, 'kragel_wager_cross_corr.png');

