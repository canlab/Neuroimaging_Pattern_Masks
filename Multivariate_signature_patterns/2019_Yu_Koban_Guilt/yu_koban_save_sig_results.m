weight_map_obj = statistic_image('Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii');
p = fmri_data('Yu_SXPO_SXPX_emotion_Mask_pVal.nii');

% The results below might be off if P-vals are not accurate/resolution
% saved in is too low...  The distribution does not look reasonable...

% use z-vals instead
z = fmri_data('Yu_SXPO_SXPX_emotion_Mask_Z.nii');
pvals = [normcdf(double(z.dat)) normcdf(double(z.dat), 'upper')];
pvals = 2 * min(pvals, [], 2); % 2-tailed

weight_map_obj.p = pvals; %p.dat;


%% Summarize results

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

figure;
montage(weight_map_obj)
saveas(gcf, 'Yu_Koban_Guilt_montage.png');

weight_map_obj = threshold(weight_map_obj, .01, 'unc'); % 2-tailed. spm = 1tailed
weight_map_obj.fullpath = fullfile(pwd, 'Yu_Koban_Guilt_01unc_thresh.nii');
write(weight_map_obj)

figure;
montage(weight_map_obj)
saveas(gcf, 'Yu_Koban_Guilt_montage_01unc_thresh.png');

%%
diary Yu_Koban_Guilt_01unc_region_table.txt  % save table output
printhdr('Yu_Koban_Guilt_01unc: Guilt-related, p<.01 unc'); % header

r = region(weight_map_obj);               % create regions

[rpos rneg] = table(r);                   % attach labels
region_obj = [rpos rneg];                 % save regions with labels
diary off

save Yu_Koban_Guilt_01unc_regions region_obj

montage(region_obj, 'regioncenters', 'colormap');
saveas(gcf, 'Yu_Koban_Guilt_01unc_region_montage.png');

