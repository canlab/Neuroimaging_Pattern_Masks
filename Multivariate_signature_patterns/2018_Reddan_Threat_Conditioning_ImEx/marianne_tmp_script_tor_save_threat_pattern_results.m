load('/Users/torwager/Google Drive/Projects/2018_Reddan_IE_Imagined_Extinction/Threat_Map/weight maps/IE_ACQmap_forpattexpapplication2015.mat')

w = stats.weight_obj;
w = fmri_data(w);

cd('/Users/torwager/Google Drive/Projects/2018_Reddan_IE_Imagined_Extinction/Threat_Map')

z = statistic_image('acq_Z_wtmap.nii');

% w obj has slightly different mask and space...
w = resample_space(w, z);


%% turn into statistic_image object with weights, z, p
weight_map_obj = z;
weight_map_obj.dat = w.dat;
weight_map_obj.p = p.dat;

% clean up by applying mask
weight_map_obj = apply_mask(weight_map_obj, which('gray_matter_mask.img'));

weight_map_obj.fullpath = fullfile(pwd, 'IE_ImEx_Acq_Threat_SVM_nothresh.nii');
write(weight_map_obj)

weight_map_obj = threshold(weight_map_obj, .01, 'unc');
weight_map_obj.fullpath = fullfile(pwd, 'IE_ImEx_Acq_Threat_SVM_01thresh.nii');
write(weight_map_obj)

weight_map_obj = threshold(weight_map_obj, .05, 'fdr');
weight_map_obj.fullpath = fullfile(pwd, 'IE_ImEx_Acq_Threat_SVM_05FDR.nii');
write(weight_map_obj)

%% Summarize results

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

figure;
montage(weight_map_obj)
saveas(gcf, 'IE_ImEx_Acq_Threat_SVM_05FDR.png');

diary IE_ImEx_Acq_Threat_SVM_05FDR_region_table.txt  % save table output
printhdr('IE_ImEx_Acq_Threat_SVM_05FDR: Threat acquisition, FDR q < .05'); % header
r = region(weight_map_obj);               % create regions
[rpos rneg] = table(r);                   % attach labels
region_obj = [rpos rneg];                 % save regions with labels
diary off

save IE_ImEx_Acq_Threat_SVM_05FDR_regions region_obj

montage(region_obj, 'regioncenters', 'colormap');
saveas(gcf, 'IE_ImEx_Acq_Threat_SVM_05FDR_region_montage.png');

%% Uncorrected summary

weight_map_obj = threshold(weight_map_obj, .01, 'unc');

figure;
montage(weight_map_obj)
saveas(gcf, 'IE_ImEx_Acq_Threat_SVM_01unc.png');

diary IE_ImEx_Acq_Threat_SVM_01unc_region_table.txt  % save table output
printhdr('IE_ImEx_Acq_Threat_SVM_01unc: Threat acquisition, uncorrected p < .01'); % header
r = region(weight_map_obj);               % create regions
r(cat(1, r.numVox) < 30) = [];            % size threshold: 30+ voxels
[rpos rneg] = table(r);                   % attach labels
region_obj = [rpos rneg];                 % save regions with labels
diary off

save IE_ImEx_Acq_Threat_SVM_01unc_regions region_obj

montage(region_obj, 'regioncenters', 'colormap');
saveas(gcf, 'IE_ImEx_Acq_Threat_SVM_01unc_region_montage.png');

%% Save the intercept

diary IE_ImEx_Acq_Threat_SVM_intercept
stats.other_output{3}
diary off

