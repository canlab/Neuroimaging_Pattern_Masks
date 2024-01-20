% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'CIT168_MNI152NLin2009cAsym_amygdala_v1.1.0';
space_description = 'MNI152NLin2009cAsym';
references = 'Tyszka, J. M. & Pauli, W. M. In vivo delineation of subdivisions of the human amygdaloid complex in a high-resolution group template. Hum. Brain Mapp. 37, 3979–3998 (2016).';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the full probabilistic atlas file:
parcellation_file = which('CIT168_pAmyNuc_1mm_MNI_to_MNI152NLin2009cAsym.nii.gz');
gunzip(parcellation_file)
parcellation_file = which('CIT168_pAmyNuc_1mm_MNI_to_MNI152NLin2009cAsym.nii');

% Get labels
% -----------------------------------------------------------------------

labels = {'AMY_BLN_La','AMY_BLN_BL_BLDI','AMY_BLN_BM','AMY_CEN','AMY_CMN',...
    'AMY_BL_BLV','AMY_ATA','AMY_ATA_ASTA','AMY_AAA','AMY'};
labels_2 = {'BL','BL','BL','CE','CM','BL','CM','CE','CE','CE'};
labels_3 = repmat({'Amygdala'},1,length(labels));
label_descriptions = {'Amygdalar basolateral complex: lateral nucleus',...
    'Amygdalar basolateral complex: basal nucleus' ...
    'Amygdalar basolateral complex: accessory basal nucleus',...
    'Amygdalar central nucleus',...
    'Amygdalar corticomedial group: medial nucleus',...
    'Amygdalar basolateral complex: paralaminar nucleus',...
    'Amygdalar corticomedial group: periamygdaloid cortex',...
    'Amygdalostriatal triansition area',...
    'Anteroir amygdaloid area',...
    'Amygdalar intercalated nuclei'};
% Create object
% -----------------------------------------------------------------------

pmap = fmri_data(parcellation_file);
% get rid of negatives, probably interpolation errors
pmap.dat(pmap.dat < 0) = 0;
total_p = sum(pmap.dat,2);
renorm = total_p > 1;
pmap.dat(renorm,:) = pmap.dat(renorm,:)./total_p(renorm);

atlas_obj = atlas(pmap, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'labels_2',labels_2,...
    'labels_3',labels_3,...
    'label_descriptions',label_descriptions,...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% Process object
% -----------------------------------------------------------------------

% Threshold at probability 0.2 or greater and k = 3 voxels or greater
atlas_obj = threshold(atlas_obj, .05, 'k', 3);

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique','overlay',which('fmriprep20_template.nii.gz'));
 
% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1,'overlay',which('fmriprep20_template.nii.gz'));
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2);
    
    savedir = fullfile(pwd, 'png_images');
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    
    scn_export_papersetup(600);
    savename = fullfile(savedir, sprintf('%s_montage.png', atlas_name));
    saveas(gcf, savename);

    
end
 
%% save object

if dosave
    
    savename = sprintf('%s_atlas_object.mat', atlas_name);
    save(savename, 'atlas_obj');
    
end

%% write - this writes only the label image

if dosave
    
    savename = sprintf('%s_atlas_regions.nii', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj, 'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    arrayfun(@(x1)set(x1,'FaceAlpha', .5), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end