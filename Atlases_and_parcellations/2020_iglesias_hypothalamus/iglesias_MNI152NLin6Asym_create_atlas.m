addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

space_description = 'MNI152NLin6Asym';
atlas_name = sprintf('iglesias_hypothal_hcp278_%s', space_description);
references = char({'Billot B, Bocchetta M, Todd E, Dalca A, Rohrer J, Iglesias J. (2020). Automated segmentation of the hypothalamus and associated subunits in brain MRI. Neuroimage, 223.'});

dosave = true;

% there are two available parcellations. One based on T1 and T2 contrast,
% the other using FA and T1 contrast. The former has more regions, but is
% not as well segmented. I particular it overlaps with the coricospinal
% tracts.
parcellation_file = 'hypothalamic_subunits_seg.v1.MNI152NLin6Asym.nii.gz';


tbl = readtable(which('atlas_labels.csv'));
labels = cellfun(@(x1)strrep(x1,'-','_'),tbl.labels,'UniformOutput',false);
labels = labels(:)';
labels_2 = cellfun(@(x1)strrep(x1,'-','_'),tbl.labels_2,'UniformOutput',false);
labels_2 = labels_2(:)';
labels_3 = cellfun(@(x1)strrep(x1,'-','_'),tbl.labels_3,'UniformOutput',false);
labels_3 = labels_3(:)';
label_descriptions = tbl.label_descriptions(:)';
ind = tbl.id;

parcellation = fmri_data(parcellation_file);

n_regions = length(unique(parcellation.dat)) - 1;
uniq_rois = unique(parcellation.dat);
uniq_rois(uniq_rois == 0) = [];

pmap = zeros(size(parcellation.dat,1), n_regions);
for i = 1:n_regions
    this_roi = uniq_rois(i);
    this_mask = (parcellation.dat == this_roi);
    pmap(:,i) = mean(this_mask,2);
end
parcellation.dat = pmap;
has_lbl = ismember(tbl.id,uniq_rois);

atlas_obj = atlas(parcellation, ...
    'labels',labels(has_lbl), ...
    'label_descriptions', label_descriptions(has_lbl), ...
    'labels_2', labels_2(has_lbl), ...
    'labels_3', labels_3(has_lbl), ...
    'atlas_name', atlas_name ,...
    'space_description', space_description, ...
    'references', references);


% Process object
% -----------------------------------------------------------------------

% Threshold at probability 0.2 or greater and k = 3 voxels or greateratlas_obj = threshold(atlas_obj, 0.2, 'k', 3);
atlas_obj = atlas_obj.threshold(0.01);


% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
atlas_obj.orthviews('unique',which([space_description '_T1_1mm.nii.gz']));


% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

cmap = scn_standard_colors(n_regions/2);
cmap = cell2mat(cat(2,cmap'));
cmap = [cmap; cmap];

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full hcp', ...
        'overlay', which(sprintf('%s_T1_1mm.nii.gz',space_description)));
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2, 'indexmap', cmap);
    
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
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj,'overwrite');
    
end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

labels = labels(has_lbl);
for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    arrayfun(@(x1)(set(x1,'FaceAlpha', .5)), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
