addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

space_description = 'MNI152NLin6Asym';
atlas_name = sprintf('harvard_aan_v2_%s', space_description);
references = char({'Edlow, BL; Kinney, HC (2023a), Harvard Ascending Arousal Network Atlas – Version 2.0, Dryad Digital Repository, https://doi.org/10.5061/dryad.zw3r228d2'});

dosave = true;

% there are two available parcellations. One based on T1 and T2 contrast,
% the other using FA and T1 contrast. The former has more regions, but is
% not as well segmented. I particular it overlaps with the coricospinal
% tracts.
parcellation_files = arrayfun(@(x1)[x1.folder, '/', x1.name], dir('src/*nii.gz'), 'UniformOutput', false);


tbl = readtable(which('harvard_aan_labels.csv'));
labels = cellfun(@(x1)strrep(x1,'-','_'),tbl.labels,'UniformOutput',false);
labels = labels(:)';
label_descriptions = tbl.label_descriptions(:)';

parcellation = fmri_data(parcellation_files);

n_regions = length(parcellation_files);
lateralized_maps = find(contains(parcellation_files, {'_L_','_R_'}));
nonlateralized_maps = find(contains(parcellation_files, {'DR','MnR','PAG','VTA'}));

pmap = zeros(size(parcellation.dat,1), n_regions);
new_labels = {};
new_descriptions = {};
for i = 1:length(lateralized_maps)
    map_ind = lateralized_maps(i);
    region = regexprep(parcellation_files{map_ind},'.*AAN_([A-Za-z]*)_([LR])_MNI.*','$1');
    h = regexprep(parcellation_files{map_ind},'.*AAN_[A-Za-z]*_([LR])_MNI.*','$1');
    this_label = [h, '_', region];
    switch h
        case 'L'
            suffix = '(left)';
        case 'R'
            suffix = '(right)';
    end
    this_description = [label_descriptions{contains(labels, region)}, ' ', suffix];

    pmap(:,i) = parcellation.dat(:,map_ind);
    new_labels{end+1} = this_label;
    new_descriptions{end+1} = this_description;
end

for i = 1:length(nonlateralized_maps)
    map_ind = nonlateralized_maps(i);
    this_label = regexprep(parcellation_files{map_ind},'.*AAN_([A-Za-z]*)_.*','$1');
    
    this_description = [label_descriptions{contains(labels, this_label)}, ' ', suffix];

    pmap(:,i+length(lateralized_maps)) = parcellation.dat(:,map_ind);
    new_labels{end+1} = this_label;
    new_descriptions{end+1} = this_description;
end
parcellation.dat = pmap;

% there's a bug with the PTg_L file, and it's actually bilateral, so let's
% subtract it out
parcellation.dat(:,contains(new_labels,'L_PTg')) = parcellation.dat(:,contains(new_labels,'L_PTg')) - parcellation.dat(:,contains(new_labels,'R_PTg'));

atlas_obj = atlas(parcellation, ...
    'labels',new_labels, ...
    'label_descriptions', new_descriptions(:), ...
    'atlas_name', atlas_name ,...
    'space_description', space_description, ...
    'references', references);


% Process object
% -----------------------------------------------------------------------


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
labels = atlas_obj.labels;
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
