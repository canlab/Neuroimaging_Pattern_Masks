addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

space_description = 'MNI152NLin2009cAsym';
atlas_name = sprintf('levinson_bari_limbic_brainstem_atlas_%s', space_description);
references = char({'Levinson S, Miller M, Iftekhar A, Justo M, Arriola D, Wei W, Hazany S, Avecillas-Chasin J, Kuhn T, Horn A, Bari A. (2023). A structural connectivity atlas of limbic brainstem nuclei. Frontiers in Neuroimaging.'});

dosave = true;

% there are two available parcellations. One based on T1 and T2 contrast,
% the other using FA and T1 contrast. The former has more regions, but is
% not as well segmented. I particular it overlaps with the coricospinal
% tracts.
parcellation_files = arrayfun(@(x1)[x1.folder, '/', x1.name], dir('src_to_MNI152NLin2009cAsym/*nii.gz'), 'UniformOutput', false);


tbl = readtable(which('levinson_bari_limbic_brainstem_atlas_labels.csv'));
labels = cellfun(@(x1)strrep(x1,'-','_'),tbl.labels,'UniformOutput',false);
labels = labels(:)';
label_descriptions = tbl.label_descriptions(:)';

parcellation = fmri_data(parcellation_files);

n_regions = length(parcellation_files);

new_labels = {};
new_descriptions = {};
for i = 1:n_regions
    this_label = regexprep(parcellation_files{i},'.*[0-9]*_([A-Za-z]*)_(ATLAS_2022a).*','$1');
    this_description = label_descriptions{contains(labels, this_label)};

    new_labels{end+1} = this_label;
    new_descriptions{end+1} = this_description;
end

% this is natively in 0.5mm space. At some future date it may make sense to
% skip downsampling (but make sure to update the underlay in the
% visualizations below or this script will break), but for now we don't
% need the storage overhead of that.
parcellation = parcellation.resample_space(which([space_description '_T1_1mm.nii.gz']));

atlas_obj = atlas(parcellation, ...
    'labels',new_labels, ...
    'label_descriptions', new_descriptions(:), ...
    'atlas_name', atlas_name ,...
    'space_description', space_description, ...
    'references', references);

bilateral_regions = atlas_obj.select_atlas_subset({'PAG','DR'});

% the next line borrows from a canlab2023 function used to subdivide
% regions around the x=0 plane
lateralized_regions = lateralize(atlas_obj.select_atlas_subset({'NTS','LC','VTA'}));

atlas_obj = [bilateral_regions, lateralized_regions];

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

cmap = scn_standard_colors(num_regions(atlas_obj));
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
