addpath('/home/bogdan/.matlab/spm/spm12');
addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore'))
addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks'))

space_description = 'MNI152NLin2009cAsym';
atlas_name = sprintf('iglesias_hcp278_%s', space_description);
references = char({'Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K, Greve DN, van der Kouwe A, Fischl B, Caballero-Gaudes C, Paz-Alonso PM. (2018). A probablistic atlas of the human thalamuc nuclei combining ex vivo MRI and histology. Neuroimage, 314-326, 183.'; ...
             'Tregidgo HFJ, Soskic S, Althonayan J, Maffei C, Van Leemput K, Golland P, Insausti R, Lerma-Usabiaga G, Caballero-Gaudes C, Paz-Alonso PM, Yendiki A, Alexander DC, Bocchetta M, Rohrer JD, Iglesias JE. (2023). Accurate Bayesian segmentation of thalamic nuclei using diffusion MRI and an improved histological atlas. Neuroimage,  274, 120129.'});

dosave = true;

% Note this file was obtained by first running the MNI152NLin6Asym script
% and then saving the probability maps, transforming them to
% MNI152NLin2009cAsym space and saving the results for use in this script.
parcellation_file = 'iglesias_MNI152NLin2009cAsym_probability_maps.nii.gz';


tbl = readtable(which('subnuclei_labels.csv'));
labels = cellfun(@(x1)strrep(x1,'-','_'),tbl.Var2,'UniformOutput',false);
labels = labels(:)';
ind = tbl.Var1;

extra_lbls = readtable('thalamic_labels.csv');

labels_2 = {};
labels_3 = {};
labels_4 = {};
label_descriptions = {};

for i = 1:length(labels)
    lbl_ind = strcmp(extra_lbls.labels, regexprep(labels(i),'[a-zA-Z]*_(.*)','$1'));
    side = regexprep(labels(i),'([a-zA-Z]*)_.*','$1');
    assert(sum(lbl_ind) == 1);
    labels{i} = extra_lbls.labels{lbl_ind};
    label_descriptions{end+1} = extra_lbls.label_descriptions{lbl_ind};
    labels_2{i} = extra_lbls.labels_2{lbl_ind};
    labels_3{i} = extra_lbls.labels_3{lbl_ind};
    labels_4{i} = extra_lbls.labels_4{lbl_ind};

    switch side{1}
        case 'Left'
            pfx = 'L_';
            sfx = ' (left)';
        case 'Right'
            pfx = 'R_';
            sfx = ' (right)';
        otherwise
            error('Could not determine laterality');
    end

    labels{i} = [pfx, labels{i}];
    labels_2{end} = strrep([pfx, labels_2{end}],' ','_');
    labels_3{end} = strrep([pfx, labels_3{end}], ' ','_');
    labels_4{end} = strrep([pfx, labels_4{end}], ' ', '_');
    label_descriptions{end} = strrep([label_descriptions{end}, sfx], ' ', '_');
end



parcellation = fmri_data(parcellation_file);

fid= fopen(which('iglesias_MNI152NLin6Asym_probability_maps.labels'));
labels2 = textscan(fid,'%s');
labels2 = labels2{1};
fclose(fid)

n_regions = size(parcellation.dat,2);
assert(length(labels2) == n_regions)

label_ind = [];
for i = 1:n_regions
    label_ind(end+1) = find(strcmp(labels2(i), labels));
end

atlas_obj = atlas(parcellation, ...
    'labels',labels(label_ind), ...
    'label_descriptions', label_descriptions(label_ind), ...
    'labels_2', labels_2(label_ind), ...
    'labels_3', labels_3(label_ind), ...
    'labels_4', labels_4(label_ind), ...
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

labels = labels(label_ind);
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
