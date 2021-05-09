% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'delaVega2017_neurosynth';
space_description = 'MNI152 space';
references = 'De La Vega et al. 2017, Cerebral Cortex';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the full probabilistic atlas file:
parcellation_file = which('wb_70.nii');  

cd(fileparts(parcellation_file))

% load and fix rounding issues
delavega = fmri_data(parcellation_file);
delavega.dat = round(delavega.dat);


% Get labels
% -----------------------------------------------------------------------
labels = {};
for i = 1:70
    labels{i} = sprintf('R%d', i);
end

% Autolabel with canlab2018 combined atlas
% Best matching region from Glasser 2016, or multiple.
% Best matching rsfMRI network from Schafer/Yeo 2018

r = region(delavega, 'unique_mask_values');
[r, region_table] = r.autolabel_regions_using_atlas;

for i = 1:70
    labels{i} = [region_table.modal_label{i} '_' region_table.modal_label_descriptions{i}];
end

% Clean up names and remove spaces
labels = format_text_letters_only(labels, 'numbers', 'cleanup', 'squeeze', 'underscore_ok');


% Create object
% -----------------------------------------------------------------------


atlas_obj = atlas(delavega, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% Process object
% -----------------------------------------------------------------------

% Note: not probabilistic, so no thresholding
% Threshold at probability 0.2 or greater and k = 3 voxels or greater
%atlas_obj = threshold(atlas_obj, .2, 'k', 3);

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique');
 
% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1);
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
    
    savename = sprintf('%s_atlas_regions.img', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj);
    
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
    
    set(han,'FaceAlpha', .5)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end