% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'Pauli2016_striatum';
space_description = 'MNI152 space';
references = 'Pauli, Wolfgang M., Randall C. O?Reilly, Tal Yarkoni, and Tor D. Wager. 2016. ?Regional Specialization within the Human Striatum for Diverse Psychological Functions.? Proceedings of the National Academy of Sciences of the United States of America 113 (7): 1907?12.';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space

% parcellation_file = which('Pauli_bg_cluster_mask_5_p.nii');
% 
% obj = fmri_data(parcellation_file);
% obj.dat = obj.dat ./ 100; % convert to p from percent
% obj.fullpath = fullfile(pwd, 'Pauli_bg_cluster_mask_5_prob.img');
% write(obj);

parcellation_file = which('Pauli_bg_cluster_mask_5_prob.img');

cd(fileparts(parcellation_file))

%% Get labels
% -----------------------------------------------------------------------

labels = {'Caudate_Cp' 'Putamen_Pa' 'Caudate_Ca' 'V_Striatum' 'Putamen_Pp'};
labels_2 = repmat({'Striatum'}, 1, 5);

% Create object
% -----------------------------------------------------------------------

atlas_obj = atlas(parcellation_file, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'labels_2', labels_2, ...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% Process object
% -----------------------------------------------------------------------

% Threshold: for probability images only
% Threshold at probability 0.2 or greater and k = 1 voxels or greater
atlas_obj = threshold(atlas_obj, .2, 'k', 1);

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique');
 
% Convert to regions
% -----------------------------------------------------------------------

 r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);


%% Subdivide into L and R hemispheres, save separate copy

atlas_obj_split_lr = split_atlas_by_hemisphere(atlas_obj); 

 
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
    save(savename, 'atlas_obj', 'atlas_obj_split_lr');
    
end

%% write - this writes only the label image

% Save: if creating from probabilty images only
%
% if dosave
%     
%     savename = sprintf('%s_atlas_regions.img', atlas_name);
%     atlas_obj.fullpath = fullfile(pwd, savename);
%     write(atlas_obj);
%     
% end

%% Turn regions into separate list of names, for canlab_load_ROI
% which loads regions by name from mat files.

clear region_names

for i = 1:length(r)
    
    eval([labels{i} ' = r(i);']);
    
    region_names{i} = r(i).shorttitle;
    
end

savename = sprintf('%s_atlas_regions.mat', atlas_name);
save(savename, 'r', 'region_names', labels{:});

%% Isosurface

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