% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'Brainnetome';
space_description = 'MNI152 space';
references = 'Fan 2016 Cerebral Cortex; doi:10.1093/cercor/bhw157';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the index image only file:
parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');

cd(fileparts(parcellation_file))

% Get labels
% -----------------------------------------------------------------------

load(fullfile(pwd, 'cluster_names.mat'))

labels = names';
labels = format_text_letters_only(labels, 'numbers', 'cleanup'); % Replace chars we don't want

% Create object
% -----------------------------------------------------------------------

atlas_obj = atlas(which(parcellation_file), ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% Process object
% -----------------------------------------------------------------------

% Threshold: for probability images only
% Threshold at probability 0.2 or greater and k = 3 voxels or greater
% atlas_obj = threshold(atlas_obj, .2, 'k', 3);

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

