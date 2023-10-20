% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'Shen_MNI152NLin2009Asym';
space_description = 'MNI152NLin2009cAsym';
references = 'Shen 2013 Neuroimage';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the index image only file:
parcellation_file = which('shen_1mm_268_parcellation_MNI152NLin2009cAsym.nii.gz');

cd(fileparts(parcellation_file))

% Get labels
% -----------------------------------------------------------------------

labeldata = importdata(which('shen_268_parcellation_networklabels.csv'))';
labels = {};
% no labels, assign some, with network values
for i = 1:size(labeldata.data, 1)
    
    labels{i} = sprintf('R_%d_Network_%d', i, labeldata.data(i, 2));
    labels_2{i} = sprintf('Network_%d', labeldata.data(i, 2));
    
end

labels = format_text_letters_only(labels, 'numbers', 'cleanup'); % Replace chars we don't want


% Create object
% -----------------------------------------------------------------------

atlas_obj = atlas(which(parcellation_file), ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'labels_2', labels_2, ...
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
orthviews(atlas_obj, 'unique','overlay',which('fmriprep20_template.nii.gz'));
 
% Convert to regions
% -----------------------------------------------------------------------

 r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full2', 'overlay', which('fmriprep20_template.nii.gz'));
    brighten(.6)
    
    o2 = montage(r, o2);
    
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

%%
if dosave
    
    figure; han = isosurface(atlas_obj);
    
    cellfun(@(x1)(set(x1,'FaceAlpha', .5)),han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end