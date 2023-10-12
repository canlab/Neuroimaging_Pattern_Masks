% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'tian_3t_s3_fmriprep20';
space_description = 'MNI152NLin2009cAsym';
references = 'Tian Y, Margulies D, Breakspear M, Zalesky A (2020). Nature Neuroscience. 23(11) 1421-1432.';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% Authors only provide 1mm parcellations in this space together with
% cortical parcellations, so we'll need to take this and truncate it.
parcellation_file = which('Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S3_3T_MNI152NLin2009cAsym_1mm.nii.gz');


% Get labels
% -----------------------------------------------------------------------
fid = fopen(which('Tian_Subcortex_S3_3T_label.txt'));
labels = textscan(fid,'%s');
labels = cellfun(@(x1)strrep(x1,'-','_'),labels{1},'UniformOutput',false);
fclose(fid);

% truncate cortical areas
parcellation = fmri_data(parcellation_file);
parcellation.dat(parcellation.dat > length(labels)) = 0;

% Create object
% -----------------------------------------------------------------------
atlas_obj = atlas(parcellation, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% I've emailed Tian and Zalesky asking about probability maps, they have
% them in figure 7c of their paper, but for now we default to a binary map.
uniq_id = unique(atlas_obj.remove_empty.dat);
n_roi = length(uniq_id);
pmap = zeros(atlas_obj.volInfo.n_inmask, n_roi);
atlas_obj = atlas_obj.replace_empty();
for i = 1:length(uniq_id)
    pmap(:,i) = atlas_obj.dat == uniq_id(i);
end
atlas_obj.probability_maps = pmap;


% Process object
% -----------------------------------------------------------------------

% Threshold at probability 0.2 or greater and k = 3 voxels or greateratlas_obj = threshold(atlas_obj, 0.2, 'k', 3);
atlas_obj = threshold(atlas_obj, .2, 'k', 3);


% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
atlas_obj.orthviews('unique',which('MNI152NLin2009cAsym_T1_1mm.nii.gz'));


% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
 %% save figure

if dosave
   
    o2 = canlab_results_fmridisplay([], 'multirow', 1, 'overlay', which('MNI152NLin2009cAsym_T1_1mm.nii.gz'));
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
    write(atlas_obj,'overwrite');
    
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
    
    cellfun(@(x1)(set(x1,'FaceAlpha', .5)), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end
