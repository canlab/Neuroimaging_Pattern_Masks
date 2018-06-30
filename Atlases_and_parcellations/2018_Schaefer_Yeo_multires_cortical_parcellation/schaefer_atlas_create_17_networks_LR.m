% Collapse 400 regions into 17 networks separated by L/R hemisphere
% Save new atlas object

atlas_name = 'Schaefer2018Cortex_17networks';

netnames = unique(atlas_obj.labels_3);  % 17 networks

new_atlas = atlas_obj;
new_atlas.atlas_name = atlas_name;

new_vals = zeros(size(new_atlas.dat));
new_labels = {};

indx = 1;

for i = 1:length(netnames)
    
    [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {['LH_' netnames{i}]}, 'flatten');
    
    new_vals(obj_subset.dat == 1) = indx;
    
    new_labels{indx, 1} = obj_subset.labels{1};
    
    indx = indx + 1;
    
    
    [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {['RH_' netnames{i}]}, 'flatten');
    
    new_vals(obj_subset.dat == 1) = indx;
    
    new_labels{indx, 1} = obj_subset.labels{1};
    
    indx = indx + 1;
    
end


new_atlas.dat = single(new_vals);
new_atlas.labels = new_labels;

new_atlas = check_properties(new_atlas);

atlas_obj = new_atlas;
r = atlas2region(atlas_obj);

atlas_obj.history{end+1} = 'Collapsed 400 regions into 17 networks separated by L/R hemisphere.';

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

