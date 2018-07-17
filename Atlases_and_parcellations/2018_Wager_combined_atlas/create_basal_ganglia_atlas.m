atlas_obj = load_atlas('pauli_bg');

% Threshold some to clean up and avoid bleed-over
% atlas_obj = threshold(atlas_obj, .6);

atlas_obj = check_properties(atlas_obj);


%% add regions from CIT168

cit168 = load_atlas('CIT168');

group_codes = {'NAC' 'GPe' 'GPi' 'STN' 'VeP'};

for i = 1:length(group_codes)
    
    if iscell(group_codes{i})
        roi_atlas = select_atlas_subset(cit168, group_codes{i}, 'flatten'); orthviews(roi_atlas)
    else
        roi_atlas = select_atlas_subset(cit168, group_codes(i), 'flatten'); orthviews(roi_atlas)
    end
    
    atlas_obj = merge_atlases(atlas_obj, roi_atlas, 'always_replace');
    
end

%% add caudate tail from CIT168 without replacing existing regions

roi_atlas = select_atlas_subset(cit168, {'Cau'}, 'flatten'); orthviews(roi_atlas)

% cut off: keep only behind caudate tail
% ******

atlas_obj = merge_atlases(atlas_obj, roi_atlas, 'noreplace');


%% Subdivide into L and R hemispheres
atlas_obj = split_atlas_by_hemisphere(atlas_obj); 

%% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique');
 
%% Convert to regions
% -----------------------------------------------------------------------

 r = atlas2region(atlas_obj);

%% Enforce some var types and compress

atlas_obj = check_properties(atlas_obj);
atlas_obj = remove_empty(atlas_obj);


%% Save

targetdir = what('2018_Wager_combined_atlas');
cd(targetdir.path);

savefile = 'Basal_ganglia_combined_atlas_object.mat';
save(savefile, 'atlas_obj');

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
