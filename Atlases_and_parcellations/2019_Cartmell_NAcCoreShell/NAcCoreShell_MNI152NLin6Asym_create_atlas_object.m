% Creates the following:
%
% ***_atlas_object.mat  atlas_obj, atlas object with labels and integer index image
%
% ***_atlas_regions.img image with integer index image
%
% ***_atlas_regions.mat r, region object, and separate region objects for
%       each region with variable names corresponding to labels. For
%       canlab_load_ROI

atlas_name = 'NAcCoreShell_MNI152NLin6Asym';
space_description = 'MNI152NLin6Asym';
references = 'Cartmell SCD, Tian Q, Thio BJ, Leuze C, Ye L, Williams NR, Yang G, Ben-Dor G, Deisseroth K, Grill WM, McNab JA, Halpern CH. Multimodal characterization of the human nucleus accumbens. Neuroimage (2019), 137-149, 198.';

% Run this from the directory containing the atlas files
% -----------------------------------------------------------------------
dosave = true;

% need to make sure we're using the one in MNI space
% This is the full probabilistic atlas file:


% Get labels
% -----------------------------------------------------------------------

labels = {'NAc_core_L','NAc_shell_L','NAc_core_R','NAc_shell_R'};
labels_2 = {'NAc_L','NAc_L','NAc_R','NAc_R'};
label_descriptions = {'Nucleus Accumbens, putative core (left)',...
    'Nucleus Accumbens, putative shell (left)',...
    'Nucleus Accumbens, putative core (right)',...
    'Nucleus Accumbens, putative shell (right)'};
% Create object
% -----------------------------------------------------------------------

pmap = fmri_data({which('lh-presumed-core.nii.gz'),...
    which('lh-presumed-shell.nii.gz'),...
    which('rh-presumed-core.nii.gz'),...
    which('rh-presumed-shell.nii.gz')});
% get rid of negatives, probably interpolation errors
pmap.dat(pmap.dat < 0) = 0;

total_p = sum(pmap.dat,2);
renorm = total_p > 1;
assert(~any(renorm));
%pmap.dat(renorm,:) = pmap.dat(renorm,:)./total_p(renorm);

atlas_obj = atlas(pmap, ...
    'atlas_name', atlas_name, ...
    'labels', labels, ...
    'labels_2',labels_2,...
    'label_descriptions',label_descriptions,...
    'space_description', space_description, ...
    'references', references, 'noverbose');

% Process object
% -----------------------------------------------------------------------

% Threshold at probability 0.2 or greater and k = 3 voxels or greater
atlas_obj = threshold(atlas_obj, .01, 'k', 3);

% Check display
% -----------------------------------------------------------------------

% Display with unique colors for each region:
orthviews(atlas_obj, 'unique','overlay',which('fsl6_hcp_template.nii.gz'));
 
% Convert to regions
% -----------------------------------------------------------------------

r = atlas2region(atlas_obj);

% Display on montage (colors may not be the same!):
% montage(r);
 
%% save figure
cmap = scn_standard_colors(num_regions(atlas_obj));
cmap = cell2mat(cmap')

if dosave
   
    o2 = canlab_results_fmridisplay([], 'full hcp','overlay',which('fsl6_hcp_template.nii.gz'));
    brighten(.6)
    
    o2 = montage(r, o2, 'wh_montages', 1:2, 'indexmap', cmap, 'interp', 'nearest');
    
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
    
    savename = sprintf('%s_atlas_regions.nii', atlas_name);
    atlas_obj.fullpath = fullfile(pwd, savename);
    write(atlas_obj, 'overwrite');
    
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
    
    arrayfun(@(x1)set(x1,'FaceAlpha', .5), han)
    view(135, 20)
    lightFollowView;
    lightRestoreSingle
    axis off
    
    savename = fullfile(savedir, sprintf('%s_isosurface.png', atlas_name));
    saveas(gcf, savename);
    
end